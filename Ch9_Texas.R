## set a working directory 
#setwd("C:/IntroESF")
setwd("D:/proj/IntroESF/Ch9")

## load libraries and functions
library(rgdal)
library(spdep)
library(psych)
source("D:/proj/IntroESF/functions.R")

## read data and transform variables
tx.pop <- read.csv("TX-space-time-population.csv")
tx.poly <- readOGR(dsn=getwd(), layer="TX_counties_spcs", GDAL1_integer64_policy=TRUE)


## socio-economic variable transformation
tx.df <- data.frame(pd=log(tx.poly$Pop_Total / tx.poly$ALAND))
tx.df$fmratio <- tx.poly$Pop_Female / tx.poly$Pop_Male
tx.df$pctw <- tx.poly$Pop_White / tx.poly$Pop_Total
tx.df$pctaa <- log(100 * tx.poly$Pop_Black / tx.poly$Pop_Total + 0.6)
tx.df$pcthisp <- log(100 * tx.poly$Pop_Hisp / tx.poly$Pop_Total)
tx.df$pctccar <- log(100 * tx.poly$Com_Car / tx.poly$Com_SumTra)
tx.df$pctcpub <- log(100 * tx.poly$Com_Public/ tx.poly$Com_SumTra + 0.002)
tx.df$pcthied <- 100 * tx.poly$Edu_HSD_hi / tx.poly$Edu_25Olde
tx.df$pctunemp <- 100 * tx.poly$Uemp_NotLa / tx.poly$Uemp_16Old
tx.df$pcthowner <- 100 * tx.poly$House_Owne / tx.poly$Housing_To
tx.df$pcthvac <- log(100 * tx.poly$House_Vaca / tx.poly$Housing_To)
tx.df$lnincome <- log(tx.poly$Inc_MedHH)


## generate eigenvectors
tx.nb <- poly2nb(tx.poly, queen=FALSE)
n <- length(tx.nb)
B <- nb2mat(tx.nb, style="B")
M <- diag(n) - matrix(1, n, n)/n
MBM <- M %*% B %*% M
eig <- eigen(MBM, symmetric=T)
EV <- as.data.frame(eig$vectors)
colnames(EV) <- paste("EV",1:n, sep="")
EV <- as.data.frame(EV[, eig$values/eig$values[1] > 0.25])


## correlation for 3 variables
# correction among attributes
tx.3 <- tx.df[, c("pd", "pctunemp", "lnincome")]
round(cor(tx.3), 4)			# Table 9.1a
tx3.pca <- principal(tx.3, 3, rotate="none")
tx3.pca$values[1]/tx3.pca$values[3]

# MESF for the 3 variables
lmf.pd <- lm(tx.df$pd ~ ., data=EV)
esf.pd <- stepwise(lmf.pd , lm(tx.df$pd ~ 1, data=EV), verbose=FALSE)
lmf.unemp <- lm(tx.df$pctunemp ~ ., data=EV)
esf.unemp <- stepwise(lmf.unemp, lm(tx.df$pctunemp ~ 1, data=EV), verbose=FALSE)
lmf.inc <- lm(tx.df$lnincome ~ ., data=EV)
esf.inc <- stepwise(lmf.inc, lm(tx.df$lnincome ~ 1, data=EV), verbose=FALSE)

# correlation among aspatial compoents
asp.3 <- cbind(resid(esf.pd) + coef(esf.pd)[1], 
               resid(esf.unemp) + coef(esf.unemp)[1], 
			   resid(esf.inc) + coef(esf.inc)[1])
round(cor(asp.3), 4)		# Table 9.1b
asp3.pca <- principal(asp.3, 3, rotate="none")
asp3.pca$values[1]/asp3.pca$values[3]

# selected EVs
ev.pd <- names(coef(esf.pd)[-1])
ev.unemp <- names(coef(esf.unemp)[-1])
ev.inc <- names(coef(esf.inc)[-1])

# find common EVs
ev.c <- Reduce(intersect, list(ev.pd, ev.unemp, ev.inc))

# run ESF without the common EVs
ev.nc <- setdiff(colnames(EV), ev.c) 
lmf.pd.nc <- lm(tx.df$pd ~ ., data=EV[,ev.nc])
esf.pd.nc <- stepwise(lmf.pd.nc, lm(tx.df$pd ~ 1, data=EV[,ev.nc]), verbose=FALSE)
lmf.unemp.nc <- lm(tx.df$pctunemp ~ ., data=EV[,ev.nc])
esf.unemp.nc <- stepwise(lmf.unemp.nc, lm(tx.df$pctunemp ~ 1, data=EV[,ev.nc]), verbose=FALSE)
lmf.inc.nc <- lm(tx.df$lnincome ~ ., data=EV[,ev.nc])
esf.inc.nc <- stepwise(lmf.inc.nc, lm(tx.df$lnincome ~ 1, data=EV[,ev.nc]), verbose=FALSE)

# selected EVs for the updated ESFs 
ev.pd.nc <- names(coef(esf.pd.nc)[-1])
ev.unemp.nc <- names(coef(esf.unemp.nc)[-1])
ev.inc.nc <- names(coef(esf.inc.nc)[-1])

# unique eigenvectors for the updated ESFs
evu.pd <- setdiff(ev.pd.nc, union(ev.unemp.nc, ev.inc.nc)) 
evu.unemp <- setdiff(ev.unemp.nc, union(ev.pd.nc, ev.inc.nc))
evu.inc <- setdiff(ev.inc.nc, union(ev.pd.nc, ev.unemp.nc))

# correlation among aspatial compoents + unique eigenvectors 
pd.u <- lm(tx.df$pd ~ ., data=EV[,evu.pd])
unemp.u <- lm(tx.df$pctunemp ~ ., data=EV[,evu.unemp])
inc.u <- lm(tx.df$lnincome ~ ., data=EV[,evu.inc])
au3 <- cbind(resid(esf.pd) + fitted(pd.u), 
			 resid(esf.unemp) + fitted(unemp.u), 
			 resid(esf.inc) + fitted(inc.u))
round(cor(au3), 4)		#Table 9.1c
au3.pca <- principal(au3, 3, rotate="none")
au3.pca$values[1]/au3.pca$values[3]

# correlation among aspatial compoents + common eigenvectors 
pd.c <- lm(tx.df$pd ~ ., data=EV[,ev.c])
unemp.c <- lm(tx.df$pctunemp ~ ., data=EV[,ev.c])
inc.c <- lm(tx.df$lnincome ~ ., data=EV[,ev.c])
ac3 <- cbind(resid(esf.pd) + fitted(pd.c), 
			 resid(esf.unemp) + fitted(unemp.c), 
			 resid(esf.inc) + fitted(inc.c))
round(cor(ac3), 4)		#Table 9.1d
ac3.pca <- principal(ac3, 3, rotate="none")
ac3.pca$values[1]/ac3.pca$values[3]


## PCA for Texas Population
library(psych)
tx.lpd <- log(tx.pop[,2:8]/tx.pop$area)
lpd.pca <- principal(tx.lpd, 7, rotate="none")
print(lpd.pca$loadings, cutoff=0.001)	# Table 9.2

n.year <- 7
lpd.esf <- matrix(0, nrow=n, ncol=n.year)
lpd.asp <- matrix(0, nrow=n, ncol=n.year)
colnames(lpd.esf) <- paste("y", 2010:2016, sep="")
colnames(lpd.asp) <- paste("y", 2010:2016, sep="")

for (i in 1:n.year) {
	lm.full <- lm(tx.lpd[,i] ~ ., data=EV)
	sf.i <- stepwise(lm.full, lm(tx.lpd[,i] ~ 1, data=EV), verbose=FALSE)
	lpd.esf[,i] <- sf.i$fitted.values
	lpd.asp[,i] <- coef(sf.i)[1] + resid(sf.i)
}

lpd.esf.pca <- principal(lpd.esf, 7, rotate="none")
print(lpd.esf.pca$loadings, cutoff=0.001)	# Table 9.3 Part I

lpd.asp.pca <- principal(lpd.asp, 7, rotate="none")
print(lpd.asp.pca$loadings, cutoff=0.001)	# Table 9.3 Part II


## MESF for the 12 socio-economic variables
n.var <- dim(tx.df)[2]
esf.mat <- matrix(0, nrow=n, ncol=n.var)
asp.mat <- matrix(0, nrow=n, ncol=n.var)
var.names <- c("pd", "fmratio", "pctw", "pctaa", "pcthisp", "pctccar", 
               "pctcpub", "pcthied", "pctunemp", "pcthowner", "pcthvac", "lnincome")
colnames(esf.mat) <- var.names
colnames(asp.mat) <- var.names

for (i in 1:n.var) {
	lm.full <- lm(tx.df[,i] ~ ., data=EV)
	sf.i <- stepwise(lm.full, lm(tx.df[,i] ~ 1, data=EV), sle=0.1, sls=0.1001, verbose=FALSE)
	esf.mat[,i] <- sf.i$fitted.values
	asp.mat[,i] <- sf.i$coefficients[1] + residuals(sf.i)
}


## PCA 
tx.pca <- principal(tx.df, 7, rotate="none")
print(tx.pca$loadings, cutoff=0.001)	# Table 9.4

esf.pca <- principal(esf.mat, 7, rotate="none")
print(esf.pca$loadings, cutoff=0.001)	# Table 9.5 Part I

asp.pca <- principal(asp.mat, 7, rotate="none")
print(asp.pca$loadings, cutoff=0.001)	# Table 9.5 Part II


## FA 
tx.fa <- principal(tx.df, 6, rotate="varimax")
print(tx.fa$loadings, cutoff=0.001)		# Table 9.6

esf.fa <- principal(esf.mat, 4, rotate="varimax")
print(esf.fa$loadings, cutoff=0.001)	# Table 9.7 Part I

asp.fa <- principal(asp.mat, 7, rotate="varimax")
print(asp.fa$loadings, cutoff=0.001)	# Table 9.7 Part II


## ANOVA & MANOVA
mas.f <- as.factor(tx.poly$MSA)

f.av <- function (x, z){
	anova(lm(x ~ z))[[5]][1]
}

f.lev <- function(x, z){
	leveneTest.square(x, z)[[3]][1]
}

tx.p <- apply(tx.df, 2, FUN=f.av, z=mas.f)
tx.lev <- apply(tx.df, 2, FUN=f.lev, z=mas.f)

esf.p <- apply(esf.mat, 2, FUN=f.av, z=mas.f)
esf.lev <- apply(esf.mat, 2, FUN=f.lev, z=mas.f)

asp.p <- apply(asp.mat, 2, FUN=f.av, z=mas.f)
asp.lev <- apply(asp.mat, 2, FUN=f.lev, z=mas.f)

round(cbind(tx.p, tx.lev, esf.p, esf.lev, asp.p, asp.lev),4)	# Table 9.8

summary(manova(as.matrix(tx.df)~mas.f))
summary(manova(esf.mat~mas.f))
summary(manova(asp.mat~mas.f))


## Discriminant Analysis
library(MASS)
library(klaR)
tx.gw <- greedy.wilks(mas.f ~ ., data = tx.df, niveau = 0.1)
tx.gw 			# Table 9.9a

esf.df <- as.data.frame(esf.mat)
esf.gw <- greedy.wilks(mas.f ~ ., data = esf.df, niveau = 0.1)
esf.gw			# Table 9.9b

asp.df <- as.data.frame(asp.mat)
asp.gw <- greedy.wilks(mas.f ~ ., data = asp.df, niveau = 0.1)
asp.gw			# Table 9.9c

fml.tx <- update(tx.gw$formula, tx.poly$MSA ~ .)
tx.lm <- lm(fml.tx, data=tx.df)
summary(tx.lm)  # Table 9.10a

fml.esf <- update(esf.gw$formula, tx.poly$MSA ~ .)
esf.lm <- lm(fml.esf, data=esf.df)
summary(esf.lm)  # Table 9.10b

fml.asp <- update(asp.gw$formula, tx.poly$MSA ~ .)
asp.lm <- lm(fml.asp, data=asp.df)
summary(asp.lm)  # Table 9.10c


## Canonical Correlation Analysis
library(CCA)
library(CCP)
tx.cc <- cc(EV, tx.df)
cc.r <- tx.cc$cor 	
print(round(cc.r,4)) 	 		# Table 9.11 correlation

cc.ev <- cc.r^2/(1-cc.r^2)
round(cc.ev/sum(cc.ev),4)		# Table 9.11 % of variance

# Calculate p-values with Wilk's lambda
p.asym(cc.r, n, dim(EV)[2], dim(tx.df)[2])		# Table 9.11 Test

round(tx.cc$scores$corr.Y.yscores[,1:6],4) # Table 9.12


## Cluster Analysis
# for the attributes (Table 9.13)
tx.dist <- dist(scale(tx.df))^2
tx.ward <- hclust(tx.dist, method="ward.D2")
tx.ward.g <- cutree(tx.ward, k=4)
table(tx.ward.g)
plot(tx.ward)
rect.hclust(tx.ward, k=4, border="red")

tx.comp <- hclust(tx.dist, method="complete")
tx.comp.g <- cutree(tx.comp, k=4)
table(tx.comp.g)

tx.mq <- hclust(tx.dist, method="mcquitty")
tx.mq.g <- cutree(tx.mq, k=4)
table(tx.mq.g)

# for ESFs (Table 9.13)
esf.dist <- dist(scale(esf.mat))^2

esf.ward <- hclust(esf.dist, method="ward.D2")
esf.ward.g <- cutree(esf.ward, k=4)
table(esf.ward.g)

esf.comp <- hclust(esf.dist, method="complete")
esf.comp.g <- cutree(esf.comp, k=4)
table(esf.comp.g)

esf.mq <- hclust(esf.dist, method="mcquitty")
esf.mq.g <- cutree(esf.mq, k=4)
table(esf.mq.g)

# for unstructured residuals (Table 9.13)
asp.dist <- dist(scale(asp.mat))^2

asp.ward <- hclust(asp.dist, method="ward.D2")
asp.ward.g <- cutree(asp.ward, k=4)
table(asp.ward.g)

asp.comp <- hclust(asp.dist, method="complete")
asp.comp.g <- cutree(asp.comp, k=4)
table(asp.comp.g)

asp.mq <- hclust(asp.dist, method="mcquitty")
asp.mq.g <- cutree(asp.mq, k=4)
table(asp.mq.g)

# with coordinates (Table 9.14)
c.x <- as.numeric(as.vector(tx.poly$INTPTLON)) 
c.y <- as.numeric(as.vector(tx.poly$INTPTLAT))

# MESF for the longitudes and latitudes 
lm.full <- lm(c.x ~ ., data=EV)
sf.x <- stepwise(lm.full, lm(c.x ~ 1, data=EV), sle=0.1, sls=0.1001, verbose=FALSE)
esf.x <- sf.x$fitted.values
asp.x <- sf.x$coefficients[1] + residuals(sf.x)

lm.full <- lm(c.y ~ ., data=EV)
sf.y <- stepwise(lm.full, lm(c.y ~ 1, data=EV), sle=0.1, sls=0.1001, verbose=FALSE)
esf.y <- sf.y$fitted.values
asp.y <- sf.y$coefficients[1] + residuals(sf.y)

txc.df <- cbind(tx.df, c.x, c.y)
txc.dist <- dist(scale(txc.df))^2

txc.ward <- hclust(txc.dist, method="ward.D2")
txc.ward.g <- cutree(txc.ward, k=4)
table(txc.ward.g)		# Table 9.14, column 1

esfc.df <- cbind(esf.df, esf.x, esf.y)
esfc.dist <- dist(scale(esfc.df))^2

esfc.ward <- hclust(esfc.dist, method="ward.D2")
esfc.ward.g <- cutree(esfc.ward, k=4)
table(esfc.ward.g)		# Table 9.14, column 2

aspc.df <-  cbind(asp.df, asp.x, asp.y)
aspc.dist <- dist(scale(aspc.df))^2

aspc.ward <- hclust(aspc.dist, method="ward.D2")
aspc.ward.g <- cutree(aspc.ward, k=4)
table(aspc.ward.g)		# Table 9.14, column 3
