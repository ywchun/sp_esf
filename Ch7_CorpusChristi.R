
## set a working directory 
setwd("C:/IntroESF") 

## load libraries and functions
library(rgdal)
library(spdep)
library(MASS)
source("functions.R")

cc.mig <- read.csv("CorpusChristi_MSA-migration_data.csv", header=TRUE)
cc.mig$lnOi <- log(cc.mig$oi)
cc.mig$lnDj <- log(cc.mig$dj)
cc.mig$offs <- cc.mig$lnOi + cc.mig$lnDj

cc.nls <- nls(flows ~ exp(k + a*lnOi + b*lnDj + g*dist), data=cc.mig, 
                      start=list(k=-1, a=1, b=1, g=-0.1))
summary(cc.nls)

cc.poi <- glm(flows ~ lnOi + lnDj + dist, family=quasipoisson, data=cc.mig)
summary(cc.poi)

cc.dc <- glm(flows ~ io1 + io2 + io3 + io4 + io5 + 
                     id1 + id2 + id3 + id4 + id5 + dist, 
					 offset=offs, family=quasipoisson, data=cc.mig)
summary(cc.dc)

fb <- grepl("io", names(coefficients(cc.dc)))
lnai <- c(coefficients(cc.dc)[fb], 0)
tb <- grepl("id", names(coefficients(cc.dc)))
lnbj <- c(coefficients(cc.dc)[tb], 0)

cc.poly <- readOGR(dsn=getwd(), layer="CorpusChristi_MSA")
cc.nb <- poly2nb(cc.poly)

## generate eigenvectors
n <- length(cc.nb)
B <- nb2mat(cc.nb, style="B")
M <- diag(n) - matrix(1, n, n)/n
MBM <- M %*% B %*% M
eig <- eigen(MBM, symmetric=T)
EV <- as.data.frame(eig$vectors)
colnames(EV) <- paste("EV",1:n, sep="")
EV <- as.data.frame(EV[, round(eig$values, 4) != 0])

## run stepwise regression for ESF
lm.ai.init <- lm(lnai ~ 1, data=EV)
lm.ai.full <- lm(lnai ~ ., data=EV)
ai.esf <- stepwise(lm.ai.full, lm.ai.init, sle=0.1, sls=0.1001) # EV5 & EV1 are selected
summary(ai.esf)$r.square 

lm.bj.init <- lm(lnbj ~ 1, data=EV)
lm.bj.full <- lm(lnbj ~ ., data=EV)
bj.esf <- stepwise(lm.bj.full, lm.bj.init, sle=0.1, sls=0.1001) # no EV is selected

cEV <- kronecker(as.matrix(EV), as.matrix(EV))
#id.ord <- as.vector(outer(seq(1,21,5), 0:4,  FUN="+"))
id.ord <- 1:NCOL(cEV)
colnames(cEV) <- paste("EV", id.ord, sep="")

cc.esf.df <- cbind(cc.mig[,c(4, 7:17, 20)], cEV)

glm.init <- glm(flows ~ io1 + io2 + io3 + io4 + io5 + 
            id1 + id2 + id3 + id4 + id5 + dist, 
			offset=offs, family=poisson, data=cc.esf.df)
glm.full <- glm(flows ~ . - offs, 
			offset=offs, family=poisson, data=cc.esf.df)
sel.scl <- sqrt(glm.init$deviance/glm.init$df.residual)
cc.dc.na <- step(glm.init, scope=list(upper=glm.full, lower=glm.init), scale=sel.scl, trace=0)
summary(cc.dc.na)

