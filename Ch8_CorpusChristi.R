
## set a working directory 
setwd("C:/IntroESF")

## load libraries and functions
library(rgdal)
library(spdep)
source("functions.R")

## read data
ccm <- read.csv("CorpusChristi_MSA-data.csv")
ccm <- ccm[order(ccm$id, ccm$year), ]
ccm$area <- ccm$area * 100000

cc.poly <- readOGR(dsn=getwd(), layer="CorpusChristi_MSA")
cc.nb <- poly2nb(cc.poly)

## run a GLMM
library(lme4)
library(optimx)
poi.mix <- glmer(population ~ 1 + (1|county), offset=log(area),
                family = poisson, data = ccm, start=list(fixef=0),
				control=glmerControl(optimizer="optimx", optCtrl=list(method="nlminb")))
summary(poi.mix)
rande <- as.numeric(ranef(poi.mix)$county[,1]) # get random effects
rande

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
lm.init <- lm(rande ~ 1, data=EV)
lm.full <- lm(rande ~ ., data=EV)
esf.ccm <- stepwise(lm.full, lm.init, sle=0.15, sls=0.1501, verbose=T) 

ssre <- fitted.values(esf.ccm)
sure <- residuals(esf.ccm)

## map the random effects
mapping.bipolar(cc.poly, ssre, 4, legend.loc="bottomright")
mapping.bipolar(cc.poly, sure, 4, legend.loc="bottomright")




