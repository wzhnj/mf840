install.packages("expm", repos="http://R-Forge.R-project.org")   
library(expm)
library(matlib)
install.packages("fGarch")
library(sandwich)	
library(lmtest)
install.packages("tseries")

data<-read.csv('sizebtm-week.csv')
data<-data/100

model<-lm((smbm3-RF) ~ Mkt.RF , data)

#OLS
coeftest(model)

# BOTH Heteroskedasticity and Autocorrelation
coeftest(model,vcov=NeweyWest) # Newey and West 
coeftest(model,vcov=vcovHAC)	  # Andrews (1991)

require(forecast)			# A  better ACF plot
par(mfrow=c(1,2),mgp=c(1.5,0.5,0))
Acf(model$residuals,main="");title(line=0.5,"Full Acf")
Acf(model$residuals,main="",type="partial");title(line=0.5,main="Partial Acf")
title(outer=T,"Figure1: ACF of OLS regression residuals",line=-2)

acf(model$residuals,plot=F)$acf[1:5]

olsres=model$residuals
sdres <- sd(olsres)*sqrt(521/522)
cormat<- toeplitz(c(1,acf(olsres,plot=F)$acf[2:4],rep(0,518)))
dmat  <- diag(rep(sdres,522))
omga  <- dmat%*% cormat %*% dmat




ev <- eigen(omga)
# extract components
(L <- ev$values)
(p <- ev$vectors)
d=diag(L)
p <- t(p)

dnzerofive=solve(diag(sqrt(L)))
a=dnzerofive %*% d %*% dnzerofive
a[1:5,1:5]
pstar <- dnzerofive %*% p
a=pstar %*% omga %*% t(pstar)
modelgls<-lm(pstar %*% (smbm3-RF) ~ pstar %*% Mkt.RF , data)
coeftest(modelgls)

record_beta=coeftest(modelgls)[2]
difference=1
coeficientalpha=c()
coeficientbeta=c()
std=c()
stdalpha=c()

while (difference>0.01){
count=count+1
olsres=modelgls$residuals
sdres <- sd(olsres)*sqrt(521/522)
cormat<- toeplitz(c(1,acf(olsres,plot=F)$acf[2:4],rep(0,518)))
dmat  <- diag(rep(sdres,522))
omga  <- dmat%*% cormat %*% dmat

ev <- eigen(omga)
# extract components
L <- ev$values
p <- ev$vectors
d=diag(L)
p <- t(p)

dnzerofive=solve(diag(sqrt(L)))
a=dnzerofive %*% d %*% dnzerofive
a[1:5,1:5]
pstar <- dnzerofive %*% p
a=pstar %*% omga %*% t(pstar)
modelgls<-lm(pstar %*% (smbm3-RF) ~ pstar %*% Mkt.RF , data)
coeftest(modelgls)

difference=abs(record_beta-coeftest(modelgls)[2])
std=append(std,coeftest(modelgls)[4])
stdalpha=append(stdalpha,coeftest(modelgls)[3])

record_beta=coeftest(modelgls)[2]
coeficientbeta=append(coeficientbeta,record_beta)
coeficientalpha=append(coeficientalpha,coeftest(modelgls)[1])

}
coeficientalpha
coeficientbeta
std
stdalpha

par(mfrow=c(2,1),mgp=c(1.5,0.5,0),mar=c(3,3,2,0.5),oma=c(0,0,2,0))
#plot(data[,4],abs(olsres),xlab="RmRf",ylab="|OLS Residual|")
#title(line=0.2,"OLS residual: absolute value vs VW-Xrm")
Acf(abs(olsres),main="")
qqnorm(olsres,main="");qqline(olsres)
title(line=0.2,"OLS Residuals: Normal Probability plot")
title(outer=T,"Figure 2: Heteroskedasticity Diagnostic",line=0)

library(fGarch)
library(tseries)
garchols<-garch(olsres)
summary(garchols)

plot(garchols)

require(fGarch)
garch2<-garchFit(formula = ~garch(1,1),data = olsres,include.mean=FALSE)
summary(garch2)

par(mfrow=c(1,1),mgp=c(1.5,0.5,0),mar=c(3,3,1,0.3),oma=c(0,0,1,0))
ts.plot(garchols$fitted.values[,1],ylab="GARCH(1,1) Std. Dev.")
mtext(outer=T,"Figure 3: GARCH(1,1) Diagnostics",line=0)


omgainv <- diag(1/garch2@h.t)
xx    <- cbind(rep(1,522),data[,8])
xomx  <- solve(t(xx)%*% omgainv %*% xx)
betgarch <- xomx %*% (t(xx) %*% omgainv %*% (data[,4]-data[,11]))
betgarch
sqrt(diag(xomx))

# Iterating

xx     <- cbind(rep(1,522),data[,8])

betiter<-matrix(0,ncol=2,nrow=100)
res0	   <- olsres
betiter[1,]<-model$coef
i      <- 1
zcond  <- 100		# Careful if you copy paste several times!
betiter<-matrix(0,ncol=2,nrow=100)
bstd<-matrix(0,ncol=2,nrow=100)

bstd[1,]=coeftest(model)[3:4]



while (zcond > 0.000001/100) { 
i 	   <- i+1 
omgainv <- diag(1/ garchFit(formula=~garch(1,1), data=res0, include.mean=FALSE,trace=FALSE)@h.t)

garch2<-garchFit(formula=~garch(1,1), data=res0, include.mean=FALSE,trace=FALSE)
summary(garch2)

xomx    <- solve(t(xx)%*% omgainv %*% xx) 
betiter[i,] <- xomx %*% (t(xx) %*% omgainv %*% (data[,4]-data[,11])) 
bstd[i,] <- sqrt(diag(xomx))

zcond  <- abs((betiter[i,2]-betiter[i-1,2])/betiter[i-1,2])
res0   <- (data[,4]-data[,11]) - xx %*% betiter[i,]

}

zcond 
betiter[1:(i+1),]

bstd[1:(i+1),]
sqrt(diag(xomx))








