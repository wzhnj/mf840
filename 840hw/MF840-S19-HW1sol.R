#
# GLS
#

rets<-read.csv("size-day-0918.csv")
rets<-rets[rets[,1]<20120000,]/100
dim(rets)

# OLS 

olsreg<-lm(rets[,3]~rets[,12])
summary(olsreg)
olsres<-olsreg$residuals

library(sandwich);library(lmtest)
coeftest(olsreg,vcov=vcovHC)

# Diagnostic of autocorrelation in the OLS errors
# Use Acf, note also the command pacf

library(forecast)
par(mfrow=c(1,2),mgp=c(1.5,0.5,0))
Acf(olsres,main="");title(line=0.5,"Full Acf")
Acf(olsres,main="",type="partial");title(line=0.5,main="Partial Acf")
title(outer=T,"Figure1: ACF of OLS regression residuals, D2 on XRm, 2009-2012",line=-2)

acf(olsres,plot=F)$acf[1:5]
# preferred model is a zip-0

####################################################
# subtracting interest rate from Ri does not change 
# anything for this period.
facs<-read.csv("FF-3fac-day.csv")
names(facs)
mean(facs$rf[facs$date>20090000&facs$date<20120000])*252/100
# yes 10bps annualized!
rd2rf<-rets[,3]-facs$rf[facs$date>20090000&facs$date<20120000]/100
reg2<-lm(rd2rf~rets[,12])
summary(reg2)
####################################################


# MA(3) GLS
# Omega will be 756x756, 
# We apply the correlations estimated to the first 3 observations as well

sdres <- sd(olsres)*sqrt(755/754)
cormat<- toeplitz(c(1,acf(olsres,plot=F)$acf[2:4],rep(0,752)))
dmat  <- diag(rep(sdres,756))
omga  <- dmat%*% cormat %*% dmat
ominv <- solve(omga)
xx    <- cbind(rep(1,756),rets[,12])
xomx  <- solve(t(xx)%*% ominv %*% xx)

betga <- xomx %*% (t(xx) %*% ominv %*% rets[,3])
betga

sqrt(diag(xomx))

#
# Diagnostics of Heteroskedasticity in the OLS errors
#
par(mfrow=c(3,1),mgp=c(1.5,0.5,0),mar=c(3,3,2,0.5),oma=c(0,0,2,0))
plot(rets[,12],abs(olsres),xlab="RmRf",ylab="|OLS Residual|")
title(line=0.2,"OLS residual: absolute value vs VW-Xrm")
Acf(abs(olsres),main="")
title(line=0.2,"OLS residuals: ACF of absolute values")
qqnorm(olsres,main="");qqline(olsres)
title(line=0.2,"OLS Residuals: Normal Probability plot")
title(outer=T,"Figure 2: Heteroskedasticity Diagnostic, Daily D2 on XRm, 2009-12",line=0)


# Estimate GARCH(1,1)

# with the basic garch command (in tseries)
garchols<-garch(olsres)
summary(garchols)

plot(garchols)	# Primitive canned graphical diagnostic, 
				# some naive plots, covers some basics

# with garchFit (in fGarch)
# garchFit estimates initial observation variances

require(fGarch)
garch2<-garchFit(formula = ~garch(1,1),data = olsres,include.mean=FALSE)
summary(garch2)


# Plot the diagnostic

par(mfrow=c(2,1),mgp=c(1.5,0.5,0),mar=c(3,3,1,0.3),oma=c(0,0,1,0))
qqnorm(garchols$residuals,ylab="GARCH(1,1) stdized. residual",main="")
title(line=-1, "Normal Probability Plot")
qqline(garchols$residuals)
ts.plot(garchols$fitted.values[,1],ylab="GARCH(1,1) Std. Dev.")
mtext(outer=T,"Figure 3: GARCH(1,1) Diagnostics, Residuals of Daily D2 on XRm, 2009-12",line=0)

#
# Do the GARCH - GLS
#

omgainv <- diag(1/garch2@h.t)
xx    <- cbind(rep(1,756),rets[,12])
xomx  <- solve(t(xx)%*% omgainv %*% xx)
betgarch <- xomx %*% (t(xx) %*% omgainv %*% rets[,3])
betgarch
sqrt(diag(xomx))


# Iterating

xx     <- cbind(rep(1,756),rets[,12])

betiter<-matrix(0,ncol=2,nrow=100)
res0	   <- olsres
betiter[1,]<-olsreg$coef
i      <- 1
zcond  <- 100		# Careful if you copy paste several times!

while (zcond > 0.000001/100) { 
i 	   <- i+1 
omgainv <- diag(1/ garchFit(formula=~garch(1,1), data=res0, include.mean=FALSE,trace=FALSE)@h.t)
xomx    <- solve(t(xx)%*% omgainv %*% xx) 
betiter[i,] <- xomx %*% (t(xx) %*% omgainv %*% rets[,3]) 
zcond  <- abs((betiter[i,2]-betiter[i-1,2])/betiter[i-1,2])
res0   <- rets[,2] - xx %*% betiter[i,]
}

zcond; betiter[1:(i+1),]; sqrt(diag(xomx))





