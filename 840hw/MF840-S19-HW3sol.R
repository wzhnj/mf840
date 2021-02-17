#
#		Homework 3 Solutions
#

# Problem 1: Odds ratios

vix<-read.csv("vix-sp-week.csv",header=T)
vixdates<-as.Date(vix$date,"%m/%d/%y")


# a) We assume a two-sided test

TT    <- 433-366+1
muhat <- mean(vix[366:433,2]);muhat
sdhat <- sd(vix[366:433,2]);sdhat
zz    <- (muhat-15)/(sdhat/sqrt(TT));zz	# 2.11

(1-pnorm(zz,0,1))/2		# 0.009<1%


# b) Posteriors with diffuse prior on sigma
#		prior p(mu|sigma)~N(mu0,sigma^2/T0)  
# The term in exponential of p(mu,sigma|D)
# is in the top of Page 14 (Bayesintro4) inside the integral

# Careful, we did not add p(\sigma)\propto 1/sigma in the
# note because section 6.5 conditions on sigma.

mu0  <- 15			# prior parameters for mu|sigma
T0   <- TT/10

# Mu | sigma, D ~ M( mupos, sig^2/(TT+T0))
mupos<- (T0*mu0+TT*muhat)/(TT+T0)
mupos
TT+T0

# OLS sum of squares, the OLS nu*s^2
ssq <-sum((vix[366:433,2]-muhat)^2) # OLS sum of squares

# posterior sum of squares, add the Do-Not-Forgets
ssqpos<-ssq + TT*muhat^2+T0*mu0^2-(TT+T0)*mupos^2

# posterior degrees of freedom 
# (from the powers of 1/sigma: T, 1, 1
#  one of them goes to mu|sigma,D~N
#  Don't forget the prior 1/sigma

nupos <- TT		

# posterior s parameter

spos <- sqrt(ssqpos/nupos);spos

# mu | D is Student-t with nupos degrees of freedom
# centered on mupo
# precision (T+T0)/spos

# c) Section 6.5 is conditional on sigma

BF01<-sqrt(1+TT/T0)*exp(-0.5*zz^2/(1+T0/TT));BF01
P0<-BF01/(1+BF01);P0

# d) Savage Density ratio method
# Need a proper prior IG for sigma
# Let's have it vague and centered on 5

nu0<- 3
s0 <- 5

#  Prior ordinate at mu=15
#  mu|sigma ~N(0,sigma^2/T0), sigma ~ ITg (nu0=3, s0=4)
#  mu ~ T(nu0=3,  mu0=15   ,  s0^2/T0 )
#  See AZ, density for the univariate T.

hh0   <- T0/s0^2
priord<- sqrt(hh0/nu0) * (1/beta(nu0/2,1/2)) / 
	(1+ hh0*(mu0-mu0)^2/nu0)^((nu0+1)/2)
priord

# Posterior ordinate
# posterior sum of squares in kernel is now:

ssqpos<-nu0*s0^2 + ssq + TT*muhat^2 + T0*mu0^2 -(TT+T0)*mupos^2 

nu1  <- TT+nu0		# DOF of the ItG p(sigma|D) and the t p(mu|D)

# Recall that ssqpos = nupos * spos	by definition. So you
# can fund spos if you want but we don't need it here

hovnupos <- (TT+T0)/ssqpos
posd <- sqrt(hovnupos) * (1/beta(nu1/2,1/2)) / 
	(1+ hovnupos*(mu0-mupos)^2)^((nu1+1)/2)

SBF01<-posd/priord;SBF01
P0   <-SBF01/(1+SBF01)

# e) sensitivity analysis
#

T0	  <-seq(1,20)
hh0   <- T0/s0^2
priord<- sqrt(hh0/nu0) * (1/beta(nu0/2,1/2)) / 
	(1+ hh0*(mu0-mu0)^2/nu0)^((nu0+1)/2)
ssqpos<-nu0*s0^2 + ssq + TT*muhat^2 + T0*mu0^2 -(TT+T0)*mupos^2 
hovnupos <- (TT+T0)/ssqpos
posd  <- sqrt(hovnupos) * (1/beta(nu1/2,1/2)) / 
	(1+ hovnupos*(mu0-mupos)^2)^((nu1+1)/2)
SBF01 <-posd/priord
P0    <-SBF01/(1+SBF01)

plot(T0,P0,type="l")
abline(h=0.5)


#
# PROBLEM 2			##################################
#

mcavg<-0.67		# MC mean
mcstd<-0.35		# MC std. err

# Want a 95% interval to not affect the 7 
# Choose 0.0045 which is less than 0.005
# 0.0045 < 1.96 * 0.35/sqrt(S)  where S is the MC sample size

(1.96*0.35/0.0045)^2		# 23240 draws

# var(sample mean) > sig^2/S
# Look at the variance under autocorrelation eqn.(1) P4, in the bracket
# [1+ 2 Sum(rho_i) - (2/S) Sum(i*rho_i)]
# With Large S, and only a few rhos non-zero or rhos go quickly to zero
# it is approximately   [1 + 2 Sum(rho_i)]

varinflation <- 1 + 2*(0.7+0.4)
1/varinflation		# is our RNE

#  We now need

varinflation*(1.96*0.35/0.0045)^2 	# 67860 draws

#################################
# Problem 3 
#
# Get the sample statistics

#106:418 Estimation data

ivec <-seq(106,418)
TT    <- length(ivec)
ydat  <- log(vix[ivec,2])
xmat  <- cbind(rep(1,TT),log(vix[ivec-1,2]))
xpxinv<- solve(t(xmat)%*%xmat)
bhat  <- xpxinv %*% t(xmat) %*% ydat;bhat
nupos <- TT-2
ssqpos<- sum((ydat-xmat%*% bhat)^2)
spos  <- sqrt(ssqpos/nupos);spos

#Drawing from the inverted gamma
# ItG ~ sqrt( nu*s^2/ chiq*=(nu)  )
# or get alpha/beta, draw the gamma, and take 1/sqrt(Gamma)
#

nsim<- 10000
sigdraw<-sqrt(ssqpos/rchisq(nsim,nupos))

# Or
alpha <- nupos/2 
scale <- 2 / ssqpos 
sigdraw<- 1 / sqrt(rgamma(nsim,alpha,scale=scale))

par(mgp=c(1.5,0.5,0))
hist(sigdraw,prob=T,nclass=100);box()
lines(density(sigdraw),col="red",lwd=2)

# Now Betas ~ N(betapos, sigma (X'X)-1)

library(mnormt)
help(rmnorm)

betadraw<-rmnorm(nsim,varcov=xpxinv)
betadraw<-betadraw * sigdraw
betadraw<- t( t(betadraw) + drop(bhat))

hist(betadraw[,2],nclass=100,prob=T);box()
lines(density(betadraw[,2]),col="red",lwd=2)

# Now the epsilons
etp1<-rnorm(nsim,0,sigdraw)
etp2<-rnorm(nsim,0,sigdraw)

lvtp1<-betadraw %*% c(1,log(vix[418,2])) + etp1
lvtp2<-betadraw[,1] + betadraw[,2]*lvtp1 + etp2 

vixtp1<-exp(lvtp1)
vixtp2<-exp(lvtp2)
mu1<-mean(vixtp1);sd1<-sd(vixtp1);qq1<-quantile(vixtp1,c(0.25,0.5,0.75))
mu2<-mean(vixtp2);sd2<-sd(vixtp2);qq2<-quantile(vixtp2,c(0.25,0.5,0.75))

par(mgp=c(1.5,0.5,0),mfrow=c(1,2),mar=c(2.5,2.5,3,0.5),oma=c(0,0,0,0))
plot(density(vixtp1),main="",xlab="VIX");abline(v=vix[419,2])
title("1 week ahead VIX",cex=0.9,line=0.2)
plot(density(vixtp2),main="",xlab="VIX");abline(v=vix[420,2])
title("2 week ahead VIX",cex=0.9,line=0.2)
mtext(outer=T,"Figure 2",line=-1,font=2)

ivecl<-seq(384,420)
ivec <-seq(384,418)
par(mgp=c(1.5,0.5,0),mfrow=c(1,1))
plot(vixdates[ivecl],vix$vix[ivecl],type="p",xaxt="n",
ylab="VIX",pch="*",xlab="") 
axis.Date(1,at=seq(vixdates[364],vixdates[420],"week"),format="%m-%d",
cex.axis=0.7)
lines(vixdates[ivec],vix$vix[ivec])
lines(vixdates[419:420],c(mean(vixtp1),mean(vixtp2)))
lines(vixdates[419:420],c(qq1[3],qq2[3]),col="red")
lines(vixdates[419:420],c(qq1[1],qq2[1]),col="red")
title("Figure 1: VIX - 2 week-ahead prediction from AR(1)",line=0.2)

