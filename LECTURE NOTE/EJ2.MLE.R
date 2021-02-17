##########################################
# How good is the Normality Assumption?
# Easy to check
# (You can also use the widget from MF793)

# Example for simple variance estimation

tt <- 1000 
kk <- 4 
truvar<-0.2^2		

# Choose a range of variances to plot
# Use the exact (small sample) distribution - we know it!
# nu * s^2 / var ~ chisq(tt-kk)

varlow <- truvar * qchisq(0.0001,tt-kk) / tt
varhigh<- truvar * qchisq(0.9999,tt-kk) / tt
varlow;varhigh

# Now pick points to compute densities
xvalues<-seq(varlow, varhigh ,length=400) 

# What would you do if you did not know the true distribution ?
# You would simulate the estimator by simulating the data and 
# estimating the parameter for simulated data.


# Exact Distribution
plot(xvalues,dchisq(xvalues*tt/truvar,tt-4)*tt/truvar, type="l",
main=paste("Sample Size:",tt,", True Var:",truvar), 
ylab="PDF",xlab="Estimate")

abline(v=truvar*(tt-kk)/tt,lwd=2,col="blue") # Exact Mean

# Approximate Normal from MLE theory

lines(xvalues,dnorm(xvalues,truvar, truvar*sqrt(2/tt)),lwd=2,col="red") 



##########################################
# What does the Likelihood function look like
# Let us compute it and plot it
#

# For our sample, the data are now numbers
# We can write the likelihod of the parameters
#  ...  for any parameter we choose

# This function computes the log-likelihood of mean,stdev
# for a given i.i.d normal sample yy
# ####################

normlogl<-function(yy,mean,stdev)
{
tt<-length(yy)
logl<- -0.5*sum((yy-mean)^2)/stdev^2-tt*log(stdev)
logl
}

######################
#
# Simulate some data
# Simple i.i.d. normal model
#

trumu  <- 0.1
trusig <- 0.2

TT	   <- 800
rdata  <-rnorm(TT,trumu,trusig)
muhat  <-mean(rdata)
sighat <-sd(rdata)*sqrt((TT-1)/TT)

normlogl(rdata,muhat,sighat)

# The likelihood can be huge or small
# We did not put the 1/sqrt(2*pi)^TT

normlogl(rdata,muhat,sighat) - (TT/2)*log(2*pi)
 
#
# Let's plot the likelihood for 
# a grid  (25x25) of values
# of the mean ...

muvec<-seq(muhat-2.5*sighat/sqrt(TT),muhat+2.5*sighat/sqrt(TT),length=25)

#  ... and the standard deviation

sigvec<-seq(sighat-2.5*sighat/(2*sqrt(TT)),sighat+2.5*sighat/(2*sqrt(TT)), length=25)

#
# put all these computed values in a matrix zz
#
zz<-matrix(0,ncol=25,nrow=25)
for (i in 1:25){	for (j in 1:25){
	zz[i,j]<-normlogl(rdata,muvec[i],sigvec[j])
}}


# And now plot the log-likelihood as a function of (mu,sigma)
# rotating a perspective requires rgl package

persp(muvec,sigvec,exp(zz),col="blue",xlab="mu",ylab="sig",theta=45,phi=35,
zlab="Lhood")
# oops what happened?

persp(muvec,sigvec,exp(zz-max(zz)),col="blue",xlab="mu",ylab="sig",
zlab="Lhood",theta=45,phi=35,ticktype="detailed")
# trick 

persp(muvec,sigvec,zz,col="blue",xlab="mu",ylab="sig",theta=45,phi=35,zlab="Logl")

# What is the shape of the likelihood vs the log-likelihood ?
# Log concavity is a nice feature of functions for optimization

library(rgl)
persp3d(muvec,sigvec,zz,col="blue",xlab="mu",ylab="sig",zlab="Logl")
muhat;sighat

contour(muvec,sigvec,zz,nlevels=25,xlab="mu",ylab="sig")
abline(h=trusig); abline(v=trumu)
points(muhat,sighat,pch=7,col="blue")

# Where is the likelihood the highest, at (muhat,sighat)
# or at (True mu, True sigma)?


