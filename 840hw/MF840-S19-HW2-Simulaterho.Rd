#
# Uncertainty in cross-correlation estimate
# Checking the asymptotics:	sd(\hat(\rho))=1/sqrt(T)
#							Asympt. approx only correct if true rho = 0

# MV normal simulation libraries
#	MASS				command mvrnorm
#	mnormt			rmnorm
#					Also has commands to compute the density for normal and t
#	mvtnorm			rmvnorm	(also has t)	

library(MASS)

# True Cov

rho <- 0.4		
Dim <- 10		# Dim*(Dim-1)/2 correlations
TT  <- 100

muvec<-rep(0,Dim)
Sigma<-diag(rep(1-rho,Dim)) + matrix(rho,ncol=Dim,nrow=Dim)

nsim <-1000
cormu<-rep(0,nsim);corsd<-cormu

rdata<-matrix(rnorm(TT*Dim*nsim,0,1),ncol=Dim,nrow=TT*nsim)	# Use rnorm if rho=0
rdata<-mvrnorm(TT*nsim,muvec,Sigma)
for (ii in 1:nsim){ 
	zdata<-rdata[(TT*(ii-1)+1):(TT*ii),]
	cors<- unique(c(cor(zdata)))
	cors<- cors[cors<1]
	cormu[ii]<-mean(cors)	
	corsd[ii]<-sd(cors)
	}

quantile(cormu)
mean(corsd)
sd(corsd)
sqrt(1/TT)
	
	
	
	
