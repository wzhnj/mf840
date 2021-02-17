#
# Problem 1 

indus<-read.csv("16indus-mon.csv") 

# Check ranges, names of variables, etc 
names(indus);dim(indus);range(indus[,1]) 

rets <-indus 
NN	 <- length(indus[1,])-3;NN 
rets[,2:(NN+1)]<-rets[,2:(NN+1)]-rets[,(NN+3)] 

# Not asked: Mean Variance efficient Set

par(mgp=c(1.5,0.5,0),mfrow=c(1,2),oma=c(0,0,1,0.5),mar=c(3,3,2,0))

ybeg<-2005; yend<-2018
zret<-rets[rets[,1]> ybeg*100 & rets[,1]< (yend+1)*100,]
zm  <-apply(zret[,2:(NN+2)],2,mean)*12
zs  <-apply(zret[,2:(NN+2)],2,sd)*sqrt(12)
zm;zs 
plot(zs,zm,xlim=range(c(0,30)),ylim=range(c(0,13))) 
title(paste(bquote(.(ybeg)),"-",bquote(.(yend))),line=0.1)
points(zs[NN+1],zm[NN+1],pch="M")
abline(0,zm[NN+1]/zs[NN+1]) 

title(outer=T,"Mean/Variance plot: Industries and Market (not asked)",line=0.1)

# Filling Table 1

rets<- rets[rets[,1]>200900&rets[,1]<201400,]
TT    <-length(rets[,1]);TT
mregs <-lsfit(rets[,NN+2],rets[,2:(NN+1)]) 
alphs <-mregs$coefficients[1,];alphs*12 
betas <-mregs$coefficients[2,];betas 
xdat  <-cbind(rep(1,length(rets[,1])),rets[,NN+2]) 
xpxinv<-solve(t(xdat)%*%xdat) 
Omga  <-cov(mregs$residuals)%x%xpxinv 		# See the Kronecker product

covbet  <-Omga[seq(2,2*NN,by=2),seq(2,2*NN,by=2)]
covalph <-Omga[seq(1,2*NN-1,by=2),seq(1,2*NN-1,by=2)]
sigbetas<-sqrt(diag(covbet)*(TT-1)/(TT-2))
sigbetas;quantile(sigbetas)

write.csv(
round(cbind(
12*alphs,  alphs/sqrt(diag(covalph)*(TT-1)/(TT-2)),
betas, (betas-1)/sigbetas
),2)
,file="table1.csv")

# t cut-off
qt(0.05/2,60-2)

# Bonferonni cutoff
0.05/NN
qt((0.05/NN)/2,TT-2)

# Hotelling cut-off
NN<-16
sqrt((TT-1)*NN*qf(0.95,NN,TT-NN)/(TT-NN))
sqrt(qchisq(0.95,NN))

#
# Problem 2
#
indus<-read.csv("30indus-mon.csv")

# Check ranges, names of variables, etc

names(indus);dim(indus);range(indus[,1])

rets <-indus
NN	<- length(indus[1,])-3;NN
rets[,2:(NN+1)]<-rets[,2:(NN+1)]-rets[,(NN+3)]


# Let's use lm this time

rets  <- rets[rets[,1]>200900&rets[,1]<201400,]
TT    <- length(rets[,1]);TT
mregs <- lm(as.matrix(rets[,2:(NN+1)])~rets[,NN+2]) # csv creates lists
													# lm wants matrices!
alphs <-coef(mregs)[1,];alphs*12 
betas <-coef(mregs)[2,];betas 

# The R-squares are one of the variables in ls.print(mregs)$summary

rsquares<-lapply(summary(mregs),'[','r.squared')
quantile(unlist(rsquares));mean(unlist(rsquares))

# The standard errors of slopes are the fourth in the coefficients list
# in summary.lm

stdbetas<-sapply(lapply(summary(mregs),coefficients),'[',4)
mean(unlist(stdbetas)) 

# Correlations of industries and residuals

zcors<-unique(sort(cor(rets[,2:(NN+1)])))
zcors<-zcors[zcors<1]
mean(zcors)


rescors<-unique(sort(cor(mregs$residuals)));rescors
rescors<-rescors[rescors<1];rescors
mean(rescors);sd(rescors);1/sqrt(TT)
mean(abs(cors))


# Hotelling T2 test
# We need the SUR framewok to get the 30x30 Covmat of the betas_MLE   

xdat  <-cbind(rep(1,length(rets[,1])),rets[,NN+2]) 
xpxinv<-solve(t(xdat)%*%xdat) 
Omga  <-cov(mregs$residuals)%x%xpxinv 		# See the Kronecker product

covbet  <-Omga[seq(2,2*NN,by=2),seq(2,2*NN,by=2)]
covalph <-Omga[seq(1,2*NN-1,by=2),seq(1,2*NN-1,by=2)]
sigbetas<-sqrt(diag(covbet)*(TT-1)/(TT-2))
sigbetas;mean(sigbetas)		# We better find the same as before with lapply!!

ls.print(mregs)

# The matrix of coefficient constraints

Rmat<-toeplitz(c(1,-1,rep(0,28)))
Rmat[lower.tri(Rmat)]<-0
Rmat<-Rmat[1:29,]

constcov<-Rmat %*% covbet %*% t(Rmat)
constvec<-Rmat %*% betas

HotTea<- t(constvec) %*% solve(constcov) %*% constvec
HotTea
qchisq(0.95,29)

((60-1)/(60-29)) *
29* qf(0.95,29,60-29)

# Hotelling T2 for the alpha=0 hypothesis

HotTeaA <- t(alphs) %*% solve(covalph) %*% alphs
HotTeaA
qchisq(0.95,30)
(60-1)*30*qf(0.95,30,60-30)/(60-30)


# Stability of Betas

indus<-read.csv("47indus-mon.csv")
dim(indus);range(indus[,1])
rets<-indus

NN	<- length(indus[1,])-3;NN
rets[,2:(NN+1)]<-rets[,2:(NN+1)]-rets[,(NN+3)]


zm<-apply(rets[,2:49],2,mean)*12
zs<-apply(rets[,2:49],2,sd)*sqrt(12)
zm;zs
plot(zs,zm,xlim=c(0,max(zs)),ylim=range(c(0,zm)),xlab="Std. Dev",ylab="Mean-Rf") 
title("47 Industries and the VW US Market, 1970-2018")
points(zs[48],zm[48],pch="M") 
abline(0,zm[48]/zs[48])


# Get beta 1 and beta 2 for the 2 5-year periods

zret1 <- rets[rets[,1]>200900&rets[,1]<201400,]
zret2 <- rets[rets[,1]>201400&rets[,1]<201900,]
TT    <- length(zret1[,1]);TT
beta1<-lsfit(zret1[,NN+2],zret1[,2:(NN+1)])$coefficients[2,]
beta2<-lsfit(zret2[,NN+2],zret2[,2:(NN+1)])$coefficients[2,]

lsfit(beta1,beta2)$coefficients

par(mgp=c(1.5,0.5,0))
zrange<-range(c(beta1,beta2))
plot(beta1,beta2, xlim=zrange,ylim=zrange,xlab="Betas 2009-14",ylab="Betas 2014-18")
abline(0,1,lty=2)
abline(lsfit(beta1,beta2),col="red",lwd=2)
title("Figure 1: Betas 2014-18 vs 2009-13",line=0.1)

#
# Do it for the 9 periods from 1974
#

zrets <-rets[rets[,1]>197400,]

betas<-matrix(0,nrow=9,ncol=NN)
for (ii in 1:9) {
	ivec<-seq(60*(ii-1)+1,60*ii)
	betas[ii,]<-lsfit(rets[ivec,NN+2],rets[ivec,2:(NN+1)])$coefficients[2,]
}

# Now put them in two vectors of beta(t) and betas(t-1)

newbets <-c((betas[2:9,]))
oldbets <-c((betas[1:8,]))

par(mgp=c(1.5,0.5,0),mar=c(2.5,2.5,2,0.1))
zrange<-range(oldbets,newbets)
plot(oldbets,newbets,xlim=zrange,ylim=zrange,xlab="Old Beta",ylab="New Beta")
abline(0,1,lty=2)
abline(lsfit(oldbets,newbets),lwd=2,col="red")
title("Figure 2: 5-year Betas vs Previous 5-year Betas, 1974-2018",line=0.1)

ls.print(lsfit(oldbets,newbets))

plot(betas[1:8,],betas[2:9,],xlab="Old Beta",ylab="New Beta")

