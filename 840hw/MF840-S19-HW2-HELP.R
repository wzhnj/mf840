#
# HELP AND TIPS FOR HW 2
#

# Problem 1 and 2

# Easiest way to output to a table as in Table 1 is to output to  csv with
# Then you open XL and copy paste the matrix on the WORD table.
# Done

# Here you would do the 2 columns for beta, the two columns for
# alpha and do the stars by hand

indus<-read.csv("16indus-mon.csv")
names(indus)		# What's in there?
dim(indus)			# How many columns, observations
range(indus[,1])	# Dates?


rets<- indus[indus[,1]>200900&indus[,1]<201400,]	# Select some dates

allregs<-lsfit(xvector,ymatrix)	# lsfit does all regressions at once

zcoefs <-allres$coefficients		# the 16 by 2 matrix of intercepts and slopes
zalphas<-zcoefs[2,]		# The intercepts

zres   <-allregs$residuals	
		# Matrix of 16 residual series use it to compute Sigma
		# as in the SUR formulas					
Sigma  <-cov(zres)		

xdat   <-cbind(rep(1,?),xrm in here)
xpxinv <-solve(t(xdat)%*%xdat)	#(X'X)-1	
Sigma %x% xpxinv					# The Kronecker is %x% with the letter x

# Now you have a 32 x 32 covariance matrix of interept/slope for 16 regressions
# we often want a subset of the matrix, very easy with seq()

Zmatrix[seq(1,31,by=2),seq(1,31,by=2)]	# would extract what?
Zmatrix[seq(2,32,by=2),seq(1,32,by=2)]	# .. and that one?


# Then diag and square root if you want the individual standard
# deviations of the alphas or the betas

zstdevs<-sqrt(diag(zcovmatrix))

# And a vector of tstats

(zbetas-1)/zstdevs

# Easy way to put correlations into a vector

zcors<-sort(unique(cor(zres)))	# The one is still there, remove it
zcors[zcors<1]					# <- like this
								# Then you can do
mean(zcors);sd(zcors)			# etc...

(zbetas-1)


# FOR TABLE 	1, Multiply the intercepts by 12 to annualize them


write.csv(zmatrix,"filename.csv")

#
# Problem 3
#

# Many ways to roll, just a loop and save
# or you can rollapply lm but you need to make it a function
# look at the rollapply help for rolling regressions
# Here a loop is really easy

# Can go by years
zyear<-round(indus[,1]/100,0) # just gives you the year of the observation

#Here we know we can grab every 60 observations
#

allbets<-matrix(0,nrow=9,ncol=47)	# Prepare a matrix of the 9 period betas

# Then loop the regression for 9 periods of 5 years
for (ii in 1:9)	{	
ivec<-seq(60*(ii-1)+1,60*ii)		# Take these 60 observations
coefs<-lsfit(xvector[ivec],ymatrix[ivec,])$coefficients	
	# Do the regressions
	# Just keep what you need, the betas, not the intercepts
allbets[ii,]<- z betas	# and put them in the estimated row ii
}

# Now you have the "second period betas" in allbets[2:9,]
# and their corresponding previous betas in allbets[1:8,]
# Just plot

# To run the regression of all pairs of new (ydata) on old (xdata) 
# betas, you need to make them a vector


newbetas<-c(allbets[2:9,])		# Make sure you understand what is where!
oldbetas<-c(allbets[1:8,])		# and that you are matching same industry
								# estimate from one 5-year period to then next.
# Then plot

# and regress !
lsfit(oldbetas,newbetas)


# Make sure to read this critically. It is a very bad
# idea to use someone's code without understanding it.





