##This code simulations the evolution of several traits (plumage similarity, song similarity, guild overlap, hybridization, mass difference, bill difference)
##And conducts logistic regressions to create a confidence interval of z values
##For more details, see description in Appendix S1

##first, load packages, data, and phylogenetic tree
require(phytools)
require(geiger)
source("var.matrix.R")
source("var.mv.matrix.R") #functions defined in var.matrix.R and var.mv.matrixR need to be run in order for the scripts to work


tree<-read.tree("warblertree.tre") #this tree is trimmed from the ultrametric r8 tree corresponding to fig5 in lovette 
N<-5000 #number of simulations



##next, generate data, first using fitContinuous to estimate sigma2, then adding a as mean and bounds set at the max and min values +-10%

#mass, centering data with z transformation (xi-xbar/sdev), bounds set as 10% of the upper and lower bounds of empirical z scores
simdata.mass<-fastBM(tree,a=0,sig2= 0.016434, bounds=c(-1.567149, 3.491771),nsim=N)
variable2.matrix<-var.matrix(simdata.mass)
variable2.matrix<-sqrt(variable2.matrix)
for (i in 1:N){
	variable2.matrix[,i]<-(variable2.matrix[,i]-mean(variable2.matrix[,i]))/sd(variable2.matrix[,i])
	} #z transform

#now using Z score transfomration to simulate bill, bounds set as 10% of the upper and lower bounds of empirical z scores
simdata.bill<-fastBM(tree,sig2= 0.022336,a= 0,bounds=c(-1.734458, 2.938706), nsim=N)
variable3.matrix<-var.matrix(simdata.bill)
variable3.matrix<-sqrt(variable3.matrix)
for (i in 1:N){
	variable3.matrix[,i]<-(variable3.matrix[,i]-mean(variable3.matrix[,i]))/sd(variable3.matrix[,i])
	} #z transform

##In this iteration of the simulation, the below measurements which are pairwise in nature are simulated 
#as two continuous variables, and then the euclidian distance is calculated between species for two values, rather than one


R<-0.044672 #rate parameter to use min=0.016434 max=0.044672
#plumage similarity 
simdata.1<-fastBM(tree,sig2= R,a= 0,nsim=N)
simdata.2<-fastBM(tree,sig2= R, a= 0,nsim=N)
simdata.plumarray<-array(dim=c(42,N,2)) #42 is number of tips, N is number of simulations, 2 is number of variables simulated
simdata.plumarray[,,1]<-simdata.1
simdata.plumarray[,,2]<-simdata.2
variable4.matrix<-var.mv.matrix(simdata.plumarray)
variable4.matrix<-log(variable4.matrix) #skews distribution to represent that "4" is ceiling value in raw data
for (i in 1:N){
	variable4.matrix[,i]<-(variable4.matrix[,i]-mean(variable4.matrix[,i]))/sd(variable4.matrix[,i])
	} #z transform


#song similarity
simdata.1<-fastBM(tree,sig2= R,a= 0,nsim=N)
simdata.2<-fastBM(tree,sig2= R, a= 0,nsim=N)
simdata.humsongarray<-array(dim=c(42,N,2)) #42 is number of tips, N is number of simulations, 2 is number of variables simulated
simdata.humsongarray[,,1]<-simdata.1
simdata.humsongarray[,,2]<-simdata.2
variable5.matrix<-var.mv.matrix(simdata.humsongarray)
variable5.matrix<-log(variable5.matrix) #skews distribution to represent that "4" is ceiling value in raw data
for (i in 1:N){
	variable5.matrix[,i]<-(variable5.matrix[,i]-mean(variable5.matrix[,i]))/sd(variable5.matrix[,i])
	} #z transform


#spcc data 
simdata.1<-fastBM(tree,sig2= R,a= 0,nsim=N)
simdata.2<-fastBM(tree,sig2= R,a= 0,nsim=N)
simdata.songarray<-array(dim=c(42,N,2)) #42 is number of tips, N is number of simulations, 2 is number of variables simulated
simdata.songarray[,,1]<-simdata.1
simdata.songarray[,,2]<-simdata.2

variable7.matrix<-var.mv.matrix(simdata.songarray) #function extracts euclidian distance between tips
for(i in 1:N){
	variable7.matrix[,i]<-max(variable7.matrix[,i])-variable7.matrix[,i]
	} #turns into a similarity index
variable7.matrix<-sqrt(variable7.matrix) #take the square root
for (i in 1:N){
	variable7.matrix[,i]<-(variable7.matrix[,i]-mean(variable7.matrix[,i]))/sd(variable7.matrix[,i])
	} #z transform


#guild overlap (0.2931034 [85] are 1's)
simdata.guild<-fastBM(tree,sig2=R,a=0,nsim=N)
variable6.matrix<-var.matrix(simdata.guild)

#this loop pulls out 85 differences for each run of the code for identifying as 1's (note: this loop is computationally expensive, and takes several minutes to run)
int.matrix<-matrix(nrow=290,ncol=N)
for(j in 1:N){
	i=1
	while(i<291){
		int.matrix[i,j]<-ifelse(variable6.matrix[i,j]<quantile(variable6.matrix[,j],c(0.2931034)),1,0) 
		i=i+1
		}
	}
variable6.matrix<-int.matrix
	
#hybridization, 24 pairs, 0.08275862
simdata.hybrid<-fastBM(tree,sig2=R,a=0,nsim=N)
variable8.matrix<-var.matrix(simdata.hybrid)
#this loop pulls out 24 differences for each run of the code for identifying as 1's (note: this loop is computationally expensive, and takes several minutes to run)
int.matrix<-matrix(nrow=290,ncol=N)
for(j in 1:N){
	i=1
	while(i<291){
		int.matrix[i,j]<-ifelse(variable7.matrix[i,j]<quantile(variable7.matrix[,j],c(.08275862)),1,0)
		i=i+1
		}
	}
variable8.matrix<-int.matrix

##run a logistic regression on each of these variable.matrices 
##NOTE this requires loading the Z transformed log(syntopy) data (Zsyn) and IT data (all.IT)

results.matrix<-matrix(ncol=N,nrow=9)
for (i in 1:N){
	m1<-glm(data$all.IT~variable2.matrix[,i]+variable3.matrix[,i]+variable4.matrix[,i]+variable6.matrix[,i]+variable8.matrix[,i]+data$Zsyn+variable7.matrix[,i]+variable5.matrix[,i],family="binomial")	
	results.matrix[1,i]<-as.numeric(summary(m1)[5]) #this pulls out the aic score
	results.matrix[2,i]<-as.numeric(coef(summary(m1))[2,3]) ##this pulls out z value for the first predictor variable in a logistic regression, for p-value pull out [2,4]
	results.matrix[3,i]<-as.numeric(coef(summary(m1))[3,3])
	results.matrix[4,i]<-as.numeric(coef(summary(m1))[4,3])
	results.matrix[5,i]<-as.numeric(coef(summary(m1))[5,3])
	results.matrix[6,i]<-as.numeric(coef(summary(m1))[6,3])
	results.matrix[7,i]<-as.numeric(coef(summary(m1))[7,3])
	results.matrix[8,i]<-as.numeric(coef(summary(m1))[8,3])
	results.matrix[9,i]<-as.numeric(coef(summary(m1))[9,3])
	}
	
##create 95% confidence intervals for each variable's z statistic, NOTE refer to glm() equation above for identifying the order of variables

aic.CI<-quantile(results.matrix[1,],c(0.05,0.95))	
var2.CI<-quantile(results.matrix[2,],c(0.05,0.95))
var3.CI<-quantile(results.matrix[3,],c(0.05,0.95))
var4.CI<-quantile(results.matrix[4,],c(0.05,0.95))
var5.CI<-quantile(results.matrix[5,],c(0.05,0.95))
var6.CI<-quantile(results.matrix[6,],c(0.05,0.95))
var7.CI<-quantile(results.matrix[7,],c(0.05,0.95))
var8.CI<-quantile(results.matrix[8,],c(0.05,0.95))
var9.CI<-quantile(results.matrix[9,],c(0.05,0.95))

var2.CI
var3.CI
var4.CI
var5.CI
var6.CI
var7.CI
var8.CI
var9.CI


##just running the best model now:

results.matrix<-matrix(ncol=N,nrow=5)
for (i in 1:N){
	m1<-glm(data$all.IT~variable4.matrix[,i]+data$Zsyn+variable7.matrix[,i],family="binomial")
	results.matrix[1,i]<-as.numeric(summary(m1)[5]) #this pulls out the aic score
	results.matrix[2,i]<-as.numeric(coef(summary(m1))[2,3]) ##this grabs the z value for the first predictor variable in a logistic regression
	results.matrix[3,i]<-as.numeric(coef(summary(m1))[3,3])
	results.matrix[4,i]<-as.numeric(coef(summary(m1))[4,3])

	}
	
##create 90% and 95% confidence intervals for each variable's z statistic

aic.CI<-quantile(results.matrix[1,],c(0.05,0.95))	
var2.CI90<-quantile(results.matrix[2,],c(0.05,0.95))
var3.CI90<-quantile(results.matrix[3,],c(0.05,0.95))
var4.CI90<-quantile(results.matrix[4,],c(0.05,0.95))

var2.CI90
var3.CI90
var4.CI90

