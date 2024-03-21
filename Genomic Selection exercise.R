###################################
### Genomic Prediction exercise ###
###################################
#for more information: doi: 10.1534/genetics.114.164442.

rm(list=ls())
#install.packages("BGLR") #install the BGLR package
library(BGLR) #import package


#Loading and preparing the input data

#############################
#### Wheat line data set #### comprised of 599 lines; 4 locations; Pedigree; Marker matrix)
#############################

data(wheat) # wheat lines data
Y<-wheat.Y # (599 lines, 4 location) (records are centered and standardized to a unit variance within environment)
X<-wheat.X # 1,279 Diversity Array Technology (DArT) markers
A<-wheat.A #pedigree relatioship matrix
y<-Y[,1] # Yeild for Environment 1





#Setting the linear predictor
ETA<-list(list(X=X, model='BRR')) #Gaussian prior
ETA<-list(list(X=X, model='BL')) #Double exponential


#Fitting the model
fm<-BGLR(y=y,ETA=ETA, nIter=5000, burnIn=1000, thin = 5) #1000 interactions in total
cor(fm$y,fm$yHat) #calculate correlation
plot(fm) #plot graph pred x observed

# Extracting results from the model
bHat <- fm$ETA[[1]]$b # Estimated Marker Effects
SD.bHat <- fm$ETA[[1]]$SD.b
plot(bHat, ylab='Estimated Marker Effect',
     type='o',cex=.5,col=4,main='Marker Effects')

# lambda (regularization parameter of the Bayesian Lasso)
lambda<-scan('ETA_1_lambda.dat')
plot(lambda,type='o',col=2,cex=.5,ylab=expression(lambda))
abline(h=fm$ETA[[1]]$lambda,col=4,lwd=2)
abline(v=fm$burnIn/fm$thin,col=4)

unlink("*.dat") #delete files created by BGLR

#EXERCISE
# increase the number of interactions, burnin, thin
# change priors
# what happened?

#___________________________________________________
rm(list=ls())

#Assessment of Prediction Accuracy using GBLUP
#Loading and preparing the input data

data(wheat) # wheat lines data (centered and standardized to a unit variance within environment)
Y<-wheat.Y # yield (599 lines, 4 location)
X<-wheat.X # 1,279 Diversity Array Technology (DArT) markers
A<-wheat.A #pedigree relatioship matrix
y<-Y[,1] # Yeild for Environment 1
n<-nrow(X) #number of lines
p<-ncol(X) #number of markers

#Creating a Testing set
yNA<-y
set.seed(123)
tst<-sample(1:n,size=100,replace=FALSE)
yNA[tst]<-NA
cbind(y,yNA) #show which phenotype was removed

#Computing G
X<-scale(X,center=TRUE,scale=TRUE)
G<-tcrossprod(X)/p #VanRaden (2008)

#Generating heatmaps
heatmap(G) #generate heatmap for genomic relationship matrix
heatmap(A) #generate heatmap for pedigree relationship matrix

#EXERCISE -> Why the heatmaps are different?


#Fits the G-BLUP model
ETA<-list(list(K=G,model='RKHS'))
fm<-BGLR(y=yNA,ETA=ETA, nIter=5000, burnIn=2000, thin = 5)
unlink("*.dat")

#generate correlation graph between pred x observed
plot(fm$yHat,y,xlab="Phenotype",
     ylab="Pred. Gen. Value" ,cex=.8,bty="L")

points(x=y[tst],y=fm$yHat[tst],col=2,cex=.8,pch=19)
legend("topleft", legend=c("training","testing"),
       bty="n",pch=c(1,19), col=c("black","red"))

cbind(fm$y,fm$yHat)# NA on the first collumn now have value (on the second collumn) -> GEBV

# Assesment of correlation in training (TRN) and testing (TST) data sets
cor(fm$yHat[tst],y[tst]) #TST
cor(fm$yHat[-tst],y[-tst]) #TRN





