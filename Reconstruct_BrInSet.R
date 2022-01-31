#################
# Try to reconstruct ESRID samples 
# Reading  ESRID CLAHE 750x750 centered images 
#The distances used are BGM_based.(edge-based)

library(jpeg)
library(splancs)
library(spatial)
library(network)
library(clue)
library(iterators)
library(foreach)
library(doSNOW)
library(sna)
library(gtools)
require(graphics)
require(utils)
library(MASS)
library(intervals) 
library(lubridate)
library("pracma")
library(factoextra)
library(matrixcalc)
library(fossil)
library(utils)
library(graphics)
library(pixmap)
library(stats)#to create the doubly-centered distance matrix B
#The reconstruction part is done on GADI
###################################

#Load BGM based distance
load("~/DV_ESRIDvsESRID_cidn_0.4_cide_0.4.Rdata")

#Select the target from one of each 46 groups
#Then take out the target and other samples from it's class from D
###########################
########################
#To take out Target and its class from the BrIn set
ClaSS <- c(1:46)

StarT <- (ClaSS-1)* 9 +1
EnD <-  ClaSS * 9
c(1:9)
StarT[1]
EnD[1]
GROUP <- matrix(0,46,9)
for(j in 1:length(ClaSS)){
  GROUP[j,] <- c( StarT[j] : EnD[j] )
}
################
########################
#Samples of the first data subject ar refrences
#Samples of data subject 2,...,46 are break-in set members
q=1 
print(q)


#########Phase 1) Embedding the break-in set##############################
#########################################################
###Take out the target and its class from BrIn set distances
D <- D[-GROUP[q,], -GROUP[q,]]
#The distances should be metric and Euclidean (If not enforce them)
is.symmetric.matrix(D)#

#D <- 1/2 *(D + t(D))

#Enforce triangle inequality
tri.ineq(D)#FALSE

Zeta<-0
for (j in 1:nrow(D)) {
  for (i in 1:nrow(D)) {
    for (k in 1:nrow(D)) {
      tmp <- (abs(D[i,j] - D[j,k] + D[i,k]))
      if(tmp > Zeta) Zeta <- tmp
    }}}#  1.882771


IdentityMatrix <- diag(nrow(D))
OneMatrix <- matrix(1, nrow=nrow(D), ncol = ncol(D))
DiagZetaMatrix <- Zeta * (OneMatrix - IdentityMatrix)

###Adding Zeta to the matrix
D <- D + DiagZetaMatrix


tri.ineq(D)#TRUE

is.symmetric.matrix(D)#TRUE
#Now the distances are Euclidean
##################

#MultiDimensional Scaling to embed break-in set members
####################MDS
TrainingDistanceMatrix = D
setwd("~/ReconMainFunctions")
source("CreateMatrixB_March.R")

PDF=0
c=0
CumC=0
while (PDF==0) {
  
  B <- CreateMatB(TrainingDistanceMatrix)$B
  c <- CreateMatB(TrainingDistanceMatrix)$c
  CumC =  CumC +c
  
  TrainingDistanceMatrix = CreateMatB(TrainingDistanceMatrix)$ModTrainingDistanceMatrix
  PDF = CreateMatB(TrainingDistanceMatrix)$PDF
}

Delta_EVD = Re(eigen(B)$value)#all positive, the last 2 are close to zero

NegativeEigenVAls <-  which(Delta_EVD < 1e-3 , arr.ind = TRUE)
Delta_EVD = Delta_EVD[-NegativeEigenVAls]
#
V_EVD = eigen(B)$vectors
V_EVD = V_EVD[,1:length(Delta_EVD)]
str(V_EVD)# num [1:405, 1:403] 

sqrtDelta_EVD = vector()
for (i in 1: length(Delta_EVD)) {
  sqrtDelta_EVD[i] = sqrt(Delta_EVD[i])
}


DeltaSQRT <- diag(sqrtDelta_EVD, nrow = length(sqrtDelta_EVD), ncol = length(sqrtDelta_EVD))

Y = t(V_EVD %*% DeltaSQRT)

str(Y)


#########Phase 1) Embedding the break-in set##############################
#########################################################


###############################


#########Phase 2) Finding affine approximation##############################
#########################################################
###############################
#Because the distances are calculated for grouped graphs, we need to import the images with the same order here.
#Load information about grouped graphs
load("~/groupings_BGM_ECGHao.Rdata")

IndeX <- as.vector(t(groupings))
########
#Read ESRID samples
setwd("~/pgm")#Directory to pgm images
listSamples <- list.files()

ImageVecSize<- 750 * 750  # Image width x Image height

BreakInMatrix<- matrix(0,ImageVecSize,length(listSamples))
for (j in 1:length(listSamples)) {
  BreakInMatrix[,j] <- cbind( as.vector( read.pnm(listSamples[IndeX[j]])@grey ))
}
#(readJPEG(Clahe(listSamples[j], nx=8, ny=8, limit=2)))))
str(BreakInMatrix)#num

BreakInMatrix <- BreakInMatrix[,-GROUP[q,]]
#Now, take out the target and it's class
##########################
######################################
####BE CAREFULE ABOUT WHAT YOU ARE TAKING OUT
str(BreakInMatrix)#

########Meangraph?! 
MeanBreakIn <- (1/ncol(BreakInMatrix))*rowSums(BreakInMatrix)
str(MeanBreakIn)# 
#MeanBreakInVec <- cbind(MeanBreakInVec, MeanBreakIn)
###########
#Generate a matrix  with columns equal to mean face image
# in order to subtract the mean face image from the Break-in set matrix
rowOFones <- rep(1,nrow(BreakInMatrix))
columnOFones <- t(rep(1,ncol(BreakInMatrix)))
as.matrix(rowOFones)
as.matrix(columnOFones)
matrixOFones <- rowOFones %*% columnOFones #Different from matrix_of_ones
str(matrixOFones)# num [1:562500, 1:405]
#matrix_of_ones<- matrix(1, nrow=nrow(DoubledD), ncol = ncol(DoubledD))
#####
MeangraphMatrix <- MeanBreakIn %*% (columnOFones)
str(MeangraphMatrix)#[1:562500, 1:405]
####
#Now, subtract mean image from the break-in set
X <- (BreakInMatrix)-(MeangraphMatrix)
str(X) #

###########AofR and XofR
#?svd
SVDresult <- svd(X)

U <- SVDresult $ u
str(U)#Ar is transpose of U
# The results of svd should be squared
DofSVD = SVDresult $ d
str(DofSVD)#num [1:405] 
NegSingVals <- which(DofSVD < 1e-3, arr.ind = TRUE)#405
DofSVD = DofSVD[-NegSingVals]
## Find the most significant value (Acconted for 99%)
NintyNinePer <- which(cumsum(DofSVD)/sum(DofSVD)>0.99, arr.ind = TRUE)#386
DofSVD = DofSVD[1:NintyNinePer[1]]
#Square the values


#Square the values

Lambda_PCA <- vector()
for (j in 1:length(DofSVD)) {
  Lambda_PCA[j] = (DofSVD[j])^2
}
###
A_r <- t(U)#Rigid affine transformation
A_r = A_r[1:length(DofSVD),]

Xr <- A_r %*% X
str(Xr)#
######################################
#Non-rigid affine approximation
#Find A-nr
Lambda_PCA_Inv = vector()
for (k in 1:length(Lambda_PCA)) {
  Lambda_PCA_Inv[k] = (Lambda_PCA[k])^-1
}

Lambda_InverseMat <- diag(Lambda_PCA_Inv, nrow = length(Lambda_PCA_Inv), ncol = length(Lambda_PCA_Inv))

A_nr <- Y %*% t(Xr) %*% Lambda_InverseMat 
str(A_nr)# 

#########Affine approximation is done ##############################
#########################################################

#######################################
#########################################

##########Reconstruct
##############Reading distances between the target and break-in set
#Distances are distances of ESRID to ESRID
#BGM-based distances
load("~/DV_ESRIDvsESRID_cidn_0.4_cide_0.4.Rdata")

#Reconstruct break-in set members
##############
dhat<- D[-GROUP[q,],-GROUP[q,]]
str(dhat)

DummyMat <- matrix(0, nrow = ImageVecSize, ncol = 1)
for(l in 1:ncol(BreakInMatrix)){
  
  dHat = dhat[,l] + rep(Zeta,length(dhat[,l]))  
  #dHat = dhat + rep(Zeta,length(dhat))  
  squaredDHAT<- vector()
  for (k in 1:length(dHat)) {
    squaredDHAT[k] = (dHat[k])^2
  }
  dHat = sqrt(squaredDHAT+ rep(CumC,length(dhat)) )
  
  str(Y)
  yNorm<- vector()
  for( i in 1:ncol(Y)){
    yNorm[i] <- Norm(Y[,i])
  }
  
  #Calculate f_i
  f_i <- vector()
  for(i in 1:ncol(Y)-1){
    f_i [i] <- 1/2 *( ( (dHat[i])^2 - (yNorm[i])^2 ) - ( (dHat[i+1])^2 - (yNorm[i+1])^2 )  )
  }
  
  
  str(f_i) 
  f_i<-as.matrix(f_i)
  MatrixF<- t(f_i)
  str(MatrixF)  
  
  E <- matrix(0, nrow(Y), ncol(Y)-1)
  for (i in 1:ncol(Y)-1) {
    E[,i] <- (Y[,i+1] -Y[,i])
  }
  
  str(E)
  
  EdaGGer <- pinv(E)
  
  YofZ <- t( MatrixF %*% EdaGGer)
  #EdaGGer %*% MatrixF
  Pinv_A_nr <- pinv(A_nr)
  
  #Pseudo-invert affine approximation to find approximation of Y_of_target(=X_of_target)
  xofZ <- t(A_r) %*% Pinv_A_nr %*% YofZ
  ####################
  DummyMat <- cbind(DummyMat, (xofZ+MeanBreakIn) )
  
}
str(DummyMat)

TargetZ <- DummyMat[,2:406] 

wd <- "/~/Results/ImageRecon"
outputfile <- paste(wd , "/CircularCropped/ImBased_ReconBrIn_ESRIDs_BGM.Rdata", sep = "")
save(TargetZ, file = outputfile)
#load(outputfile)
str(TargetZ)
#################