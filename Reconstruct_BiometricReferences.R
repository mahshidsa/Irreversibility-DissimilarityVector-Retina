
#################
# Try to reconstruct ESRID samples  using samples from other  
#data subjects in ESRID
# Reading  ESRID CLAHE 750x750 centered images 
#The distances used are BGM_based.
#The didtances are calculated using Edge-based distances and ci_e, cid_n= 0.4.
library(pixmap)
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
library(stats)#to create the doubly-centered distance matrix B
#The reconstruction part is done on GADI(NCI's high performance super computer)
###################################

#Load D, ESRID vs ESRID distances:
#The distances are calculated for grouped graphs
load("~/DV_ESRIDvsESRID_cidn_0.4_cide_0.4.Rdata")#D

#Select the target from one of each 46 gropups
#Then take out the target and other samples from its class from D
###########################

#str(D)# num [1:414, 1:414] 
########################

########################
#Generate GROUPED list
#To take out Target and its class from the BrIN set
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
IMG_Size <- 750
DummyMat <- matrix(0, nrow = IMG_Size*IMG_Size , ncol = 1)
MeanBreakInVec <- matrix(0, nrow =  IMG_Size*IMG_Size, ncol = 1)
for(q in 1:nrow(GROUP)){
  print(q)

  #########Phase 1) Embedding the break-in set##############################
  #########################################################
  ###Take out the target and its class from BrIn distances
  D <- D[-GROUP[q,], -GROUP[q,]]
  #The distances should be metric and Euclidean (If not enforce them)
  is.symmetric.matrix(D)#Our distances are symmetric
  
  #If they were not, to make them symmetric use the following line
  # D <- 1/2 *(D + t(D))
  
  tri.ineq(D)#FALSE
  
  #Enforce triangle inequality
  Zeta<-0
  for (j in 1:nrow(D)) {
    for (i in 1:nrow(D)) {
      for (k in 1:nrow(D)) {
        tmp <- (abs(D[i,j] - D[j,k] + D[i,k]))
        if(tmp > Zeta) Zeta <- tmp
      }}}
  
  
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
  #Load the function to generate the doubly centered matrix
  setwd("~/ReconMainFunctions")
  source("CreateMatrixB_March.R")
  ##Positive semidefinite flag
  PDF=0
  c=0
  CumC=0#cummulativeC
  while (PDF==0) {
    
    B <- CreateMatB(TrainingDistanceMatrix)$B
    c <- CreateMatB(TrainingDistanceMatrix)$c
    CumC =  CumC +c
    
    TrainingDistanceMatrix = CreateMatB(TrainingDistanceMatrix)$ModTrainingDistanceMatrix
    PDF = CreateMatB(TrainingDistanceMatrix)$PDF
  }
  
  Delta_EVD = Re(eigen(B)$value)#We don't need negative (or close to zero) values
  
  NegativeEigenVAls <-  which(Delta_EVD < 1e-3 , arr.ind = TRUE)
  Delta_EVD = Delta_EVD[-NegativeEigenVAls]
  #eigenvectors
  V_EVD = eigen(B)$vectors
  V_EVD = V_EVD[,1:length(Delta_EVD)]
  str(V_EVD)# 
  
  sqrtDelta_EVD = vector()
  for (i in 1: length(Delta_EVD)) {
    sqrtDelta_EVD[i] = sqrt(Delta_EVD[i])
  }
  
  
  DeltaSQRT <- diag(sqrtDelta_EVD, nrow = length(sqrtDelta_EVD), ncol = length(sqrtDelta_EVD))
  #Embeded points
  Y = t(V_EVD %*% DeltaSQRT)
  
  str(Y)#Ceck dimensionality reduction
  #########Embedding of break-in set is done ##############################
  #########################################################
  
  
  
  
  #########Phase 2) Finding affine approximation##############################
  #########################################################
  ###############################
  #Because the distances are calculated for grouped graphs, we need to import the images with the same order here.
  #Load information about grouped graphs
  load("~/groupings_BGM_ECGHao.Rdata")#groupings
  
  IndeX <- as.vector(t(groupings))
  ########
  #Read ESRID samples
  setwd("~/pgm")#Directory to pgm images
  listSamples <- list.files()
  #clahe(GreenChannel_I1, nx= 8, ny= 8, limit = 2)
  ImageVecSize<- 750 * 750  # Image width x Image height
  
  BreakInMatrix<- matrix(0,ImageVecSize,length(listSamples))
  for (j in 1:length(listSamples)) {
    BreakInMatrix[,j] <- cbind( as.vector( read.pnm(listSamples[IndeX[j]])@grey ))
  }
  #(readJPEG(Clahe(listSamples[j], nx=8, ny=8, limit=2)))))
  str(BreakInMatrix)#num
  #Take out reference images
  BreakInMatrix <- BreakInMatrix[,-GROUP[q,]]
  #Now, take out the target and it's class
  ##########################
  ########Mean_Image (center images by subtracting their mean from each of them)
  MeanBreakIn <- (1/ncol(BreakInMatrix))*rowSums(BreakInMatrix)
  str(MeanBreakIn)# 
  MeanBreakInVec <- cbind(MeanBreakInVec, MeanBreakIn)
  ###########
  #Generate a matrix  with columns equal to mean face image
  # in order to subtract the mean face image from the Break-in set matrix
  #rowOFones <- rep(1,nrow(BreakInMatrix))
  columnOFones <- t(rep(1,ncol(BreakInMatrix)))
  #as.matrix(rowOFones)
  as.matrix(columnOFones)
  #matrixOFones <- rowOFones %*% columnOFones #Different from matrix_of_ones
  #str(matrixOFones)# num [1:562500, 1:405]
  #matrix_of_ones<- matrix(1, nrow=nrow(DoubledD), ncol = ncol(DoubledD))
  #####
  MeangraphMatrix <- MeanBreakIn %*% (columnOFones)
  str(MeangraphMatrix)#[1:562500, 1:405]
  ####
  #Now, we can subtract mean image from the break-in set
  X <- (BreakInMatrix)-(MeangraphMatrix)
  #str(X) # num [1:562500, 1:405]
  
  
  #Find rigid affine transformation and principal coordinates
  ########### AofR and XofR
  #?svd
  SVDresult <- svd(X)
  
  U <- SVDresult $ u
  str(U)# num[1:562500, 1:405]Ar is transpose of U
  # The results of svd should be squared
  DofSVD = SVDresult $ d
  str(DofSVD)#num [1:405] 
  NegSingVals <- which(DofSVD < 1e-3, arr.ind = TRUE)#405
  DofSVD = DofSVD[-NegSingVals]
  ## Find the most significant value (Accounted for 99%)
  NintyNinePer <- which(cumsum(DofSVD)/sum(DofSVD)>0.99, arr.ind = TRUE)#386
  DofSVD = DofSVD[1:NintyNinePer[1]]

  #Square the values
  Lambda_PCA <- vector()
  for (j in 1:length(DofSVD)) {
    Lambda_PCA[j] = (DofSVD[j])^2
  }
  ###
  A_r <- t(U)
  A_r = A_r[1:length(DofSVD),]#Rigid affine trasformation
  
  Xr <- A_r %*% X 
  str(Xr)##num [1:386, 1:405] 
  ######################################
  #Non-rigid affine approximation
  #Find A-nr
  Lambda_PCA_Inv = vector()
  for (k in 1:length(Lambda_PCA)) {
    Lambda_PCA_Inv[k] = (Lambda_PCA[k])^-1
  }
  
  Lambda_InverseMat <- diag(Lambda_PCA_Inv, nrow = length(Lambda_PCA_Inv), ncol = length(Lambda_PCA_Inv))
  
  A_nr <- Y %*% t(Xr) %*% Lambda_InverseMat 

  
  #########Affine approximation is done ##############################
  #########################################################
  
  #######################################
  #########################################
  
  ##########Reconstruct
  ##############Reading distances between the target and break-in set
  #Distances are distances of ESRID to ESRID samples
  #BGM was used to calculate them.
  load("~/DV_ESRIDvsESRID_cidn_0.4_cide_0.4.Rdata")
  
  #Select the target from one of each 46 gropups
  #Then take out the target and other samples from its class from D
  ###########################
  
  for (s in 1:ncol(GROUP)) {
    ##############
    dhat<- D[GROUP[q,s],-GROUP[q,]]# 
    str(dhat)#num 
    
    dHat = dhat + rep(Zeta,length(dhat))  
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
    
    #Embed the target
    YofZ <- t( MatrixF %*% EdaGGer)
    #EdaGGer %*% MatrixF
    
    #Pseudo-invert affine approximation to find approximation of Y_of_target(=X_of_target)
    Pinv_A_nr <- pinv(A_nr)
    #Find X_of_Z
    xofZ <- t(A_r) %*% Pinv_A_nr %*% YofZ
    ####################
    #Should add the meanImage before reconstructing finishes
    #xofZ+MeanBreakIn is the reconstructed vector. Plot it as an image to observe it.
    DummyMat <- cbind(DummyMat, (xofZ+MeanBreakIn) )
    
  }
}
str(DummyMat)

TargetZ <- DummyMat[,2:415]

wd <- "~/Results/ImageRecon"
outputfile <- paste(wd , "/CircularCropped/ImBased_ReconBR_ESRIDs_BGM.Rdata", sep = "")
save(TargetZ, file = outputfile)
###The reconstruction took a long time.
#We performed it on GADI HPC.
#load(outputfile)
str(TargetZ)
#################
