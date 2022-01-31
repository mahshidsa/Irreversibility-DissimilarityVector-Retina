
CreateMatB<-function(TrainingDistanceMatrix){
  
SQRDTrainingDistanceMatrix <- matrix(0, nrow(TrainingDistanceMatrix), ncol(TrainingDistanceMatrix))
for(i in 1:nrow(TrainingDistanceMatrix)){
  for (j in 1:ncol(TrainingDistanceMatrix)) {
    SQRDTrainingDistanceMatrix[i,j] <- (TrainingDistanceMatrix[i,j])^2
  }
}

A = - 0.5 * SQRDTrainingDistanceMatrix

OneMat = matrix(1, nrow = nrow(TrainingDistanceMatrix), ncol = ncol(TrainingDistanceMatrix))


I = diag(1, nrow = nrow(TrainingDistanceMatrix), ncol = ncol(TrainingDistanceMatrix))

H = ( I - ( (1/nrow(TrainingDistanceMatrix)) * (OneMat) ) )

B = H %*% A %*% H

EigenValues <- eigen(B)$value

c = -2 * min(Re(EigenValues))#

cMatriX = c *( OneMat -I )

ModTrainingDistanceMatrix = sqrt((SQRDTrainingDistanceMatrix)+ cMatriX )

if(c > 1e-9)#TRUE
{
 PDF=0
}

if(c <= 1e-13 ) {
  PDF=1
}

list(B=B, ModTrainingDistanceMatrix= ModTrainingDistanceMatrix ,c= -2 * min(Re(EigenValues)), PDF=PDF)

}
#########################

