   #@author Matej Mihelcic
  #@institution Rudjer Boskovic Institute, Zagreb, Croatia
  #@mail matmih1@gmail.com
  #@description function to compute AUPRC values from the .oob prediction file

require(PRROC)

computeAUPRCOOB<- function(MatrixPred ,numFeatures){
  
  auprcs<-c();
  
  for (i in 1:numFeatures){
    class <- as.numeric(MatrixPred[,2*numFeatures+i+1]);
    #prediction <- as.numeric(Matrix[,numFeatures*2+(i)*2]);
    prediction1 <- as.numeric(MatrixPred[,(2*i-1+1)]);
    prediction0 <- as.numeric(MatrixPred[,2*i+1]);
    prediction<-(prediction1/(prediction1+prediction0));
    
    prediction[which(is.na(prediction))]<-0;
    
    fg<-prediction[class == 1];
    bg<-prediction[class == 0];
    
    pr<-pr.curve(scores.class0 = fg, scores.class1 = bg);
    auprcs<-c(auprcs,pr$auc.integral);
  }
  
  return(auprcs);
  
}
