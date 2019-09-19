  #@author Matej Mihelcic
  #@institution Rudjer Boskovic Institute, Zagreb, Croatia
  #@mail matmih1@gmail.com
  #@description function to compute the AUPRC values for the ensemble classifier used for phenotype prediction

require(PRROC)

computeAUPRCEnsemble<- function(MatrixPred, MatrixPred2 ,numFeatures,w1,w2){
  
  auprcs<-c();
  
  for (i in 1:numFeatures){
    class <- as.numeric(MatrixPred[,2*numFeatures+i+1]);
    #prediction <- as.numeric(Matrix[,numFeatures*2+(i)*2]);
    prediction1 <- as.numeric(MatrixPred[,(2*i-1+1)]);
    prediction0 <- as.numeric(MatrixPred[,2*i+1]);
    prediction<-(prediction1/(prediction1+prediction0));
    maxP<-max(prediction);
    prediction<-prediction/maxP;
    predictionC2<-as.numeric(MatrixPred2[i,]);
    maxP1<-max(predictionC2,na.rm=TRUE);
    if(maxP1>0)
      predictionC2<-predictionC2/maxP1;
    if((w1+w2)!=1){
      w1<-w1/(w1+w2);
      w2<-w2/(w1+w2);
    }
    predFin<-(w1*prediction+w2*predictionC2);
    
    predFin[which(is.na(predFin))]<-0;
    
    fg<-predFin[class == 1];
    bg<-predFin[class == 0];
    
    pr<-pr.curve(scores.class0 = fg, scores.class1 = bg);
    auprcs<-c(auprcs,pr$auc.integral);
  }
  
  return(auprcs);
  
}