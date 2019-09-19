   #@author Matej Mihelcic
  #@institution Rudjer Boskovic Institute, Zagreb, Croatia
  #@mail matmih1@gmail.com
  #@description function to compute AUC/AUPRC for each target variable

#require(AUC)
require(PRROC)

computeAUPRCB<- function(Predictions,Targets, GOS ,numFeatures){
  
  auprcs<-c();
  feat<-c();
  
  for (i in 1:numFeatures){
    class <- as.numeric(Targets[,i+1]);
  
    #prediction <- as.numeric(Matrix[,numFeatures*2+(i)*2]);
    prediction <- as.numeric(Predictions[i,1:dim(Predictions)[2]]);
    prediction[which(is.na(prediction))]<-0;
    #print(length(prediction))
    #print(length(class))
    fg<-prediction[class == 1];
    bg<-prediction[class == 0];
   
    #print("NAs");
    #print(sum(is.na(fg))); print(sum(is.na(bg)));
    #print(which(is.na(bg)))
    #print(bg);
    #print(prediction[which(is.na(bg))]);
    
    if(sum(class) == 0)
      next;
    print(sum(class))
    #roc_obj <- roc(class, prediction);
    pr<-pr.curve(scores.class0 = fg, scores.class1 = bg);
    auprcs<-c(auprcs,pr$auc.integral);
    feat<-c(feat,GOS[i,1]);
    #aucs<-c(aucs,auc(class,prediction));
    #aucs<-c(aucs,auc(class,prediction));
    #aucs<-c(aucs,auc(roc(prediction,class)));
    print(auprcs)
  }
  
  return(data.frame(auprcs,feat));
  
}
