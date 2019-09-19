   #@author Matej Mihelcic
  #@institution Rudjer Boskovic Institute, Zagreb, Croatia
  #@mail matmih1@gmail.com
  #@description function to compute AUC/AUPRC values for CAFFA experiments

computeAUPRROCCaffa<-function(ScoreMatrix, LabelMatrix){
  require(PRROC)
  require(pROC)
  aucs<-c();
  gos<-c();
  auprc<-c();
  
  for (i in 1:dim(ScoreMatrix)[1]){
    print(paste('i: ',i));
    vTrue<-data.matrix(LabelMatrix[i,2:dim(LabelMatrix)[2]]);
    vPred<-data.matrix(ScoreMatrix[i,2:dim(ScoreMatrix)[2]]);
    fg <- vPred[vTrue == 1];
    bg <- vPred[vTrue == 0];
    if(sum(vTrue)==0)
      next();
    roc_obj <- roc(vTrue, vPred);
    a<-auc(roc_obj);
    aucs<-c(aucs,a);
    #print(LabelMatrix[i,1]);
    gos<-c(gos,as.character(LabelMatrix[i,1]));
    pr <- pr.curve(scores.class0 = fg, scores.class1 = bg);
    auprc<-c(auprc,pr$auc.integral);
    #print(dim(data.matrix(pr[2])));
  }
  
  out<-data.frame(gos,aucs,auprc);
  
  write.table(out, file = "AUCCaffa.txt",
        append = FALSE, sep = " ",row.names = FALSE,
        col.names = FALSE);
  
  return (data.frame(gos,aucs,auprc));
  
}