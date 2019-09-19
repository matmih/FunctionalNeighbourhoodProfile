require(pROC)

 #@author Matej Mihelcic
 #@institution Rudjer Boskovic Institute, Zagreb, Croatia
 #@mail matmih1@gmail.com
 #@description function that plots comparative AUROC curves between a baseline and ensemble method for phenotype prediction
 
plotComparativeROC<- function(MatrixPred, MatrixPred2,numFeatures,w1,w2, FeatureNames, index){
  
  aucs<-c();
  recals<-c();
  recalsBase<-c();
  

    class <- as.numeric(MatrixPred[,2*numFeatures+index+1]);
    prediction1 <- as.numeric(MatrixPred[,(2*index-1+1)]);
    prediction0 <- as.numeric(MatrixPred[,2*index+1]);
    prediction<-(prediction1/(prediction1+prediction0));
    maxP<-max(prediction);
    prediction<-prediction/maxP;
    predictionC2<-as.numeric(MatrixPred2[index,]);
    maxP1<-max(predictionC2,na.rm=TRUE);
    if(maxP1>0)
      predictionC2<-predictionC2/maxP1;
    if((w1+w2)!=1){
      w1<-w1/(w1+w2);
      w2<-w2/(w1+w2);
    }
    predFin<-(w1*prediction+w2*predictionC2);
    
    roc_obj1 <- plot(roc(class, predFin), grid=TRUE,
                    print.auc=TRUE, show.thres=TRUE, col = "blue", main = paste("Method comparison for phenotype: ",FeatureNames[[1]][index]), legacy.axes = TRUE, asp = NA, cex.axis = 2.0, cex.lab = 2.0);
    roc_obj2 <- plot(roc(class, predictionC2), grid=TRUE,
                    print.auc=TRUE, print.auc.y = .4, show.thres=TRUE, add = TRUE, col = "darkgreen");
    legend(x="bottom",legend = c("Baseline","Ensemble"), col = c("darkgreen","blue"),lty = 1, cex = 1.0, bty = 'n');
  
}