  #@author Matej Mihelcic
  #@institution Rudjer Boskovic Institute, Zagreb, Croatia
  #@mail matmih1@gmail.com
  #@description function to draw OR heatmap. 
  #mode = 0 -> cluster by LOR, mode = 1, cluster by semantic similarity
  

drawORHeatmap<-function(ORMatrix, SimilarityMatrix,cuttof, colorPalete, mode, orders){
  library(gplots)
  library('seriation')
  library('viridis')
  
  tmpSimMat<-data.matrix(SimilarityMatrix,rownames.force=NA);
  tmpORMat<-data.matrix(ORMatrix,rownames.force=NA);
  tmpORMat[which(tmpORMat>6)] = 6;
  tmpORMat[which(tmpORMat< (-6))] = -6;
  tmpSimMat[is.na(tmpSimMat)] = 0.0;
  
  order<-as.vector(orders[,1]);
  orderT<-as.vector(orders[,2]);
  
  if(mode==0){
    heatmap.2(tmpORMat,col=colorPalete, key.xlab = "Log OR", key.ylab ="Frequeny" , Rowv=order,Colv=orderT,scale="none",trace="none",breaks = c(-6,-2,0,0.001,2,4,6));#col=redgreen
  }
  else{
    tmpSimMat<-data.matrix(SimilarityMatrix,rownames.force=NA);
    tmpSimMat[which(tmpSimMat>6)] = 6;
    tmpSimMat[which(tmpSimMat< (-6))] = -6;
    t<-heatmap.2(tmpSimMat,col=colorPalete, key.xlab = "Semantic similarity", key.ylab ="Frequeny" , Rowv=order,Colv=orderT,scale="none",trace="none");#col=redgreen
    tmpORMat[which(is.na(tmpSimMat))]<-NaN;
    heatmap.2(tmpORMat,col=colorPalete, key.xlab = "Semantic similarity", key.ylab ="Frequeny" , Rowv=t$rowDendrogram,Colv=t$colDendrogram,scale="none",trace="none", breaks = c(0,2,3.99,4,6));#col=redgreen
  }
  return (data.frame(order,orderT));
}