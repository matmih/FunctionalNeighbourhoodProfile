  #@author Matej Mihelcic
  #@institution Rudjer Boskovic Institute, Zagreb, Croatia
  #@mail matmih1@gmail.com
  #@description function to draw scatter plot showing the average percentage of gene functional Neighborhood of different COGs that contain GO functions GO1 and GO2
  

drawScatterPlotN<-function(cuttof1,cuttof2,drawEB,writeText,GO1Lab,GO2Lab, COGDesc){
  
  GO1<-read.delim(paste("~/NetBeansProjects/Recursive file search/FilesForScatterPlot/COGNeighCompPerc",GO1Lab,".txt",sep=""));
  GO2<- read.delim(paste("~/NetBeansProjects/Recursive file search/FilesForScatterPlot/COGNeighCompPerc",GO2Lab,".txt",sep=""));
  GO1Rand<- read.delim(paste("~/NetBeansProjects/Recursive file search/FilesForScatterPlot/COGNeighCompPercRand",GO1Lab,".txt",sep=""));
  GO2Rand<- read.delim(paste("~/NetBeansProjects/Recursive file search/FilesForScatterPlot/COGNeighCompPercRand",GO2Lab,".txt",sep=""));
    
  plot(GO1$mean,GO2$mean,xlim=c(0,max(GO1$mean[which((GO1$mean>=cuttof1) & (GO2$mean>=cuttof2))]+GO1$std[which((GO1$mean>=cuttof1) & (GO2$mean>=cuttof2))])),ylim=c(0,max(GO2$mean[which((GO1$mean>=cuttof1) & (GO2$mean>=cuttof2))]+ GO2$std[which((GO1$mean>=cuttof1) & (GO2$mean>=cuttof2))])),xlab=GO1Lab,ylab=GO2Lab,col=addTrans("black",100));
  points(GO1Rand$mean,GO2Rand$mean,col=addTrans("red",100));
  title(main=paste("COG neighbourhood comparissons ","cut1: ", cuttof1, ", cut2: ",cuttof2));
  
  if(writeText==TRUE){
    labs<-GO1$COG[which((GO1$mean>=cuttof1) & (GO2$mean>=cuttof2))];
    
    labNew<-c();
    
    for(i in 1:length(labs)){
      desc<-findCOGLabel(labs[i],COGDesc);
      print("desc")
     # print(desc)
      lab<-paste(desc," \n ",labs[i]);
      #print("lab")
      #print(lab)
      labNew<-c(labNew,lab);
    }
    
    #print(labNew)
    text(GO1$mean[which((GO1$mean>=cuttof1) & (GO2$mean>=cuttof2))],GO2$mean[which((GO1$mean>= cuttof1) & (GO2$mean>=cuttof2))], labels=labNew, cex= 0.7, pos=3);
    
    #text(GO1$mean[which((GO1$mean>=cuttof1) & (GO2$mean>=cuttof2))],GO2$mean[which((GO1$mean>= cuttof1) & (GO2$mean>=cuttof2))], labels=GO1$COG[which((GO1$mean>=cuttof1) & (GO2$mean>=cuttof2))], cex= 0.7, pos=3)
  }
  
  if(drawEB==TRUE){
    arrows(GO1$mean[which((GO1$mean>=cuttof1) & (GO2$mean>=cuttof2))],max(0,GO2$mean[which((GO1$mean>= cuttof1) & (GO2$mean>=cuttof2))]-GO2$std[which((GO1$mean>= cuttof1) & (GO2$mean>=cuttof2))]), GO1$mean[which((GO1$mean>=cuttof1) & (GO2$mean>=cuttof2))],GO2$mean[which((GO1$mean>= cuttof1) & (GO2$mean>=cuttof2))]+GO2$std[which((GO1$mean>= cuttof1) & (GO2$mean>=cuttof2))],length=0.05, angle=90, code=3);
    arrows(max(0,GO1$mean[which((GO1$mean>=cuttof1) & (GO2$mean>=cuttof2))]-GO1$std[which((GO1$mean>=cuttof1) & (GO2$mean>=cuttof2))]),GO2$mean[which((GO1$mean>= cuttof1) & (GO2$mean>=cuttof2))], GO1$mean[which((GO1$mean>=cuttof1) & (GO2$mean>=cuttof2))]+GO1$std[which((GO1$mean>=cuttof1) & (GO2$mean>=cuttof2))],GO2$mean[which((GO1$mean>= cuttof1) & (GO2$mean>=cuttof2))],length=0.05, angle=90, code=3);
  }
}