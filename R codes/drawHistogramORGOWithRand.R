  #@author Matej Mihelcic
  #@institution Rudjer Boskovic Institute, Zagreb, Croatia
  #@mail matmih1@gmail.com
  #@description function to draw comparative histograms (random vs original data) or divided by semantic similarity (CLPar vs CLMed vs Dist)
  

drawHistogramORGOWRand<-function(RealSamp, RandomSamp){
  
  min<-min(min(RealSamp),min(RandomSamp));
  max<-max(max(RealSamp),max(RandomSamp));
  
  print(min)
  print(max)
  
  breaks<-seq((min-1),(max+1),0.1);
  
  a<-hist(data.matrix(RealSamp),col=rgb(0,0,1,0.5),breaks=breaks);
  b<-hist(data.matrix(RandomSamp),col=rgb(0,0,1,0.5),breaks=breaks);
  
  mc<-max(max(a$counts),max(b$counts));
  
  #hist(CLPar,col=rgb(1,0,0,0.5),breaks=breaks,main = "LogOR frequency", freq=TRUE,plot=TRUE,xlab="logOR",cex.axis=2.0, cex.lab=1.5,cex.main=1.8,ylim=c(0,200));
  hist(data.matrix(RealSamp),col=rgb(0,0,1,0.5),breaks=breaks,main ="LogOR distribution", freq=TRUE,plot=TRUE,xlab="logOR",cex.axis=2.0, cex.lab=2.0,cex.main=2.0,ylim=c(0,mc));
  hist(data.matrix(RandomSamp),add=T,col=rgb(0,1,0,0.5),breaks=breaks,main ="Dist", freq=TRUE,plot=TRUE,cex.axis=1.2);
  
  #hist(CLPar,col="red",breaks=breaks,main = "OR frequency", freq=TRUE,plot=TRUE,xlab="logOR",cex.axis=2.0, cex.lab=1.5,cex.main=1.8,ylim=c(0,200));29600
  #hist(CLMed,add=T,col=rgb(0.1,0.1,0.1,0.5),breaks=breaks,main ="CLMed", freq=TRUE,plot=TRUE,cex.axis=2.0);
  #hist(Dist,add=T,col=rgb(0.8,0.8,0.8,0.5),breaks=breaks,main ="Dist", freq=TRUE,plot=TRUE,cex.axis=2.0);
  legend("topright",c("Original", "Random"),fill=c(rgb(0,0,1,0.5),rgb(0,1,0,0.5)),cex=1.8)
  
  #legend("topright",c("ClPar", "ClMed", "Dist"),fill=c(rgb(1,0,0,0.5),rgb(0,0,1,0.5),rgb(0,1,0,0.5)),cex=2.0)
  
}