   #@author Matej Mihelcic
  #@institution Rudjer Boskovic Institute, Zagreb, Croatia
  #@mail matmih1@gmail.com
  #@description function to draw a LOR histogram of a given function divided by semantic similarity. Bins containing LORs for semanticaly close pairs to a given
  #GO are denoted by arrow. 
  #Inputs to this function are three arrays: CLPar, containing all LORs for pairs GOt - GOx (where GOx equals GOt or is its parent in GO ontology)
  # CLMed, containing all LORs for pairs GOt - GOx (where GOx is a function that is close or medium close with respect to semantical distance to GOt)
  # Dist, containing all LORs for pairs GOt - GOx (where GOx is a function that is semantically distant to GOt)
  

drawHistogramORGOFunction<-function(funcName){
  
   CLPar<- read.csv(paste("~/NetBeansProjects/Recursive file search/HistDivided",funcName,"CLPar.txt",sep=""), header=FALSE);
   CLMed<- read.csv(paste("~/NetBeansProjects/Recursive file search/HistDivided",funcName,"MedCl.txt",sep=""), header=FALSE);
   Dist<- read.csv(paste("~/NetBeansProjects/Recursive file search/HistDivided",funcName,"Dist.txt",sep=""), header=FALSE);
  
   CLPar<-data.matrix(CLPar);
   CLMed<-data.matrix(CLMed);
   Dist<-data.matrix(Dist);
   
  breaks<-seq(-3,3,0.1);
  Dist[which(Dist<(-3))]<-(-3);
  Dist[which(Dist>3)]<-(3);
  CLPar[which(CLPar>3)]<-3;
  CLPar[which(CLPar<(-3))]<-(-3);
  CLMed[which(CLMed<(-3))]<-(-3);
  CLMed[which(CLMed>3)]<-3;
  #hist(CLPar,col=rgb(1,0,0,0.5),breaks=breaks,main = "OR frequency", freq=TRUE,plot=TRUE,xlab="logOR",cex.axis=2.0, cex.lab=1.5,cex.main=1.8,ylim=c(0,200));
  hist(CLMed,col=rgb(0,0,1,0.5),breaks=breaks,main ="Log OR frequency", freq=TRUE,plot=TRUE,xlab="logOR",cex.axis=2.0, cex.lab=1.5,cex.main=1.8,ylim=c(0,70));
  hist(Dist,add=T,col=rgb(0,1,0,0.5),breaks=breaks,main ="Dist", freq=TRUE,plot=TRUE,cex.axis=2.0);
  
  #hist(CLPar,col="red",breaks=breaks,main = "log OR frequency", freq=TRUE,plot=TRUE,xlab="logOR",cex.axis=2.0, cex.lab=1.5,cex.main=1.8,ylim=c(0,200));
  #hist(CLMed,add=T,col=rgb(0.1,0.1,0.1,0.5),breaks=breaks,main ="CLMed", freq=TRUE,plot=TRUE,cex.axis=2.0);
  #hist(Dist,add=T,col=rgb(0.8,0.8,0.8,0.5),breaks=breaks,main ="Dist", freq=TRUE,plot=TRUE,cex.axis=2.0);
  legend("topright",c(expression(CLPar %->% " "),"ClMed", "Dist"),fill=c(NA,rgb(0,0,1,0.5),rgb(0,1,0,0.5)),cex=1.0)
  
  #legend("topright",c("ClPar", "ClMed", "Dist"),fill=c(rgb(1,0,0,0.5),rgb(0,0,1,0.5),rgb(0,1,0,0.5)),cex=2.0)
  
  for (i in 1:length(CLPar)){
    arrows(CLPar[i],10,CLPar[i],1)
  }
  
}