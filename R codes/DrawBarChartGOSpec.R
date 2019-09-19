  #@author Matej Mihelcic
  #@institution Rudjer Boskovic Institute, Zagreb, Croatia
  #@mail matmih1@gmail.com
  #@description function to draw a barchart showing the frequency of GO functions with some level of semantical similarity (it also displays the fraction of 
  #signficantly enriched pairs.


drawBarchartGOSpec<-function(input){
  
  data <- read.delim(paste("~/NetBeansProjects/Recursive file search/DensityOverlayBySignif",input,".txt",sep=""));
  
  tmpOR<-data$logOR;
  h<-hist(tmpOR,plot=FALSE);
  intVec<-c();
  print(h$breaks)
  print(max(data$logOR))
  print(min(data$logOR))
  print(length(tmpOR))
  orderVec<-c();
  levels1<-c();
  
  for(j in 1:(length(h$breaks)-1)){
    levels1<-c(levels1,paste("[ ",h$breaks[j]," - ",h$breaks[j+1]," ]"));
  }
    
  
  for(i in 1:dim(data)[1]){
    count<-0;
      for(j in 1:(length(h$breaks)-1)){
          if(tmpOR[i]>=h$breaks[j] && tmpOR[i]<h$breaks[j+1]){
            intVec<-c(intVec,paste("[ ",h$breaks[j]," - ",h$breaks[j+1]," ]"));
            count<-1;
            orderVec<-c(orderVec,j);
            break;
          }
      }
      if(tmpOR[i]==h$breaks[length(h$breaks)]){
        intVec<-c(intVec,paste("[ ",h$breaks[length(h$breaks)-1]," - ",h$breaks[length(h$breaks)]," ]"));
        orderVec<-c(orderVec,length(h$breaks));
      }
        
    if(count==0)
      print(tmpOR[i])
  }
  
  TmpData<-data;
  TmpData$intervals <- intVec;
  TmpData$order<-orderVec;
  library(ggplot2)
  
  ggplot(transform(TmpData,
                   intervals=factor(intervals,levels=levels1)), aes(x = hist, fill = category)) + geom_bar(stat = 'count', position = 'stack') + facet_grid(~ intervals);
  
}