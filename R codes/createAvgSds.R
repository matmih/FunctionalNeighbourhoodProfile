createAvgSds <- function(data){
  
  avgs<-c();
  sds<-c();
  
  for(i in 1:dim(data)[1]){
    avg<-sum(data[i,])/length(data[i,]);
    sd<-sd(data[i,]);
    avgs<-c(avgs,avg);
    sds<-c(sds,sd);
  }
    
  
  res<-rbind(avgs,sds);
  
  return(res);
  
}