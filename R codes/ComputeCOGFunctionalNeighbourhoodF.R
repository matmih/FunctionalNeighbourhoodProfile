   #@author Matej Mihelcic
  #@institution Rudjer Boskovic Institute, Zagreb, Croatia
  #@mail matmih1@gmail.com
  #@description function to compute COG functional neighbourhood

computeCOGFN<-function(input){
  
  vecFreq<-c();
  vecSd<-c();
  
  for(i in 1:dim(input)[1]){
    tmp<-as.numeric(input[i,2:length(input[i,])]);
    vecFreq<-c(vecFreq,sum(tmp[1:(length(tmp)-1)])/tmp[length(tmp)]);
    vecSd<-c(vecSd,sd(tmp[1:(length(tmp)-1)]));
  }
  
  outSummary<-data.frame(input[,1],vecFreq,vecSd);
  return(outSummary);
  
}