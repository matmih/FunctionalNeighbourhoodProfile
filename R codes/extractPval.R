   #@author Matej Mihelcic
  #@institution Rudjer Boskovic Institute, Zagreb, Croatia
  #@mail matmih1@gmail.com
  #@description helper function needed in FDR adjust

extractPVal<-function(GOMap, ORinput){
  
  out<-c();
  
  for(i in 1:dim(ORinput)[1])
      out<-c(out,ORinput[i,4]);
  
  return (out);
  
}