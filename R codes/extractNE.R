  #@author Matej Mihelcic
  #@institution Rudjer Boskovic Institute, Zagreb, Croatia
  #@mail matmih1@gmail.com
  #@description helper function to extract enrichments of a given GO


extractNE<-function(GOMap, ORinput, GO){
  
  out<-c();
  
  index<--1;
  
  for(i in 1:dim(GOMap)[1]){
    if(GOMap[i,1]==GO){
      index<-GOMap[i,2];
      break;
    }
  }
  
  print(index)
  for(i in 1:dim(ORinput)[1])
    if(ORinput[i,3]<1){
      #print("row found...")
      out<-rbind(out,ORinput[i,]);
    }
  
  return (out);
  
}