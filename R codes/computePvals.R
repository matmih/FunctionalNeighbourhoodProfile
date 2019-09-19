   #@author Matej Mihelcic
  #@institution Rudjer Boskovic Institute, Zagreb, Croatia
  #@mail matmih1@gmail.com
  #@description function to compute OR and pval from contingency tables

ComputePvals<-function(GOMap, ContingencyOutput){
  
  out<-c();
  
  for(i in 1:(dim(ContingencyOutput)[1]/3)){
    
    if(i==1){
      index<-1;
    }
    else{
      index<-4*(i-1)-(i-2);
    }
    
    name<-strsplit(as.character(ContingencyOutput[index,2]),"-");
    print(index)
    mapIndex<-c();
    for(j in 1:dim(GOMap)[1])
      if(GOMap[j,1]==name[[1]][[2]] || GOMap[j,1]==name[[1]][[1]]){
        mapIndex<-c(mapIndex,GOMap[j,2]);
      }
    
    tab<-extractCont(ContingencyOutput,index);
    res<-fisher.test(tab);
    OR<-res[[3]][[1]];
    pval<-res[[1]][[1]];
    print(OR)
    t<-c(mapIndex,OR,pval);
    print(t)
    out<-rbind(out,t);
  }
  
  return(out);
}