   #@author Matej Mihelcic
  #@institution Rudjer Boskovic Institute, Zagreb, Croatia
  #@mail matmih1@gmail.com
  #@description function to compute the number of enriched functions for a given GO

CheckEnrichment<-function(GOMap, ContingencyOutput,GO){
  
  out<-c();

  rpval<-ComputePvals(GOMap, ContingencyOutput);
  rpvalNE<-extractNE(GOMap, rpval, GO);
  Pval<-extractPVal(GOMap, rpvalNE);
  
  outNEPvalAdj<-p.adjust(Pval, method = "fdr", n = length(Pval));
  
  out<-c(out,sum(outNEPvalAdj<0.01),sum(outNEPvalAdj>=0.01));
  
  return(out);
}