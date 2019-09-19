  #@author Matej Mihelcic
  #@institution Rudjer Boskovic Institute, Zagreb, Croatia
  #@mail matmih1@gmail.com
  #@description function to FDR adjust p-values
  
computeFDR<-function(PvalsLORAll,output){
  LorPvalVec<-PvalsLORAll$V1;
  padjLor<-p.adjust(LorPvalVec, method = "fdr", n = length(LorPvalVec));
  write(padjLor, file =output ,append=FALSE,sep="\n");
  #"C:\\Users\\matej\\Documents\\NetBeansProjects\\Recursive file search\\LORcorrectedPvalsAll.txt"
  }