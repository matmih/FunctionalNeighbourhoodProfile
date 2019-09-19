  #@author Matej Mihelcic
  #@institution Rudjer Boskovic Institute, Zagreb, Croatia
  #@mail matmih1@gmail.com
  #@description helper function to join tables obtained by different models trained on CAFFA dataset into one single table

res<-computePvalsCaffa(ClassifierScoresCaffa,TrueLabelsCaffa);
OOBScores <- read.table("~/NetBeansProjects/Recursive file search/OOBScores.txt", quote="\"")
res<-left_join(OOBScores, AUCCaffa, by= c("V1"));
library(dplyr)
colnames(res)<-c("GO","AUCCaffa","AUPRCCaffa","AUCOOB","AUPRCOOB");
drs<-rowSums(TrueLabelsCaffa[,2:dim(TrueLabelsCaffa)[2]]);
resCount<-data.frame(d,drs); #load countsCaffa.txt
resF<-left_join(res,CountsCaffa, by= c("GO"));
write.table(resF,file="CaffaEvaluation.txt",col.names=TRUE,row.names=FALSE)



res<-left_join(LOOBScores, res, by= c("V1"));