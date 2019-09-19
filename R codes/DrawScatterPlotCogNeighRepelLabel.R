 #@author Matej Mihelcic
  #@institution Rudjer Boskovic Institute, Zagreb, Croatia
  #@mail matmih1@gmail.com
  #@description function to draw scatter plot showing the average percentage of gene functional Neighborhood of different COGs that contain GO functions GO1 and GO2
  #repel label used to obtain better spatial arrangement of COG descriptions on the figure 

drawScatterRepelLabel<-function(cuttof1, cuttof2, GO1Lab,GO2Lab,CogDescription){
  
  library(ggrepel)
  set.seed(42)
  
  GO1<-read.delim(paste("~/NetBeansProjects/Recursive file search/FilesForScatterPlot/COGNeighCompPerc",GO1Lab,".txt",sep=""),strip.white = TRUE);
  GO2<- read.delim(paste("~/NetBeansProjects/Recursive file search/FilesForScatterPlot/COGNeighCompPerc",GO2Lab,".txt",sep=""),strip.white = TRUE);
  GO1Rand<- read.delim(paste("~/NetBeansProjects/Recursive file search/FilesForScatterPlot/COGNeighCompPercRand",GO1Lab,".txt",sep=""),strip.white = TRUE);
  GO2Rand<- read.delim(paste("~/NetBeansProjects/Recursive file search/FilesForScatterPlot/COGNeighCompPercRand",GO2Lab,".txt",sep=""),strip.white = TRUE);
  
  join <- data.frame(GO1$COG,GO1$mean);
  join$meanC2<-GO2$mean[match(join$GO1.COG,GO2$COG)];
  join$Desc<-CogDescription$V4[match(join$GO1.COG,CogDescription$V2)];
  facts<-sapply(1:nrow(join),function(x){if(is.na(join$Desc[x])) return(paste(join$GO1.COG[x])); desc<-as.character(join$Desc[x]);
                                         descR<-strsplit(desc,"\\|");
                                         retStr<-descR[[1]][1];
                                         # retMod<-gsub(" ","\n",retStr);
                                         retMod<-retStr;  return(paste(join$GO1.COG[x],"\n",retMod));});
  #View(facts[which(duplicated(facts))])
  join$Desc<-facts;
  join$Color<-as.factor(rep("Original",dim(join)[1]));
  rownames(join) <- join$Desc;
  #View(join);
  
  joinRand <- data.frame(GO1Rand$COG,GO1Rand$mean);
 # View(joinRand)
  joinRand$meanC2<-GO2Rand$mean[match(joinRand$GO1Rand.COG,GO2Rand$COG)];
  joinRand$Desc<-CogDescription$V4[match(joinRand$GO1Rand.COG,CogDescription$V2)];
  facts<-sapply(1:nrow(joinRand),function(x){if(is.na(joinRand$Desc[x])) return(paste(joinRand$GO1Rand.COG[x],"Rand")); desc<-as.character(join$Desc[x]);
                                             descR<-strsplit(desc,"\\|");
                                             retStr<-descR[[1]][1];
                                             # retMod<-gsub(" ","\n",retStr);
                                             retMod<-retStr; return(paste(joinRand$GO1Rand.COG[x],"Rand","\n",retMod));});
  #View(facts[which(duplicated(facts))])
  joinRand$Desc<-facts;
  joinRand$Color<-as.factor(rep("Random",dim(joinRand)[1]));
  rownames(joinRand) <- joinRand$Desc;
  #View(joinRand);
  colnames(joinRand)<-colnames(join);
  total <- rbind(join, joinRand);
 # View(total[4000:dim(total)[1],])
  
  p <- ggplot(total, aes(x=GO1.mean, y=meanC2, label = Desc, color=Color)) +
    geom_point(size=2) + scale_color_manual(values=c('#0c0504','#cc250c')) +labs(x = GO1Lab, y = GO2Lab)+theme(title = element_text(size=16, face = "bold"),axis.text=element_text(size=12),axis.title=element_text(size=14,face="bold"), axis.text.x = element_text(size = 12), axis.text.y = element_text(size = 12) ,legend.text=element_text(size=12));
    
  p2 <- p+ geom_text_repel(segment.size = 0.0,data=subset(total[1:dim(join)[1],], (GO1.mean>=cuttof1 & meanC2>=cuttof2))) + labs(title = paste(GO1Lab,"-",GO2Lab,";(",cuttof1,",",cuttof2,")"));
  
  gridExtra::grid.arrange(p2, ncol = 1);
}