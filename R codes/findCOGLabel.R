  #@author Matej Mihelcic
  #@institution Rudjer Boskovic Institute, Zagreb, Croatia
  #@mail matmih1@gmail.com
  #@description helper function to find COG description by using COG code. Function returns modified description string

findCOGLabel <- function(COGCode,COGSet)
{ 
  cog<-as.character(COGCode);
  cog<-gsub(" ","",cog);
  
  for(i in 1:dim(COGSet)[1]){
    
    cog1<-as.character(COGSet[i,2]);
    cog1<-gsub(" ","",cog1);
    cog1<-as.character(cog1);
    
    resp<-grep(cog,cog1);
    
    if(length(grep(cog,cog1))!=0){
        desc<-as.character(COGSet[i,4]);
        descR<-strsplit(desc,"\\|");
        retStr<-descR[[1]][1];
       retMod<-retStr;
        return(retMod);
    }   
  }
  
  rm<-"";
  
  return(rm);
}