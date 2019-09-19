  #@author Matej Mihelcic
  #@institution Rudjer Boskovic Institute, Zagreb, Croatia
  #@mail matmih1@gmail.com
  #@description function to visualize the neighbourhood of an input COG

visualizeNeighbourhoods<-function(inputCOG, cogName){

library(scales)


pdf(paste(cogName,".pdf"));

# Configure graphic device
par(mai = c(0.1, 0, 0.1, 0), bg = "grey85")
# Plot creation
gc.grid <- layout(matrix(1:20, nrow = 20))
for(i in 1:dim(inputCOG)[1]){
  
  print(paste('ind:  ',i));
  inval<-grep('taxId:',inputCOG[i,],fixed = TRUE);
  inval<-as.numeric(inval);
  if(identical(inval, numeric(0)))
      inval<-0;
  
  print(paste('inval ',inval));
  
  if(inval){
    plot(c(0, 10), c(0,1), 
         bg = "grey90", 
         type = "n", 
         bty="n", 
         xaxt="n", 
         yaxt="n", xlab="", ylab="")
    text(1, 0.55, labels = inputCOG[i,], cex = 1.2, font = 4);
    
  }
  else{
  
  COGS<-unlist(strsplit(as.character(inputCOG[i,]),"\\|"));  
    print(COGS);
  print(length(COGS));
  gc.ramp <- hue_pal()(i)
  plot(c(0, length(COGS)), c(0,1), 
       bg = "grey90", 
       type = "n", 
       bty="n", 
       xaxt="n", 
       yaxt="n", xlab="", ylab="")
  for(j in 1:length(COGS)){
    inval<-grep('_C',COGS[j],fixed = TRUE);
    inval1<-grep('GOx',COGS[j], fixed = TRUE);
    inval2<-grep('GOy', COGS[j], fixed = TRUE);
    
    if(identical(inval, integer(0)))
      inval<-0;
    if(identical(inval1, integer(0)))
      inval1<-0;
    if(identical(inval2, integer(0)))
      inval2<-0;
    print(inval)
    print(paste(inval," ",inval1," ",inval2));
    
    cprint<-c();
    
    if(inval1>0 && inval>0){
      gc.ramp[j]<-'tan1';
      cprint<-strsplit(as.character(COGS[j]),"\\_C")[[1]][1];
      #print("After split");
      #print(cprint)
    }
    else if(inval2>0 && inval>0){
      gc.ramp[j]<-'thistle2';
      cprint<-strsplit(as.character(COGS[j]),"\\_C")[[1]][1];
    }
    else if(inval1>0 && inval2>0){
      gc.ramp[j]<-'green';
      cprint<-strsplit(as.character(COGS[j]),":GOx:GOy")[[1]][1]; 
    }
    else if(inval2>0){
      gc.ramp[j]<-'aquamarine';
      cprint<-strsplit(as.character(COGS[j]),":GOy")[[1]][1];
    }
    else if(inval1>0){
      gc.ramp[j]<-'slategray3';
      cprint<-strsplit(as.character(COGS[j]),":GOx")[[1]][1];
    }
    else{ 
      gc.ramp[j]<-"grey90";
      cprint<-COGS[j];
    }
    
    print(gc.ramp[j]);
    
    rect(j - 1.3, 0, j - 0.3, 1, col = gc.ramp[j])
    text(j - 0.8, 0.55, labels = cprint, cex = 0.6, font = 2)
  }
 }
}
# Reset graphical device to default values
#par(op)
dev.off();
}