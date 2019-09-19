   #@author Matej Mihelcic
  #@institution Rudjer Boskovic Institute, Zagreb, Croatia
  #@mail matmih1@gmail.com
  #@description a series of commands used to plot different boxplots, comparative violin plots, RUMI curves and other

boxplot(as.vector(supp1), as.vector(supp2), as.vector(supp3),as.vector(supp4),horizontal=TRUE, names=c("[5,10]", "[11,39]","[40,99]","[100,300]"),xlab="Entropy",ylab="Support interval",col=c("red","orange","cyan","blue"),main="Entropy of level of dementia in described patients")

pl <- barchart(newsuppSmall4,scales=list(cex=1.5, cex.axis=3),
               +                stack=TRUE, xlab = "Count",auto.key = list(adj = 1,title="Redescriptions with support in [100,300]",cex=2),box.width = 0.1,cex.label=2)
> plot(pl)
> pl <- barchart(newsuppSmall4,scales=list(cex=1.5, cex.axis=3),
                 +                stack=TRUE, xlab = list("Count",cex=2),auto.key = list(adj = 1,title="Redescriptions with support in [100,300]",cex=2),box.width = 0.1)
> plot(pl)


#barPlotForMultipleFUnctions
#vectore rucno napuniti
funcTmpDataN<-rbind(NN1Func,NN1FuncW,NN10FuncW,GFP1,RFDist,RFFA);
barplot(funcTmpData,beside=T,col=c("red","red","red","blue","pink","green"),xlab="Functions",names=c("GO0046352","GO0043447","GO0006694","GO0006631","GO0000003","GO0043093","GO0034622","GO0042493","GO0007610","GO0034655","GO0008610","GO0060255","GO0019725","GO0003333","GO0032502"),ylab="AUPRC",legend=c("1-NN","5-NN","10-NN","GFP","RFDist","RFFN"),args.legend= list(title = "Methods", x = "topright", cex = .7))

#newCode
BarPlotComparativeData <- read.table("~/NetBeansProjects/Recursive file search/BarPlotComparativeData.txt", quote="\"")
barplot(data.matrix(BarPlotComparativeData),beside=T,col=c("red","red3","red4","blue","green"),xlab="Functions",names=c("5975","6260","6974","16051","6457","46700","6310","46903","8610","6865","6508"),ylab="AUPRC",legend=c("1-NN","3-NN","10-NN","GFP","FGN"),args.legend= list(title = "Methods", x = "topright", cex = .7))
barplot(data.matrix(BarPlotComparativeDataAUC),beside=T,col=c("red","red3","red4","blue","green"),xlab="Functions",names=c("5975","6260","6974","16051","6457","46700","6310","46903","8610","6865","6508"),ylab="AUC",legend=c("1-NN","3-NN","10-NN","GFP","FGN"),args.legend= list(title = "Methods", x = "top", cex = .5,horiz = TRUE,inset = -0.15))

#violin plot code
my.vioplot1(MethodViolinComparisonR[,1],MethodViolinComparisonR[,2],MethodViolinComparisonR[,3],MethodViolinComparisonR[,4],MethodViolinComparisonR[,5],MethodViolinComparisonR[,6],col=c("red","red","red","blue","pink","green"),names=c("1-NN","5-NN","10-NN","GFP","RFL","RFFN"))
my.violplot1(`genome5LogK=1`[,3],`genome5LogK=3`[,3],`genome5LogK=10`[,3],GFPResultsLab[,3],`BaselineOGsk=4_600_`[,3],col=c("red","red","red","blue","green"),names=c("1-NN","3-NN","10-NN","GFP","NFP"),cax=2)
title(main="Method comparisson", xlab="", ylab="AUPRC", cex.lab = 1.5,cex.main = 2.0)

#compute RUMI
source('C:/Users/matej/Downloads/SemDist-master/SemDist-master/R/findRuMi.R')
source('C:/Users/matej/Downloads/SemDist-master/SemDist-master/R/myviolplot1.R')#load last
source('C:/Users/matej/Downloads/SemDist-master/SemDist-master/R/computeIA.R')
source('C:/Users/matej/Downloads/SemDist-master/SemDist-master/R/utilities.R')
source('C:/Users/matej/Downloads/SemDist-master/SemDist-master/R/gene2GO.R')
LabelsCurve <- read.delim("~/NetBeansProjects/Recursive file search/LabelsCurve.txt", header=FALSE)
PredsCurve <- read.delim("~/NetBeansProjects/Recursive file search/PredsCurve.txt", header=FALSE)
Iacr<-load("C:\\Users\\matej\\Downloads\\SemDist-master\\SemDist-master\\inst\\extdata\\myIA.rda")
avgRUMIvals <- RUMIcurve("MF", "human", 0.05, LabelsCurve, PredsCurve,IAccr = IA)
avgRUMIvals4 <- RUMIcurve("MF", "human", 0.05, "~/NetBeansProjects/Recursive file search/LabelsCurveFungiGFP.txt", "~/NetBeansProjects/Recursive file search/PredsCurveFungiGFP.txt",IAccr = IA)
firstset<- avgRUMIvals[[1]]
plot/points(firstsetx$RU,firstsetx$MI,type="o",col="purple")

#RUMICurve
plot(firstset8$RU, firstset8$MI,type="o",col="green")
points(firstset6$RU,firstset6$MI,type="o",col="blue")
points(firstset$RU,firstset$MI,type="o",col="red")
points(firstset3$RU,firstset3$MI,type="o",col="black")
points(firstset5$RU,firstset5$MI,type="o",col="orange")
points(firstset4$RU,firstset4$MI,type="o",col="purple")
legend("topright", legend = c("1-NN","5-NN","10-NN","GFP", "RFL","RFFN"), fill=c("red","blue","green","orange","purple","black"))

#New RUMICurve
plot(firstset2$RU,firstset2$MI,type="b",col="green",xlim = c(3,14),ylab="Misinformation",xlab="Remaining uncertainty",cex = 1.0, lwd = 1.6,lty=1,pch=19,cex.axis = 1.4, cex.lab = 1.7,main = "Ru-Mi curve",cex.main = 2.0)
> points(firstset1$RU,firstset1$MI,type="b",col="blue",cex=1.0,lwd=1.6,lty=2,pch=17)
> points(firstset$RU,firstset$MI,type="b",col="red",cex=1.0,lwd=1.6,lty=3,pch=18)
> legend("topright",legend=c("10-NN", "GFP","FGN"),
         +        col=c("green","blue","red"),lty=1:3 ,cex=0.8)

#compute accretions eukaryots
ontFile<-"~/NetBeansProjects/Recursive file search/ontologyIA.txt";
annotations<-"~/NetBeansProjects/Recursive file search/labelsFungiIA.txt";
IAT<-computeIA("my", "values", specify.ont=TRUE, myont=ontFile, specify.annotations=TRUE, annotfile=annotations)

#compute accretions eukaryots

IACFungi <- read.table("~/NetBeansProjects/Recursive file search/IACFungi.txt", quote="\"")
IACFungi1<-data.matrix(IACFungi);
names(IACFungi2)<-rownames(IACFungi1);
IACFungi2<-as.numeric(IACFungi1);

avgRUMIvals3F <- RUMIcurve("human", "BP", 0.05, "~/NetBeansProjects/Recursive file search/LabelsCurveFungiFNP.txt", "~/NetBeansProjects/Recursive file search/PredsCurveFungiFNP.txt",IAccr = IACFungi2)

#Different k plot
plot(c(1,2,3,4,5,6,7,8,9,10),knnauprc,type="o",ylim=c(0.1,0.4),ylab="AUPRC",xlab="Number of neighbours",cex = 1.0, lwd = 1.6,lty=1,pch=19,cex.axis = 1.4, cex.lab = 1.5,main = "FGN and k-NN for different neighbourhood size",cex.main = 2.0)
lines(c(1,2,3,4,5,6,7,8,9,10),gfpauprc, type="o", pch=8,lty=3,col="green")
lines(c(1,2,3,4,5,6,7,8,9,10),fnauprc, type="o", col="red",cex=1.0,lwd=1.6,lty=2,pch=17)
legend("topright",legend=c("GFN","GFP","k-NN"),col=c("red","green","black"),lty=c(2,3,1),pch=c(17,8,19) ,cex=0.8)

plot(c(1,2,3,4,5,6,7,8,9,10),knnauc,type="o",ylim=c(0.5,1.0),ylab="AUC",xlab="Number of neighbours",cex = 1.0, lwd = 1.6,lty=1,pch=19,cex.axis = 1.4, cex.lab = 1.5,main = "FGN and k-NN for different neighbourhood size",cex.main = 2.0)
lines(c(1,2,3,4,5,6,7,8,9,10),gfpauc, type="o", pch=8,lty=3,col="green")
lines(c(1,2,3,4,5,6,7,8,9,10),fnaauc, type="o", col="red",cex=1.0,lwd=1.6,lty=2,pch=17)
legend("topright",legend=c("GFN","GFP","k-NN"),col=c("red","green","black"),lty=c(2,3,1),pch=c(17,8,19) ,cex=0.8)

plot(c(1,2,3,4,5,6,7,8,9,10),knnauprc,type="o",pch=0,ylim=c(0.1,0.4),ylab="AUPRC",xlab="Number of neighbours")
lines(c(1,2,3,4,5,6,7,8,9,10),gfpauprc, type="o", pch=8,lty=3,col="green")
lines(c(1,2,3,4,5,6,7,8,9,10),fnauprc, type="o", pch=1,lty=2,col="red")
legend("topright",lty=c(2,3,1),pch=c(1,8,0),c("GFN","GFP","k-NN"),col=c("red","green","black"))

#GFP different k - plot
gfpauc<-c(0.539715,0.619221,0.64985,0.66269,0.677192,0.686052,0.6906167,0.691925,0.696567,0.700804);
gfpauprc<-c(0.176541,0.202501,0.201745,0.199670,0.199103,0.197422,0.197731188,0.195185,0.1946597,0.194162);
plot(c(1,2,3,4,5,6,7,8,9,10),gfpauprc,type="o",ylim=c(0.1,0.4),ylab="AUPRC",xlab="Number of neighbours",cex = 1.0, lwd = 1.6,lty=1,pch=19,cex.axis = 1.4, cex.lab = 1.5,main = "GFP for different neighbourhood size",cex.main = 2.0)
plot(c(1,2,3,4,5,6,7,8,9,10),gfpauc,type="o",ylim=c(0.1,1.0),ylab="AUC",xlab="Number of neighbours",cex = 1.0, lwd = 1.6,lty=1,pch=19,cex.axis = 1.4, cex.lab = 1.5,main = "GFP for different neighbourhood size",cex.main = 2.0)


#Plot enrichment histograms
bp1<-barplot(vec1,beside=T,col=rainbow(10),ylim=c(0,0.18),legend=c("CLPar","CL","MED","DIST","AllBP","ORAll","ORNCLPar","ORNCL","DIFFNCLPar","DIFFNCL"),args.legend= list(title = "Divisions", x = "topleft", cex = .7),names=c("46352"),ylab="AUPRC",xlab="Functions")
error.bar(bp1,vec1,vecsdSE)

sdss<-createavgSds(d3); #d3 is a input file
y=matrix(sdss[1,],c(11,19));
bp1<-barplot(y,beside=T,col=rainbow(11),ylim=c(0,1.0),legend=c("CLPar","CL/CLPar","CL","MED","DIST","AllBP","ORAll","ORNCLPar","ORNCL","DIFFNCLPar","DIFFNCL"),args.legend= list(title = "Methods", x = "topright", cex = .7),names=c("46352","43447","6694","6631","3","43093","34622","42493","7610","34655","8610","60255","19725","3333","32502","9266","51128","60255","272"),ylab="AUPRC",xlab="Functions")
error.bar(bp1,sdss[1,],sdss[2,]/sqrt(200))

#new plot (real)
#SumaryNN -> najnovija verzijaž
sdss<-createAvgSds(SummaryN); #d3 is a input file
y=matrix(sdss[1,],c(11,19));
bp1<-barplot(y,beside=T,col=c('red2','sienna1','sienna3','sienna4','seagreen1','seagreen3','seagreen4','royalblue1','royalblue3','royalblue4','snow3'),ylim=c(0,0.3),legend=c("BPP","CL","CL/CLPar","CL/CLPar/Enr","Med","Med/Par","Med/Par/Enr","Dist","Dist/Par","Dist/Par/Enr","ClPar"),args.legend= list(title = "Methods", x = "topright", cex = .7),names=c("5975","6260","6974","16051","6457","46700","6310","46903","8610","6865","6508"),ylab="AUPRC",xlab="Functions")
error.bar(bp1,sdss[1,],sdss[2,]/sqrt(200))

#new plot short
sdss1<- sdss[,-c(2,5,8,13,16,19,24,27,30,35,38,41,46,49,52,57,60,63,68,71,74,79,82,85,90,93,96,101,104,107,112,115,118)];
y=matrix(sdss1[1,],c(8,19));
bp1<-barplot(y,beside=T,col=c('red2','sienna3','sienna4','seagreen3','seagreen4','royalblue3','royalblue4','snow3'),ylim=c(0,0.3),legend=c("BPP","CL/CLPar","CL/CLPar/Enr","Med/Par","Med/Par/Enr","Dist/Par","Dist/Par/Enr","ClPar"),args.legend= list(title = "Methods", x = "topright", cex = .7),names=c("5975","6260","6974","16051","6457","46700","6310","46903","8610","6865","6508"),ylab="AUPRC",xlab="Functions")
error.bar(bp1,sdss1[1,],sdss1[2,]/sqrt(200))

#less criteria
sdss1<-createAvgSds(SummaryNewSCComop);
y=matrix(sdss1[1,],c(5,20));
bp1<-barplot(y,beside=T,col=rainbow(5),ylim=c(0,1.0),legend=c("AllBP","CL","Med","Dist","CLPar"),args.legend= list(title = "Methods", x = "topright", cex = .7),names=c("44255","45184","6418","6400","5975","8360","6281","6865","6260","6457","271","910","9306","902","6310","6508","19439","16051","6974","46365"),ylab="AUPRC",xlab="Functions")
error.bar(bp1,sdss1[1,],1.96*sdss1[2,]/sqrt(200))


sdss1<-createAvgSds(SummaryNewSCComop);
y=matrix(sdss1[1,],c(5,20));
bp1<-barplot(y,beside=T,col=rainbow(5),ylim=c(0,1.0),legend=c("AllBP","CL","Med","Dist","CLPar"),args.legend= list(title = "Methods", x = "topright", cex = .7),names=c("44255","45184","6418","6400","5975","8360","6281","6865","6260","6457","271","910","9306","902","6310","6508","19439","16051","6974","46365"),ylab="AUPRC",xlab="Functions")
error.bar(bp1,sdss1[1,],sdss1[2,]/sqrt(200))

sdss2<-createAvgSds(SummaryExtended2);
y=matrix(sdss2[1,],c(11,20));
bp2<-barplot(y,beside=T,col=rainbow(11),ylim=c(0,max(SummaryExtended2)),legend=c("AllBP","CL","CL/Par","CL/Par/Enr","Med","Med/Par","Med/Par/Enr","Dist","Dist/Par","Dist/Par/Enr","CLPar"),args.legend= list(title = "Methods", x = "topright", cex = .7),names=c("44255","45184","6418","6400","5975","8360","6281","6865","6260","6457","271","910","9306","902","6310","6508","19439","16051","6974","46365"),ylab="AUPRC",xlab="Functions")
error.bar(bp2,sdss2[1,],sdss2[2,]/sqrt(200))
error.bar(bp2,sdss2[1,],1.96*sdss2[2,]/sqrt(200))

#plot plots for different ontologies
#Biological process
originalAUPRC<-c(0.2986163937759632,0.2979935679825981,0.3061824336921997,0.3082397001826956);
randomizedAUPRC<-c(0.2986163937759632,0.296215501983633,0.282492118385087,0.2797381423111656);
plot(c(0.0,0.1,0.6,1.0),originalAUPRC,type="o",pch=0,ylim=c(0.28,0.32),ylab="AUPRC",xlab="Jaccard threshold",lwd = 2, cex.lab=1.6, cex.axis=2.0)
lines(c(0.0,0.1,0.6,1.0),randomizedAUPRC, type="o", pch=1,lty=2,col="red", lwd = 2)
legend("topright",lty=c(1,2),pch=c(0,1),c("Original","Randomized"),col=c("black","red"))

originalAUC<-c(0.8414434726969463,0.8425402832620659,0.8416545761045998,0.8423779225891278);
randomizedAUC<-c(0.8414434726969463,0.8404630094468212,0.8374890883590139,0.8255733789055387);
plot(c(0.0,0.1,0.6,1.0),originalAUC,type="o",pch=0,ylim=c(0.82,0.85),ylab="AUC",xlab="Jaccard threshold", lwd = 2, cex.lab=1.6, cex.axis=2.0)
lines(c(0.0,0.1,0.6,1.0),randomizedAUC, type="o", pch=1,lty=2,col="red", lwd = 2)
legend("bottomleft",lty=c(1,2),pch=c(0,1),c("Original","Randomized"),col=c("black","red"))

#Molecular function

originalAUPRC<-c(0.21222013978253765,0.2133542413073231,0.2226658718168784,0.21659687371905162);
randomizedAUPRC<-c(0.2080786103153154,0.19053910263537005,0.1791671141377158,0.1727316673912004);
plot(c(0.0,0.1,0.6,1.0),originalAUPRC,type="o",pch=0,ylim=c(0.17,0.24),ylab="AUPRC",xlab="Jaccard threshold",lwd = 2,cex.lab=1.6, cex.axis=2.0)
lines(c(0.0,0.1,0.6,1.0),randomizedAUPRC, type="o", pch=1,lty=2,col="red",lwd = 2)
legend("topright",lty=c(1,2),pch=c(0,1),c("Original","Randomized"),col=c("black","red"))

originalAUC<-c(0.7901059605391002,0.7837015961911604,0.7820266603797988,0.7856218844704365);
randomizedAUC<-c(0.7860481823799709,0.7687507926211204,0.7486124699836294,0.745466335053871);
plot(c(0.0,0.1,0.6,1.0),originalAUC,type="o",pch=0,ylim=c(0.74,0.79),ylab="AUC",xlab="Jaccard threshold",cex.lab=1.6, cex.axis=2.0,lwd = 2)
lines(c(0.0,0.1,0.6,1.0),randomizedAUC, type="o", pch=1,lty=2,col="red",lwd = 2)
legend("bottomleft",lty=c(1,2),pch=c(0,1),c("Original","Randomized"),col=c("black","red"))


#Cellular component

originalAUPRC<-c(0.4468688291647663,0.4615138487227128,0.4869685298132838,0.48795224513392854);
randomizedAUPRC<-c(0.4252213900342395,0.4008802678065932,0.37709735130832167,0.38677508082042983);
plot(c(0.0,0.1,0.6,1.0),originalAUPRC,type="o",pch=0,ylim=c(0.37,0.49),ylab="AUPRC",xlab="Jaccard threshold",cex.lab=1.6, cex.axis=2.0,lwd = 2)
lines(c(0.0,0.1,0.6,1.0),randomizedAUPRC, type="o", pch=1,lty=2,col="red",lwd = 2)
legend("center",lty=c(1,2),pch=c(0,1),c("Original","Randomized"),col=c("black","red"))

originalAUC<-c(0.8719995775137932,0.8725252990465598,0.8732698975655642,0.8791359915031423);
randomizedAUC<-c(0.8612858671127437,0.8322433201833318,0.8413773473112495,0.8326676619328516);
plot(c(0.0,0.1,0.6,1.0),originalAUC,type="o",pch=0,ylim=c(0.83,0.88),ylab="AUC",xlab="Jaccard threshold",cex.lab=1.6, cex.axis=2.0,lwd = 2)
lines(c(0.0,0.1,0.6,1.0),randomizedAUC, type="o", pch=1,lty=2,col="red",lwd = 2)
legend("center",lty=c(1,2),pch=c(0,1),c("Original","Randomized"),col=c("black","red"))

#test correlation between JI and log OR

SignificantPairsLORG0AllS <- read.delim("~/NetBeansProjects/Recursive file search/SignificantPairsLORG0AllS.txt")
logORs<-as.numeric(SignificantPairsLORG0AllS$LOGOR.GO1.GO2);
Jaccards <- as.numeric(SignificantPairsLORG0AllS$JaccI_30);
cor.test(logORs,Jaccards,method = "pearson")

plot(logORs, Jaccards,ylab = "Jaccard index",xlab = "Log Odds Ratio", cex = 1.5, main = "Correlation between JI and LOR for significantly enriched GO pairs FDR<0.2", pch=20,col=rgb(0,0,0,alpha=0.1) )

#Plot histograms of enriched disimilar functions
a<-table(SignificantPairsDifferentPDist$V1);
v<-rep(0,131);
ab<-c(a,v);
hist(a, breaks = 75);
aMetazoa<-table(SignificantPairsMetazoaDiffPairsDist$V1)
hist(aMetazoa,75, main = "GO function enrichment with dissimilar functions from Biological process Gene ontology graph", xlim = c(100, max(aMetazoa)), xlab = "Number of enriched dissimilar functions", ylab = "Number of GO terms", cex = 2.0, col = "green", cex.lab = 1.5, cex.axis = 1.5);

#scatter
SignificantPairsLORG0All <- read.delim("~/NetBeansProjects/Recursive file search/SignificantPairsLORG0All.txt")
logORs<-as.numeric(SignificantPairsLORG0All$LOGOR.GO1.GO2);
plot(logORs, Jaccards,ylab = "Jaccard index",xlab = "Log Odds Ratio", cex = 1.5, main = "Correlation between JI and LOR for significantly enriched GO pairs FDR<0.2", pch=20,col=rgb(0,0,0,alpha=0.1) )

#draw venn diagram
library("VennDiagram")
draw.pairwise.venn(23.19, 1.301, 1.156, category = c("Known annotations", "Predicted annotations"), lty = rep("blank",2), fill = c("light blue", "pink"), alpha = rep(0.5,2), cat.pos = c(0,0), cat.dist = rep(0.025,2))

#plot cdf for qs
library(lattice)
library(latticeExtra)
ecdfplot(ltmp)

#drawing redescription set table
library(lattice)
supSmall1F<-createRedescriptionTableF(RowLabels,DementiaLevels,`5to10shortDist`);
newsuppSmall1<-supSmall1F[1:9,];
pl <- barchart(supSmall1F[2:10,],scales=list(cex=1.2),
nwsP1<-newsuppSmall1;
pl <- barchart(nwsP1,scales=list(cex=1.3),stack=TRUE, xlab = list("Count",cex=1.5),auto.key = list(adj = 1,title="Redescriptions with support in [100,300]",cex=1.5),box.width = 0.1)
rownames(nwsP1)[4]<-"MidTemp >= 11204.0 <= 15483.0 AND HMT7 >= 4.14 <= 6.58 AND Entorhinal >= 1443.0 <= 2170.0 \n EcogPtOrgan >= 0.28 <= 1.63 AND attentionMoca >= 0.0 <= 0.62 \n JS: 1.0, p-value: 2.7e-13"
