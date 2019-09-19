  #@author Matej Mihelcic
  #@institution Rudjer Boskovic Institute, Zagreb, Croatia
  #@mail matmih1@gmail.com
  #@description experimentation workflow to compute AUROC and AUPRC of models for phenotype prediction

library(readr)
conditional_scores <- read_delim("C:/Users/matej/Downloads/elife-31035-supp4-v2/conditional_scores.tsv","\t",escape_double = FALSE)
a<-colnames(conditional_scores)
Targets <- read_csv("~/PhenotypeFiles/R functions/Targets.txt", col_names = FALSE)
b<-Targets[,1]
c<-setdiff(a,b)
cs1<-conditional_scores[ ,-which(names(conditional_scores) %in% c)]
rows<-c(2,4,6,7,8,9,10,11,12,14,17,21,22,23,24,26,27,28,29,30,33,35,36,37,38,42,43,45,46,47,48,49,50,51,52,53,54,56,57,59,61,62,63,64,65,66,68,69,70,71,72,74,76,77,78,79,80,81,82,83,84,87,88,89,90,91,92,93,94,95,96,97,98,99,100,101,102,103,104,105,106,107,109,110,111,112,113,114,115,116,117,118,119,120,121,122,123,124,125,126,127,129,130,131,132,133,134,135,136,137,138,139,142,143,144,145,146,147,148,149,150,152,153,154,155,156,157,158,159,160,161,162,163,164,166,167,169,170,171,173,174,176,177,178,179,180,181,182,183,184,186)
cs2<-cs1[as.numeric(row.names(cs1)) %in% rows,]
BaselinePhenotype_oob <- read_csv("~/PhenotypeFiles/Phenotype prediction/BaselinePhenotype_oob.preds")
aucBaseline<-computeAUCB(cs1,Targets,conditional_scores,187)#outputs 151x2 matrix (only usable phenotypes)
AUCBaseline <- read_csv("~/PhenotypeFiles/R functions/AUCBaseline.txt",col_names = FALSE)#AUCBaseline contains RF on conditional_scores
AUCPCABaseline <- read_csv("~/PhenotypeFiles/R functions/AUCPCABaseline.txt",col_names = FALSE)
PCABaselinePhenotypeR_oob <- read_csv("~/PhenotypeFiles/Phenotype prediction/PCABaselinePhenotypeR_oob.preds")
aucEnsembleW2<-computeAUCEnsemble(PCABaselinePhenotypeR_oob,cs2,151,0.9,0.1)
#FunctionalPCAR_oob - PCA neighbourhoods
library(vioplot)
vioplot(as.numeric(aucBaseline$aucs),as.numeric(AUCBaseline$X1),as.numeric(AUCPCABaseline$X1), as.numeric(aucEnsembleW2),names=c("Baseline","BaselineRF","PCABaselineRF","Ensemble"))


#GreedSearch Ensemble phenotypes
aucEnsemble91<-computeAUCEnsemble(PCABaselinePhenotypeR_oob,cs2,151,0.9,0.1);
View(aucEnsemble91)
aucEnsemble82<-computeAUCEnsemble(PCABaselinePhenotypeR_oob,cs2,151,0.8,0.2);
aucEnsemble73<-computeAUCEnsemble(PCABaselinePhenotypeR_oob,cs2,151,0.7,0.3);
aucEnsemble64<-computeAUCEnsemble(PCABaselinePhenotypeR_oob,cs2,151,0.6,0.4);
aucEnsemble55<-computeAUCEnsemble(PCABaselinePhenotypeR_oob,cs2,151,0.5,0.5);
aucEnsemble46<-computeAUCEnsemble(PCABaselinePhenotypeR_oob,cs2,151,0.4,0.6);
aucEnsemble37<-computeAUCEnsemble(PCABaselinePhenotypeR_oob,cs2,151,0.3,0.7);
aucEnsemble28<-computeAUCEnsemble(PCABaselinePhenotypeR_oob,cs2,151,0.2,0.8);
aucEnsemble19<-computeAUCEnsemble(PCABaselinePhenotypeR_oob,cs2,151,0.1,0.9);
par(cex = 1.4);
vioplot(as.numeric(aucEnsemble19),as.numeric(aucEnsemble28),as.numeric(aucEnsemble37),as.numeric(aucEnsemble46), as.numeric(aucEnsemble55), as.numeric(aucEnsemble64), as.numeric(aucEnsemble73), as.numeric(aucEnsemble82), as.numeric(aucEnsemble91),names=c("(0.1,0.9)","(0.2,0.8)","(0.3,0.7)","(0.4,0.6)","(0.5,0.5)","(0.6,0.4)","(0.7,0.3)","(0.8,0.2)","(0.9,0.1)"));
title(main = "Ensemble performance for different weights of base models" , ylab="AUC", xlab="Weights used to construct ensemble classifier");
source('~/ComputeAUPRCEnsemble.R')
auprcEnsemble91<-computeAUPRCEnsemble(PCABaselinePhenotypeR_oob,cs2,151,0.9,0.1);
auprcEnsemble82<-computeAUPRCEnsemble(PCABaselinePhenotypeR_oob,cs2,151,0.8,0.2);
auprcEnsemble73<-computeAUPRCEnsemble(PCABaselinePhenotypeR_oob,cs2,151,0.7,0.3);
auprcEnsemble64<-computeAUPRCEnsemble(PCABaselinePhenotypeR_oob,cs2,151,0.6,0.4);
auprcEnsemble55<-computeAUPRCEnsemble(PCABaselinePhenotypeR_oob,cs2,151,0.5,0.5);
auprcEnsemble46<-computeAUPRCEnsemble(PCABaselinePhenotypeR_oob,cs2,151,0.4,0.6);
auprcEnsemble37<-computeAUPRCEnsemble(PCABaselinePhenotypeR_oob,cs2,151,0.3,0.7);
auprcEnsemble28<-computeAUPRCEnsemble(PCABaselinePhenotypeR_oob,cs2,151,0.2,0.8);
auprcEnsemble19<-computeAUPRCEnsemble(PCABaselinePhenotypeR_oob,cs2,151,0.1,0.9);
par(cex = 1.4);
vioplot(as.numeric(auprcEnsemble19),as.numeric(auprcEnsemble28),as.numeric(auprcEnsemble37),as.numeric(auprcEnsemble46), as.numeric(auprcEnsemble55), as.numeric(auprcEnsemble64), as.numeric(auprcEnsemble73), as.numeric(auprcEnsemble82), as.numeric(auprcEnsemble91),names=c("(0.1,0.9)","(0.2,0.8)","(0.3,0.7)","(0.4,0.6)","(0.5,0.5)","(0.6,0.4)","(0.7,0.3)","(0.8,0.2)","(0.9,0.1)"));
title(main = "Ensemble performance for different weights of base models" , ylab="AUPRC", xlab="Weights used to construct ensemble classifier");
