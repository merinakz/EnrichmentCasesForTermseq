#!/usr/bin/env Rscript
rm(list = ls())
library(dplyr) 
library(tibble)

setwd("~/Desktop/Stuff/AMeyer_Lab/Experiment_results/3seq/QuinnEnrichment/Ksg/")
RNAfold <- read.csv("~/Desktop/Stuff/AMeyer_Lab/Experiment_results/3seq/QuinnEnrichment/Ksg/RNAfoldCalculated-Ksg.csv")
RNAfold <- RNAfold[,-c(2,3)]
top <- RNAfold %>% filter(strand == "+")
comp <- RNAfold %>% filter(strand == "-")
comp <- comp %>% add_column(TSSUTR ="", .after = "PeakCoord") %>% 
  add_column(Locus ="", .after = "TSSUTR") %>% 
  add_column(termUTR ="", .after = "Locus")
top <- top %>% add_column(TSSUTR ="", .after = "PeakCoord")%>% 
  add_column(Locus ="", .after = "TSSUTR") %>% 
  add_column(termUTR ="", .after = "Locus")
####################################
#Read your annotation file
#################################### 
genes <- read.delim("~/Desktop/Stuff/AMeyer_Lab/Experiment_results/3seq/QuinnEnrichment/NC_003028.v3.17.ncrna.genes", header = F)
genes <- genes %>% add_column(UTR="", .before = "V2") %>% add_column(termUTR = "", .after = "V3")

for(i in 1:nrow(genes)){
  genes$UTR[i] <- genes$V2[i] - 150
  genes$termUTR[i] <- genes$V3[i] + 150
}

genes <- genes[complete.cases(genes),]
genes.top <- genes %>% filter(V4 == "+")
colnames(genes.top) <- c("genome", "UTR", "from", "to","termUTR", "strand", "name", "old.name", "new.name", "WP", "rfam", "other" )
genes.comp <- genes %>% filter(V4 == "-")
colnames(genes.comp) <- c("genome", "termUTR", "to", "from","UTR", "strand", "name", "old.name", "new.name", "WP", "rfam", "other" )

genes.top$UTR <- as.numeric(genes.top$UTR)
genes.top$from <- as.numeric(genes.top$from)
genes.top$to <- as.numeric(genes.top$to)
genes.top$termUTR <- as.numeric(genes.top$termUTR)
genes.top <- genes.top[complete.cases(genes.top),]
genes.comp$UTR <- as.numeric(genes.comp$UTR)
genes.comp$from <- as.numeric(genes.comp$from)
genes.comp$to <- as.numeric(genes.comp$to)
genes.comp$termUTR <- as.numeric(genes.comp$termUTR)
genes.comp <- genes.comp[complete.cases(genes.comp),]
comp$PeakCoord <- as.numeric(comp$PeakCoord)
top$PeakCoord <- as.numeric(top$PeakCoord)

for(i in 2:nrow(genes.top)){
  if(genes.top$UTR[i] < genes.top$to[i-1]){
    genes.top$UTR[i] <- genes.top$to[i-1] + 50
  }
  if(genes.top$termUTR[i-1] > genes.top$from[i]){
    genes.top$termUTR[i-1] <- genes.top$from[i] - 50
  }
}

for(i in 2:nrow(genes.comp)){
  if(genes.comp$UTR[i-1] > genes.comp$to[i]){
    genes.comp$UTR[i-1] <- genes.comp$to[i] - 50
  }
  if(genes.comp$termUTR[i] < genes.comp$from[i-1]){
    genes.comp$termUTR[i] <- genes.comp$from[i-1] + 50
  }
}

rownames(genes.top) <- 1:nrow(genes.top)
rownames(genes.comp) <- 1:nrow(genes.comp)
####################################
#Annotate your file
####################################
for(i in 1:nrow(top)){
  for(j in 1:nrow(genes.top)){
    if((top$PeakCoord[i] >= genes.top$from[j]) & (top$PeakCoord[i] <= genes.top$to[j])){
      top$Locus[i] <- genes.top$old.name[j]
    }
    if((top$PeakCoord[i] >= genes.top$to[j]) & (top$PeakCoord[i] <= genes.top$termUTR[j])){
      top$termUTR[i] <- genes.top$old.name[j]
    }
    if((top$PeakCoord[i] >= genes.top$UTR[j]) & (top$PeakCoord[i] <= genes.top$from[j])){
      top$TSSUTR[i] <- genes.top$old.name[j]
    }
  }
}
for(i in 1:nrow(comp)){
  for(j in 1:nrow(genes.comp)){
    if((comp$PeakCoord[i] <= genes.comp$from[j]) & (comp$PeakCoord[i] >= genes.comp$to[j])){
      comp$Locus[i] <- genes.comp$old.name[j]
    }
    if((comp$PeakCoord[i] <= genes.comp$to[j]) & (comp$PeakCoord[i] >= genes.comp$termUTR[j])){
      comp$termUTR[i] <- genes.comp$old.name[j]
    }
    if((comp$PeakCoord[i] <= genes.comp$UTR[j]) & (comp$PeakCoord[i] >= genes.comp$from[j])){
      comp$TSSUTR[i] <- genes.comp$old.name[j]
    }
  }
}

#############################################################
############### Gather Cases - TOP STRAND ###################
############################################################# 
case1.top <-  top %>% filter(Sig_Control == 1 & Sig_Cond1 == 1 & Sig_Cond2 == 1)
case2.top <- top %>% filter(Sig_Control == 0 & Sig_Cond1 == 0 & Sig_Cond2 == 0)
case1.top <- rbind(case1.top, case2.top)
rm(case2.top)
case2.top <- top %>% filter(Sig_Control == 1 & Sig_Cond1 == 0 & Sig_Cond2 == 0)
case3.top <- top %>% filter(Sig_Control == 0 & Sig_Cond1 == 1 & Sig_Cond2 == 0)
case4.top <- top %>% filter(Sig_Control == 0 & Sig_Cond1 == 0 & Sig_Cond2 == 1)
case5.top <-  top %>% filter(Sig_Control == 0 & Sig_Cond1 == 1 & Sig_Cond2 == 1)
#############################################################
############### Gather Cases - COMP STRAND ##################
############################################################# 
case1.comp <-  comp %>% filter(Sig_Control == 1 & Sig_Cond1 == 1 & Sig_Cond2 == 1)
case2.comp <- comp %>% filter(Sig_Control == 0 & Sig_Cond1 == 0 & Sig_Cond2 == 0)
case1.comp <- rbind(case1.comp, case2.comp)
rm(case2.comp)
case2.comp <- comp %>% filter(Sig_Control == 1 & Sig_Cond1 == 0 & Sig_Cond2 == 0)
case3.comp <- comp %>% filter(Sig_Control == 0 & Sig_Cond1 == 1 & Sig_Cond2 == 0)
case4.comp <- comp %>% filter(Sig_Control == 0 & Sig_Cond1 == 0 & Sig_Cond2 == 1)
case5.comp <-  comp %>% filter(Sig_Control == 0 & Sig_Cond1 == 1 & Sig_Cond2 == 1)

#########################################################
##########  Case 1 - Transcriptome Distribution #########
#########################################################
############  Top

genic.case1.top <- case1.top[((!case1.top$Locus == "") & (case1.top$TSSUTR == "") & (case1.top$termUTR == "")),]
TssUTR.case1.top <- case1.top[((!case1.top$TSSUTR == "") & (case1.top$Locus == "") & (case1.top$termUTR == "")),]
intergenic.case1.top <- case1.top[((case1.top$Locus == "") & (case1.top$TSSUTR == "") & (case1.top$termUTR == "")),]
termUTR.case1.top <- case1.top[((case1.top$Locus == "") & (case1.top$TSSUTR == "") & (!case1.top$termUTR == "")),]
other1 <- case1.top[((!case1.top$Locus == "") & (!case1.top$TSSUTR == "") & (!case1.top$termUTR == "")),]
other2 <- case1.top[((case1.top$Locus == "") & (!case1.top$TSSUTR == "") & (!case1.top$termUTR == "")),]
other3 <- case1.top[((!case1.top$Locus == "") & (case1.top$TSSUTR == "") & (!case1.top$termUTR == "")),]
other4 <- case1.top[((!case1.top$Locus == "") & (!case1.top$TSSUTR == "") & (case1.top$termUTR == "")),]
other.case1.top <- rbind(other1, other2, other3, other4)

############  Comp

genic.case1.comp <- case1.comp[((!case1.comp$Locus == "") & (case1.comp$TSSUTR == "") & (case1.comp$termUTR == "")),]
TssUTR.case1.comp <- case1.comp[((!case1.comp$TSSUTR == "") & (case1.comp$Locus == "") & (case1.comp$termUTR == "")),]
intergenic.case1.comp <- case1.comp[((case1.comp$Locus == "") & (case1.comp$TSSUTR == "") & (case1.comp$termUTR == "")),]
termUTR.case1.comp <- case1.comp[((case1.comp$Locus == "") & (case1.comp$TSSUTR == "") & (!case1.comp$termUTR == "")),]
other1 <- case1.comp[((!case1.comp$Locus == "") & (!case1.comp$TSSUTR == "") & (!case1.comp$termUTR == "")),]
other2 <- case1.comp[((case1.comp$Locus == "") & (!case1.comp$TSSUTR == "") & (!case1.comp$termUTR == "")),]
other3 <- case1.comp[((!case1.comp$Locus == "") & (case1.comp$TSSUTR == "") & (!case1.comp$termUTR == "")),]
other4 <- case1.comp[((!case1.comp$Locus == "") & (!case1.comp$TSSUTR == "") & (case1.comp$termUTR == "")),]
other.case1.comp <- rbind(other1, other2, other3, other4)

#########################################################
##########  Case 2 - Transcriptome Distribution #########
#########################################################
############  Top

genic.case2.top <- case2.top[((!case2.top$Locus == "") & (case2.top$TSSUTR == "") & (case2.top$termUTR == "")),]
TssUTR.case2.top <- case2.top[((!case2.top$TSSUTR == "") & (case2.top$Locus == "") & (case2.top$termUTR == "")),]
intergenic.case2.top <- case2.top[((case2.top$Locus == "") & (case2.top$TSSUTR == "") & (case2.top$termUTR == "")),]
termUTR.case2.top <- case2.top[((case2.top$Locus == "") & (case2.top$TSSUTR == "") & (!case2.top$termUTR == "")),]
other1 <- case2.top[((!case2.top$Locus == "") & (!case2.top$TSSUTR == "") & (!case2.top$termUTR == "")),]
other2 <- case2.top[((case2.top$Locus == "") & (!case2.top$TSSUTR == "") & (!case2.top$termUTR == "")),]
other3 <- case2.top[((!case2.top$Locus == "") & (case2.top$TSSUTR == "") & (!case2.top$termUTR == "")),]
other4 <- case2.top[((!case2.top$Locus == "") & (!case2.top$TSSUTR == "") & (case2.top$termUTR == "")),]
other.case2.top <- rbind(other1, other2, other3, other4)

############  Comp

genic.case2.comp <- case2.comp[((!case2.comp$Locus == "") & (case2.comp$TSSUTR == "") & (case2.comp$termUTR == "")),]
TssUTR.case2.comp <- case2.comp[((!case2.comp$TSSUTR == "") & (case2.comp$Locus == "") & (case2.comp$termUTR == "")),]
intergenic.case2.comp <- case2.comp[((case2.comp$Locus == "") & (case2.comp$TSSUTR == "") & (case2.comp$termUTR == "")),]
termUTR.case2.comp <- case2.comp[((case2.comp$Locus == "") & (case2.comp$TSSUTR == "") & (!case2.comp$termUTR == "")),]
other1 <- case2.comp[((!case2.comp$Locus == "") & (!case2.comp$TSSUTR == "") & (!case2.comp$termUTR == "")),]
other2 <- case2.comp[((case2.comp$Locus == "") & (!case2.comp$TSSUTR == "") & (!case2.comp$termUTR == "")),]
other3 <- case2.comp[((!case2.comp$Locus == "") & (case2.comp$TSSUTR == "") & (!case2.comp$termUTR == "")),]
other4 <- case2.comp[((!case2.comp$Locus == "") & (!case2.comp$TSSUTR == "") & (case2.comp$termUTR == "")),]
other.case2.comp <- rbind(other1, other2, other3, other4)

#########################################################
##########  Case 3 - Transcriptome Distribution #########
#########################################################
############  Top

genic.case3.top <- case3.top[((!case3.top$Locus == "") & (case3.top$TSSUTR == "") & (case3.top$termUTR == "")),]
TssUTR.case3.top <- case3.top[((!case3.top$TSSUTR == "") & (case3.top$Locus == "") & (case3.top$termUTR == "")),]
intergenic.case3.top <- case3.top[((case3.top$Locus == "") & (case3.top$TSSUTR == "") & (case3.top$termUTR == "")),]
termUTR.case3.top <- case3.top[((case3.top$Locus == "") & (case3.top$TSSUTR == "") & (!case3.top$termUTR == "")),]
other1 <- case3.top[((!case3.top$Locus == "") & (!case3.top$TSSUTR == "") & (!case3.top$termUTR == "")),]
other2 <- case3.top[((case3.top$Locus == "") & (!case3.top$TSSUTR == "") & (!case3.top$termUTR == "")),]
other3 <- case3.top[((!case3.top$Locus == "") & (case3.top$TSSUTR == "") & (!case3.top$termUTR == "")),]
other4 <- case3.top[((!case3.top$Locus == "") & (!case3.top$TSSUTR == "") & (case3.top$termUTR == "")),]
other.case3.top <- rbind(other1, other2, other3, other4)

############  Comp

genic.case3.comp <- case3.comp[((!case3.comp$Locus == "") & (case3.comp$TSSUTR == "") & (case3.comp$termUTR == "")),]
TssUTR.case3.comp <- case3.comp[((!case3.comp$TSSUTR == "") & (case3.comp$Locus == "") & (case3.comp$termUTR == "")),]
intergenic.case3.comp <- case3.comp[((case3.comp$Locus == "") & (case3.comp$TSSUTR == "") & (case3.comp$termUTR == "")),]
termUTR.case3.comp <- case3.comp[((case3.comp$Locus == "") & (case3.comp$TSSUTR == "") & (!case3.comp$termUTR == "")),]
other1 <- case3.comp[((!case3.comp$Locus == "") & (!case3.comp$TSSUTR == "") & (!case3.comp$termUTR == "")),]
other2 <- case3.comp[((case3.comp$Locus == "") & (!case3.comp$TSSUTR == "") & (!case3.comp$termUTR == "")),]
other3 <- case3.comp[((!case3.comp$Locus == "") & (case3.comp$TSSUTR == "") & (!case3.comp$termUTR == "")),]
other4 <- case3.comp[((!case3.comp$Locus == "") & (!case3.comp$TSSUTR == "") & (case3.comp$termUTR == "")),]
other.case3.comp <- rbind(other1, other2, other3, other4)

#########################################################
##########  Case 4 - Transcriptome Distribution #########
#########################################################
############  Top

genic.case4.top <- case4.top[((!case4.top$Locus == "") & (case4.top$TSSUTR == "") & (case4.top$termUTR == "")),]
TssUTR.case4.top <- case4.top[((!case4.top$TSSUTR == "") & (case4.top$Locus == "") & (case4.top$termUTR == "")),]
intergenic.case4.top <- case4.top[((case4.top$Locus == "") & (case4.top$TSSUTR == "") & (case4.top$termUTR == "")),]
termUTR.case4.top <- case4.top[((case4.top$Locus == "") & (case4.top$TSSUTR == "") & (!case4.top$termUTR == "")),]
other1 <- case4.top[((!case4.top$Locus == "") & (!case4.top$TSSUTR == "") & (!case4.top$termUTR == "")),]
other2 <- case4.top[((case4.top$Locus == "") & (!case4.top$TSSUTR == "") & (!case4.top$termUTR == "")),]
other3 <- case4.top[((!case4.top$Locus == "") & (case4.top$TSSUTR == "") & (!case4.top$termUTR == "")),]
other4 <- case4.top[((!case4.top$Locus == "") & (!case4.top$TSSUTR == "") & (case4.top$termUTR == "")),]
other.case4.top <- rbind(other1, other2, other3, other4)

############  Comp

genic.case4.comp <- case4.comp[((!case4.comp$Locus == "") & (case4.comp$TSSUTR == "") & (case4.comp$termUTR == "")),]
TssUTR.case4.comp <- case4.comp[((!case4.comp$TSSUTR == "") & (case4.comp$Locus == "") & (case4.comp$termUTR == "")),]
intergenic.case4.comp <- case4.comp[((case4.comp$Locus == "") & (case4.comp$TSSUTR == "") & (case4.comp$termUTR == "")),]
termUTR.case4.comp <- case4.comp[((case4.comp$Locus == "") & (case4.comp$TSSUTR == "") & (!case4.comp$termUTR == "")),]
other1 <- case4.comp[((!case4.comp$Locus == "") & (!case4.comp$TSSUTR == "") & (!case4.comp$termUTR == "")),]
other2 <- case4.comp[((case4.comp$Locus == "") & (!case4.comp$TSSUTR == "") & (!case4.comp$termUTR == "")),]
other3 <- case4.comp[((!case4.comp$Locus == "") & (case4.comp$TSSUTR == "") & (!case4.comp$termUTR == "")),]
other4 <- case4.comp[((!case4.comp$Locus == "") & (!case4.comp$TSSUTR == "") & (case4.comp$termUTR == "")),]
other.case4.comp <- rbind(other1, other2, other3, other4)

#########################################################
##########  Case 5 - Transcriptome Distribution #########
#########################################################
############  Top

genic.case5.top <- case5.top[((!case5.top$Locus == "") & (case5.top$TSSUTR == "") & (case5.top$termUTR == "")),]
TssUTR.case5.top <- case5.top[((!case5.top$TSSUTR == "") & (case5.top$Locus == "") & (case5.top$termUTR == "")),]
intergenic.case5.top <- case5.top[((case5.top$Locus == "") & (case5.top$TSSUTR == "") & (case5.top$termUTR == "")),]
termUTR.case5.top <- case5.top[((case5.top$Locus == "") & (case5.top$TSSUTR == "") & (!case5.top$termUTR == "")),]
other1 <- case5.top[((!case5.top$Locus == "") & (!case5.top$TSSUTR == "") & (!case5.top$termUTR == "")),]
other2 <- case5.top[((case5.top$Locus == "") & (!case5.top$TSSUTR == "") & (!case5.top$termUTR == "")),]
other3 <- case5.top[((!case5.top$Locus == "") & (case5.top$TSSUTR == "") & (!case5.top$termUTR == "")),]
other4 <- case5.top[((!case5.top$Locus == "") & (!case5.top$TSSUTR == "") & (case5.top$termUTR == "")),]
other.case5.top <- rbind(other1, other2, other3, other4)

############  Comp

genic.case5.comp <- case5.comp[((!case5.comp$Locus == "") & (case5.comp$TSSUTR == "") & (case5.comp$termUTR == "")),]
TssUTR.case5.comp <- case5.comp[((!case5.comp$TSSUTR == "") & (case5.comp$Locus == "") & (case5.comp$termUTR == "")),]
intergenic.case5.comp <- case5.comp[((case5.comp$Locus == "") & (case5.comp$TSSUTR == "") & (case5.comp$termUTR == "")),]
termUTR.case5.comp <- case5.comp[((case5.comp$Locus == "") & (case5.comp$TSSUTR == "") & (!case5.comp$termUTR == "")),]
other1 <- case5.comp[((!case5.comp$Locus == "") & (!case5.comp$TSSUTR == "") & (!case5.comp$termUTR == "")),]
other2 <- case5.comp[((case5.comp$Locus == "") & (!case5.comp$TSSUTR == "") & (!case5.comp$termUTR == "")),]
other3 <- case5.comp[((!case5.comp$Locus == "") & (case5.comp$TSSUTR == "") & (!case5.comp$termUTR == "")),]
other4 <- case5.comp[((!case5.comp$Locus == "") & (!case5.comp$TSSUTR == "") & (case5.comp$termUTR == "")),]
other.case5.comp <- rbind(other1, other2, other3, other4)

EnrichedOutput <- as.data.frame(matrix(nrow = 10,ncol = 13))
colnames(EnrichedOutput) <- c("cases", "totalnumber", "genic.number","intergenic.number", "UTR.number","termUTR.number","other","genicpercentage","intergenic.percentage","UTR.percentage","termUTR.percentage",  "other.percentage", "total.percentage")
EnrichedOutput$cases <- c("case1.top", "case2.top", "case3.top", "case4.top", "case5.top", "case1.comp", "case2.comp", "case3.comp", "case4.comp", "case5.comp")

EnrichedOutput$totalnumber <- c(nrow(case1.top), nrow(case2.top), nrow(case3.top), nrow(case4.top), nrow(case5.top), nrow(case1.comp), nrow(case2.comp), nrow(case3.comp), nrow(case4.comp), nrow(case5.comp))

EnrichedOutput$genic.number <- c(nrow(genic.case1.top), nrow(genic.case2.top), nrow(genic.case3.top), nrow(genic.case4.top), nrow(genic.case5.top), nrow(genic.case1.comp), nrow(genic.case2.comp), nrow(genic.case3.comp), nrow(genic.case4.comp), nrow(genic.case5.comp))

EnrichedOutput$intergenic.number <- c(nrow(intergenic.case1.top), nrow(intergenic.case2.top), nrow(intergenic.case3.top), nrow(intergenic.case4.top), nrow(intergenic.case5.top), nrow(intergenic.case1.comp), nrow(intergenic.case2.comp), nrow(intergenic.case3.comp), nrow(intergenic.case4.comp), nrow(intergenic.case5.comp))

EnrichedOutput$UTR.number <- c(nrow(TssUTR.case1.top), nrow(TssUTR.case2.top), nrow(TssUTR.case3.top), nrow(TssUTR.case4.top), nrow(TssUTR.case5.top), nrow(TssUTR.case1.comp), nrow(TssUTR.case2.comp), nrow(TssUTR.case3.comp), nrow(TssUTR.case4.comp), nrow(TssUTR.case5.comp))

EnrichedOutput$termUTR.number <- c(nrow(termUTR.case1.top), nrow(termUTR.case2.top), nrow(termUTR.case3.top), nrow(termUTR.case4.top), nrow(termUTR.case5.top), nrow(termUTR.case1.comp), nrow(termUTR.case2.comp), nrow(termUTR.case3.comp), nrow(termUTR.case4.comp), nrow(termUTR.case5.comp))

EnrichedOutput$other <- c(nrow(other.case1.top), nrow(other.case2.top), nrow(other.case3.top), nrow(other.case4.top), nrow(other.case5.top), nrow(other.case1.comp), nrow(other.case2.comp), nrow(other.case3.comp), nrow(other.case4.comp), nrow(other.case5.comp))

for(i in 1:nrow(EnrichedOutput)){
  EnrichedOutput$genicpercentage[i] <- 100*EnrichedOutput$genic.number[i]/EnrichedOutput$totalnumber[i]
  EnrichedOutput$intergenic.percentage[i] <- 100*EnrichedOutput$intergenic.number[i]/EnrichedOutput$totalnumber[i]
  EnrichedOutput$UTR.percentage[i] <- 100*EnrichedOutput$UTR.number[i]/EnrichedOutput$totalnumber[i]
  EnrichedOutput$termUTR.percentage[i] <- 100*EnrichedOutput$termUTR.number[i]/EnrichedOutput$totalnumber[i]
  EnrichedOutput$other.percentage[i] <- 100*EnrichedOutput$other[i]/EnrichedOutput$totalnumber[i]
  EnrichedOutput$total.percentage[i] <- EnrichedOutput$genicpercentage[i] + EnrichedOutput$intergenic.percentage[i] + EnrichedOutput$UTR.percentage[i] + EnrichedOutput$other.percentage[i] + EnrichedOutput$termUTR.percentage[i]
}

write.csv(EnrichedOutput, file = "Termseq-EnrichedOutput-Ksg-TranscriptomicDistributionAnalysis.csv", sep = "", col.names = F, row.names = F)
