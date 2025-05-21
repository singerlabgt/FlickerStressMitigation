#--------------------------------------------------------Sara Bitarafan_WOOD Lab 2023 --------------------------------------------------------------#
#Setting Up directory
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))  # setting up working directory
#load Library 
library(DESeq2)
library(tidyverse)
library(biomaRt)
source("WoodLabFunctions.R")
source("summarySE.R")
if (!require("pacman")) install.packages("pacman")
if (!require("BiocManager")) install.packages("BiocManager")



pacman::p_load(pacman,rio,tidyverse,readxl,matrixStats,ggpubr,heatmap3,Rtsne,limma,
               gplots,GSVA,ropls,RColorBrewer,ggsci, EnhancedVolcano,DESeq2,readxl,xlsx,ggvenn,grid,forcats,viridis,org.Mm.eg.db,ggbeeswarm,WGCNA,impute,preprocessCore) # install/load these packages
#---Set heatmap3 color bar parameters
breakBarColors=c(-200,seq(-1.5, 1.5, 0.01),200) 
barColors = colorpanel(length(breakBarColors)-1, "blue", "white", "red") # create spectrum of blue to red for heatmaps




#Neurons Male
#Modifiable Gen

NeuronMaleModGenes <- read_excel("Stress DEGs modifiable by AV Flicker_Across CNS cells_WorkingSheet.xlsx",sheet= "ModifiableNeuronDEGs_Males")

NeuronMaleModGenes<-NeuronMaleModGenes[-c(1) , c(2,8)]

colnames(NeuronMaleModGenes)<-c("10Hz","40Hz")

#reading DEGs

Male40zDEGs<-read_csv("Neurons/061924/Deseq2/Neuronal DEGs_Male 40Hz Stress vs Male Dark Stress DESeq2 ThreshholdPassing FC2 v2.csv")


#40Hz

indForty=unique(match(NeuronMaleModGenes$'40Hz', Male40zDEGs$rownames.sigOE_ordered.))
indForty

Forty<-Male40zDEGs[indForty ,]
FortyExport<-data.frame(Forty,NeuronMaleModGenes$'40Hz')

export(data.frame(FortyExport),file="Data/ModifiableGenes_Male_Forty_Neurons_Export.csv")

#10Hz

Male10zDEGs<-read_csv("Neurons/061924/Deseq2/Neuronal DEGs_Male 10Hz Stress vs Male Dark Stress DESeq2 ThreshholdPassing FC2 v2.csv")
Male10zDEGs<-drop_na(Male10zDEGs)

indTen=unique(match(NeuronMaleModGenes$'10Hz', Male10zDEGs$rownames.sigOE_ordered.))
indTen

ModTen<-NeuronMaleModGenes$'10Hz'
ModTen<-na.omit(ModTen)

Ten<-Male10zDEGs[indTen ,]
MISSING <- is.na(Ten$rownames.sigOE_ordered.) |
  is.na(Ten$baseMean) |
  is.na(Ten$log2FoldChange) |
  is.na(Ten$lfcSE) |
  is.na(Ten$stat) |
  is.na(Ten$pvalue)

TenV2<-subset(Ten, 
              subset = !MISSING)

TenExport<-data.frame(TenV2,ModTen)

export(data.frame(TenExport),file="Data/ModifiableGenes_Male_Neurons_Ten_Export.csv")



#Neurons Female

#Modifiable Gen

NeuronFemaleModGenes <- read_excel("Stress DEGs modifiable by AV Flicker_Across CNS cells_WorkingSheet.xlsx",sheet= "ModifiableNeuronDEGs_Females")

NeuronFemaleModGenes<-NeuronFemaleModGenes[-c(1) , c(2,9)]

colnames(NeuronFemaleModGenes)<-c("10Hz","40Hz")

#reading DEGs

Female40zDEGs<-read_csv("Neurons/061924/Deseq2/Neuronal DEGs_Female 40Hz Stress vs Female Dark Stress DESeq2 ThreshholdPassing FC2 v2.csv")


#40Hz

indForty=unique(match(NeuronFemaleModGenes$'40Hz', Female40zDEGs$rownames.sigOE_ordered.))
indForty

Forty<-Female40zDEGs[indForty ,]
FortyExport<-data.frame(Forty,NeuronFemaleModGenes$'40Hz')

export(data.frame(FortyExport),file="Data/ModifiableGenes_female_Forty_Neurons_Export.csv")



#10Hz

Female10zDEGs<-read_csv("Neurons/061924/Deseq2/Neuronal DEGs_Female 10Hz Stress vs Female Dark Stress DESeq2 ThreshholdPassing FC2 v2.csv")
Female10zDEGs<-drop_na(Female10zDEGs)

indTen=unique(match(NeuronFemaleModGenes$'10Hz', Female10zDEGs$rownames.sigOE_ordered.))
indTen

ModTen<-NeuronFemaleModGenes$'10Hz'
ModTen<-na.omit(ModTen)

Ten<-Female10zDEGs[indTen ,]
MISSING <- is.na(Ten$rownames.sigOE_ordered.) |
  is.na(Ten$baseMean) |
  is.na(Ten$log2FoldChange) |
  is.na(Ten$lfcSE) |
  is.na(Ten$stat) |
  is.na(Ten$pvalue)

TenV2<-subset(Ten, 
              subset = !MISSING)

TenExport<-data.frame(TenV2,ModTen)

export(data.frame(TenExport),file="Data/ModifiableGenes_Female_Neurons_Ten_Export.csv")












#Astrocyte Male
#Modifiable Gen

AstrocyteMaleModGenes <- read_excel("Stress DEGs modifiable by AV Flicker_Across CNS cells_WorkingSheet.xlsx",sheet= "ModifiableAstrocyteDEGs_Males")

AstrocyteMaleModGenes<-AstrocyteMaleModGenes[-c(1) , c(2,9)]

colnames(AstrocyteMaleModGenes)<-c("10Hz","40Hz")

#reading DEGs

Male40zDEGs<-read_csv("Astrocyte/062524/Deseq2/Astrocyte DEGs_Stress Male 40Hz vs Stress Male Dark_DESeq2 ThreshholdPassing FC2 v2.csv")


#40Hz

indForty=unique(match(AstrocyteMaleModGenes$'40Hz', Male40zDEGs$rownames.sigOE_ordered.))
indForty

Forty<-Male40zDEGs[indForty ,]
FortyExport<-data.frame(Forty,AstrocyteMaleModGenes$'40Hz')

export(data.frame(FortyExport),file="Data/ModifiableGenes_Male_Forty_Astrocyte_Export.csv")

#10Hz

Male10zDEGs<-read_csv("Astrocyte/062524/Deseq2/Astrocyte DEGs_Stress Male 10Hz vs Stress Male Dark_DESeq2 ThreshholdPassing FC2 v2.csv")
Male10zDEGs<-drop_na(Male10zDEGs)

indTen=unique(match(AstrocyteMaleModGenes$'10Hz', Male10zDEGs$rownames.sigOE_ordered.))
indTen

ModTen<-AstrocyteMaleModGenes$'10Hz'
ModTen<-na.omit(ModTen)

Ten<-Male10zDEGs[indTen ,]
MISSING <- is.na(Ten$rownames.sigOE_ordered.) |
  is.na(Ten$baseMean) |
  is.na(Ten$log2FoldChange) |
  is.na(Ten$lfcSE) |
  is.na(Ten$stat) |
  is.na(Ten$pvalue)

TenV2<-subset(Ten, 
              subset = !MISSING)

TenExport<-data.frame(TenV2,ModTen)

export(data.frame(TenExport),file="Data/ModifiableGenes_Male_Astrocyte_Ten_Export.csv")



#Astrocyte Female

#Modifiable Gen

AstrocyteFemaleModGenes <- read_excel("Stress DEGs modifiable by AV Flicker_Across CNS cells_WorkingSheet.xlsx",sheet= "ModifiableAstrocyteDEGs_Females")

AstrocyteFemaleModGenes<-AstrocyteFemaleModGenes[-c(1) , c(2,8)]

colnames(AstrocyteFemaleModGenes)<-c("10Hz","40Hz")

#reading DEGs

Female40zDEGs<-read_csv("Astrocyte/062524/Deseq2/Astrocyte DEGs_Stress Female 40Hz vs Stress Female Dark_DESeq2 ThreshholdPassing FC2 v2.csv")


#40Hz

indForty=unique(match(AstrocyteFemaleModGenes$'40Hz', Female40zDEGs$rownames.sigOE_ordered.))
indForty

Forty<-Female40zDEGs[indForty ,]
FortyExport<-data.frame(Forty,AstrocyteFemaleModGenes$'40Hz')

export(data.frame(FortyExport),file="Data/ModifiableGenes_female_Forty_Astrocyte_Export.csv")


# 
# #10Hz
# 
# Female10zDEGs<-read_csv("Astrocyte/062524/Deseq2/Astrocyte DEGs_Stress Female 10Hz vs Stress Female Dark_DESeq2 ThreshholdPassing FC2 v2.csv")
# Female10zDEGs<-drop_na(Female10zDEGs)
# 
# indTen=unique(match(NeuronFemaleModGenes$'10Hz', Female10zDEGs$rownames.sigOE_ordered.))
# indTen
# 
# ModTen<-NeuronFemaleModGenes$'10Hz'
# ModTen<-na.omit(ModTen)
# 
# Ten<-Female10zDEGs[indTen ,]
# MISSING <- is.na(Ten$rownames.sigOE_ordered.) |
#   is.na(Ten$baseMean) |
#   is.na(Ten$log2FoldChange) |
#   is.na(Ten$lfcSE) |
#   is.na(Ten$stat) |
#   is.na(Ten$pvalue)
# 
# TenV2<-subset(Ten, 
#               subset = !MISSING)
# 
# TenExport<-data.frame(TenV2,ModTen)
# 
# export(data.frame(TenExport),file="Data/ModifiableGenes_Female_Astrocyte_Ten_Export.csv")
# 
# 
# 






#Microglia Male
#Modifiable Gen

MicrogliaMaleModGenes <- read_excel("Stress DEGs modifiable by AV Flicker_Across CNS cells_WorkingSheet.xlsx",sheet= "ModifiableMicrogliaDEGs_Males")

MicrogliaMaleModGenes<-MicrogliaMaleModGenes[-c(1) , c(2,9)]

colnames(MicrogliaMaleModGenes)<-c("10Hz","40Hz")

#reading DEGs

Male40zDEGs<-read_csv("Microglia/062524/Deseq2/Microglia DEGs_Male 40Hz Stress vs Male Dark Stress DESeq2 ThreshholdPassing FC2 V2.csv")


#40Hz

indForty=unique(match(MicrogliaMaleModGenes$'40Hz', Male40zDEGs$rownames.sigOE_ordered.))
indForty

ModForty<-MicrogliaMaleModGenes$'40Hz'
ModForty<-na.omit(ModForty)

Forty<-Male40zDEGs[indForty ,]

MISSING <- is.na(Forty$rownames.sigOE_ordered.) |
  is.na(Forty$baseMean) |
  is.na(Forty$log2FoldChange) |
  is.na(Forty$lfcSE) |
  is.na(Forty$stat) |
  is.na(Forty$pvalue)
FortyV2<-subset(Forty, 
              subset = !MISSING)


FortyExport<-data.frame(FortyV2,ModForty)

export(data.frame(FortyExport),file="Data/ModifiableGenes_Male_Forty_Microglia_Export.csv")

#10Hz

Male10zDEGs<-read_csv("Microglia/062524/Deseq2/Microglia DEGs_Male 10Hz Stress vs Male Dark Stress DESeq2 ThreshholdPassing FC2 V2.csv")
Male10zDEGs<-drop_na(Male10zDEGs)

indTen=unique(match(MicrogliaMaleModGenes$'10Hz', Male10zDEGs$rownames.sigOE_ordered.))
indTen

ModTen<-MicrogliaMaleModGenes$'10Hz'
ModTen<-na.omit(ModTen)

Ten<-Male10zDEGs[indTen ,]
MISSING <- is.na(Ten$rownames.sigOE_ordered.) |
  is.na(Ten$baseMean) |
  is.na(Ten$log2FoldChange) |
  is.na(Ten$lfcSE) |
  is.na(Ten$stat) |
  is.na(Ten$pvalue)

TenV2<-subset(Ten, 
              subset = !MISSING)

TenExport<-data.frame(TenV2)

export(data.frame(TenExport),file="Data/ModifiableGenes_Male_Microglia_Ten_Export.csv")



#Microglia Female

#Modifiable Gene

MicrogliaFemaleModGenes <- read_excel("Stress DEGs modifiable by AV Flicker_Across CNS cells_WorkingSheet.xlsx",sheet= "ModifiableMicrogliaDEGs_Females")

MicrogliaFemaleModGenes<-MicrogliaFemaleModGenes[-c(1:2) , c(1,8)]

colnames(MicrogliaFemaleModGenes)<-c("10Hz","40Hz")

#reading DEGs

Female40zDEGs<-read_csv("Microglia/062524/Deseq2/Microglia DEGs_Female 40Hz Stress vs Female Dark Stress DESeq2 ThreshholdPassing FC2 V2.csv")


#40Hz

indForty=unique(match(MicrogliaFemaleModGenes$'40Hz', Female40zDEGs$rownames.sigOE_ordered.))
indForty

Forty<-Female40zDEGs[indForty ,]
FortyExport<-data.frame(Forty,MicrogliaFemaleModGenes$'40Hz')

export(data.frame(FortyExport),file="Data/ModifiableGenes_female_Forty_Microglia_Export.csv")



#10Hz

Female10zDEGs<-read_csv("Microglia/062524/Deseq2/Microglia DEGs_Female 10Hz Stress vs Female Dark Stress DESeq2 ThreshholdPassing FC2 V2.csv")
Female10zDEGs<-drop_na(Female10zDEGs)

indTen=unique(match(MicrogliaFemaleModGenes$'10Hz', Female10zDEGs$rownames.sigOE_ordered.))
indTen

ModTen<-MicrogliaFemaleModGenes$'10Hz'
ModTen<-na.omit(ModTen)

Ten<-Female10zDEGs[indTen ,]
MISSING <- is.na(Ten$rownames.sigOE_ordered.) |
  is.na(Ten$baseMean) |
  is.na(Ten$log2FoldChange) |
  is.na(Ten$lfcSE) |
  is.na(Ten$stat) |
  is.na(Ten$pvalue)

TenV2<-subset(Ten, 
              subset = !MISSING)

TenExport<-data.frame(TenV2,ModTen)

export(data.frame(TenExport),file="Data/ModifiableGenes_Female_Microglia_Ten_Export.csv")






# Behavioral data normalization

#Male
Behvarial <- read_excel("Suplementary Table_Behavioral Data.xlsx",sheet= "Male")

Behvarial=data.frame(Behvarial)
Group=Behvarial$Group
Score=Behvarial$Score


datagroup=data.frame(Score,Group)


dataSum=summarySE(data.frame(datagroup), measurevar="Score", groupvars = c("Group"))

write.xlsx(data.frame(dataSum), file=paste("./Male Behavioral Stat Summary.xlsx"))

Zscore=(apply(t(Behvarial[ , c(2)]),1,scale)); #z-score the data

NormtoCtrl=(Score-(-0.1473946))/0.4556077 

NormtoNoStim=(Score-(-3.9795767))/0.6357543

Dataout=data.frame(datagroup,Zscore,NormtoCtrl,NormtoNoStim)

write.xlsx(data.frame(Dataout), file=paste("./Male Behavioral normalizing different ways.xlsx"))

# Female
Behvarial <- read_excel("Suplementary Table_Behavioral Data.xlsx",sheet= "female")

Behvarial=data.frame(Behvarial)
Group=Behvarial$Group
Score=Behvarial$Score


datagroup=data.frame(Score,Group)


dataSum=summarySE(data.frame(datagroup), measurevar="Score", groupvars = c("Group"))

write.xlsx(data.frame(dataSum), file=paste("./Female Behavioral Stat Summary.xlsx"))

Zscore=(apply(t(Behvarial[ , c(2)]),1,scale)); #z-score the data
dataSum
NormtoCtrl=(Score-(0.0000020 ))/0.4318171  

NormtoNoStim=(Score-(-1.0532250))/0.4893728 

Dataout=data.frame(datagroup,Zscore,NormtoCtrl,NormtoNoStim)

write.xlsx(data.frame(Dataout), file=paste("./Female Behavioral normalizing different ways.xlsx"))

# 
# Data <- read_excel("Suplementary Table_Behavioral Data.xlsx",sheet= "Female composite score values")
# Data=data.frame(Data)
# 
# Ctrl=Data$Control
# Ctrl=na.omit(Ctrl)
# NoStim=Data$No.stim
# NoStim=na.omit(NoStim)
# Ten=Data$X10Hz
# Ten=na.omit(Ten)
# twenty=Data$X20Hz
# twenty=na.omit(twenty)
# Forty=Data$X40Hz
# Fory=na.omit(Forty)
# 
# 
# FC=log2(mean(NoStim)-mean(Ctrl))







#----GO terms FDR adjustment----#
 
# Neuron Male No Stim UP DEGs

GoTerms <- read_excel("Data/Neuron Analysis_Biological Processes_Shared Modifiable Mark.xlsx",sheet= "Sheet2")
GoTerms=data.frame(GoTerms)

GOP=GoTerms$GO.biological.process.complete
pVal=GoTerms$Client.Text.Box.Input..raw.P.value.
FE=GoTerms$Client.Text.Box.Input..fold.Enrichment.
Count=GoTerms$Count

DataSet=data.frame(GOP,pVal,FE,Count)
DataSet$FDRAjusted_SB=p.adjust(DataSet$pVal, method="fdr",n=nrow(DataSet))
pvalue.cutoff <- 5*10e-2
lfc.cutoff <- 1
threshold <- DataSet$FDRAjusted_SB < pvalue.cutoff & abs(DataSet$FE) > lfc.cutoff
length(which(threshold))

write.xlsx(data.frame(DataSet), file=paste("./Data/Raw GO terms Feb 25 update/FDR adjusted/Neuron Analysis_Biological Processes_Shared Modifiable Mark_Count 0.xlsx"))



#----GO terms FDR adjustment----#

# Neuron Male No Stim Down DEGs

GoTerms <- read_excel("Data/Raw GO terms/male neurons No Stim_downregulated_No Correction.xlsx",sheet= "Sheet1")
GoTerms=data.frame(GoTerms)

GOP=GoTerms$GO.biological.process.complete
pVal=GoTerms$Client.Text.Box.Input..raw.P.value.
FE=GoTerms$Client.Text.Box.Input..fold.Enrichment.
Count=GoTerms$Count

DataSet=data.frame(GOP,pVal,FE,Count)
DataSet$FDRAjusted_SB=p.adjust(DataSet$pVal, method="fdr",n=nrow(DataSet))
pvalue.cutoff <- 5*10e-2
lfc.cutoff <- 1
threshold <- DataSet$FDRAjusted_SB < pvalue.cutoff & abs(DataSet$FE) > lfc.cutoff
length(which(threshold))

write.xlsx(data.frame(DataSet), file=paste("./Data/Raw GO terms/FDR adjusted/GoTerms_male neurons No Stim_downregulated_Count 0.xlsx"))



#----GO terms FDR adjustment----#

# Neuron Female No Stim Up DEGs

GoTerms <- read_excel("Data/Raw GO terms/female neurons No Stim_upregulated_No Correction.xlsx",sheet= "Sheet1")
GoTerms=data.frame(GoTerms)

GOP=GoTerms$GO.biological.process.complete
pVal=GoTerms$Client.Text.Box.Input..raw.P.value.
FE=GoTerms$Client.Text.Box.Input..fold.Enrichment.
Count=GoTerms$Count

DataSet=data.frame(GOP,pVal,FE,Count)
DataSet$FDRAjusted_SB=p.adjust(DataSet$pVal, method="fdr",n=nrow(DataSet))
pvalue.cutoff <- 5*10e-2
lfc.cutoff <- 1
threshold <- DataSet$FDRAjusted_SB < pvalue.cutoff & abs(DataSet$FE) > lfc.cutoff
length(which(threshold))

write.xlsx(data.frame(DataSet), file=paste("./Data/Raw GO terms/FDR adjusted/GoTerms_female neurons No Stim_upregulated_Count 0.xlsx"))




#----GO terms FDR adjustment----#

# Neuron Female No Stim vs 10Hz  UP DEGs

GoTerms <- read_excel("Data/Raw GO terms/Biological Processes_10Hz_F_Up_No correction.xlsx",sheet= "Sheet1")
GoTerms=data.frame(GoTerms)

GOP=GoTerms$GO.biological.process.complete
pVal=GoTerms$Client.Text.Box.Input..raw.P.value.
FE=GoTerms$Client.Text.Box.Input..fold.Enrichment.
Count=GoTerms$Count

DataSet=data.frame(GOP,pVal,FE,Count)
DataSet$FDRAjusted_SB=p.adjust(DataSet$pVal, method="fdr",n=nrow(DataSet))

write.xlsx(data.frame(DataSet), file=paste("./Data/Raw GO terms/FDR adjusted/GoTerms_10Hz_F_Up_Count 0.xlsx"))


# Neuron male No Stim vs 10Hz  UP DEGs

GoTerms <- read_excel("Data/Raw GO terms/Biological Processes_10Hz_M_UP_No correction.xlsx",sheet= "Sheet2")
GoTerms=data.frame(GoTerms)

GOP=GoTerms$GO.biological.process.complete
pVal=GoTerms$Client.Text.Box.Input..raw.P.value.
FE=GoTerms$Client.Text.Box.Input..fold.Enrichment.
Count=GoTerms$Count

DataSet=data.frame(GOP,pVal,FE,Count)
DataSet$FDRAjusted_SB=p.adjust(DataSet$pVal, method="fdr",n=nrow(DataSet))

write.xlsx(data.frame(DataSet), file=paste("./Data/Raw GO terms/FDR adjusted/GoTerms_10Hz_M_UP_Count 0.xlsx"))


# Neuron female No Stim vs 40Hz  UP DEGs

GoTerms <- read_excel("Data/Raw GO terms/Biological Processes_40Hz_F_UP_No correction.xlsx",sheet= "Sheet2")
GoTerms=data.frame(GoTerms)

GOP=GoTerms$GO.biological.process.complete
pVal=GoTerms$Client.Text.Box.Input..raw.P.value.
FE=GoTerms$Client.Text.Box.Input..fold.Enrichment.
Count=GoTerms$Count

DataSet=data.frame(GOP,pVal,FE,Count)
DataSet$FDRAjusted_SB=p.adjust(DataSet$pVal, method="fdr",n=nrow(DataSet))

write.xlsx(data.frame(DataSet), file=paste("./Data/Raw GO terms/FDR adjusted/GoTerms_40Hz_F_UP_Count 0.xlsx"))



# Neuron Male No Stim vs 40Hz Down DEGs

GoTerms <- read_excel("Data/Raw GO terms/Biological Processes_40Hz_M_Down_No correction.xlsx",sheet= "Sheet2")
GoTerms=data.frame(GoTerms)

GOP=GoTerms$GO.biological.process.complete
pVal=GoTerms$Client.Text.Box.Input..raw.P.value.
FE=GoTerms$Client.Text.Box.Input..fold.Enrichment.
Count=GoTerms$Count

DataSet=data.frame(GOP,pVal,FE,Count)
DataSet$FDRAjusted_SB=p.adjust(DataSet$pVal, method="fdr",n=nrow(DataSet))

write.xlsx(data.frame(DataSet), file=paste("./Data/Raw GO terms/FDR adjusted/GoTerms_40Hz_M_Down_Count 0.xlsx"))



# Neuron Male No Stim vs 40Hz Up DEGs

GoTerms <- read_excel("Data/Raw GO terms/Biological Processes_40Hz_M_UP_No correction.xlsx",sheet= "Sheet2")
GoTerms=data.frame(GoTerms)

GOP=GoTerms$GO.biological.process.complete
pVal=GoTerms$Client.Text.Box.Input..raw.P.value.
FE=GoTerms$Client.Text.Box.Input..fold.Enrichment.
Count=GoTerms$Count

DataSet=data.frame(GOP,pVal,FE,Count)
DataSet$FDRAjusted_SB=p.adjust(DataSet$pVal, method="fdr",n=nrow(DataSet))

write.xlsx(data.frame(DataSet), file=paste("./Data/Raw GO terms/FDR adjusted/GoTerms_40Hz_M_UP_Count 0.xlsx"))




# Microglia Male Control vs Stress Up DEGs

GoTerms <- read_excel("Data/Raw GO terms/male microglia No Stim_upregulated_No Correction.xlsx",sheet= "Sheet2")
GoTerms=data.frame(GoTerms)

GOP=GoTerms$GO.biological.process.complete
pVal=GoTerms$Client.Text.Box.Input..raw.P.value.
FE=GoTerms$Client.Text.Box.Input..fold.Enrichment.
Count=GoTerms$Count

DataSet=data.frame(GOP,pVal,FE,Count)
DataSet$FDRAjusted_SB=p.adjust(DataSet$pVal, method="fdr",n=nrow(DataSet))

write.xlsx(data.frame(DataSet), file=paste("./Data/Raw GO terms/FDR adjusted/GoTerms_male microglia No Stim_upregulated_Count 0.xlsx"))


# Microglia Male Control vs Stress Down DEGs

GoTerms <- read_excel("Data/Raw GO terms/male microglia No Stim_downregulated_No Correction.xlsx",sheet= "Sheet2")
GoTerms=data.frame(GoTerms)

GOP=GoTerms$GO.biological.process.complete
pVal=GoTerms$Client.Text.Box.Input..raw.P.value.
FE=GoTerms$Client.Text.Box.Input..fold.Enrichment.
Count=GoTerms$Count

DataSet=data.frame(GOP,pVal,FE,Count)
DataSet$FDRAjusted_SB=p.adjust(DataSet$pVal, method="fdr",n=nrow(DataSet))

write.xlsx(data.frame(DataSet), file=paste("./Data/Raw GO terms/FDR adjusted/GoTerms_male microglia No Stim_downregulated_Count 0.xlsx"))




# Microglia Female Control vs Stress Up DEGs

GoTerms <- read_excel("Data/Raw GO terms/female microglia No Stim_upregulated_No Correction.xlsx",sheet= "Sheet2")
GoTerms=data.frame(GoTerms)

GOP=GoTerms$GO.biological.process.complete
pVal=GoTerms$Client.Text.Box.Input..raw.P.value.
FE=GoTerms$Client.Text.Box.Input..fold.Enrichment.
Count=GoTerms$Count

DataSet=data.frame(GOP,pVal,FE,Count)
DataSet$FDRAjusted_SB=p.adjust(DataSet$pVal, method="fdr",n=nrow(DataSet))

write.xlsx(data.frame(DataSet), file=paste("./Data/Raw GO terms/FDR adjusted/GoTerms_female microglia No Stim_upregulated_Count 0.xlsx"))






# Microglia male No Stim vs 40Hz Up DEGs

GoTerms <- read_excel("Data/Raw GO terms/Microglia_Biological Processes_40Hz_M_UP_No correction.xlsx",sheet= "Sheet2")
GoTerms=data.frame(GoTerms)

GOP=GoTerms$GO.biological.process.complete
pVal=GoTerms$Client.Text.Box.Input..raw.P.value.
FE=GoTerms$Client.Text.Box.Input..fold.Enrichment.
Count=GoTerms$Count

DataSet=data.frame(GOP,pVal,FE,Count)
DataSet$FDRAjusted_SB=p.adjust(DataSet$pVal, method="fdr",n=nrow(DataSet))

write.xlsx(data.frame(DataSet), file=paste("./Data/Raw GO terms/FDR adjusted/GoTerms_40Hz_M_UP_Count 0.xlsx"))



# Microglia male No Stim vs 40Hz Down DEGs

GoTerms <- read_excel("Data/Raw GO terms/Microglia_Biological Processes_40Hz_M_down_No correction.xlsx",sheet= "Sheet2")
GoTerms=data.frame(GoTerms)

GOP=GoTerms$GO.biological.process.complete
pVal=GoTerms$Client.Text.Box.Input..raw.P.value.
FE=GoTerms$Client.Text.Box.Input..fold.Enrichment.
Count=GoTerms$Count

DataSet=data.frame(GOP,pVal,FE,Count)
DataSet$FDRAjusted_SB=p.adjust(DataSet$pVal, method="fdr",n=nrow(DataSet))

write.xlsx(data.frame(DataSet), file=paste("./Data/Raw GO terms/FDR adjusted/GoTerms_40Hz_M_down_Count 0.xlsx"))




# Microglia Female No Stim vs 40Hz up DEGs

GoTerms <- read_excel("Data/Raw GO terms/Microglia_Biological Processes_40Hz_F_UP_No correction.xlsx",sheet= "Sheet2")
GoTerms=data.frame(GoTerms)

GOP=GoTerms$GO.biological.process.complete
pVal=GoTerms$Client.Text.Box.Input..raw.P.value.
FE=GoTerms$Client.Text.Box.Input..fold.Enrichment.
Count=GoTerms$Count

DataSet=data.frame(GOP,pVal,FE,Count)
DataSet$FDRAjusted_SB=p.adjust(DataSet$pVal, method="fdr",n=nrow(DataSet))

write.xlsx(data.frame(DataSet), file=paste("./Data/Raw GO terms/FDR adjusted/GoTerms_40Hz_F_UP_Count 0.xlsx"))


# Microglia Male No Stim vs 10Hz up DEGs

GoTerms <- read_excel("Data/Raw GO terms/Microglia_Biological Processes_10Hz_M_UP_No correction.xlsx",sheet= "Sheet2")
GoTerms=data.frame(GoTerms)

GOP=GoTerms$GO.biological.process.complete
pVal=GoTerms$Client.Text.Box.Input..raw.P.value.
FE=GoTerms$Client.Text.Box.Input..fold.Enrichment.
Count=GoTerms$Count

DataSet=data.frame(GOP,pVal,FE,Count)
DataSet$FDRAjusted_SB=p.adjust(DataSet$pVal, method="fdr",n=nrow(DataSet))

write.xlsx(data.frame(DataSet), file=paste("./Data/Raw GO terms/FDR adjusted/GoTerms_10Hz_M_UP_Count 0.xlsx"))



# Microglia Male No Stim vs 10Hz Down DEGs

GoTerms <- read_excel("Data/Raw GO terms/Microglia_Biological Processes_10Hz_M_down_No correction.xlsx",sheet= "Sheet2")
GoTerms=data.frame(GoTerms)

GOP=GoTerms$GO.biological.process.complete
pVal=GoTerms$Client.Text.Box.Input..raw.P.value.
FE=GoTerms$Client.Text.Box.Input..fold.Enrichment.
Count=GoTerms$Count

DataSet=data.frame(GOP,pVal,FE,Count)
DataSet$FDRAjusted_SB=p.adjust(DataSet$pVal, method="fdr",n=nrow(DataSet))

write.xlsx(data.frame(DataSet), file=paste("./Data/Raw GO terms/FDR adjusted/GoTerms_10Hz_M_down_Count 0.xlsx"))




# Microglia female No Stim vs 10Hz Up DEGs

GoTerms <- read_excel("Data/Raw GO terms/Microglia_Biological Processes_10Hz_F_UP_No correction.xlsx",sheet= "Sheet2")
GoTerms=data.frame(GoTerms)

GOP=GoTerms$GO.biological.process.complete
pVal=GoTerms$Client.Text.Box.Input..raw.P.value.
FE=GoTerms$Client.Text.Box.Input..fold.Enrichment.
Count=GoTerms$Count

DataSet=data.frame(GOP,pVal,FE,Count)
DataSet$FDRAjusted_SB=p.adjust(DataSet$pVal, method="fdr",n=nrow(DataSet))

write.xlsx(data.frame(DataSet), file=paste("./Data/Raw GO terms/FDR adjusted/GoTerms_10Hz_F_UP_Count 0.xlsx"))





# Microglia Female No Stim vs 10Hz Down DEGs

GoTerms <- read_excel("Data/Raw GO terms/Microglia_Biological Processes_10Hz_F_down_No correction.xlsx",sheet= "Sheet2")
GoTerms=data.frame(GoTerms)

GOP=GoTerms$GO.biological.process.complete
pVal=GoTerms$Client.Text.Box.Input..raw.P.value.
FE=GoTerms$Client.Text.Box.Input..fold.Enrichment.
Count=GoTerms$Count

DataSet=data.frame(GOP,pVal,FE,Count)
DataSet$FDRAjusted_SB=p.adjust(DataSet$pVal, method="fdr",n=nrow(DataSet))

write.xlsx(data.frame(DataSet), file=paste("./Data/Raw GO terms/FDR adjusted/GoTerms_10Hz_F_down_Count 0.xlsx"))



# Astrocyte male Stress vs Control Up DEGs

GoTerms <- read_excel("Data/Raw GO terms/male astrocyte No Stim_upregulated_No Correction.xlsx",sheet= "Sheet2")
GoTerms=data.frame(GoTerms)

GOP=GoTerms$GO.biological.process.complete
pVal=GoTerms$Client.Text.Box.Input..raw.P.value.
FE=GoTerms$Client.Text.Box.Input..fold.Enrichment.
Count=GoTerms$Count

DataSet=data.frame(GOP,pVal,FE,Count)
DataSet$FDRAjusted_SB=p.adjust(DataSet$pVal, method="fdr",n=nrow(DataSet))

write.xlsx(data.frame(DataSet), file=paste("./Data/Raw GO terms/FDR adjusted/GoTerms_male astrocyte No Stim_upregulated_Count 0.xlsx"))




# Astrocyte male Stress vs Control Down DEGs

GoTerms <- read_excel("Data/Raw GO terms/male astrocyte No Stim_downregulated_No Correction.xlsx",sheet= "Sheet2")
GoTerms=data.frame(GoTerms)

GOP=GoTerms$GO.biological.process.complete
pVal=GoTerms$Client.Text.Box.Input..raw.P.value.
FE=GoTerms$Client.Text.Box.Input..fold.Enrichment.
Count=GoTerms$Count

DataSet=data.frame(GOP,pVal,FE,Count)
DataSet$FDRAjusted_SB=p.adjust(DataSet$pVal, method="fdr",n=nrow(DataSet))

write.xlsx(data.frame(DataSet), file=paste("./Data/Raw GO terms/FDR adjusted/GoTerms_male astrocyte No Stim_downregulated_Count 0.xlsx"))


# Astrocyte Female Stress vs Control Up DEGs

GoTerms <- read_excel("Data/Raw GO terms/female astrocyte No Stim_upregulated_No Correction.xlsx",sheet= "Sheet2")
GoTerms=data.frame(GoTerms)

GOP=GoTerms$GO.biological.process.complete
pVal=GoTerms$Client.Text.Box.Input..raw.P.value.
FE=GoTerms$Client.Text.Box.Input..fold.Enrichment.
Count=GoTerms$Count

DataSet=data.frame(GOP,pVal,FE,Count)
DataSet$FDRAjusted_SB=p.adjust(DataSet$pVal, method="fdr",n=nrow(DataSet))

write.xlsx(data.frame(DataSet), file=paste("./Data/Raw GO terms/FDR adjusted/GoTerms_female astrocyte No Stim_upregulated_Count 0.xlsx"))


# Astrocyte Female No Stim vs 10Hz Down DEGs

GoTerms <- read_excel("Data/Raw GO terms/Astrocytes_Biological Processes_10Hz_F_down_No correction.xlsx",sheet= "Sheet2")
GoTerms=data.frame(GoTerms)

GOP=GoTerms$GO.biological.process.complete
pVal=GoTerms$Client.Text.Box.Input..raw.P.value.
FE=GoTerms$Client.Text.Box.Input..fold.Enrichment.
Count=GoTerms$Count

DataSet=data.frame(GOP,pVal,FE,Count)
DataSet$FDRAjusted_SB=p.adjust(DataSet$pVal, method="fdr",n=nrow(DataSet))

write.xlsx(data.frame(DataSet), file=paste("./Data/Raw GO terms/FDR adjusted/Astrocytes_Biological Processes_10Hz_F_down_Count 0.xlsx"))


# Astrocyte Female No Stim vs 10Hz Up DEGs

GoTerms <- read_excel("Data/Raw GO terms/Astrocytes_Biological Processes_10Hz_F_UP_No correction.xlsx",sheet= "Sheet2")
GoTerms=data.frame(GoTerms)

GOP=GoTerms$GO.biological.process.complete
pVal=GoTerms$Client.Text.Box.Input..raw.P.value.
FE=GoTerms$Client.Text.Box.Input..fold.Enrichment.
Count=GoTerms$Count

DataSet=data.frame(GOP,pVal,FE,Count)
DataSet$FDRAjusted_SB=p.adjust(DataSet$pVal, method="fdr",n=nrow(DataSet))

write.xlsx(data.frame(DataSet), file=paste("./Data/Raw GO terms/FDR adjusted/Astrocytes_Biological Processes_10Hz_F_UP_Count 0.xlsx"))

# Astrocyte Male No Stim vs 10Hz Down DEGs

GoTerms <- read_excel("Data/Raw GO terms/Astrocytes_Biological Processes_10Hz_M_down_No correction.xlsx",sheet= "Sheet2")
GoTerms=data.frame(GoTerms)

GOP=GoTerms$GO.biological.process.complete
pVal=GoTerms$Client.Text.Box.Input..raw.P.value.
FE=GoTerms$Client.Text.Box.Input..fold.Enrichment.
Count=GoTerms$Count

DataSet=data.frame(GOP,pVal,FE,Count)
DataSet$FDRAjusted_SB=p.adjust(DataSet$pVal, method="fdr",n=nrow(DataSet))

write.xlsx(data.frame(DataSet), file=paste("./Data/Raw GO terms/FDR adjusted/Astrocytes_Biological Processes_10Hz_M_down_Count 0.xlsx"))


# Astrocyte Male No Stim vs 10Hz Up DEGs

GoTerms <- read_excel("Data/Raw GO terms/Astrocytes_Biological Processes_10Hz_M_UP_No correction.xlsx",sheet= "Sheet2")
GoTerms=data.frame(GoTerms)

GOP=GoTerms$GO.biological.process.complete
pVal=GoTerms$Client.Text.Box.Input..raw.P.value.
FE=GoTerms$Client.Text.Box.Input..fold.Enrichment.
Count=GoTerms$Count

DataSet=data.frame(GOP,pVal,FE,Count)
DataSet$FDRAjusted_SB=p.adjust(DataSet$pVal, method="fdr",n=nrow(DataSet))

write.xlsx(data.frame(DataSet), file=paste("./Data/Raw GO terms/FDR adjusted/Astrocytes_Biological Processes_10Hz_M_UP_Count 0.xlsx"))


# Astrocyte Female No Stim vs 40Hz Down DEGs

GoTerms <- read_excel("Data/Raw GO terms/Astrocytes_Biological Processes_40Hz_F_DOWN_No correction.xlsx",sheet= "Sheet2")
GoTerms=data.frame(GoTerms)

GOP=GoTerms$GO.biological.process.complete
pVal=GoTerms$Client.Text.Box.Input..raw.P.value.
FE=GoTerms$Client.Text.Box.Input..fold.Enrichment.
Count=GoTerms$Count

DataSet=data.frame(GOP,pVal,FE,Count)
DataSet$FDRAjusted_SB=p.adjust(DataSet$pVal, method="fdr",n=nrow(DataSet))

write.xlsx(data.frame(DataSet), file=paste("./Data/Raw GO terms/FDR adjusted/Astrocytes_Biological Processes_40Hz_F_DOWN_Count 0.xlsx"))


# Astrocyte Female No Stim vs 40Hz Up DEGs

GoTerms <- read_excel("Data/Raw GO terms/Astrocytes_Biological Processes_40Hz_F_UP_No correction.xlsx",sheet= "Sheet2")
GoTerms=data.frame(GoTerms)

GOP=GoTerms$GO.biological.process.complete
pVal=GoTerms$Client.Text.Box.Input..raw.P.value.
FE=GoTerms$Client.Text.Box.Input..fold.Enrichment.
Count=GoTerms$Count

DataSet=data.frame(GOP,pVal,FE,Count)
DataSet$FDRAjusted_SB=p.adjust(DataSet$pVal, method="fdr",n=nrow(DataSet))

write.xlsx(data.frame(DataSet), file=paste("./Data/Raw GO terms/FDR adjusted/Astrocytes_Biological Processes_40Hz_F_UP_Count 0.xlsx"))



# Astrocyte male No Stim vs 40Hz Down DEGs

GoTerms <- read_excel("Data/Raw GO terms/Astrocytes_Biological Processes_40Hz_M_down_No correction.xlsx",sheet= "Sheet2")
GoTerms=data.frame(GoTerms)

GOP=GoTerms$GO.biological.process.complete
pVal=GoTerms$Client.Text.Box.Input..raw.P.value.
FE=GoTerms$Client.Text.Box.Input..fold.Enrichment.
Count=GoTerms$Count

DataSet=data.frame(GOP,pVal,FE,Count)
DataSet$FDRAjusted_SB=p.adjust(DataSet$pVal, method="fdr",n=nrow(DataSet))

write.xlsx(data.frame(DataSet), file=paste("./Data/Raw GO terms/FDR adjusted/Astrocytes_Biological Processes_40Hz_M_down_Count 0.xlsx"))



# Astrocyte male No Stim vs 40Hz Up DEGs

GoTerms <- read_excel("Data/Raw GO terms/Astrocytes_Biological Processes_40Hz_M_UP_No correction.xlsx",sheet= "Sheet2")
GoTerms=data.frame(GoTerms)

GOP=GoTerms$GO.biological.process.complete
pVal=GoTerms$Client.Text.Box.Input..raw.P.value.
FE=GoTerms$Client.Text.Box.Input..fold.Enrichment.
Count=GoTerms$Count

DataSet=data.frame(GOP,pVal,FE,Count)
DataSet$FDRAjusted_SB=p.adjust(DataSet$pVal, method="fdr",n=nrow(DataSet))

write.xlsx(data.frame(DataSet), file=paste("./Data/Raw GO terms/FDR adjusted/Astrocytes_Biological Processes_40Hz_M_UP_Count 0.xlsx"))




#------------------ Differential Expression Analysis--------#

#Selecting Condition of interest for volcano Plots

rawcounts <- read.csv("Data/NeuronalCounts_22173R.csv", header = TRUE, row.names = 1)


#filter criteria as recommended on WGCNA FAQ 
# https://horvath.genetics.ucla.edu/html/CoexpressionNetwork/Rpackages/WGCNA/faq.html
#(genes where 50% of samples have count >10)
n=ncol(rawcounts)*0.75 # 75% of the sample
indKeepNew=which(rowSums(rawcounts>50)> n)
datafliter=rawcounts[indKeepNew,]
genes=datafliter
x = rowSums(rawcounts)
numberofFilteredGenes=count(x<50) # run this to see how many genes filterefd




# Reading experimental Group
Expt_data_in <- read_excel("Data/numericMeta Neuron.xlsx" , sheet="numericMeta") %>%  # imports the metadata
  mutate(Sex = case_when(
    Sex== "Males" ~ "M",
    Sex== "Females" ~ "F"
  )) %>%
  mutate(sampleID = str_c(Sample,Sex,sep=" ")) # creates new variable sampleID that combines sample and region

colnames(genes)=Expt_data_in$sampleID


Expt_data <- tibble(.rows = ncol(genes)) # create new tibble data frame with one row for each sample
Expt_data$sampleID <- colnames(genes) # create a column with sample names in same order as genes data
Expt_data <- left_join(Expt_data,Expt_data_in,by="sampleID") # add the experiment metadata in same order as genes data

Expt_data <- Expt_data %>%
  mutate(Condition = paste0(Sex," ",Frequency,"Hz ",Stressed," "))


#Sorting
ideal_sorted = arrange(Expt_data,Sex,Frequency,Stressed)
sort_ind = match(ideal_sorted$sampleID,Expt_data$sampleID)





#-------------DEGs------------#
#----- Female Stress vs Control _ Dark -----#

sortby <- which(Expt_data$Condition =="F DarkHz Stress " | Expt_data$Condition  =="F DarkHz Control ")#This  can be changed to any group or condition
Expt_data_sub=Expt_data[sortby,]
genes_sub=genes[,sortby]

dds <- DESeqDataSetFromMatrix(countData = genes_sub,
                              colData = Expt_data_sub,
                              design = ~ Stressed)

dds$Stressed <- relevel( dds$Stressed, ref = "Control")


dds= DESeq(dds)
summary(results(dds))
DesqData <- results(dds)



#Simplifying Gene Names
GeneName=rownames(genes_sub)
GeneName2=sub("\\|.*", "", GeneName)

rownames(DesqData)=GeneName2


#Choosing Threshhold for volcanos
pvalue.cutoff <- 5*10e-3
lfc.cutoff <- 1
threshold <- DesqData$pvalue < pvalue.cutoff & abs(DesqData$log2FoldChange) > lfc.cutoff
length(which(threshold))
DesqData$threshold<-threshold
sigOE <- data.frame(subset(DesqData, threshold==TRUE))
sigOE_ordered <- sigOE[order(sigOE$pvalue), ]

outMale=data.frame(rownames(sigOE_ordered),sigOE_ordered)
export(outMale,"Neurons/061924/Deseq2/Neuronal DEGs_Female Dark Stress vs Female Dark Control DESeq2 ThreshholdPassing FC2 v2.csv")

#####Export Data#####

DesqDataOut=as.matrix(data.frame(DesqData))
DseqExport <-data.frame(rownames(DesqDataOut),DesqDataOut)
export(DseqExport,"Neurons/061924/Deseq2/NeuronalDEGs_Female Dark Stress vs Female Dark Control DESeq2 All v2.csv")

# Setting up output folder for the figure
output_folder=output_folder = paste0("Neurons/061924/Deseq2/Volcano/")
dir.create(output_folder,recursive = TRUE, showWarnings=FALSE)


keyvals <- ifelse(
  DesqData$log2FoldChange < -1 & DesqData$pvalue < 5*10e-3, '#0070c0',
  ifelse(DesqData$log2FoldChange > 1 & DesqData$pvalue < 5*10e-3, '#FA0000',
         'gray'))
keyvals[is.na(keyvals)] <- 'gray'
names(keyvals)[keyvals == '#FA0000'] <- 'Up'
names(keyvals)[keyvals == 'gray'] <- 'NS'
names(keyvals)[keyvals == '#0070c0'] <- 'Down'




#---- Generating Volcano Plot---#
# Male 10Hz vs Dark selectLab=c('Syde1','Abi1','Dnajc6','Nlgn3','Ngdn','Mctp1','Nae1','Adra2c','Atp2b3','Bcas1','Eef1d','Hapln1','Map2','Necab2','Ngdn','Rgs9','S1pr5','Sema4b','Slc30a1')
#Male 40Hz vs Dark   

# LabNames <- read_excel("Neurons/061924/Deseq2/Copy of SynGO annotations_synaptic markers_Neuronal.xlsx" , sheet="40Fvs10M_all v2")
# 
# LabNames$gene=str_to_title(LabNames$Genes)
# LabNames$gene=as.character(LabNames$gene)
# GeneListVolcano=unique(LabNames$gene)


LabNames <- read_excel("Data/Genes to label for updating volcano plots_12Feb2025.xlsx" , sheet="Female Neurons")
GeneList<-LabNames$`Stress No Stim`
GeneListVolcano<-unique(str_to_title(GeneList))



#pdf(paste0(output_folder,"Neuronal DEGs Female Dark Stress vs Female Dark Control FC2.pdf"), width=7,height=8,pointsize = 14,useDingbats = FALSE)
png(paste0(output_folder,"Neuronal DEGs Female Dark Stress vs Female Dark Control FC2.png"), width=7,height=8,pointsize = 14,units="in",res=600)
EnhancedVolcano(DesqData,
                lab = rownames(DesqData),
                selectLab=GeneListVolcano,
               # selectLab =rownames(DesqData)[which(names(keyvals) %in% c('Up', 'Down'))],
                labCol = 'black',
                labFace = 'italic',
                #boxedLabels = TRUE,
                colCustom = keyvals,
                colAlpha = 1,
                legendPosition = 'right',
                legendLabSize = 15,
                legendIconSize = 5.0,
                labSize=3,
                x = 'log2FoldChange',
                y = 'pvalue',
                pCutoff = 5*10e-3,
                FCcutoff=1,
                # col=c('grey', 'grey', 'grey', 'magenta'),
                gridlines.major = FALSE,
                gridlines.minor = FALSE,
                # xlim = c(min(DesqData$log2FoldChange, na.rm = TRUE) - 1.5, max(DesqData$log2FoldChange, na.rm = TRUE) +
                #1.5),
                xlim=c(-6,6) , ylim=c(0,6),
                # ylim = c(0, max(-log10(DesqData$pvalue), na.rm = TRUE) + 5),
                border = 'full')
dev.off()



#----- Female 10Hz vs Dark _ Stress -----#

sortby <- which(Expt_data$Condition =="F 10HzHz Stress " | Expt_data$Condition  =="F DarkHz Stress ")#This  can be changed to any group or condition
Expt_data_sub=Expt_data[sortby,]
genes_sub=genes[,sortby]

 
 dds <- DESeqDataSetFromMatrix(countData = genes_sub,
                             colData = Expt_data_sub,
                              design = ~ Frequency)
 dds$Frequency <- relevel( dds$Frequency, ref = "Dark")



dds= DESeq(dds)
summary(results(dds))
DesqData <- results(dds)



#Simplifying Gene Names
GeneName=rownames(genes_sub)
GeneName2=sub("\\|.*", "", GeneName)

rownames(DesqData)=GeneName2


#Choosing Threshhold for volcanos
pvalue.cutoff <- 5*10e-3
lfc.cutoff <- 1
threshold <- DesqData$pvalue < pvalue.cutoff & abs(DesqData$log2FoldChange) > lfc.cutoff
length(which(threshold))
DesqData$threshold<-threshold
sigOE <- data.frame(subset(DesqData, threshold==TRUE))
sigOE_ordered <- sigOE[order(sigOE$pvalue), ]

outMale=data.frame(rownames(sigOE_ordered),sigOE_ordered)
export(outMale,"Neurons/061924/Deseq2/Neuronal DEGs_Female 10Hz Stress vs Female Dark Stress DESeq2 ThreshholdPassing FC2 v2.csv")

#####Export Data#####

DesqDataOut=as.matrix(data.frame(DesqData))
DseqExport <-data.frame(rownames(DesqDataOut),DesqDataOut)
export(DseqExport,"Neurons/061924/Deseq2/NeuronalDEGs_Female 10Hz Stress vs Female Dark Stress DESeq2 All v2.csv")

# Setting up output folder for the figure
output_folder=output_folder = paste0("Neurons/061924/Deseq2/Volcano/")
dir.create(output_folder,recursive = TRUE, showWarnings=FALSE)


keyvals <- ifelse(
  DesqData$log2FoldChange < -1 & DesqData$pvalue < 5*10e-3, '#0070c0',
  ifelse(DesqData$log2FoldChange > 1 & DesqData$pvalue < 5*10e-3, '#FA0000',
         'gray'))
keyvals[is.na(keyvals)] <- 'gray'
names(keyvals)[keyvals == '#FA0000'] <- 'Up'
names(keyvals)[keyvals == 'gray'] <- 'NS'
names(keyvals)[keyvals == '#0070c0'] <- 'Down'




#---- Generating Volcano Plot---#


LabNames <- read_excel("Data/Genes to label for updating volcano plots_12Feb2025.xlsx" , sheet="Female Neurons")
GeneList<-LabNames$'10Hz'
GeneListVolcano<-unique(str_to_title(GeneList))



pdf(paste0(output_folder,"Neuronal DEGs Female 10Hz Stress vs Female Dark Stress FC2.pdf"), width=7,height=8,pointsize = 14,useDingbats = FALSE)
#png(paste0(output_folder,"Neuronal DEGs Female 10Hz Stress vs Female Dark Stress FC2.png"), width=7,height=8,pointsize = 14,units="in",res=600)
EnhancedVolcano(DesqData,
                lab = rownames(DesqData),
                selectLab=GeneListVolcano,
                # selectLab =rownames(DesqData)[which(names(keyvals) %in% c('Up', 'Down'))],
                labCol = 'black',
                labFace = 'italic',
                #boxedLabels = TRUE,
                colCustom = keyvals,
                colAlpha = 1,
                legendPosition = 'right',
                legendLabSize = 15,
                legendIconSize = 5.0,
                labSize=3,
                x = 'log2FoldChange',
                y = 'pvalue',
                pCutoff = 5*10e-3,
                FCcutoff=1,
                # col=c('grey', 'grey', 'grey', 'magenta'),
                gridlines.major = FALSE,
                gridlines.minor = FALSE,
                # xlim = c(min(DesqData$log2FoldChange, na.rm = TRUE) - 1.5, max(DesqData$log2FoldChange, na.rm = TRUE) +
                #1.5),
                xlim=c(-6,6) , ylim=c(0,6),
                # ylim = c(0, max(-log10(DesqData$pvalue), na.rm = TRUE) + 5),
                border = 'full')
dev.off()



#----- Female 40Hz vs Dark _ Stress -----#

sortby <- which(Expt_data$Condition =="F 40HzHz Stress " | Expt_data$Condition  =="F DarkHz Stress ")#This  can be changed to any group or condition
Expt_data_sub=Expt_data[sortby,]
genes_sub=genes[,sortby]


dds <- DESeqDataSetFromMatrix(countData = genes_sub,
                              colData = Expt_data_sub,
                              design = ~ Frequency)
dds$Frequency <- relevel( dds$Frequency, ref = "Dark")



dds= DESeq(dds)
summary(results(dds))
DesqData <- results(dds)



#Simplifying Gene Names
GeneName=rownames(genes_sub)
GeneName2=sub("\\|.*", "", GeneName)

rownames(DesqData)=GeneName2


#Choosing Threshhold for volcanos
pvalue.cutoff <- 5*10e-3
lfc.cutoff <- 1
threshold <- DesqData$pvalue < pvalue.cutoff & abs(DesqData$log2FoldChange) > lfc.cutoff
length(which(threshold))
DesqData$threshold<-threshold
sigOE <- data.frame(subset(DesqData, threshold==TRUE))
sigOE_ordered <- sigOE[order(sigOE$pvalue), ]

outMale=data.frame(rownames(sigOE_ordered),sigOE_ordered)
export(outMale,"Neurons/061924/Deseq2/Neuronal DEGs_Female 40Hz Stress vs Female Dark Stress DESeq2 ThreshholdPassing FC2 v2.csv")

#####Export Data#####

DesqDataOut=as.matrix(data.frame(DesqData))
DseqExport <-data.frame(rownames(DesqDataOut),DesqDataOut)
export(DseqExport,"Neurons/061924/Deseq2/NeuronalDEGs_Female 40Hz Stress vs Female Dark Stress DESeq2 All v2.csv")

# Setting up output folder for the figure
output_folder=output_folder = paste0("Neurons/061924/Deseq2/Volcano/")
dir.create(output_folder,recursive = TRUE, showWarnings=FALSE)


keyvals <- ifelse(
  DesqData$log2FoldChange < -1 & DesqData$pvalue < 5*10e-3, '#0070c0',
  ifelse(DesqData$log2FoldChange > 1 & DesqData$pvalue < 5*10e-3, '#FA0000',
         'gray'))
keyvals[is.na(keyvals)] <- 'gray'
names(keyvals)[keyvals == '#FA0000'] <- 'Up'
names(keyvals)[keyvals == 'gray'] <- 'NS'
names(keyvals)[keyvals == '#0070c0'] <- 'Down'




#---- Generating Volcano Plot---#


LabNames <- read_excel("Data/Genes to label for updating volcano plots_12Feb2025.xlsx" , sheet="Female Neurons")
GeneList<-LabNames$'40Hz'
GeneListVolcano<-unique(str_to_title(GeneList))



#pdf(paste0(output_folder,"Neuronal DEGs Female 40Hz Stress vs Female Dark Stress FC2.pdf"), width=7,height=8,pointsize = 14,useDingbats = FALSE)
png(paste0(output_folder,"Neuronal DEGs Female 40Hz Stress vs Female Dark Stress FC2.png"), width=7,height=8,pointsize = 14,units="in",res=600)
EnhancedVolcano(DesqData,
                lab = rownames(DesqData),
                selectLab=GeneListVolcano,
                # selectLab =rownames(DesqData)[which(names(keyvals) %in% c('Up', 'Down'))],
                labCol = 'black',
                labFace = 'italic',
                #boxedLabels = TRUE,
                colCustom = keyvals,
                colAlpha = 1,
                legendPosition = 'right',
                legendLabSize = 15,
                legendIconSize = 5.0,
                labSize=3,
                x = 'log2FoldChange',
                y = 'pvalue',
                pCutoff = 5*10e-3,
                FCcutoff=1,
                # col=c('grey', 'grey', 'grey', 'magenta'),
                gridlines.major = FALSE,
                gridlines.minor = FALSE,
                # xlim = c(min(DesqData$log2FoldChange, na.rm = TRUE) - 1.5, max(DesqData$log2FoldChange, na.rm = TRUE) +
                #1.5),
                xlim=c(-6,6) , ylim=c(0,6),
                # ylim = c(0, max(-log10(DesqData$pvalue), na.rm = TRUE) + 5),
                border = 'full')
dev.off()



#-----------Male----------#


#----- Male Stress vs Control _ Dark -----#

sortby <- which(Expt_data$Condition =="M DarkHz Stress " | Expt_data$Condition  =="M DarkHz Control ")#This  can be changed to any group or condition
Expt_data_sub=Expt_data[sortby,]
genes_sub=genes[,sortby]

dds <- DESeqDataSetFromMatrix(countData = genes_sub,
                              colData = Expt_data_sub,
                              design = ~ Stressed)

dds$Stressed <- relevel( dds$Stressed, ref = "Control")


dds= DESeq(dds)
summary(results(dds))
DesqData <- results(dds)



#Simplifying Gene Names
GeneName=rownames(genes_sub)
GeneName2=sub("\\|.*", "", GeneName)

rownames(DesqData)=GeneName2


#Choosing Threshhold for volcanos
pvalue.cutoff <- 5*10e-3
lfc.cutoff <- 1
threshold <- DesqData$pvalue < pvalue.cutoff & abs(DesqData$log2FoldChange) > lfc.cutoff
length(which(threshold))
DesqData$threshold<-threshold
sigOE <- data.frame(subset(DesqData, threshold==TRUE))
sigOE_ordered <- sigOE[order(sigOE$pvalue), ]

outMale=data.frame(rownames(sigOE_ordered),sigOE_ordered)
export(outMale,"Neurons/061924/Deseq2/Neuronal DEGs_Male Dark Stress vs Male Dark Control DESeq2 ThreshholdPassing FC2 v2.csv")

#####Export Data#####

DesqDataOut=as.matrix(data.frame(DesqData))
DseqExport <-data.frame(rownames(DesqDataOut),DesqDataOut)
export(DseqExport,"Neurons/061924/Deseq2/Male Dark Stress vs Male Dark Control DESeq2 All v2.csv")

# Setting up output folder for the figure
output_folder=output_folder = paste0("Neurons/061924/Deseq2/Volcano/")
dir.create(output_folder,recursive = TRUE, showWarnings=FALSE)


keyvals <- ifelse(
  DesqData$log2FoldChange < -1 & DesqData$pvalue < 5*10e-3, '#0070c0',
  ifelse(DesqData$log2FoldChange > 1 & DesqData$pvalue < 5*10e-3, '#FA0000',
         'gray'))
keyvals[is.na(keyvals)] <- 'gray'
names(keyvals)[keyvals == '#FA0000'] <- 'Up'
names(keyvals)[keyvals == 'gray'] <- 'NS'
names(keyvals)[keyvals == '#0070c0'] <- 'Down'




#---- Generating Volcano Plot---#

LabNames <- read_excel("Data/Genes to label for updating volcano plots_12Feb2025.xlsx" , sheet="Male Neurons")
GeneList<-LabNames$`Stress No Stim`
GeneListVolcano<-unique(str_to_title(GeneList))



pdf(paste0(output_folder,"Neuronal DEGs Male Dark Stress vs Male Dark Control FC2.pdf"), width=7,height=8,pointsize = 14,useDingbats = FALSE)
#png(paste0(output_folder,"Neuronal DEGs Male Dark Stress vs Male Dark Control FC2.png"), width=7,height=8,pointsize = 14,units="in",res=600)
EnhancedVolcano(DesqData,
                lab = rownames(DesqData),
                selectLab=GeneListVolcano,
                # selectLab =rownames(DesqData)[which(names(keyvals) %in% c('Up', 'Down'))],
                labCol = 'black',
                labFace = 'italic',
                #boxedLabels = TRUE,
                colCustom = keyvals,
                colAlpha = 1,
                legendPosition = 'right',
                legendLabSize = 15,
                legendIconSize = 5.0,
                labSize=3,
                x = 'log2FoldChange',
                y = 'pvalue',
                pCutoff = 5*10e-3,
                FCcutoff=1,
                # col=c('grey', 'grey', 'grey', 'magenta'),
                gridlines.major = FALSE,
                gridlines.minor = FALSE,
                # xlim = c(min(DesqData$log2FoldChange, na.rm = TRUE) - 1.5, max(DesqData$log2FoldChange, na.rm = TRUE) +
                #1.5),
                xlim=c(-6,6) , ylim=c(0,6),
                # ylim = c(0, max(-log10(DesqData$pvalue), na.rm = TRUE) + 5),
                border = 'full')
dev.off()



#----- Male 10Hz vs Dark _ Stress -----#

sortby <- which(Expt_data$Condition =="M 10HzHz Stress " | Expt_data$Condition  =="M DarkHz Stress ")#This  can be changed to any group or condition
Expt_data_sub=Expt_data[sortby,]
genes_sub=genes[,sortby]


dds <- DESeqDataSetFromMatrix(countData = genes_sub,
                              colData = Expt_data_sub,
                              design = ~ Frequency)
dds$Frequency <- relevel( dds$Frequency, ref = "Dark")



dds= DESeq(dds)
summary(results(dds))
DesqData <- results(dds)



#Simplifying Gene Names
GeneName=rownames(genes_sub)
GeneName2=sub("\\|.*", "", GeneName)

rownames(DesqData)=GeneName2


#Choosing Threshhold for volcanos
pvalue.cutoff <- 5*10e-3
lfc.cutoff <- 1
threshold <- DesqData$pvalue < pvalue.cutoff & abs(DesqData$log2FoldChange) > lfc.cutoff
length(which(threshold))
DesqData$threshold<-threshold
sigOE <- data.frame(subset(DesqData, threshold==TRUE))
sigOE_ordered <- sigOE[order(sigOE$pvalue), ]

outMale=data.frame(rownames(sigOE_ordered),sigOE_ordered)
export(outMale,"Neurons/061924/Deseq2/Neuronal DEGs_Male 10Hz Stress vs Male Dark Stress DESeq2 ThreshholdPassing FC2 v2.csv")

#####Export Data#####

DesqDataOut=as.matrix(data.frame(DesqData))
DseqExport <-data.frame(rownames(DesqDataOut),DesqDataOut)
export(DseqExport,"Neurons/061924/Deseq2/NeuronalDEGs_Male 10Hz Stress vs Male Dark Stress DESeq2 All v2.csv")

# Setting up output folder for the figure
output_folder=output_folder = paste0("Neurons/061924/Deseq2/Volcano/")
dir.create(output_folder,recursive = TRUE, showWarnings=FALSE)


keyvals <- ifelse(
  DesqData$log2FoldChange < -1 & DesqData$pvalue < 5*10e-3, '#0070c0',
  ifelse(DesqData$log2FoldChange > 1 & DesqData$pvalue < 5*10e-3, '#FA0000',
         'gray'))
keyvals[is.na(keyvals)] <- 'gray'
names(keyvals)[keyvals == '#FA0000'] <- 'Up'
names(keyvals)[keyvals == 'gray'] <- 'NS'
names(keyvals)[keyvals == '#0070c0'] <- 'Down'




#---- Generating Volcano Plot---#


LabNames <- read_excel("Data/Genes to label for updating volcano plots_12Feb2025.xlsx" , sheet="Male Neurons")
GeneList<-LabNames$'10Hz'
GeneListVolcano<-unique(str_to_title(GeneList))



pdf(paste0(output_folder,"Neuronal DEGs Male 10Hz Stress vs Male Dark Stress FC2.pdf"), width=7,height=8,pointsize = 14,useDingbats = FALSE)
#png(paste0(output_folder,"Neuronal DEGs Male 10Hz Stress vs Male Dark Stress FC2.png"), width=7,height=8,pointsize = 14,units="in",res=600)
EnhancedVolcano(DesqData,
                lab = rownames(DesqData),
                selectLab=GeneListVolcano,
                # selectLab =rownames(DesqData)[which(names(keyvals) %in% c('Up', 'Down'))],
                labCol = 'black',
                labFace = 'italic',
                #boxedLabels = TRUE,
                colCustom = keyvals,
                colAlpha = 1,
                legendPosition = 'right',
                legendLabSize = 15,
                legendIconSize = 5.0,
                labSize=3,
                x = 'log2FoldChange',
                y = 'pvalue',
                pCutoff = 5*10e-3,
                FCcutoff=1,
                # col=c('grey', 'grey', 'grey', 'magenta'),
                gridlines.major = FALSE,
                gridlines.minor = FALSE,
                # xlim = c(min(DesqData$log2FoldChange, na.rm = TRUE) - 1.5, max(DesqData$log2FoldChange, na.rm = TRUE) +
                #1.5),
                xlim=c(-6,6) , ylim=c(0,6),
                # ylim = c(0, max(-log10(DesqData$pvalue), na.rm = TRUE) + 5),
                border = 'full')
dev.off()



#----- Male 40Hz vs Dark _ Stress -----#

sortby <- which(Expt_data$Condition =="M 40HzHz Stress " | Expt_data$Condition  =="M DarkHz Stress ")#This  can be changed to any group or condition
Expt_data_sub=Expt_data[sortby,]
genes_sub=genes[,sortby]


dds <- DESeqDataSetFromMatrix(countData = genes_sub,
                              colData = Expt_data_sub,
                              design = ~ Frequency)
dds$Frequency <- relevel( dds$Frequency, ref = "Dark")



dds= DESeq(dds)
summary(results(dds))
DesqData <- results(dds)



#Simplifying Gene Names
GeneName=rownames(genes_sub)
GeneName2=sub("\\|.*", "", GeneName)

rownames(DesqData)=GeneName2


#Choosing Threshhold for volcanos
pvalue.cutoff <- 5*10e-3
lfc.cutoff <- 1
threshold <- DesqData$pvalue < pvalue.cutoff & abs(DesqData$log2FoldChange) > lfc.cutoff
length(which(threshold))
DesqData$threshold<-threshold
sigOE <- data.frame(subset(DesqData, threshold==TRUE))
sigOE_ordered <- sigOE[order(sigOE$pvalue), ]

outMale=data.frame(rownames(sigOE_ordered),sigOE_ordered)
export(outMale,"Neurons/061924/Deseq2/Neuronal DEGs_Male 40Hz Stress vs Male Dark Stress DESeq2 ThreshholdPassing FC2 v2.csv")

#####Export Data#####

DesqDataOut=as.matrix(data.frame(DesqData))
DseqExport <-data.frame(rownames(DesqDataOut),DesqDataOut)
export(DseqExport,"Neurons/061924/Deseq2/NeuronalDEGs_Male 40Hz Stress vs Male Dark Stress DESeq2 All v2.csv")

# Setting up output folder for the figure
output_folder=output_folder = paste0("Neurons/061924/Deseq2/Volcano/")
dir.create(output_folder,recursive = TRUE, showWarnings=FALSE)


keyvals <- ifelse(
  DesqData$log2FoldChange < -1 & DesqData$pvalue < 5*10e-3, '#0070c0',
  ifelse(DesqData$log2FoldChange > 1 & DesqData$pvalue < 5*10e-3, '#FA0000',
         'gray'))
keyvals[is.na(keyvals)] <- 'gray'
names(keyvals)[keyvals == '#FA0000'] <- 'Up'
names(keyvals)[keyvals == 'gray'] <- 'NS'
names(keyvals)[keyvals == '#0070c0'] <- 'Down'




#---- Generating Volcano Plot---#


LabNames <- read_excel("Data/Genes to label for updating volcano plots_12Feb2025.xlsx" , sheet="Male Neurons")
GeneList<-LabNames$'40Hz'
GeneListVolcano<-unique(str_to_title(GeneList))



pdf(paste0(output_folder,"Neuronal DEGs Male 40Hz Stress vs Male Dark Stress FC2.pdf"), width=7,height=8,pointsize = 14,useDingbats = FALSE)
#png(paste0(output_folder,"Neuronal DEGs Male 40Hz Stress vs Male Dark Stress FC2.png"), width=7,height=8,pointsize = 14,units="in",res=600)
EnhancedVolcano(DesqData,
                lab = rownames(DesqData),
                selectLab=GeneListVolcano,
                # selectLab =rownames(DesqData)[which(names(keyvals) %in% c('Up', 'Down'))],
                labCol = 'black',
                labFace = 'italic',
                #boxedLabels = TRUE,
                colCustom = keyvals,
                colAlpha = 1,
                legendPosition = 'right',
                legendLabSize = 15,
                legendIconSize = 5.0,
                labSize=3,
                x = 'log2FoldChange',
                y = 'pvalue',
                pCutoff = 5*10e-3,
                FCcutoff=1,
                # col=c('grey', 'grey', 'grey', 'magenta'),
                gridlines.major = FALSE,
                gridlines.minor = FALSE,
                # xlim = c(min(DesqData$log2FoldChange, na.rm = TRUE) - 1.5, max(DesqData$log2FoldChange, na.rm = TRUE) +
                #1.5),
                xlim=c(-6,6) , ylim=c(0,6),
                # ylim = c(0, max(-log10(DesqData$pvalue), na.rm = TRUE) + 5),
                border = 'full')
dev.off()





#-------------------------------------Venn Diaram-------------------------------#

#Female

FemaleNeuron <- read_excel("Data/Neuronal DEGs_Female Dark Stress vs Female Dark Control.xlsx",sheet="Up")
FemalesigOENeuron <- data.frame(subset( FemaleNeuron, threshold==TRUE))

FemaleMicroglia <- read_excel("Data/Microglia DEGs_Female Dark Stress vs Female Dark Control.xlsx",sheet="Up")
FemalesigOEMicroglia <- data.frame(subset( FemaleMicroglia, threshold==TRUE))

FemaleAstrocyte <- read_excel("Data/Astrocyte DEGs_Stress Female Stress vs Female Control.xlsx",sheet="Up")
FemalesigOEAstrocyte <- data.frame(subset( FemaleAstrocyte, threshold==TRUE))

#Male

MaleNeuron <- read_excel("Data/Neuronal DEGs_Male Dark Stress vs Male Dark Control.xlsx",sheet="Up")
MalesigOENeuron <- data.frame(subset(MaleNeuron, threshold==TRUE))


MaleMicroglia <- read_excel("Data/Microglia DEGs_Male Dark Stress vs Male Dark Control.xlsx",sheet="Up")
MalesigOEMicroglia <- data.frame(subset( MaleMicroglia, threshold==TRUE))

MaleAstrocyte <- read_excel("Data/Astrocyte DEGs_Stress Male Stress vs Male Control.xlsx",sheet="Up")
MalesigOEAstrocyte <- data.frame(subset( MaleAstrocyte, threshold==TRUE))




#Cutoff Threshhold
pvalue.cutoff <- 5*10e-3
lfc.cutoff <- 1



MaleNeuro=MalesigOENeuron$rownames.sigOE_ordered.
MaleAstro=MalesigOEAstrocyte$rownames.sigOE_ordered.
MaleMig=MalesigOEMicroglia$rownames.sigOE_ordered.

FemaleNeuro=FemalesigOENeuron$rownames.sigOE_ordered.
FemaleAstro=FemalesigOEAstrocyte$rownames.sigOE_ordered.
FemaleMig=FemalesigOEMicroglia$rownames.sigOE_ordered.




VennFreq1=list(MaleNeuro=MalesigOENeuron$rownames.sigOE_ordered.,MaleAstro=MalesigOEAstrocyte$rownames.sigOE_ordered.,MaleMig=MalesigOEMicroglia$rownames.sigOE_ordered.)

VennFreq2=list(FemaleNeuro=FemalesigOENeuron$rownames.sigOE_ordered.,FemaleAstro=FemalesigOEAstrocyte$rownames.sigOE_ordered.,FemaleMig=FemalesigOEMicroglia$rownames.sigOE_ordered.)

output_folder = paste0("All Combined/")



#png(paste0(output_folder,"Venn Diagram DEGs Stress vs Control Male Up Only.png"), width=7,height=8,pointsize = 14,units="in",res=600)
pdf(paste0(output_folder,"Venn Diagram DEGs Stress vs Control Male Up Only.pdf"), width=8,height=8,pointsize = 12)
ggvenn(VennFreq1,
       #fill_color = c("#E41A1C","#377EB8"),
       stroke_size = 0.5, set_name_size = 6, text_size = 5)
dev.off()



#-----------------------------------------------------------------------------------------------------------#
# 
# # Reading experimental Group
# Expt_data_in <- read_excel("Data/numericMeta Neuron.xlsx" , sheet="numericMeta") %>%  # imports the metadata
#   mutate(Sex = case_when(
#     Sex== "Males" ~ "M",
#     Sex== "Females" ~ "F"
#   )) %>%
#   mutate(sampleID = str_c(Sample,Sex,sep=" ")) # creates new variable sampleID that combines sample and region
# 
# colnames(genes)=Expt_data_in$sampleID
# 
# 
# Expt_data <- tibble(.rows = ncol(genes)) # create new tibble data frame with one row for each sample
# Expt_data$sampleID <- colnames(genes) # create a column with sample names in same order as genes data
# Expt_data <- left_join(Expt_data,Expt_data_in,by="sampleID") # add the experiment metadata in same order as genes data
# 
# Expt_data <- Expt_data %>%
#   mutate(Condition = paste0(Sex," ",Frequency,"Hz ",Stressed," "))
# 
# 
# #Sorting
# ideal_sorted = arrange(Expt_data,Sex,Frequency,Stressed)
# sort_ind = match(ideal_sorted$sampleID,Expt_data$sampleID)
# 
# 
# #Female
# 
# sortby1 <- which(Expt_data$Condition =="F 10HzHz Stress " | Expt_data$Condition  =="F 40HzHz Stress " | Expt_data$Condition  =="F LightHz Stress " )
# #This  can be changed to any group or condition
# Expt_data_sub_F=Expt_data[sortby1,]
# genes_sub_F=genes[,sortby1]
# GeneName_F=rownames(genes_sub_F)
# rownames(genes_sub_F)=sub("\\|.*", "", GeneName_F)
# genes_sub_Z_F <- genes_sub_F %>%
#   apply(1,"scale") %>%
#   t() %>%
#   as_tibble()
# colnames(genes_sub_Z_F) = colnames(genes_sub_F)
# genes_sub_Z_F$geneID = rownames(genes_sub_F)
# genes_sub_Z_F[is.na(genes_sub_Z_F)] <- 0   # set any NA values to 0
# 
# # Create a matrix of Z-scored subset data
# genes_subZ_mat_F <- as.matrix(genes_sub_Z_F[,1:{ncol(genes_sub_Z_F)-1}])
# rownames(genes_subZ_mat_F) <- genes_sub_Z_F$geneID
# 
# 
# #Male
# 
# sortby2 <- which(Expt_data$Condition =="M 10HzHz Stress " | Expt_data$Condition  =="M 40HzHz Stress " | Expt_data$Condition  =="M LightHz Stress ")
# #This  can be changed to any group or condition
# Expt_data_sub_M=Expt_data[sortby2,]
# genes_sub_M=genes[,sortby2]
# 
# GeneName_M=rownames(genes_sub_M)
# GeneName2_M=sub("\\|.*", "", GeneName_M)
# 
# 
# rownames(genes_sub_M)=GeneName2_M
# genes_sub_Z_M <- genes_sub_M %>%
#   apply(1,"scale") %>%
#   t() %>%
#   as_tibble()
# colnames(genes_sub_Z_M) = colnames(genes_sub_M)
# genes_sub_Z_M$geneID = rownames(genes_sub_M)
# genes_sub_Z_M[is.na(genes_sub_Z_M)] <- 0   # set any NA values to 0
# 
# # Create a matrix of Z-scored subset data
# 
# genes_subZ_mat_M <- as.matrix(genes_sub_Z_M[,1:{ncol(genes_sub_Z_M)-1}])
# rownames(genes_subZ_mat_M) <- genes_sub_Z_M$geneID
# 
# 
# 
# 
# #Venn Diagram of all DEGs
# #FC =2 , p<0.05
# 
# 
# #Female
# 
# Female40 <- read_excel("Neurons/Deseq2/Neuronal DEGs_Stress 40Hz vs Light_Female DESeq2 ThreshholdPassing FC2.xlsx",sheet="Up") 
# FemalesigOE40 <- data.frame(subset( Female40, threshold==TRUE))
# 
# Female10 <- read_excel("Neurons/Deseq2/Neuronal DEGs_Stress 10Hz vs Light_Female DESeq2 ThreshholdPassing FC2.xlsx",sheet="Up") 
# FemalesigOE10 <- data.frame(subset( Female10, threshold==TRUE))
# # 
# # Female1040 <- read_excel("Neurons/Deseq2/Tina Data Stress 40Hz vs 10 Female DESeq2 All FC2 111323SB.xlsx",sheet="Up") 
# # FemalesigOE1040 <- data.frame(subset( Female1040, threshold==TRUE))
# # 
# 
# #Male
# 
# 
# Male40 <- read_excel("Neurons/Deseq2/Neuronal DEGs_Stress 40Hz vs Light_Male DESeq2 ThreshholdPassing FC2.xlsx",sheet="Up") 
# MalesigOE40 <- data.frame(subset(Male40, threshold==TRUE))
# 
# 
# Male10 <- read_excel("Neurons/Deseq2/Neuronal DEGs_Stress 10Hz vs Light_Male DESeq2 ThreshholdPassing FC2.xlsx",sheet="Up") 
# MalesigOE10 <- data.frame(subset( Male10, threshold==TRUE))
# 
# #Cutoff Threshhold
# pvalue.cutoff <- 5*10e-3
# lfc.cutoff <- 1
# 
# 
# 
# Male40L=MalesigOE40$rownames.sigOE_ordered.
# Male10L=MalesigOE10$rownames.sigOE_ordered.
# Female10L=FemalesigOE10$rownames.sigOE_ordered.
# Female40L=FemalesigOE40$rownames.sigOE_ordered.
# 
# #Male down 40Hz specific 
# MaleData40=genes_subZ_mat_M[match(MalesigOE40$rownames.sigOE_ordered. , sub("\\|.*", "",rownames(genes_subZ_mat_M))) ,]
# MaleDataOthers<-genes_subZ_mat_M[which(Male40L %in% c(Male10L,Female10L,Female40L)),]
# 
# indRemove=match(rownames(MaleDataOthers), rownames(genes_subZ_mat_M))
# 
# MaleData40only<- MaleData40[-c(indRemove) ,]
# 
# export(data.frame(rownames(MaleData40only),MaleData40only),file="Neurons/Deseq2/Male Neurons 40Hz specific DEGs on venn diagram Up.csv")
# 
# #Male down 10Hz specific 
# 
# MaleData10=genes_subZ_mat_M[match(MalesigOE10$rownames.sigOE_ordered. , sub("\\|.*", "",rownames(genes_subZ_mat_M))) ,]
# MaleDataOthers<-genes_subZ_mat_M[which(Male10L %in% c(Male40L,Female10L,Female40L)),]
# 
# 
# indRemove=match(rownames(MaleDataOthers), rownames(genes_subZ_mat_M))
# 
# MaleData10only<- MaleData10[-c(indRemove) ,]
# 
# export(data.frame(rownames(MaleData10only),MaleData10only),file="Neurons/Deseq2/Male Neurons 10Hz specific DEGs on venn diagram Up.csv")
# 
# 
# 
# #Female down 10Hz specific 
# 
# FemaleData10=genes_subZ_mat_F[match(FemalesigOE10$rownames.sigOE_ordered. , sub("\\|.*", "",rownames(genes_subZ_mat_F))) ,]
# FemaleDataOthers<-genes_subZ_mat_F[which(Female10L %in% c(Male40L,Male10L,Female40L)),]
# 
# 
# indRemove=match(rownames(FemaleDataOthers), rownames(genes_subZ_mat_F))
# 
# FemaleData10only<- FemaleData10[-c(indRemove) ,]
# 
# export(data.frame(rownames(FemaleData10only),FemaleData10only),file="Neurons/Deseq2/Female Neurons 10Hz specific DEGs on venn diagram Up.csv")
# 
# 
# 
# 
# #Female down 40Hz specific 
# 
# FemaleData40=genes_subZ_mat_F[match(FemalesigOE40$rownames.sigOE_ordered. , sub("\\|.*", "",rownames(genes_subZ_mat_F))) ,]
# FemaleDataOthers<-genes_subZ_mat_F[which(Female40L %in% c(Male40L,Male10L,Female10L)),]
# 
# 
# indRemove=match(rownames(FemaleDataOthers), rownames(genes_subZ_mat_F))
# 
# FemaleData40only<- FemaleData40[-c(indRemove) ,]
# 
# export(data.frame(rownames(FemaleData40only),FemaleData40only),file="Neurons/Deseq2/Female Neurons 40Hz specific DEGs on venn diagram Up.csv")
# 
# 
# output_folder = paste0("Neurons/Deseq2/")
# VennFreq=list(Male40L=MalesigOE40$rownames.sigOE_ordered.,Male10L=MalesigOE10$rownames.sigOE_ordered.,Female10L=FemalesigOE10$rownames.sigOE_ordered., Female40L=FemalesigOE40$rownames.sigOE_ordered.)
# 
# png(paste0(output_folder,"Venn Diagram DEGs Neurons Up Only.png"), width=7,height=8,pointsize = 14,units="in",res=600)
# #pdf(paste0(output_folder,"Venn Diagram DEGs Neurons Up Only.pdf"), width=8,height=8,pointsize = 12)
# ggvenn(VennFreq, 
#        #fill_color = c("#E41A1C","#377EB8"),
#        stroke_size = 0.5, set_name_size = 6, text_size = 5)
# dev.off()
# 
# 
# 
# 
# 
# 
# 
# #Female
# 
# Female40 <- read_excel("Neurons/Deseq2/Neuronal DEGs_Stress 40Hz vs Light_Female DESeq2 ThreshholdPassing FC2.xlsx",sheet="Up") 
# FemalesigOE40 <- data.frame(subset( Female40, threshold==TRUE))
# 
# Female10 <- read_excel("Neurons/Deseq2/Neuronal DEGs_Stress 10Hz vs Light_Female DESeq2 ThreshholdPassing FC2.xlsx",sheet="Up") 
# FemalesigOE10 <- data.frame(subset( Female10, threshold==TRUE))
# # 
# # Female1040 <- read_excel("Neurons/Deseq2/Tina Data Stress 40Hz vs 10 Female DESeq2 All FC2 111323SB.xlsx",sheet="Up") 
# # FemalesigOE1040 <- data.frame(subset( Female1040, threshold==TRUE))
# # 
# 
# #Male
# 
# 
# Male40 <- read_excel("Neurons/Deseq2/Neuronal DEGs_Stress 40Hz vs Light_Male DESeq2 ThreshholdPassing FC2.xlsx",sheet="Up") 
# MalesigOE40 <- data.frame(subset(Male40, threshold==TRUE))
# 
# 
# Male10 <- read_excel("Neurons/Deseq2/Neuronal DEGs_Stress 10Hz vs Light_Male DESeq2 ThreshholdPassing FC2.xlsx",sheet="Up") 
# MalesigOE10 <- data.frame(subset( Male10, threshold==TRUE))
# 
# 
# #Cutoff Threshhold
# pvalue.cutoff <- 5*10e-3
# lfc.cutoff <- 1
# output_folder = paste0("Neurons/Deseq2/")
# VennFreq=list(Male10L=MalesigOE10$rownames.sigOE_ordered., Female40L=FemalesigOE40$rownames.sigOE_ordered.)
# 
# 
# 
# #png(paste0(output_folder,"Venn Diagram Male10 Female 40 Down.png"), width=7,height=8,pointsize = 14,units="in",res=600)
# pdf(paste0(output_folder,"Venn Diagram Male10 Female 40 Down.pdf"), width=8,height=8,pointsize = 12)
# ggvenn(VennFreq, 
#        fill_color = c("dodgerblue2","hotpink"),
#        stroke_size = 0.5, set_name_size = 6, text_size = 5)
# dev.off()
# 
# 
# 
# 
# ind=match(MalesigOE10$rownames.sigOE_ordered.,FemalesigOE40$rownames.sigOE_ordered.)
# 
# Shared=MalesigOE10[ind,]
# Shared=drop_na(Shared)
# 
# export(Shared,"Neurons/Deseq2/Shared Gene Venn Diagram Male 10Hz vs Light Female 40Hz vs light stressed.csv",row.names=T)
# 
# 


#######################################################################################################################################
#load  RAW Counts Data 

rawcounts <- read.csv("Data/NeuronalCounts_22173R.csv", header = TRUE, row.names = 1)

genes=rawcounts

#filter criteria as recommended on WGCNA FAQ 
# https://horvath.genetics.ucla.edu/html/CoexpressionNetwork/Rpackages/WGCNA/faq.html
#(genes where 50% of samples have count >10)
n=ncol(rawcounts)*0.75 # 50% of the sample
indKeepNew=which(rowSums(rawcounts>50)> n)
datafliter=rawcounts[indKeepNew,]
genes=datafliter
x = rowSums(rawcounts)
numberofFilteredGenes=count(x<50) # run this to see how many genes filterefd

# Reading experimental Group
Expt_data_in <- read_excel("Data/numericMeta Neuron.xlsx" , sheet="numericMeta") %>%  # imports the metadata
  mutate(Sex = case_when(
    Sex== "Males" ~ "M",
    Sex== "Females" ~ "F"
  )) %>%
  mutate(sampleID = str_c(Sample,Sex,sep=" ")) # creates new variable sampleID that combines sample and region

colnames(genes)=Expt_data_in$sampleID


Expt_data <- tibble(.rows = ncol(genes)) # create new tibble data frame with one row for each sample
Expt_data$sampleID <- colnames(genes) # create a column with sample names in same order as genes data
Expt_data <- left_join(Expt_data,Expt_data_in,by="sampleID") # add the experiment metadata in same order as genes data

Expt_data <- Expt_data %>%
  mutate(Condition = paste0(Sex," ",Frequency,"Hz ",Stressed," "))








dds <- DESeqDataSetFromMatrix(countData = genes,
                              colData = Expt_data,
                              design = ~  Sex + Stressed + Frequency)

#dds$Frequency <- relevel( dds$Frequency, ref = "Dark")
dds <- DESeq(dds)
results(dds)
#
#
dds <- estimateSizeFactors(dds)
normcounts <- counts(dds, normalized=TRUE)
saveRDS(normcounts,"NeuronNormCount_NewFilter_75PercentAbove50.RDS")



#Reading Norm Count for further analysis


normcounts = readRDS("NeuronNormCount_NewFilter_75PercentAbove50.RDS")
genes=normcounts
# colnames(normcounts)=colnames(genes)
# SAMPLE=colnames(genes)
# GeneID=rownames(normcounts)
# normcounts$geneID=GeneID
# dataall=data.frame(t(genes),SAMPLE)
# export(dataall, rownames =T, "normcountsTinaDataaaaa.xlsx")

#Sorting
ideal_sorted = arrange(Expt_data,Sex,Frequency,Stressed)
sort_ind = match(ideal_sorted$sampleID,Expt_data$sampleID)

genes <- genes[,sort_ind]
Expt_data <- Expt_data[sort_ind,]

#Z-scoring

genes_Z <- genes %>%
   apply(1,"scale") %>%
    t() %>%
   as_tibble()
 colnames(genes_Z) = colnames(genes)
 genes_Z$geneID = rownames(genes)
  genes_Z[is.na(genes_Z)] <- 0   # set any NA values to 0

  # Create a matrix of Z-scored subset data
  genes_Z_mat <- as.matrix(genes_Z[,1:{ncol(genes_Z)-1}])
  rownames(genes_Z_mat) <- str_to_upper(genes_Z$geneID)
  

pheno=Expt_data$Condition
pheno=as.factor(pheno)



#Stresspal <- brewer.pal(2,"YlGnBu")

StressColorBar <- case_when(
  Expt_data$Stressed == "Control" ~ "#FBE3D6",
  Expt_data$Stressed == "Stress" ~ "#E98F2B"
)

#Freqpal<-brewer.pal(4,"Dark2")
freqColorBar <- case_when(
  Expt_data$Frequency == "Dark" ~ "#861C33",
  Expt_data$Frequency == "10Hz" ~ "#2D567C",
  Expt_data$Frequency == "40Hz" ~ "#19B799"
)


#Sexpal<-brewer.pal(5,"Set1")
SexColorBar <- case_when(
  Expt_data$Sex == "M" ~ "#872EFE",
  Expt_data$Sex == "F" ~ "#D98DFF" 
)


#Clustering Genes 
hrGenes= hclust(dist((genes_Z_mat),method = "euclidean"), method="ward.D2")


#Generating Heatmap 

# Genes Clustered Heatmap
#pdf(paste0("./Neurons/061924/All Genes Clustered Heatmap Neurons 75percen Above 50_061924.pdf"), width=6,height=8,pointsize = 14, useDingbats = FALSE)
png(paste0("Neurons/061924/All Genes Clustered Heatmap Neurons 75percen Above 50_061924.png"), width=7,height=8,pointsize = 14, units="in",res=600)
heatmap3(genes_Z_mat, ColSideColors =cbind(c(SexColorBar),c(freqColorBar),c(StressColorBar)), ColSideWidth=1, ColSideLabs=NA,  
         col=barColors, breaks=breakBarColors, legendfun=function()showLegend(legend=c(NA),col=c(NA)),
         Rowv=as.dendrogram(hrGenes), 
         Colv=NA,  scale="none",
         margins = c(3,10),
         labRow = NA, labCol=NA,
         cexCol=1.3)
title("Genes Cluster Neurons", adj = 0.5, line = 2)
dev.off()



new_meanSEM <- function(x,y,colors=x,x_str="",y_str="",title=y_str,axis.text.x.size=20,axis.text.y.size=20,axis.text.x.angle=90,title.size=20,text.size=20){
  dataPlot <- data.frame(factor(x),y,factor(colors))
  colnames(dataPlot)=c("x_dot","y_dot","colors_dot")
  plotOut <- ggplot(dataPlot,aes(x=x_dot,y=y_dot)) +
    geom_boxplot(aes(color=colors))+
     geom_point(size=6,shape=21,aes(fill=colors))+
    geom_beeswarm(aes(fill = colors),size =6,shape=21, stroke = 1, cex = 2,groupOnX = TRUE)+
    #stat_summary(geom="errorbar",color="black",size=1.5,aes(width=0.3),)+
    #stat_summary(fun = "mean", geom="point", size=4,show.legend = F, color=colors) +
    theme(panel.background = element_rect(fill = 'transparent'),
          panel.grid.major=element_blank(),
          axis.line=element_line(color="black"),
          text = element_text(size=text.size),
          legend.key = element_blank(),
          axis.text.x = element_text(color = "black",
                                     size = axis.text.x.size,angle= axis.text.x.angle),
          axis.text.y = element_text( color = "black",
                                      size = axis.text.y.size),
          plot.title=element_text(hjust=0.5,face = "bold",size=title.size))+
    xlab(x_str)+
    ylab(y_str)+
    #ylim(c(-0, 10000))+
    ggtitle(paste(title))
  return(plotOut)
}

GeneName=rownames(genes)
GeneName2=sub("\\|.*", "", GeneName)

rownames(genes)=GeneName2
GeneList=c("Csf1r","Aif1","Vim","P2ry12")
titles=c("Csf1r","Aif1","Vim","P2ry12")
Condition=Expt_data$Condition
for(i in 1:length(GeneList)){
  ind <- which(str_detect(rownames(genes),GeneList[i]))
  plotOut <- new_meanSEM(x=Expt_data$Condition,y=genes[ind,],colors =Condition,
                         title = titles[i],y_str = "Expression",
                         axis.text.x.size = 12,axis.text.y.size = 12,title.size = 15,text.size = 15)
 pdf(paste0("./Neurons/061924/SanityCheck",titles[i],".pdf"),height=5,width=6,useDingbats = FALSE,useKerning = FALSE)
  #png(paste0("./Neurons/061924/SanityCheck",titles[i],".png"),height=5,width=6,units="in",res=600)
  print(plotOut)
  dev.off()
}





















#-------------------------------------------------Male Data------------------------------------------------

sortby <- which(Expt_data$Condition =="M 10HzHz Stress " | Expt_data$Condition  =="M 40HzHz Stress " | Expt_data$Condition  =="M DarkHz Stress " | Expt_data$Condition=="M DarkHz Control " )
                  #This  can be changed to any group or condition
Expt_data_sub=Expt_data[sortby,]
genes_sub=genes[,sortby]


GeneName=rownames(genes_sub)
GeneName2=sub("\\|.*", "", GeneName)



genes_sub_Z <- genes_sub %>%
  apply(1,"scale") %>%
  t() %>%
  as_tibble()
colnames(genes_sub_Z) = colnames(genes_sub)
genes_sub_Z$geneID = rownames(genes_sub)
genes_sub_Z[is.na(genes_sub_Z)] <- 0   # set any NA values to 0
rownames(genes_sub)=GeneName2
# Create a matrix of Z-scored subset data
genes_sub_Z_mat <- as.matrix(genes_sub_Z[,1:{ncol(genes_sub_Z)-1}])
rownames(genes_sub_Z_mat) <- str_to_upper(genes_sub_Z$geneID)



StressColorBar <- case_when(
  Expt_data_sub$Stressed == "Control" ~ "#FBE3D6",
  Expt_data_sub$Stressed == "Stress" ~ "#E98F2B"
)


freqColorBar <- case_when(
  Expt_data_sub$Frequency == "Dark" ~ "#861C33",
  Expt_data_sub$Frequency == "10Hz" ~ "#2D567C",
  Expt_data_sub$Frequency == "40Hz" ~ "#19B799"
)



SexColorBar <- case_when(
  Expt_data_sub$Sex == "M" ~  "#872EFE",
  Expt_data_sub$Sex == "F" ~ "#D98DFF" 
)





#Clustering Genes 
hrGenes= hclust(dist((genes_sub_Z_mat),method = "euclidean"), method="ward.D2")


#Generating Heatmap 

# Genes Clustered Heatmap
pdf(paste0("./Neurons/061924/Neuron All Genes Clustered Heatmap All Subjects Male Only.pdf"), width=6,height=8,pointsize = 14, useDingbats = FALSE)
#png(paste0("Neurons/061924/Neuron All Genes Clustered Heatmap All Subjects Male Only.png"), width=7,height=8,pointsize = 14, units="in",res=600)
heatmap3(genes_sub_Z_mat, ColSideColors =cbind(c(SexColorBar),c(freqColorBar),c(StressColorBar)), ColSideWidth=1, ColSideLabs=NA,  
         col=barColors, breaks=breakBarColors, legendfun=function()showLegend(legend=c(NA),col=c(NA)),
         Rowv=as.dendrogram(hrGenes), 
         Colv=NA,  scale="none",
         margins = c(3,10),
         labRow = NA, labCol=NA,
         cexCol=1.3)
title("Genes Cluster Neuron", adj = 0.5, line = 2)
dev.off()




#Cluster Analysis



cutThreshCTD=3.76 #use with the euclidian distance method, average
myclCTD2 = cutree(hrGenes, k=8);
colorBrewerInterp=colorRampPalette(brewer.pal(8,"Spectral"))
mycolhcCTD2 = sample(colorBrewerInterp(length(unique(myclCTD2))))
mycolhcCTD2 = mycolhcCTD2[as.vector(myclCTD2)]
t1CTD2=as.dendrogram(hrGenes)
rowOrder2=order.dendrogram(t1CTD2);
GeneName=rownames(genes_sub_Z_mat)
GeneName2=sub("\\|.*", "", GeneName)

GroupsOut2=tibble(mycolhcCTD2[rowOrder2],rownames(genes_sub_Z_mat)[rowOrder2], str_to_title(GeneName2[rowOrder2]))


colnames(GroupsOut2)=c("ClusterColor","Genes")
GroupsOut2 <- mutate(GroupsOut2, Cluster=case_when(
  ClusterColor == unique(GroupsOut2$ClusterColor)[1] ~ 1,
  ClusterColor == unique(GroupsOut2$ClusterColor)[2] ~ 2,
  ClusterColor == unique(GroupsOut2$ClusterColor)[3] ~ 3,
  ClusterColor == unique(GroupsOut2$ClusterColor)[4] ~ 4,
  ClusterColor == unique(GroupsOut2$ClusterColor)[5] ~ 5,
  ClusterColor == unique(GroupsOut2$ClusterColor)[6] ~ 6,
  ClusterColor == unique(GroupsOut2$ClusterColor)[7] ~ 7,
  ClusterColor == unique(GroupsOut2$ClusterColor)[8] ~ 8
))
saveRDS(GroupsOut2,file="./Neurons/061924/Cluster Analysis Genes Male Only.RDS")

export(data.frame(GroupsOut2),file="Neurons/061924/Cluster Analysis Male Only.csv")
Clusters=readRDS("./Neurons/061924/Cluster Analysis Genes Stressed Male Only.RDS")

png(paste0("./Neurons/061924/Genes Clustered Heatmap Male Only side color.png"), width=7,height=8,pointsize = 14, units="in",res=600)
heatmap3(genes_sub_Z_mat, RowSideColors=mycolhcCTD2,ColSideColors =cbind(c(SexColorBar),c(freqColorBar),c(StressColorBar)), ColSideWidth=1, ColSideLabs=NA,  
         col=barColors, breaks=breakBarColors, legendfun=function()showLegend(legend=c(NA),col=c(NA)),
         Rowv=as.dendrogram(hrGenes), 
         Colv=NA,  scale="none",
         margins = c(3,10),
         labRow = NA, labCol=NA,
         cexCol=1.3)
title("Genes Cluster", adj = 0.5, line = 2)
dev.off()


###Neuron  Specific Gene Set (Custome)
sortby <- which(Expt_data$Condition =="M 10HzHz Stress " | Expt_data$Condition  =="M 40HzHz Stress " | Expt_data$Condition  =="M DarkHz Stress " | Expt_data$Condition=="M DarkHz Control " )
#This  can be changed to any group or condition
Expt_data_sub=Expt_data[sortby,]
genes_sub=genes[,sortby]



GeneName=rownames(genes_sub)
GeneName2=sub("\\|.*", "", GeneName)



genes_sub_Z <- genes_sub %>%
  apply(1,"scale") %>%
  t() %>%
  as_tibble()
colnames(genes_sub_Z) = colnames(genes_sub)
genes_sub_Z$geneID = rownames(genes_sub)
genes_sub_Z[is.na(genes_sub_Z)] <- 0   # set any NA values to 0
rownames(genes_sub)=GeneName2
# Create a matrix of Z-scored subset data
genes_sub_Z_mat <- as.matrix(genes_sub_Z[,1:{ncol(genes_sub_Z)-1}])
rownames(genes_sub_Z_mat) <- str_to_upper(sub("\\|.*", "", (genes_sub_Z$geneID)))




# GSVA using Neuron Annotations
geneSetsMat=read_excel("Reference Files/GeneSets_Neurons_FXN1.xlsx")
geneSets=as.list(geneSetsMat)
geneSets=lapply(geneSets, function(x) x[!is.na(x)])


geneSetEnrich=gsva(genes_sub_Z_mat, geneSets, mx.diff=TRUE, kcdf="Gaussian")
geneSetEnrichSortZ=t(apply(geneSetEnrich,1,scale)); #z-score the data
geneSetEnrichSortZ[is.na(geneSetEnrichSortZ)]=0
colnames(geneSetEnrichSortZ)=colnames(genes_sub_Z_mat)


Condition=Expt_data_sub$Condition
Condition <- as.character(Condition)
Condition <- factor(Condition, levels=c("M DarkHz Control ","M DarkHz Stress ","M 10HzHz Stress ","M 40HzHz Stress "))
datagroupVars=data.frame(t(geneSetEnrichSortZ), Condition)
for (ind in 1:16)
{
  
 # png(paste0("Neurons/061924/NeuronCustom GS/GS boxplot All Male",ind,".png"),height=4,width=5,units="in",res=1000)
   pdf(paste0("Neurons/061924/NeuronCustom GS/GS boxplot All Male",ind,".pdf"),height=4,width=5)
  print({
    ggplot(datagroupVars, aes(x=Condition, y=geneSetEnrichSortZ[ind,], color=Condition)) + 
      geom_boxplot()+
      geom_beeswarm(aes(fill = Condition),size =6,shape=21, stroke = 1, cex = 2,groupOnX = TRUE)+
      theme(panel.background = element_rect(fill = 'white', colour = 'black'),
            text = element_text(size=20),
            axis.text.x = element_text(face = "bold", color = "black",
                                       size = 8 ,angle=45),
            axis.text.y = element_text(face = "bold", color = "black",
                                       size = 12),
            plot.title = element_text(hjust = 0.5))+
      scale_color_manual(values=c("#747474","#861C33","#2D567C","#19B799"))+
      scale_fill_manual(values=c("#747474","#861C33","#2D567C","#19B799"))+
      
      #stat_compare_means(method="wilcox", comparisons=list(c(3,4),c(2,4),c(1,4)),label="p.format",tip.length = 0.03) +
      xlab("")+
      ylab("Enrichment Score [a.u.]")+
      labs(color="Condition")+
      ylim(c(-2, 2)) +
      ggtitle(colnames(datagroupVars)[ind])
  })
  dev.off()
}


hrGSVA= hclust(dist((geneSetEnrichSortZ),method = "euclidean"), method="ward.D2")


png(paste0("./Neurons/061924/NeuronCustom GS/Custom GS Neuron Male All Heatmap.png"), width=6,height=8,units="in",res=600, pointsize = 14)
#pdf(paste0("./Neurons/061924/NeuronCustom GS/Custom GS Neuron Male All Heatmap.pdf"), width=6,height=8,pointsize = 14, useDingbats = FALSE)
heatmap3(geneSetEnrichSortZ, ColSideColors =cbind(c(SexColorBar),c(freqColorBar),c(StressColorBar)), ColSideWidth=1, ColSideLabs=NA,
         col=barColors, breaks=breakBarColors, legendfun=function()showLegend(legend=c(NA),col=c(NA)),
         Rowv=as.dendrogram(hrGSVA), Colv=NA,  scale="none",
         labRow = rownames(geneSetEnrich), labCol = NA, cexRow=1, margins=c(5,10))
title("Custom GS Neurons", adj = 0.5, line = 2)
dev.off()

# 
# #Stats
# # Light vs XX
# ind1 = unique(which(Expt_data_sub2$ConditionOrder==1))
# ind2 = unique(which(Expt_data_sub2$ConditionOrder==10))
# resultMat = tibble(geneSetEnrichSortZ[,ind1],geneSetEnrichSortZ[,ind2])
# colnames(resultMat)=c("Light 1hr","40Hz 4hr")
# resultMat$pvalsWilcox=NA
# resultMat$pvalsttest=NA
# resultMat$pvalsadjustBF=NA
# resultMat$pvalsadjustfdr=NA
# 
# for(i in 1:nrow(geneSetEnrichSortZ)){
#   resultMat$pvalsWilcox[i] = wilcox.test(geneSetEnrichSortZ[i,ind1],geneSetEnrichSortZ[i,ind2])$p.value
#   resultMat$pvalsttest[i] = t.test(geneSetEnrichSortZ[i,ind1],geneSetEnrichSortZ[i,ind2])$p.value
# }
# resultMat$pvalsadjustBF = p.adjust(resultMat$pvalsttest , method="bonferroni",n=nrow(resultMat))
# resultMat$pvalsadjustfdr = p.adjust(resultMat$pvalsttest , method="fdr",n=nrow(resultMat))
# 
# resultMat$GS = rownames(geneSetEnrichSortZ)
# rownames(resultMat)=rownames(geneSetEnrichSortZ)
# write.xlsx(data.frame(resultMat), file=paste("./OCT 23/WGNCA/5xFAD VC/Neuron Custome GS/021224_Stats 40Hz 4hr vs Light.xlsx"))
# 



#-------------------------------------------------Female Data------------------------------------------------

sortby <- which(Expt_data$Condition =="F 10HzHz Stress " | Expt_data$Condition  =="F 40HzHz Stress " | Expt_data$Condition  =="F DarkHz Stress " | Expt_data$Condition=="F DarkHz Control " )
#This  can be changed to any group or condition
Expt_data_sub=Expt_data[sortby,]
genes_sub=genes[,sortby]


GeneName=rownames(genes_sub)
GeneName2=sub("\\|.*", "", GeneName)



genes_sub_Z <- genes_sub %>%
  apply(1,"scale") %>%
  t() %>%
  as_tibble()
colnames(genes_sub_Z) = colnames(genes_sub)
genes_sub_Z$geneID = rownames(genes_sub)
genes_sub_Z[is.na(genes_sub_Z)] <- 0   # set any NA values to 0
rownames(genes_sub)=GeneName2
# Create a matrix of Z-scored subset data
genes_sub_Z_mat <- as.matrix(genes_sub_Z[,1:{ncol(genes_sub_Z)-1}])
rownames(genes_sub_Z_mat) <- str_to_upper(genes_sub_Z$geneID)



StressColorBar <- case_when(
  Expt_data_sub$Stressed == "Control" ~ "#FBE3D6",
  Expt_data_sub$Stressed == "Stress" ~ "#E98F2B"
)


freqColorBar <- case_when(
  Expt_data_sub$Frequency == "Dark" ~ "#861C33",
  Expt_data_sub$Frequency == "10Hz" ~ "#2D567C",
  Expt_data_sub$Frequency == "40Hz" ~ "#19B799"
)



SexColorBar <- case_when(
  Expt_data_sub$Sex == "M" ~  "#872EFE",
  Expt_data_sub$Sex == "F" ~ "#D98DFF" 
)





#Clustering Genes 
hrGenes= hclust(dist((genes_sub_Z_mat),method = "euclidean"), method="ward.D2")


#Generating Heatmap 

# Genes Clustered Heatmap
pdf(paste0("./Neurons/061924/Neuron All Genes Clustered Heatmap All Subjects Female Only.pdf"), width=6,height=8,pointsize = 14, useDingbats = FALSE)
#png(paste0("Neurons/061924/Neuron All Genes Clustered Heatmap All Subjects Female Only.png"), width=7,height=8,pointsize = 14, units="in",res=600)
heatmap3(genes_sub_Z_mat, ColSideColors =cbind(c(SexColorBar),c(freqColorBar),c(StressColorBar)), ColSideWidth=1, ColSideLabs=NA,  
         col=barColors, breaks=breakBarColors, legendfun=function()showLegend(legend=c(NA),col=c(NA)),
         Rowv=as.dendrogram(hrGenes), 
         Colv=NA,  scale="none",
         margins = c(3,10),
         labRow = NA, labCol=NA,
         cexCol=1.3)
title("Genes Cluster Neuron", adj = 0.5, line = 2)
dev.off()




#Cluster Analysis



cutThreshCTD=3.76 #use with the euclidian distance method, average
myclCTD2 = cutree(hrGenes, k=8);
colorBrewerInterp=colorRampPalette(brewer.pal(8,"Spectral"))
mycolhcCTD2 = sample(colorBrewerInterp(length(unique(myclCTD2))))
mycolhcCTD2 = mycolhcCTD2[as.vector(myclCTD2)]
t1CTD2=as.dendrogram(hrGenes)
rowOrder2=order.dendrogram(t1CTD2);
GeneName=rownames(genes_sub_Z_mat)
GeneName2=sub("\\|.*", "", GeneName)

GroupsOut2=tibble(mycolhcCTD2[rowOrder2],rownames(genes_sub_Z_mat)[rowOrder2], str_to_title(GeneName2[rowOrder2]))


colnames(GroupsOut2)=c("ClusterColor","Genes")
GroupsOut2 <- mutate(GroupsOut2, Cluster=case_when(
  ClusterColor == unique(GroupsOut2$ClusterColor)[1] ~ 1,
  ClusterColor == unique(GroupsOut2$ClusterColor)[2] ~ 2,
  ClusterColor == unique(GroupsOut2$ClusterColor)[3] ~ 3,
  ClusterColor == unique(GroupsOut2$ClusterColor)[4] ~ 4,
  ClusterColor == unique(GroupsOut2$ClusterColor)[5] ~ 5,
  ClusterColor == unique(GroupsOut2$ClusterColor)[6] ~ 6,
  ClusterColor == unique(GroupsOut2$ClusterColor)[7] ~ 7,
  ClusterColor == unique(GroupsOut2$ClusterColor)[8] ~ 8
))
saveRDS(GroupsOut2,file="./Neurons/061924/Cluster Analysis Genes Female Only.RDS")

export(data.frame(GroupsOut2),file="Neurons/061924/Cluster Analysis Female Only.csv")
Clusters=readRDS("./Neurons/061924/Cluster Analysis Genes Female Only.RDS")

png(paste0("./Neurons/061924/Genes Clustered Heatmap Female Only side color.png"), width=7,height=8,pointsize = 14, units="in",res=600)
heatmap3(genes_sub_Z_mat, RowSideColors=mycolhcCTD2,ColSideColors =cbind(c(SexColorBar),c(freqColorBar),c(StressColorBar)), ColSideWidth=1, ColSideLabs=NA,  
         col=barColors, breaks=breakBarColors, legendfun=function()showLegend(legend=c(NA),col=c(NA)),
         Rowv=as.dendrogram(hrGenes), 
         Colv=NA,  scale="none",
         margins = c(3,10),
         labRow = NA, labCol=NA,
         cexCol=1.3)
title("Genes Cluster", adj = 0.5, line = 2)
dev.off()


###Neuron  Specific Gene Set (Custome)
sortby <- which(Expt_data$Condition =="F 10HzHz Stress " | Expt_data$Condition  =="F 40HzHz Stress " | Expt_data$Condition  =="F DarkHz Stress " | Expt_data$Condition=="F DarkHz Control " )
#This  can be changed to any group or condition
Expt_data_sub=Expt_data[sortby,]
genes_sub=genes[,sortby]



GeneName=rownames(genes_sub)
GeneName2=sub("\\|.*", "", GeneName)



genes_sub_Z <- genes_sub %>%
  apply(1,"scale") %>%
  t() %>%
  as_tibble()
colnames(genes_sub_Z) = colnames(genes_sub)
genes_sub_Z$geneID = rownames(genes_sub)
genes_sub_Z[is.na(genes_sub_Z)] <- 0   # set any NA values to 0
rownames(genes_sub)=GeneName2
# Create a matrix of Z-scored subset data
genes_sub_Z_mat <- as.matrix(genes_sub_Z[,1:{ncol(genes_sub_Z)-1}])
rownames(genes_sub_Z_mat) <- str_to_upper(sub("\\|.*", "", (genes_sub_Z$geneID)))




# GSVA using Neuron Annotations
geneSetsMat=read_excel("Reference Files/GeneSets_Neurons_FXN1.xlsx")
geneSets=as.list(geneSetsMat)
geneSets=lapply(geneSets, function(x) x[!is.na(x)])


geneSetEnrich=gsva(genes_sub_Z_mat, geneSets, mx.diff=TRUE, kcdf="Gaussian")
geneSetEnrichSortZ=t(apply(geneSetEnrich,1,scale)); #z-score the data
geneSetEnrichSortZ[is.na(geneSetEnrichSortZ)]=0
colnames(geneSetEnrichSortZ)=colnames(genes_sub_Z_mat)


Condition=Expt_data_sub$Condition
Condition <- as.character(Condition)
Condition <- factor(Condition, levels=c("F DarkHz Control ","F DarkHz Stress ","F 10HzHz Stress ","F 40HzHz Stress "))
datagroupVars=data.frame(t(geneSetEnrichSortZ), Condition)
for (ind in 1:16)
{
  
  #png(paste0("Neurons/061924/NeuronCustom GS/GS boxplot All Female",ind,".png"),height=4,width=5,units="in",res=1000)
  pdf(paste0("Neurons/061924/NeuronCustom GS/GS boxplot All Female",ind,".pdf"),height=4,width=5)
  print({
    ggplot(datagroupVars, aes(x=Condition, y=geneSetEnrichSortZ[ind,], color=Condition)) + 
      geom_boxplot()+
      geom_beeswarm(aes(fill = Condition),size =6,shape=21, stroke = 1, cex = 2,groupOnX = TRUE)+
      theme(panel.background = element_rect(fill = 'white', colour = 'black'),
            text = element_text(size=20),
            axis.text.x = element_text(face = "bold", color = "black",
                                       size = 8 ,angle=45),
            axis.text.y = element_text(face = "bold", color = "black",
                                       size = 12),
            plot.title = element_text(hjust = 0.5))+
      scale_color_manual(values=c("#747474","#861C33","#2D567C","#19B799"))+
      scale_fill_manual(values=c("#747474","#861C33","#2D567C","#19B799"))+
      
      #stat_compare_means(method="wilcox", comparisons=list(c(3,4),c(2,4),c(1,4)),label="p.format",tip.length = 0.03) +
      xlab("")+
      ylab("Enrichment Score [a.u.]")+
      labs(color="Condition")+
      ylim(c(-2, 2)) +
      ggtitle(colnames(datagroupVars)[ind])
  })
  dev.off()
}


hrGSVA= hclust(dist((geneSetEnrichSortZ),method = "euclidean"), method="ward.D2")


#png(paste0("./Neurons/061924/NeuronCustom GS/Custom GS Neuron Female All Heatmap.png"), width=6,height=8,units="in",res=600, pointsize = 14)
pdf(paste0("./Neurons/061924/NeuronCustom GS/Custom GS Neuron Female All Heatmap.pdf"), width=6,height=8,pointsize = 14, useDingbats = FALSE)
heatmap3(geneSetEnrichSortZ, ColSideColors =cbind(c(SexColorBar),c(freqColorBar),c(StressColorBar)), ColSideWidth=1, ColSideLabs=NA,
         col=barColors, breaks=breakBarColors, legendfun=function()showLegend(legend=c(NA),col=c(NA)),
         Rowv=as.dendrogram(hrGSVA), Colv=NA,  scale="none",
         labRow = rownames(geneSetEnrich), labCol = NA, cexRow=1, margins=c(5,10))
title("Custom GS Neurons", adj = 0.5, line = 2)
dev.off()

#-------------------------------------------------------------------------------------------------------------------#

#Exporting Gene Names All
GeneName=rownames(genes_sub_Z_mat)
GeneName2=sub("\\|.*", "", GeneName)

export(data.frame(str_to_title(GeneName2)),file="Neurons/061924/All Genes Neurons 061224.csv")

#-----Exporting Legend colors
png(("legend_Sex 061924.png"), width=6,height=8,units="in",res=600, pointsize = 14)
#pdf("legend_Sex 061924.pdf",height=5,width=4,useDingbats = FALSE)
plot(NULL,xaxt="n",yaxt="n",bty="n",ylab="",xlab="",xlim=0:1,ylim=0:1)
legend("topleft",title="Sex",legend=c("M","F"),lty=1,lwd=10,cex=1.25,bty="n",col=c("#872EFE","#D98DFF" ))
dev.off()

png(("legend_freq 061924.png"), width=6,height=8,units="in",res=600, pointsize = 14)
#pdf("legend_freq 061924.pdf",height=5,width=4,useDingbats = FALSE)
plot(NULL,xaxt="n",yaxt="n",bty="n",ylab="",xlab="",xlim=0:1,ylim=0:1)
legend("topleft",title="Frequency",legend=c("Dark","10Hz","40Hz"),lty=1,lwd=10,cex=1.25,bty="n",col=c("#861C33","#2D567C","#19B799"))
dev.off()

png(("legend_Stress 061924.png"), width=6,height=8,units="in",res=600, pointsize = 14)
#pdf("legend_Stress 061924.pdf",height=5,width=4,useDingbats = FALSE)
plot(NULL,xaxt="n",yaxt="n",bty="n",ylab="",xlab="",xlim=0:1,ylim=0:1)
legend("topleft",title="Stress",legend=c("Control","Stress"),lty=1,lwd=10,cex=1.25,bty="n",col=c("#FBE3D6","#E98F2B"))
dev.off()




#-----------------------------------------------------------------------------------------------------------------------#


# Count Plots for GO terms


Data <- read_excel("./Neurons/061924/Neuronal Biological Processes_6 GO terms selection_June2024 v2.xlsx",sheet= "BP_10Hz_M_Up_NC v2")
Data=data.frame(Data)

Name=Data$GO.biological.process.complete
pValFDR=Data$Client.Text.Box.Input..FDR.
FE=Data$Client.Text.Box.Input..fold.Enrichment.
Count=Data$Count
pval=Data$Client.Text.Box.Input..raw.P.value.
FELog=log2(FE)
data2=data.frame(Name,log2(FE),Count,pval)

#png(paste0("./Neurons/061924/GOterms/BP_10Hz_M_Up_NC v2.png"), width=10,height=8,pointsize = 20, units="in",res=600)
pdf("./Neurons/061924/GOterms/BP_10Hz_M_Up_NC v2.pdf", width=10,height=10,pointsize = 6)
data2 %>% 
  as.tibble() %>% 
  ggplot(aes(reorder(Name,FELog), FELog)) + 
  geom_count(aes(size = Count, color= -log10(pval)) , position = position_dodge2(width = 0.25, preserve = "single", padding = -0.25)) +
  scale_colour_gradient(low = "violetred",
                        high = "royalblue")+
  #ylim(0.1,6)+
  coord_flip()+
  labs(y= "Log2 fold enrichment", x = "")+
  scale_size(range=c(5, 10))+
  theme_bw()
dev.off()




#-----------------Venn Diagram of Neuronal DEGs vs publihsed study markers-------------------#

Data_GeneName= read_excel("Neurons/061924/062124_Venn Diagram_Flicker DEGs vs Published Genes.xlsx",sheet="Res Sus") 

DEG_M_10Hz_All=Data_GeneName$`10Hz_M DEG all`
DEG_M_40Hz_All=Data_GeneName$`40Hz_M DEG all`
DEG_M_40Hzvs10Hz_All=Data_GeneName$`40v10Hz_M DEG all`
DEG_M_10Hzvs_F_40Hz_All=Data_GeneName$`40Fv10M DEG all`
DEG_F_10Hz_All=Data_GeneName$`10Hz_F DEG all`
DEG_F_40Hz_All=Data_GeneName$`40Hz_F DEG all`
DEG_F_40Hzvs10Hz_All=Data_GeneName$`40v10Hz_F DEG all`
Resiliant=Data_GeneName$Resilient
Susceptible=Data_GeneName$Susceptible
AdultMale=Data_GeneName$`Adult male DEG`
ELSMale=Data_GeneName$`ELS male DEG`
AdultFemale=Data_GeneName$`Adult female DEG`
ELSFemale=Data_GeneName$`ELS female DEG`
PublishedGenes=Data_GeneName$`Genes previously associated with stress in either humans_animal models or both`

indKeep1=match(DEG_M_10Hz_All,AdultMale)
DEG_M_10Hz_AllSharedWithAdultMale<- DEG_M_10Hz_All[c(indKeep1)]
export(data.frame(DEG_M_10Hz_AllSharedWithAdultMale),file="Neurons/061924/Deseq2/Shared Neuronal Genes_Male Adult vs Male 10Hz all .csv")

indKeep2=match(DEG_M_40Hz_All,AdultMale)
DEG_M_40Hz_AllSharedWithAdultMale<- DEG_M_40Hz_All[c(indKeep2)]
export(data.frame(DEG_M_40Hz_AllSharedWithAdultMale),file="Neurons/061924/Deseq2/Shared Neuronal Genes_Male Adult vs Male 40Hz all .csv")


indKeep3=match(DEG_M_40Hz_All,ELSMale)
DEG_M_40Hz_AllSharedWithAdultMale<- DEG_M_40Hz_All[c(indKeep3)]
export(data.frame(DEG_M_40Hz_AllSharedWithAdultMale),file="Neurons/061924/Deseq2/Shared Neuronal Genes_Male ESL vs Male 40Hz all .csv")


indKeep4=match(DEG_M_10Hz_All,ELSMale)
DEG_M_10Hz_AllSharedWithAdultMale<- DEG_M_10Hz_All[c(indKeep4)]
export(data.frame(DEG_M_10Hz_AllSharedWithAdultMale),file="Neurons/061924/Deseq2/Shared Neuronal Genes_Male ESL vs Male 10Hz all .csv")





indKeep5=match(DEG_F_10Hz_All,AdultFemale)
DEG_F_10Hz_AllSharedWithAdultMale<- DEG_F_10Hz_All[c(indKeep5)]
export(data.frame(DEG_F_10Hz_AllSharedWithAdultMale),file="Neurons/061924/Deseq2/Shared Neuronal Genes_Female Adult vs Female 10Hz all .csv")

indKeep6=match(DEG_F_40Hz_All,AdultFemale)
DEG_F_40Hz_AllSharedWithAdultFemale<- DEG_F_40Hz_All[c(indKeep6)]
export(data.frame(DEG_F_40Hz_AllSharedWithAdultFemale),file="Neurons/061924/Deseq2/Shared Neuronal Genes_Female Adult vs Female 40Hz all .csv")

indKeep7=match(DEG_F_40Hz_All,ELSFemale)
DEG_F_40Hz_AllSharedWithAdultFemale<- DEG_F_40Hz_All[c(indKeep7)]
export(data.frame(DEG_F_40Hz_AllSharedWithAdultFemale),file="Neurons/061924/Deseq2/Shared Neuronal Genes_Female ESL vs Female 40Hz all .csv")


indKeep8=match(DEG_F_10Hz_All,ELSFemale)
DEG_F_10Hz_AllSharedWithAdultFemale<- DEG_F_10Hz_All[c(indKeep8)]
export(data.frame(DEG_F_10Hz_AllSharedWithAdultFemale),file="Neurons/061924/Deseq2/Shared Neuronal Genes_Female ESL vs Female 10Hz all .csv")


indKeep9=match(DEG_F_40Hz_All,PublishedGenes)
DEG_F_40Hz_AllShared<- DEG_F_40Hz_All[c(indKeep9)]
export(data.frame(DEG_F_40Hz_AllShared),file="Neurons/061924/Deseq2/Shared Neuronal Genes_Previously Published vs Female 40Hz all .csv")


indKeep10=match(DEG_F_10Hz_All,PublishedGenes)
DEG_F_10Hz_AllShared<- DEG_F_10Hz_All[c(indKeep10)]
export(data.frame(DEG_F_10Hz_AllShared),file="Neurons/061924/Deseq2/Shared Neuronal Genes_Previously PublishedL vs Female 10Hz all .csv")


indKeep11=match(DEG_M_40Hz_All,PublishedGenes)
DEG_M_40Hz_AllShared<- DEG_M_40Hz_All[c(indKeep11)]
export(data.frame(DEG_M_40Hz_AllShared),file="Neurons/061924/Deseq2/Shared Neuronal Genes_Previously Published vs Male 40Hz all .csv")


indKeep12=match(DEG_M_10Hz_All,PublishedGenes)
DEG_M_10Hz_AllShared<- DEG_M_10Hz_All[c(indKeep12)]
export(data.frame(DEG_M_10Hz_AllShared),file="Neurons/061924/Deseq2/Shared Neuronal Genes_Previously Published vs Male 10Hz all .csv")





output_folder = paste0("Neurons/061924/Deseq2/")
#VennFreq1=list(Male10L=DEG_M_10Hz_All, Female40L=DEG_F_40Hz_All)
#VennFreq2=list(Male10L=DEG_M_10Hz_All, Male40L= DEG_M_40Hz_All,Female10L=DEG_F_10Hz_All , Female40L=DEG_F_40Hz_All )

VennFreq1=list(Male10L=DEG_M_10Hz_All, Male40HzL=DEG_M_40Hz_All,AdultM= AdultMale , ELSM=ELSMale)
#png(paste0(output_folder,"Venn Diagram Male 10Hz vs Male 40Hz vs AdultMale VS ELSMale.png"), width=7,height=8,pointsize = 14,units="in",res=600)
pdf(paste0(output_folder,"Venn Diagram Male 10Hz vs Male 40Hz vs AdultMale VS ELSMale.pdf"), width=8,height=8,pointsize = 12)
ggvenn(VennFreq1, 
       #fill_color = c("#872EFE","#D98DFF", ),
       stroke_size = 0.5, set_name_size = 6, text_size = 5)
dev.off()

VennFreq2=list(Male10L=DEG_M_10Hz_All,Male40HzL=DEG_M_40Hz_All, Resil= Resiliant , Supc=Susceptible)
#png(paste0(output_folder,"Venn Diagram Male 10Hz vs Male 40Hz vs Susceptible vs Resiliant.png"), width=7,height=8,pointsize = 14,units="in",res=600)
pdf(paste0(output_folder,"Venn Diagram Male 10Hz vs Male 40Hz vs Susceptible vs Resiliant.pdf"), width=8,height=8,pointsize = 12)
ggvenn(VennFreq2, 
       #fill_color = c("#872EFE","#D98DFF", ),
       stroke_size = 0.5, set_name_size = 6, text_size = 5)
dev.off()

VennFreq3=list(Female10L=DEG_F_10Hz_All,Female40HzL=DEG_F_40Hz_All, Resil= Resiliant , Supc=Susceptible)
#png(paste0(output_folder,"Venn Diagram Female 10Hz vs Female 40Hz vs Susceptible vs Resiliant.png"), width=7,height=8,pointsize = 14,units="in",res=600)
pdf(paste0(output_folder,"Venn Diagram Female 10Hz vs Female 40Hz vs Susceptible vs Resiliant.pdf"), width=8,height=8,pointsize = 12)
ggvenn(VennFreq3, 
       #fill_color = c("#872EFE","#D98DFF", ),
       stroke_size = 0.5, set_name_size = 6, text_size = 5)
dev.off()


VennFreq4=list(Female10L=DEG_F_10Hz_All, Female40HzL=DEG_F_40Hz_All,AdultF= AdultFemale , ELSF=ELSFemale)
#png(paste0(output_folder,"Venn Diagram Female 10Hz vs Female 40Hz vs AdultFemale VS ELSFemale.png"), width=7,height=8,pointsize = 14,units="in",res=600)
pdf(paste0(output_folder,"Venn Diagram Female 10Hz vs Female 40Hz vs AdultFemale VS ELSFemale.pdf"), width=8,height=8,pointsize = 12)
ggvenn(VennFreq4, 
       #fill_color = c("#872EFE","#D98DFF", ),
       stroke_size = 0.5, set_name_size = 6, text_size = 5)
dev.off()


VennFreq5=list(Female10L=DEG_F_10Hz_All, Female40HzL=DEG_F_40Hz_All,StressGene=PublishedGenes)
png(paste0(output_folder,"Venn Diagram Female 10Hz vs Female 40Hz vs previosuly Published Genes.png"), width=7,height=8,pointsize = 14,units="in",res=600)
#pdf(paste0(output_folder,"Venn Diagram Female 10Hz vs Female 40Hz vs previosuly Published Genes.pdf"), width=8,height=8,pointsize = 12)
ggvenn(VennFreq5, 
       #fill_color = c("#872EFE","#D98DFF", ),
       stroke_size = 0.5, set_name_size = 6, text_size = 5)
dev.off()

VennFreq6=list(Male10L=DEG_M_10Hz_All, Male40HzL=DEG_M_40Hz_All,StressGene=PublishedGenes)
png(paste0(output_folder,"Venn Diagram Male 10Hz vs Male 40Hz vs previosuly Published Genes.png"), width=7,height=8,pointsize = 14,units="in",res=600)
#pdf(paste0(output_folder,"Venn Diagram Male 10Hz vs Male 40Hz vs previosuly Published Genes.pdf"), width=8,height=8,pointsize = 12)
ggvenn(VennFreq6, 
       #fill_color = c("#872EFE","#D98DFF", ),
       stroke_size = 0.5, set_name_size = 6, text_size = 5)
dev.off()

# #Male 40Hz specific 
# MaleData40=genes_subZ_mat_M[match(MalesigOE40$rownames.sigOE_ordered. , sub("\\|.*", "",rownames(genes_subZ_mat_M))) ,]
# MaleDataOthers<-genes_subZ_mat_M[which(Male40L %in% c(Male10L,Female10L,Female40L)),]
# 
# indRemove=match(rownames(MaleDataOthers), rownames(genes_subZ_mat_M))
# 
# MaleData40only<- MaleData40[-c(indRemove) ,]
# 
# export(data.frame(rownames(MaleData40only),MaleData40only),file="Neurons/Deseq2/Male Neurons 40Hz specific DEGs on venn diagram Up.csv")
# 
# #Male down 10Hz specific 
# 
# MaleData10=genes_subZ_mat_M[match(MalesigOE10$rownames.sigOE_ordered. , sub("\\|.*", "",rownames(genes_subZ_mat_M))) ,]
# MaleDataOthers<-genes_subZ_mat_M[which(Male10L %in% c(Male40L,Female10L,Female40L)),]
# 
# 
# indRemove=match(rownames(MaleDataOthers), rownames(genes_subZ_mat_M))
# 
# MaleData10only<- MaleData10[-c(indRemove) ,]
# 
# export(data.frame(rownames(MaleData10only),MaleData10only),file="Neurons/Deseq2/Male Neurons 10Hz specific DEGs on venn diagram Up.csv")
# 
# 
# 
# #Female down 10Hz specific 
# 
# FemaleData10=genes_subZ_mat_F[match(FemalesigOE10$rownames.sigOE_ordered. , sub("\\|.*", "",rownames(genes_subZ_mat_F))) ,]
# FemaleDataOthers<-genes_subZ_mat_F[which(Female10L %in% c(Male40L,Male10L,Female40L)),]
# 
# 
# indRemove=match(rownames(FemaleDataOthers), rownames(genes_subZ_mat_F))
# 
# FemaleData10only<- FemaleData10[-c(indRemove) ,]
# 
# export(data.frame(rownames(FemaleData10only),FemaleData10only),file="Neurons/Deseq2/Female Neurons 10Hz specific DEGs on venn diagram Up.csv")
# 
# 
# 
# 
# #Female down 40Hz specific 
# 
# FemaleData40=genes_subZ_mat_F[match(FemalesigOE40$rownames.sigOE_ordered. , sub("\\|.*", "",rownames(genes_subZ_mat_F))) ,]
# FemaleDataOthers<-genes_subZ_mat_F[which(Female40L %in% c(Male40L,Male10L,Female10L)),]
# 
# 
# indRemove=match(rownames(FemaleDataOthers), rownames(genes_subZ_mat_F))
# 
# FemaleData40only<- FemaleData40[-c(indRemove) ,]
# 
# export(data.frame(rownames(FemaleData40only),FemaleData40only),file="Neurons/Deseq2/Female Neurons 40Hz specific DEGs on venn diagram Up.csv")
# 
# 
# output_folder = paste0("Neurons/Deseq2/")
# VennFreq=list(Male40L=MalesigOE40$rownames.sigOE_ordered.,Male10L=MalesigOE10$rownames.sigOE_ordered.,Female10L=FemalesigOE10$rownames.sigOE_ordered., Female40L=FemalesigOE40$rownames.sigOE_ordered.)
# 
# png(paste0(output_folder,"Venn Diagram DEGs Neurons Up Only.png"), width=7,height=8,pointsize = 14,units="in",res=600)
# #pdf(paste0(output_folder,"Venn Diagram DEGs Neurons Up Only.pdf"), width=8,height=8,pointsize = 12)
# ggvenn(VennFreq, 
#        #fill_color = c("#E41A1C","#377EB8"),
#        stroke_size = 0.5, set_name_size = 6, text_size = 5)
# dev.off()
# 
# 
# 
# 
# 
# 
# 
# #Female
# 
# Female40 <- read_excel("Neurons/Deseq2/Neuronal DEGs_Stress 40Hz vs Light_Female DESeq2 ThreshholdPassing FC2.xlsx",sheet="Up") 
# FemalesigOE40 <- data.frame(subset( Female40, threshold==TRUE))
# 
# Female10 <- read_excel("Neurons/Deseq2/Neuronal DEGs_Stress 10Hz vs Light_Female DESeq2 ThreshholdPassing FC2.xlsx",sheet="Up") 
# FemalesigOE10 <- data.frame(subset( Female10, threshold==TRUE))
# # 
# # Female1040 <- read_excel("Neurons/Deseq2/Tina Data Stress 40Hz vs 10 Female DESeq2 All FC2 111323SB.xlsx",sheet="Up") 
# # FemalesigOE1040 <- data.frame(subset( Female1040, threshold==TRUE))
# # 
# 
# #Male
# 
# 
# Male40 <- read_excel("Neurons/Deseq2/Neuronal DEGs_Stress 40Hz vs Light_Male DESeq2 ThreshholdPassing FC2.xlsx",sheet="Up") 
# MalesigOE40 <- data.frame(subset(Male40, threshold==TRUE))
# 
# 
# Male10 <- read_excel("Neurons/Deseq2/Neuronal DEGs_Stress 10Hz vs Light_Male DESeq2 ThreshholdPassing FC2.xlsx",sheet="Up") 
# MalesigOE10 <- data.frame(subset( Male10, threshold==TRUE))



output_folder = paste0("Neurons/Deseq2/")
VennFreq=list(Male10L=MalesigOE10$rownames.sigOE_ordered., Female40L=FemalesigOE40$rownames.sigOE_ordered.)



#png(paste0(output_folder,"Venn Diagram Male10 Female 40 Down.png"), width=7,height=8,pointsize = 14,units="in",res=600)
pdf(paste0(output_folder,"Venn Diagram Male10 Female 40 Down.pdf"), width=8,height=8,pointsize = 12)
ggvenn(VennFreq, 
       fill_color = c("dodgerblue2","hotpink"),
       stroke_size = 0.5, set_name_size = 6, text_size = 5)
dev.off()




ind=match(MalesigOE10$rownames.sigOE_ordered.,FemalesigOE40$rownames.sigOE_ordered.)

Shared=MalesigOE10[ind,]
Shared=drop_na(Shared)

export(Shared,"Neurons/Deseq2/Shared Gene Venn Diagram Male 10Hz vs Light Female 40Hz vs light stressed.csv",row.names=T)



















#----------------------------------------------------Male Only---------------------------------------------#

sortby <- which(Expt_data$Condition =="F 10HzHz Stress " | Expt_data$Condition  =="F 40HzHz Stress " | Expt_data$Condition  =="F LightHz Stress ")
#This  can be changed to any group or condition
Expt_data_sub=Expt_data[sortby,]
genes_sub=genes[,sortby]

GeneName=rownames(genes_sub)
GeneName2=sub("\\|.*", "", GeneName)


rownames(genes_sub)=GeneName2
genes_sub_Z <- genes_sub %>%
  apply(1,"scale") %>%
  t() %>%
  as_tibble()
colnames(genes_sub_Z) = colnames(genes_sub)
genes_sub_Z$geneID = rownames(genes_sub)
genes_sub_Z[is.na(genes_sub_Z)] <- 0   # set any NA values to 0

# Create a matrix of Z-scored subset data

genes_sub_Z_mat <- as.matrix(genes_sub_Z[,1:{ncol(genes_sub_Z)-1}])
rownames(genes_sub_Z_mat) <- str_to_upper(genes_sub_Z$geneID)






Stresspal <- brewer.pal(2,"YlGnBu")

StressColorBar <- case_when(
  Expt_data_sub$Stressed == "Control" ~ Stresspal[1],
  Expt_data_sub$Stressed == "Stress" ~ Stresspal[2]
)

Freqpal<-brewer.pal(4,"Dark2")
freqColorBar <- case_when(
  Expt_data_sub$Frequency == "Light" ~ Freqpal[1],
  Expt_data_sub$Frequency == "10Hz" ~ Freqpal[2],
  Expt_data_sub$Frequency == "40Hz" ~ Freqpal[3]
)


#Sexpal<-brewer.pal(5,"Set1")
SexColorBar <- case_when(
  Expt_data_sub$Sex == "M" ~ "dodgerblue2",
  Expt_data_sub$Sex == "F" ~ "hotpink" 
)





#Clustering Genes 
hrGenes= hclust(dist((genes_sub_Z_mat),method = "euclidean"), method="ward.D2")


#Generating Heatmap 

# Genes Clustered Heatmap
pdf(paste0("Neurons/Neuron Genes Clustered Heatmap Stressed Subjects Male Only.pdf"),height=4,width=5)
#png(paste0("Neurons/Neuron Genes Clustered Heatmap Stressed Subjects Male Only.png"), width=7,height=8,pointsize = 14, units="in",res=600)
heatmap3(genes_sub_Z_mat, ColSideColors =cbind(c(SexColorBar),c(freqColorBar)), ColSideWidth=1, ColSideLabs=NA,  
         col=barColors, breaks=breakBarColors, legendfun=function()showLegend(legend=c(NA),col=c(NA)),
         Rowv=as.dendrogram(hrGenes), 
         Colv=NA,  scale="none",
         margins = c(3,10),
         labRow = NA, labCol=NA,
         cexCol=1.3)
title("Genes Cluster", adj = 0.5, line = 2)
dev.off()




#Cluster Analysis



cutThreshCTD=3.76 #use with the euclidian distance method, average
myclCTD2 = cutree(hrGenes, k=8);
colorBrewerInterp=colorRampPalette(brewer.pal(8,"Spectral"))
mycolhcCTD2 = sample(colorBrewerInterp(length(unique(myclCTD2))))
mycolhcCTD2 = mycolhcCTD2[as.vector(myclCTD2)]
t1CTD2=as.dendrogram(hrGenes)
rowOrder2=order.dendrogram(t1CTD2);
GeneName=rownames(genes_sub_Z_mat)
GeneName2=sub("\\|.*", "", GeneName)

GroupsOut2=tibble(mycolhcCTD2[rowOrder2],rownames(genes_sub_Z_mat)[rowOrder2], str_to_title(GeneName2[rowOrder2]))


colnames(GroupsOut2)=c("ClusterColor","Genes")
GroupsOut2 <- mutate(GroupsOut2, Cluster=case_when(
  ClusterColor == unique(GroupsOut2$ClusterColor)[1] ~ 1,
  ClusterColor == unique(GroupsOut2$ClusterColor)[2] ~ 2,
  ClusterColor == unique(GroupsOut2$ClusterColor)[3] ~ 3,
  ClusterColor == unique(GroupsOut2$ClusterColor)[4] ~ 4,
  ClusterColor == unique(GroupsOut2$ClusterColor)[5] ~ 5,
  ClusterColor == unique(GroupsOut2$ClusterColor)[6] ~ 6,
  ClusterColor == unique(GroupsOut2$ClusterColor)[7] ~ 7,
  ClusterColor == unique(GroupsOut2$ClusterColor)[8] ~ 8
))
saveRDS(GroupsOut2,file="Cluster Analysis Genes Stressed Male Only.RDS")

export(data.frame(GroupsOut2),file="Cluster Analysis Stressed Male Only.csv")
Clusters=readRDS("Cluster Analysis Genes Stressed Male Only.RDS")

png(paste0("Genes Clustered Heatmap Stressed Subjects Male Only side color.png"), width=7,height=8,pointsize = 14, units="in",res=600)
heatmap3(genes_sub_Z_mat, RowSideColors=mycolhcCTD2,ColSideColors =cbind(c(SexColorBar),c(freqColorBar),c(StressColorBar)), ColSideWidth=1, ColSideLabs=NA,  
         col=barColors, breaks=breakBarColors, legendfun=function()showLegend(legend=c(NA),col=c(NA)),
         Rowv=as.dendrogram(hrGenes), 
         Colv=NA,  scale="none",
         margins = c(3,10),
         labRow = NA, labCol=NA,
         cexCol=1.3)
title("Genes Cluster", adj = 0.5, line = 2)
dev.off()




#Cell Type Specific 

geneSets=read_excel("Reference Files/MyGene-Mouse-SharmaZhangUnion.xlsx") %>%
  as.list() %>%
  lapply(function(x) x[!is.na(x)])

genes_subZ <- genes_sub %>%
  apply(1,"scale") %>%
  t() %>%
  as_tibble()
colnames(genes_subZ) = colnames(genes_sub)
genes_subZ$geneID = rownames(genes_sub)
genes_subZ <- drop_na(genes_subZ)   # remove ANY rows with NA values

# Create a matrix of Z-scored subset data
genes_subZ_mat <- as.matrix(genes_subZ[,1:{ncol(genes_subZ)-1}])



rownames(genes_subZ_mat) <- genes_subZ$geneID
geneSetEnrich=gsva(genes_subZ_mat, geneSets, mx.diff=TRUE, kcdf="Gaussian")
geneSetEnrichSortZ=t(apply(geneSetEnrich,1,scale)); #z-score the data


hrGSVA= hclust(dist((geneSetEnrichSortZ),method = "euclidean"), method="ward.D2")



png(paste0("./Neurons/GSVA_CTE_Male Stressed.png"), width=6,height=8,units="in",res=600, pointsize = 14)
#pdf(paste0("./Neurons/GSVA_CTE_Male Stressed.pdf"), width=6,height=8,pointsize = 14, useDingbats = FALSE)
heatmap3(geneSetEnrichSortZ, ColSideColors =cbind(c(SexColorBar),c(freqColorBar)), ColSideWidth=1, ColSideLabs=NA,
         col=barColors, breaks=breakBarColors, legendfun=function()showLegend(legend=c(NA),col=c(NA)),
         Rowv=as.dendrogram(hrGSVA), Colv=NA,  scale="none",
         labRow = rownames(geneSetEnrich), labCol = NA, cexRow=1, margins=c(5,10))
title("GSVA CTE ", adj = 0.5, line = 2)
dev.off()



# Specific Cell type Gene Sets
Neuron_gene_set <- geneSets$Neuron 
Ast_gene_set <- geneSets$Astrocytes 
Mig_gene_set <- geneSets$Microglia 
olig_gene_set<-geneSets$Oligodendrocytes
Endoth_gen_set<-geneSets$Endothelia 



#Neuron
Neuron_gene <- genes_subZ_mat[which(rownames(genes_subZ_mat) %in% Neuron_gene_set),] 
indNeu=match(rownames(Neuron_gene), rownames(genes_subZ_mat))
NeuGenes=genes_subZ_mat[indNeu,]

#Astrocyte

Astro_gene <- genes_subZ_mat[which(rownames(genes_subZ_mat) %in% Ast_gene_set),] 
indAstr=match(rownames(Astro_gene), rownames(genes_subZ_mat))
astroGenes=genes_subZ_mat[indAstr,]

#Microglia
Mig_gene <- genes_subZ_mat[which(rownames(genes_subZ_mat) %in% Mig_gene_set),]
indMig=match(rownames(Mig_gene), rownames(genes_subZ_mat))
MigGenes=genes_subZ_mat[indMig,]


#Oligodendrocytes

Olig_gene <- genes_subZ_mat[which(rownames(genes_subZ_mat) %in% olig_gene_set),] 
indOlig=match(rownames(Olig_gene), rownames(genes_subZ_mat))
OligGenes=genes_subZ_mat[indOlig,]

#Endothelial

Endo_gene <- genes_subZ_mat[which(rownames(genes_subZ_mat) %in% Endoth_gen_set),] 
indEndo=match(rownames(Endo_gene), rownames(genes_subZ_mat))
EndoGenes=genes_subZ_mat[indEndo,]


AllGenes=rownames(genes_subZ_mat)
MicrogliaGenes=rownames(MigGenes)
EnodthelialGenes=rownames(EndoGenes)
OligodenrocyteGenes=rownames(OligGenes)
NeuronalGenes=rownames(NeuGenes)
AstrocyteGenes=rownames(astroGenes)


#Cell type Distribuation within neurons

others=length(AllGenes)-sum (length(AstrocyteGenes),length(NeuronalGenes),length(MicrogliaGenes),length(EnodthelialGenes),length(OligodenrocyteGenes))
value=c(length(AstrocyteGenes),length(NeuronalGenes),length(MicrogliaGenes),length(EnodthelialGenes),length(OligodenrocyteGenes),others)
group=c("Astrocyte","Neuron","Microglia","Endothelial","Oligodendrocyte","Unclassified")
labels= c(length(AstrocyteGenes)/length(AllGenes) *100 ,length(NeuronalGenes)/length(AllGenes) *100 ,length(MicrogliaGenes)/length(AllGenes) *100 ,length(EnodthelialGenes)/length(AllGenes) *100 ,length(OligodenrocyteGenes)/length(AllGenes) *100 ,others/length(AllGenes) *100)
labels=round(labels, digits=2)
MaleStressedpie=data.frame(value,group,labels)

png(paste0("./Neurons/Neuron Male Stressed Cell Type Distribution.png"), width=6,height=8,units="in",res=600, pointsize = 14)
#pdf(paste0("./Neurons/Neuron Male Stressed Cell Type Distribution.pdf"), width=6,height=8,pointsize = 14, useDingbats = FALSE)
ggplot(MaleStressedpie, aes(x = "", y = value, fill = group)) +
  geom_col() +
  geom_label(aes(label = labels),label.size = 0.25,
             position = position_stack(vjust = 0.5),
             show.legend = FALSE) +
  coord_polar(theta = "y") +
  scale_fill_brewer(palette="Set2")+
  theme_void()
dev.off()



#Bar Plot

png(paste0("./Neurons/Neuron Male Stressed Cell Type DistributionBar plot all.png"), width=6,height=8,units="in",res=600, pointsize = 14)
#pdf(paste0("./Neurons/Neuron Male Stressed Cell Type DistributionBar plot all.pdf"), width=6,height=8,pointsize = 14, useDingbats = FALSE)
ggplot(data=MaleStressedpie, aes(x=1, y=labels, fill=group)) +
  geom_bar(stat="identity") + 
  geom_label(aes(label = labels),label.size = 0.25,
             position = position_stack(vjust = 0.5),
             show.legend = FALSE) +
  scale_fill_manual(values=c("#1a476f","#90353b","#55752f", "#e37e00","#ffd200","#d9e6eb"))+
  theme_classic()
dev.off()




#Heatmaps for each Cell type
#Clustering Genes 

# Specific Cell type Gene Sets
Neuron_gene_set <- geneSets$Neuron 
Ast_gene_set <- geneSets$Astrocytes 
Mig_gene_set <- geneSets$Microglia 
olig_gene_set<-geneSets$Oligodendrocytes
Endoth_gen_set<-geneSets$Endothelia 



#Neuron
Neuron_gene <- genes_sub[which(rownames(genes_sub) %in% Neuron_gene_set),] #1384
indNeu=match(rownames(Neuron_gene), rownames(genes_sub))
NeuGenesRaw=genes_sub[indNeu,]

#Astrocyte

Astro_gene <- genes_sub[which(rownames(genes_sub) %in% Ast_gene_set),] #1066
indAstr=match(rownames(Astro_gene), rownames(genes_sub))
astroGenesRaw=genes_sub[indAstr,]

#Microglia
Mig_gene <- genes_sub[which(rownames(genes_sub) %in% Mig_gene_set),] #1125
indMig=match(rownames(Mig_gene), rownames(genes_sub))
MigGenesRaw=genes_sub[indMig,]


#Oligodendrocytes

Olig_gene <- genes_sub[which(rownames(genes_sub) %in% olig_gene_set),] #895
indOlig=match(rownames(Olig_gene), rownames(genes_sub))
OligGenesRaw=genes_sub[indOlig,]

#Endothelial

Endo_gene <- genes_sub[which(rownames(genes_sub) %in% Endoth_gen_set),] #718
indEndo=match(rownames(Endo_gene), rownames(genes_sub))
EndoGenesRaw=genes_sub[indEndo,]



genes_subZ <- EndoGenesRaw %>%
  apply(1,"scale") %>%
  t() %>%
  as_tibble()
colnames(genes_subZ) = colnames(EndoGenesRaw)
genes_subZ$geneID = rownames(EndoGenesRaw)
genes_subZ <- drop_na(genes_subZ)   # remove ANY rows with NA values

# Create a matrix of Z-scored subset data
genes_subZ_mat <- as.matrix(genes_subZ[,1:{ncol(genes_subZ)-1}])



rownames(genes_subZ_mat) <- genes_subZ$geneID

hrGenes= hclust(dist((genes_subZ_mat),method = "euclidean"), method="ward.D2")


#pdf(paste0("./Neurons/CTE_Endothelial Clustered Heatmap Stressed Subjects Male.pdf"), width=6,height=8,pointsize = 14, useDingbats = FALSE)
png(paste0("Neurons/CTE_Endothelial Clustered Heatmap Stressed Subjects Male.png"), width=7,height=8,pointsize = 14, units="in",res=600)
heatmap3(genes_subZ_mat, ColSideColors =cbind(c(SexColorBar),c(freqColorBar),c(StressColorBar)), ColSideWidth=1, ColSideLabs=NA,  
         col=barColors, breaks=breakBarColors, legendfun=function()showLegend(legend=c(NA),col=c(NA)),
         Rowv=as.dendrogram(hrGenes), 
         Colv=NA,  scale="none",
         margins = c(3,10),
         labRow = NA, labCol=nrow(genes_subZ_mat),
         cexCol=1.3)
title("Genes Cluster Endothelial", adj = 0.5, line = 2)
dev.off()





sortby <- which(Expt_data$Condition =="M 10HzHz Stress " | Expt_data$Condition  =="M 40HzHz Stress " | Expt_data$Condition  =="M LightHz Stress ")
#This  can be changed to any group or condition
Expt_data_sub=Expt_data[sortby,]
genes_sub=genes[,sortby]

GeneName=rownames(genes_sub)
GeneName2=sub("\\|.*", "", GeneName)


rownames(genes_sub)=GeneName2
genes_sub_Z <- genes_sub %>%
  apply(1,"scale") %>%
  t() %>%
  as_tibble()
colnames(genes_sub_Z) = colnames(genes_sub)
genes_sub_Z$geneID = rownames(genes_sub)
genes_sub_Z[is.na(genes_sub_Z)] <- 0   # set any NA values to 0

# Create a matrix of Z-scored subset data

genes_subZ_mat <- as.matrix(genes_sub_Z[,1:{ncol(genes_sub_Z)-1}])
rownames(genes_subZ_mat) <- genes_sub_Z$geneID


#Heatmap of All genes with cell type marker
# Specific Cell type Gene Sets
Neuron_gene_set <- geneSets$Neuron 
Ast_gene_set <- geneSets$Astrocytes 
Mig_gene_set <- geneSets$Microglia 
olig_gene_set<-geneSets$Oligodendrocytes
Endoth_gen_set<-geneSets$Endothelia 



#Neuron
Neuron_gene <- genes_subZ_mat[which(rownames(genes_subZ_mat) %in% Neuron_gene_set),]
indNeu=match(rownames(Neuron_gene), rownames(genes_subZ_mat))
NeuGenes=genes_subZ_mat[indNeu,]
NeuGenes=data.frame(NeuGenes)
NeuGenes <- NeuGenes %>%
  mutate(Celltype = paste0("Neurons"))


#Astrocyte

Astro_gene <- genes_subZ_mat[which(rownames(genes_subZ_mat) %in% Ast_gene_set),] 
indAstr=match(rownames(Astro_gene), rownames(genes_subZ_mat))
astroGenes=genes_subZ_mat[indAstr,]

astroGenes=data.frame(astroGenes)
astroGenes <- astroGenes %>%
  mutate(Celltype = paste0("Astrocyte"))


#Microglia
Mig_gene <- genes_subZ_mat[which(rownames(genes_subZ_mat) %in% Mig_gene_set),] 
indMig=match(rownames(Mig_gene), rownames(genes_subZ_mat))
MigGenes=genes_subZ_mat[indMig,]

MigGenes=data.frame(MigGenes)
MigGenes <- MigGenes %>%
  mutate(Celltype = paste0("Microglia"))

#Oligodendrocytes

Olig_gene <- genes_subZ_mat[which(rownames(genes_subZ_mat) %in% olig_gene_set),] 
indOlig=match(rownames(Olig_gene), rownames(genes_subZ_mat))
OligGenes=genes_subZ_mat[indOlig,]
OligGenes=data.frame(OligGenes)
OligGenes <- OligGenes %>%
  mutate(Celltype = paste0("Oligodendrocyte"))
#Endothelial

Endo_gene <- genes_subZ_mat[which(rownames(genes_subZ_mat) %in% Endoth_gen_set),] 
indEndo=match(rownames(Endo_gene), rownames(genes_subZ_mat))
EndoGenes=genes_subZ_mat[indEndo,]

EndoGenes=data.frame(EndoGenes)
EndoGenes <- EndoGenes %>%
  mutate(Celltype = paste0("Endothelial"))


classified<-genes_subZ_mat[which(rownames(genes_subZ_mat) %in% c(Endoth_gen_set,olig_gene_set,Mig_gene_set,Ast_gene_set,Neuron_gene_set)),]

indclassified=match(rownames(classified), rownames(genes_subZ_mat))
unclassified<- genes_subZ_mat[-c(indclassified) ,]
unclassified=data.frame(unclassified)
unclassified <- unclassified %>%
  mutate(Celltype = paste0("Unclassified"))



Data <- unclassified %>%
  bind_rows(NeuGenes,astroGenes,MigGenes,OligGenes,EndoGenes) 


CellTypecolorBar <- case_when(
  Data$Celltype == "Neurons" ~ "#e37e00",
  Data$Celltype == "Microglia" ~ "#55752f",
  Data$Celltype == "Astrocyte" ~ "#1a476f",
  Data$Celltype == "Endothelial" ~ "#90353b",
  Data$Celltype == "Oligodendrocyte" ~ "#ffd200",
  Data$Celltype == "Unclassified" ~ "#d9e6eb"
)

export(data.frame(rownames(Data),Data),file="./Neurons/Male Data with cell type annotation.csv")

hrGenes= hclust(dist((Data[, c(1:9)]),method = "euclidean"), method="ward.D2")

#pdf(paste0("./Neurons/Gene Heatmap Male Stressed CTE Color Bar with cluster.pdf"), width=6,height=8,pointsize = 14, useDingbats = FALSE)
png(paste0("./Neurons/Gene Heatmap Male Stressed CTE Color Bar with cluster.png"), width=7,height=8,pointsize = 14, units="in",res=600)
heatmap3(Data[, c(1:9)], RowSideColors=CellTypecolorBar,ColSideColors =cbind(c(SexColorBar),c(freqColorBar),c(StressColorBar)), ColSideWidth=1, ColSideLabs=NA,  
         col=barColors, breaks=breakBarColors, legendfun=function()showLegend(legend=c(NA),col=c(NA)),
         Rowv=as.dendrogram(hrGenes), 
         Colv=NA,  scale="none",
         margins = c(3,10),
         labRow = NA, labCol=NA,
         cexCol=1.3)
title("Genes Cluster", adj = 0.5, line = 2)
dev.off()










###Neuron  Specific Gene Set (Custome)

sortby <- which(Expt_data$Condition =="M 10HzHz Stress " | Expt_data$Condition  =="M 40HzHz Stress " | Expt_data$Condition  =="M LightHz Stress ")
#This  can be changed to any group or condition
Expt_data_sub=Expt_data[sortby,]
genes_sub=genes[,sortby]


GeneName=rownames(genes_sub)
GeneName2=sub("\\|.*", "", GeneName)



genes_sub_Z <- genes_sub %>%
  apply(1,"scale") %>%
  t() %>%
  as_tibble()
colnames(genes_sub_Z) = colnames(genes_sub)
genes_sub_Z$geneID = rownames(genes_sub)
genes_sub_Z[is.na(genes_sub_Z)] <- 0   # set any NA values to 0
rownames(genes_sub)=GeneName2
# Create a matrix of Z-scored subset data
genes_sub_Z_mat <- as.matrix(genes_sub_Z[,1:{ncol(genes_sub_Z)-1}])
rownames(genes_sub_Z_mat) <- str_to_upper(genes_sub_Z$geneID)




# GSVA using Neuron Annotations
geneSetsMat=read_excel("Reference Files/GeneSets_Neurons_FXN1.xlsx")
geneSets=as.list(geneSetsMat)
geneSets=lapply(geneSets, function(x) x[!is.na(x)])


geneSetEnrich=gsva(genes_sub_Z_mat, geneSets, mx.diff=TRUE, kcdf="Gaussian")
geneSetEnrichSortZ=t(apply(geneSetEnrich,1,scale)); #z-score the data
geneSetEnrichSortZ[is.na(geneSetEnrichSortZ)]=0
colnames(geneSetEnrichSortZ)=colnames(genes_sub_Z_mat)





Freq=Expt_data_sub$Frequency
datagroupVars=data.frame(t(geneSetEnrichSortZ), Freq)


for (ind in 1:16)
{
  
 # png(paste0("Neurons/NeuronCustom GS/Male/GS Male boxplot",ind,".png"),height=4,width=5,units="in",res=1000)
  pdf(paste0("Neurons/NeuronCustom GS/Male/GS Male boxplot",ind,".pdf"),height=4,width=5)
  print({
    ggplot(datagroupVars, aes(x=Freq, y=geneSetEnrichSortZ[ind,], color=Freq)) + 
      geom_boxplot()+
      geom_beeswarm(aes(fill = Freq),size =6,shape=21, stroke = 1, cex = 2,groupOnX = TRUE)+
      theme(panel.background = element_rect(fill = 'white', colour = 'black'),
            text = element_text(size=20),
            axis.text.x = element_text(face = "bold", color = "black",
                                       size = 16),
            axis.text.y = element_text(face = "bold", color = "black",
                                       size = 16),
            plot.title = element_text(hjust = 0.5))+
      scale_color_manual(values=c("#D95F02","#7570B3","#1B9E77"))+
      scale_fill_manual(values=c("#D95F02","#7570B3","#1B9E77"))+
      
      #stat_compare_means(method="wilcox", comparisons=list(c(3,4),c(2,4),c(1,4)),label="p.format",tip.length = 0.03) +
      xlab("")+
      ylab(rownames(geneSetEnrichSortZ)[ind])+
      labs(color="Freq")+
      ylim(c(-2, 2)) +
      ggtitle(colnames(datagroupVars)[ind])
  })
  dev.off()
}



#Cluster Analysis



cutThreshCTD=3.76 #use with the euclidian distance method, average
myclCTD2 = cutree(hrGenes, k=8);
colorBrewerInterp=colorRampPalette(brewer.pal(8,"Spectral"))
mycolhcCTD2 = sample(colorBrewerInterp(length(unique(myclCTD2))))
mycolhcCTD2 = mycolhcCTD2[as.vector(myclCTD2)]
t1CTD2=as.dendrogram(hrGenes)
rowOrder2=order.dendrogram(t1CTD2);
GeneName=rownames(genes_sub)
GeneName2=sub("\\|.*", "", GeneName)

GroupsOut2=tibble(mycolhcCTD2[rowOrder2],rownames(genes_Z_mat)[rowOrder2], GeneName[rowOrder2])

colnames(GroupsOut2)=c("ClusterColor","Genes")
GroupsOut2 <- mutate(GroupsOut2, Cluster=case_when(
  ClusterColor == unique(GroupsOut2$ClusterColor)[1] ~ 1,
  ClusterColor == unique(GroupsOut2$ClusterColor)[2] ~ 2,
  ClusterColor == unique(GroupsOut2$ClusterColor)[3] ~ 3,
  ClusterColor == unique(GroupsOut2$ClusterColor)[4] ~ 4,
  ClusterColor == unique(GroupsOut2$ClusterColor)[5] ~ 5,
  ClusterColor == unique(GroupsOut2$ClusterColor)[6] ~ 6,
  ClusterColor == unique(GroupsOut2$ClusterColor)[7] ~ 7,
  ClusterColor == unique(GroupsOut2$ClusterColor)[8] ~ 8
))
saveRDS(GroupsOut2,file="Cluster Analysis Genes.RDS")

export(data.frame(GroupsOut2),file="Cluster Analysis.csv")
Clusters=readRDS("Cluster Analysis Genes.RDS")

png(paste0("Genes Clustered Heatmap side color.png"), width=7,height=8,pointsize = 14, units="in",res=600)
heatmap3(genes_Z_mat, RowSideColors=mycolhcCTD2,ColSideColors =cbind(c(SexColorBar),c(freqColorBar),c(StressColorBar)), ColSideWidth=1, ColSideLabs=NA,  
         col=barColors, breaks=breakBarColors, legendfun=function()showLegend(legend=c(NA),col=c(NA)),
         Rowv=as.dendrogram(hrGenes), 
         Colv=NA,  scale="none",
         margins = c(3,10),
         labRow = NA, labCol=NA,
         cexCol=1.3)
title("Genes Cluster", adj = 0.5, line = 2)
dev.off()






