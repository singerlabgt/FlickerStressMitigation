#-------------------------------------------------------- Sara Bitarafan _ WOOD Lab 2023 --------------------------------------------------------------#
#Generating Volcano Plots for Neuronal Data
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
barColors = colorpanel(length(breakBarColors)-1, "#84206b", "#e55c30", "#f6d746") # create spectrum of blue to red for heatmaps



barColors = viridis_pal(199,-200,200, "inferno") # create spectrum of blue to red for heatmaps

#---Making a legend figure for the colorbar---
pdf(paste0("Colorbar_bluered.pdf"), width=2,height=5,pointsize = 14)
color.bar(barColors,-1.5,1.5,nticks=5)
dev.off()

#--- Behavioral Data Batch Correction---

DataBehavior <- read_excel("./Behavioral Data Batch Correction.xlsx",sheet= "Male Norm to CTRL")
DataBehavior=data.frame(DataBehavior)

Data=data.frame(DataBehavior[, c(4:12)])
batch=as.vector(DataBehavior$Batch)
df=t(Data)

DataNorm=removeBatchEffect(df,batch)
DataNormT=t(DataNorm)
rownames(DataNormT)=rownames(Data)
DatExp=data.frame(DataNormT,DataBehavior$Batch,DataBehavior$Condition,DataBehavior$Sample)
rownames(DatExp)=DataBehavior$Sample

write.xlsx(data.frame(DatExp), file=paste("./Behavior post Batch Adjust/Jan25_sb_Behavioral Data_Male_Norm to CTRL_BatchAjusted.xlsx"))


#------------ Figure 4 _Panel D------------#


# Count Plots for GO terms


Data <- read_excel("./Figure 4_ Panel D_ biological process _Male and Female.xlsx",sheet= "Sheet1")
Data=data.frame(Data)

Name=Data$GO.biological.process.complete
FE=data.frame(Data$Neurons.FE,Data$Astrocyte.FE,Data$Microglia.FE)
pval=data.frame(Data$Neurons.pValue,Data$Astrocyte.pValue,Data$Microglia.pValue)
Group=Data$Group
data2=data.frame(Name,FE,pval,Group)

# #png(paste0("./Stress vs no stress_no stim_ GO terms_Neurons_Male 010724.png"), width=10,height=8,pointsize = 20, units="in",res=600)
# pdf("./Stress vs no stress_no stim_ GO terms_Neurons_Male 010724.pdf", width=10,height=10,pointsize = 6)
# data2 %>% 
#   as.tibble() %>% 
#   ggplot(aes(reorder(Name,Group), Name)) + 
#   geom_count(aes(size = log2(FE), color= -log10(pval)) , position = position_dodge2(width = 0.25, preserve = "single", padding = -0.25)) +
#   scale_colour_gradient(low =  "royalblue",
#                         high = "red")+
#   #ylim(0.1,6)+
#   coord_flip()+
#   labs(y= "", x = "")+
#   scale_size(range=c(5, 10))+
#   theme_bw()
# dev.off()



png(paste0("./Stress vs no stress_no stim_ GO terms_Neuron_Female010724.png"), width=10,height=8,pointsize = 20, units="in",res=600)
#pdf("./Stress vs no stress_no stim_ GO terms_Neuron_Female 010724.pdf", width=10,height=10,pointsize = 6)
data2 %>%
  mutate(Name = fct_reorder(Name, Group)) %>%
  ggplot(aes(x=Name, y=FE[,1])) +
  geom_count(aes(size = FE[,1], color= -log10(pval[,1])) , position = position_dodge2(width = 0.25, preserve = "single", padding = -0.25)) +
  scale_colour_gradient(low =  "royalblue",
                        high = "red",limits=c(0,15))+
  
  coord_flip() +
  xlab("") +
 scale_size(breaks =c(0, 15, 1))+
  theme_bw()
dev.off()










Data <- read_excel("./Figure 4_ Panel D_ biological process _Male and Female.xlsx",sheet= "Sheet1")
Data=data.frame(Data)

Name=Data$GO.biological.process.complete
FE=Data$FE
pval=Data$pValue
Group=Data$Group
data2=data.frame(Name,FE,pval,Group)

# #png(paste0("./Stress vs no stress_no stim_ GO terms_Neurons_Male 010724.png"), width=10,height=8,pointsize = 20, units="in",res=600)
# pdf("./Stress vs no stress_no stim_ GO terms_Neurons_Male 010724.pdf", width=10,height=10,pointsize = 6)
# data2 %>% 
#   as.tibble() %>% 
#   ggplot(aes(reorder(Name,Group), Name)) + 
#   geom_count(aes(size = log2(FE), color= -log10(pval)) , position = position_dodge2(width = 0.25, preserve = "single", padding = -0.25)) +
#   scale_colour_gradient(low =  "royalblue",
#                         high = "red")+
#   #ylim(0.1,6)+
#   coord_flip()+
#   labs(y= "", x = "")+
#   scale_size(range=c(5, 10))+
#   theme_bw()
# dev.off()



#png(paste0("./Stress vs no stress_no stim_ GO terms_Neuron_Female 011025.png"), width=10,height=8,pointsize = 20, units="in",res=600)
pdf("./Stress vs no stress_no stim_ GO terms_Neuron_Female 011025.pdf", width=10,height=10,pointsize = 6)
data2 %>%
  mutate(Name = fct_reorder(Name, Group)) %>%
  ggplot(aes(x=Name, y=FE)) +
  geom_count(aes(size = FE, color= -log10(pval)) , position = position_dodge2(width = 0.25, preserve = "single", padding = -0.25)) +
  scale_colour_gradient(low =  "royalblue",
                        high = "red",limits=c(0,15))+
  
  coord_flip() +
  xlab("") +
  scale_size(c(5, 10))+
  theme_bw()
dev.off()

#------------------ Differential Expression Analysis--------#

#Selecting Condition of interest for volcano Plots

rawcountsNeu <- read.csv("Data/NeuronalCounts_22173R.csv", header = TRUE, row.names = 1)
colnames(rawcountsNeu)=paste("Neuron",1:ncol(rawcountsNeu))


rawcountsMG <- read.csv("Data/RawCountMIG.csv", header = TRUE)
geneNameMG=rawcountsMG[,1]
geneIDMG=paste0(geneNameMG,"|", rawcountsMG$Gene.ID)
rawcountsMG=rawcountsMG[, -c(1,2)]
rownames(rawcountsMG)=geneIDMG
genesMG=rawcountsMG
colnames(rawcountsMG)=paste("Microglia",1:ncol(rawcountsMG))


rawcountsAst <- read.csv("Data/Astrocye Raw Counts V2.csv", header = TRUE)
geneNameAst=rawcountsAst[,1]
geneIDAst=paste0(geneNameAst,"|", rawcountsAst$Gene.ID)
rawcountsAst=rawcountsAst[, -c(1,2)]
rownames(rawcountsAst)=geneIDAst
genesAst=rawcountsAst
colnames(rawcountsAst)=paste("Astrocyte",1:ncol(rawcountsAst))

Datacom=data.frame(rawcountsMG,rawcountsAst)

ind=match(rownames(rawcountsNeu),rownames(Datacom))

Datacom2=Datacom[ind,]
DataCombined=data.frame(Datacom2,rawcountsNeu)

rawcounts=DataCombined


n=ncol(rawcounts)*0.75 # 50% of the sample
indKeepNew=which(rowSums(rawcounts>50)> n)
datafliter=rawcounts[indKeepNew,]
genes=datafliter
x = rowSums(rawcounts)

rownames(genes)=str_to_title(sub("\\|.*", "", rownames(genes)))


Data_Ref= read_excel("All Combined/Validate_cell specific markers.xlsx",sheet="Sheet1") 
GeneNameRef=str_to_title(Data_Ref$Gene)
Groups=Data_Ref$Group
Data_Ref=data.frame(GeneNameRef,Groups)
rownames(Data_Ref)=GeneNameRef
GeneSelected=genes[which(rownames(genes) %in% GeneNameRef ) ,]
Data_Ref2=Data_Ref[which(rownames(GeneSelected) %in% rownames(Data_Ref) ) ,]

Df=merge(GeneSelected,Data_Ref2,by=0)

Df=Df[order(Df$Groups, decreasing = TRUE), ] 

rownames(Df)=Df$Row.names

Df=Df[, c(2:73)]



genes_Z <- Df %>%
  apply(1,"scale") %>%
  t() %>%
  as_tibble()
colnames(genes_Z) = colnames(Df)
genes_Z$geneID = rownames(Df)
genes_Z[is.na(genes_Z)] <- 0   # set any NA values to 0

# Create a matrix of Z-scored subset data
genes_Z_mat <- as.matrix(genes_Z[,1:{ncol(genes_Z)-1}])
rownames(genes_Z_mat) <- str_to_upper(genes_Z$geneID)

hrGenes= hclust(dist((genes_Z_mat),method = "euclidean"), method="ward.D2")




# Genes Clustered Heatmap
pdf(paste0("./All Combined/All Genes all cell type validation Ordered.pdf"), width=8,height=12,pointsize = 14, useDingbats = FALSE)
#png(paste0("All Combined/All Genes all cell type validation.png"), width=7,height=8,pointsize = 14, units="in",res=600)
heatmap3(genes_Z_mat,  
         col=viridis(length(breakBarColors)-1,option="B"), breaks=breakBarColors, legendfun=function()showLegend(legend=c(NA),col=c(NA)),
         Rowv=NA, 
         Colv=NA,  scale="none",
         margins = c(3,10),
         labRow = str_to_title(rownames(genes_Z_mat)), labCol=colnames(genes_Z_mat),
         cexCol=0.2)
title("Genes Cluster", adj = 0.2, line = 2)
dev.off()

