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



# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# 
# BiocManager::install("WGCNA")
# 
# BiocManager::install("htmlwidgets")
#if (!require("WGCNA")) install.packages("WGCNA")
# Load contributed packages with pacman
pacman::p_load(pacman,rio,tidyverse,readxl,matrixStats,ggpubr,heatmap3,Rtsne,limma,
               gplots,GSVA,ropls,RColorBrewer,ggsci, EnhancedVolcano,DESeq2,readxl,xlsx,ggvenn,grid,forcats,viridis,org.Mm.eg.db,ggbeeswarm,WGCNA,impute,preprocessCore) # install/load these packages
# pacman: for package handling
#---Set heatmap3 color bar parameters
breakBarColors=c(-200,seq(-1.5, 1.5, 0.01),200) 
barColors = colorpanel(length(breakBarColors)-1, "blue", "white", "red") # create spectrum of blue to red for heatmaps









#######################################################################################################################################
#load your RAW Counts Data 

rawcounts <- read.csv("Data/RawCountMIG.csv", header = TRUE)
geneName=rawcounts[,1]
geneID=paste0(geneName,"|", rawcounts$Gene.ID)
rawcounts=rawcounts[, -c(1,2)]
rownames(rawcounts)=geneID
genes=rawcounts

#filter criteria as recommended on WGCNA FAQ 
# https://horvath.genetics.ucla.edu/html/CoexpressionNetwork/Rpackages/WGCNA/faq.html
#(genes where 50% of samples have count >10)
n=ncol(rawcounts)*0.70 # 70% of the sample
indKeepNew=which(rowSums(rawcounts>50)> n)
datafliter=rawcounts[indKeepNew,]
genes=datafliter
x = rowSums(rawcounts)
numberofFilteredGenes=count(x<50) # run this to see how many genes filterefd

# Reading experimental Group
Expt_data_in <- read_excel("Data/MetaDataMIG.xlsx" , sheet="Sheet1") %>%  # imports the metadata
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

#dds$Frequency <- relevel( dds$Frequency, ref = "Light")
dds <- DESeq(dds)
results(dds)
#
#
dds <- estimateSizeFactors(dds)
normcounts <- counts(dds, normalized=TRUE)
saveRDS(normcounts,"MicrogliaNormCount_70PercentAbove50.RDS")



#Reading Norm Count for further analysis


normcounts = readRDS("MicrogliaNormCount_70PercentAbove50.RDS")
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

output_folder=paste0("Microglia/062524/")

#Generating Heatmap 

# Genes Clustered Heatmap
#png(paste0(output_folder,"All Genes Clustered Heatmap Microglia.png"),height=6,width=7.8,units="in",res=600)
pdf(paste0("./Microglia/062524/All Genes 70 percent above 50 Clustered Heatmap Microglia.pdf"), width=6,height=8,pointsize = 14, useDingbats = FALSE)
heatmap3(genes_Z_mat, ColSideColors =cbind(c(SexColorBar),c(freqColorBar),c(StressColorBar)), ColSideWidth=1, ColSideLabs=NA,  
         col=barColors, breaks=breakBarColors, legendfun=function()showLegend(legend=c(NA),col=c(NA)),
         Rowv=as.dendrogram(hrGenes), 
         Colv=NA,  scale="none",
         margins = c(3,10),
         labRow = NA, labCol=NA,
         cexCol=1.3)
title("Genes Cluster Microglia", adj = 0.5, line = 2)
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
GeneList=c("Csf1r","Aif1","Vim")
titles=c("Csf1r","Aif1","Vim")
Condition=Expt_data$Condition
for(i in 1:length(GeneList)){
  ind <- which(str_detect(rownames(genes),GeneList[i]))
  plotOut <- new_meanSEM(x=Expt_data$Condition,y=genes[ind,],colors =Condition,
                         title = titles[i],y_str = "Expression",
                         axis.text.x.size = 12,axis.text.y.size = 12,title.size = 15,text.size = 15)
  # pdf(paste0("./Figures/Gene Set ",titles[i],".pdf"),height=5,width=6,useDingbats = FALSE,useKerning = FALSE)
  png(paste0("./Microglia/062524/SanityCheck V2 ",titles[i],".png"),height=5,width=6,units="in",res=600)
  print(plotOut)
  dev.off()
}




#Heatmap for only stressed Data


#Female Data

sortby <- which(Expt_data$Condition =="F 10HzHz Stress " | Expt_data$Condition  =="F 40HzHz Stress " | Expt_data$Condition  =="F DarkHz Stress "|Expt_data$Condition  =="F DarkHz Control " )
#This  can be changed to any group or condition
Expt_data_sub=Expt_data[sortby,]
genes_sub=genes[,sortby]



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
pdf(paste0("./Microglia/062524/Microglia Genes Clustered Heatmap  Females Only.pdf"), width=6,height=8,pointsize = 14, useDingbats = FALSE)
#png(paste0("./Microglia/062524/Microglia Genes Clustered Heatmap Females Only.png"), width=7,height=8,pointsize = 14, units="in",res=600)
heatmap3(genes_sub_Z_mat, ColSideColors =cbind(c(SexColorBar),c(StressColorBar),c(freqColorBar)), ColSideWidth=1, ColSideLabs=NA,  
         col=barColors, breaks=breakBarColors, legendfun=function()showLegend(legend=c(NA),col=c(NA)),
         Rowv=as.dendrogram(hrGenes), 
         Colv=NA,  scale="none",
         margins = c(3,10),
         labRow = NA, labCol=NA,
         cexCol=1.3)
title("Genes Cluster", adj = 0.5, line = 2)
dev.off()






# GSVA using Microglia Annotations
geneSetsMat=read_excel("Reference Files/Microglia markers for gene set categorization_Jan 25.xlsx")
geneSets=as.list(geneSetsMat)
geneSets=lapply(geneSets, function(x) x[!is.na(x)])

rownames(genes_sub_Z_mat)=sub("\\|.*", "", rownames(genes_sub_Z_mat))
rownames(genes_sub_Z_mat)=tolower(rownames(genes_sub_Z_mat))

geneSetEnrich=gsva(genes_sub_Z_mat, geneSets, mx.diff=TRUE, kcdf="Gaussian")
geneSetEnrichSortZ=t(apply(geneSetEnrich,1,scale)); #z-score the data
geneSetEnrichSortZ[is.na(geneSetEnrichSortZ)]=0
colnames(geneSetEnrichSortZ)=colnames(genes_sub_Z_mat)


Condition=Expt_data_sub$Condition
Condition <- as.character(Condition)
Condition <- factor(Condition, levels=c("F DarkHz Control ","F DarkHz Stress ","F 10HzHz Stress ","F 40HzHz Stress "))
datagroupVars=data.frame(t(geneSetEnrichSortZ), Condition)


for (ind in 1:6)
{
  
 png(paste0("Microglia/062524/Microglia Custom GS updated/Female Jan 25 ",ind,".png"),height=4,width=5,units="in",res=1000)
# pdf(paste0("Microglia/062524/Microglia Custom GS updated/ Female Jan 25 ",ind,".pdf"),height=4,width=5)
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


#png(paste0("./Microglia/062524/Microglia Custom GS updated/Custom GS Microglia Female All Heatmap Jan25.png"), width=6,height=8,units="in",res=600, pointsize = 14)
pdf(paste0("./Microglia/062524/Microglia Custom GS updated/Custom GS Microglia Female All Heatmap Jan 25.pdf"), width=6,height=8,pointsize = 14, useDingbats = FALSE)
heatmap3(geneSetEnrichSortZ, ColSideColors =cbind(c(SexColorBar),c(freqColorBar),c(StressColorBar)), ColSideWidth=1, ColSideLabs=NA,
         col=barColors, breaks=breakBarColors, legendfun=function()showLegend(legend=c(NA),col=c(NA)),
         Rowv=as.dendrogram(hrGSVA), Colv=NA,  scale="none",
         labRow = rownames(geneSetEnrich), labCol = NA, cexRow=1, margins=c(5,10))
title("Custom GS Microglia", adj = 0.5, line = 2)
dev.off()




#----Specific Microglia Gene Sets heatmaps
ClassicActivation_GeneSet <- geneSetsMat$`Classic Activation`

rownames(genes_sub)=sub("\\|.*", "", rownames(genes_sub))
rownames(genes_sub)=tolower(rownames(genes_sub))
MIG_ClassicActivation <- genes_sub[which(rownames(genes_sub) %in% ClassicActivation_GeneSet),]

MIG_ClassicActivation <- genes_sub[which(ClassicActivation_GeneSet %in% rownames(genes_sub)),]
genes_sub_Z <- MIG_ClassicActivation %>%
  apply(1,"scale") %>%
  t() %>%
  as_tibble()
colnames(genes_sub_Z) = colnames(MIG_ClassicActivation)
genes_sub_Z$geneID = rownames(MIG_ClassicActivation)
genes_sub_Z[is.na(genes_sub_Z)] <- 0   # set any NA values to 0

# Create a matrix of Z-scored subset data
genes_sub_Z_mat <- as.matrix(genes_sub_Z[,1:{ncol(genes_sub_Z)-1}])
rownames(genes_sub_Z_mat) <- str_to_upper(genes_sub_Z$geneID)


png(paste0("./Microglia/062524/Microglia Custom GS updated/Heatmap of MIG genes in classic activation geneSet.png"), width=7,height=8,units="in",res=600, pointsize = 14)
#pdf(paste0("./Microglia/062524/Microglia Custom GS updated/Heatmap of MIG genes in classic activation geneSet.pdf"), width=7,height=8,pointsize = 14,useDingbats = FALSE)
heatmap3(genes_sub_Z_mat, ColSideColors =cbind(c(SexColorBar),c(freqColorBar),c(StressColorBar)), ColSideWidth=1, ColSideLabs=NA,
         col=barColors, breaks=breakBarColors, legendfun=function()showLegend(legend=c(NA),col=c(NA)),
         #Rowv=NA, 
         Colv=NA,  scale="none",
         #labRow = NA, 
         labCol = NA, cexRow=1, margins=c(1,5))
title("Classic Activation Microglia Gene Set", adj = 0.5, line = 2)
dev.off()





#----Specific Microglia Gene Sets heatmaps
HM_GeneSet <- geneSetsMat$'Homeostatic Markers'
rownames(genes_sub)=sub("\\|.*", "", rownames(genes_sub))
rownames(genes_sub)=tolower(rownames(genes_sub))
MIG_HM <- genes_sub[which(rownames(genes_sub) %in% HM_GeneSet),]


genes_sub_Z <- MIG_HM %>%
  apply(1,"scale") %>%
  t() %>%
  as_tibble()
colnames(genes_sub_Z) = colnames(MIG_HM)
genes_sub_Z$geneID = rownames(MIG_HM)
genes_sub_Z[is.na(genes_sub_Z)] <- 0   # set any NA values to 0

# Create a matrix of Z-scored subset data
genes_sub_Z_mat <- as.matrix(genes_sub_Z[,1:{ncol(genes_sub_Z)-1}])
rownames(genes_sub_Z_mat) <- str_to_title(genes_sub_Z$geneID)


#png(paste0("./Microglia/062524/Microglia Custom GS updated/Heatmap of MIG genes in HM geneSet.png"), width=7,height=8,units="in",res=600, pointsize = 14)
pdf(paste0("./Microglia/062524/Microglia Custom GS updated/Heatmap of MIG genes in HM geneSet.pdf"), width=7,height=8,pointsize = 14,useDingbats = FALSE)
heatmap3(genes_sub_Z_mat, ColSideColors =cbind(c(SexColorBar),c(freqColorBar),c(StressColorBar)), ColSideWidth=1, ColSideLabs=NA,
         col=barColors, breaks=breakBarColors, legendfun=function()showLegend(legend=c(NA),col=c(NA)),
         #Rowv=NA, 
         Colv=NA,  scale="none",
         #labRow = NA, 
         labCol = NA, cexRow=1, margins=c(1,5))
title("Homeostatic Markers Microglia Gene Set", adj = 0.5, line = 2)
dev.off()













#Male Data

sortby <- which(Expt_data$Condition =="M 10HzHz Stress " | Expt_data$Condition  =="M 40HzHz Stress " | Expt_data$Condition  =="M DarkHz Stress "|Expt_data$Condition  =="M DarkHz Control " )
#This  can be changed to any group or condition
Expt_data_sub=Expt_data[sortby,]
genes_sub=genes[,sortby]



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
#pdf(paste0("./Microglia/062524/Microglia Genes Clustered Heatmap  Male Only.pdf"), width=6,height=8,pointsize = 14, useDingbats = FALSE)
png(paste0("./Microglia/062524/Microglia Genes Clustered Heatmap Male Only.png"), width=7,height=8,pointsize = 14, units="in",res=600)
heatmap3(genes_sub_Z_mat, ColSideColors =cbind(c(SexColorBar),c(StressColorBar),c(freqColorBar)), ColSideWidth=1, ColSideLabs=NA,  
         col=barColors, breaks=breakBarColors, legendfun=function()showLegend(legend=c(NA),col=c(NA)),
         Rowv=as.dendrogram(hrGenes), 
         Colv=NA,  scale="none",
         margins = c(3,10),
         labRow = NA, labCol=NA,
         cexCol=1.3)
title("Genes Cluster", adj = 0.5, line = 2)
dev.off()



# GSVA using Neuron Annotations
geneSetsMat=read_excel("Reference Files/Microglia markers for gene set categorization_Jan 25.xlsx")
geneSets=as.list(geneSetsMat)
geneSets=lapply(geneSets, function(x) x[!is.na(x)])

rownames(genes_sub_Z_mat)=sub("\\|.*", "", rownames(genes_sub_Z_mat))
rownames(genes_sub_Z_mat)=tolower(rownames(genes_sub_Z_mat))


geneSetEnrich=gsva(genes_sub_Z_mat, geneSets, mx.diff=TRUE, kcdf="Gaussian")
geneSetEnrichSortZ=t(apply(geneSetEnrich,1,scale)); #z-score the data
geneSetEnrichSortZ[is.na(geneSetEnrichSortZ)]=0
colnames(geneSetEnrichSortZ)=colnames(genes_sub_Z_mat)


Condition=Expt_data_sub$Condition
Condition <- as.character(Condition)
Condition <- factor(Condition, levels=c("M DarkHz Control ","M DarkHz Stress ","M 10HzHz Stress ","M 40HzHz Stress "))
datagroupVars=data.frame(t(geneSetEnrichSortZ), Condition)


for (ind in 1:6)
{
  
   png(paste0("Microglia/062524/Microglia Custom GS updated/Male Jan25   ",ind,".png"),height=4,width=5,units="in",res=1000)
  #pdf(paste0("Microglia/062524/Microglia Custom GS updated/Male Jan25  ",ind,".pdf"),height=4,width=5)
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


png(paste0("./Microglia/062524/Microglia Custom GS updated/Custom GS Microglia Male All Heatmap Jan25.png"), width=6,height=8,units="in",res=600, pointsize = 14)
#pdf(paste0("./Microglia/062524/Microglia Custom GS updated/Custom GS Microglia Male All Heatmap Jan25.pdf"), width=6,height=8,pointsize = 14, useDingbats = FALSE)
heatmap3(geneSetEnrichSortZ, ColSideColors =cbind(c(SexColorBar),c(freqColorBar),c(StressColorBar)), ColSideWidth=1, ColSideLabs=NA,
         col=barColors, breaks=breakBarColors, legendfun=function()showLegend(legend=c(NA),col=c(NA)),
         Rowv=as.dendrogram(hrGSVA), Colv=NA,  scale="none",
         labRow = rownames(geneSetEnrich), labCol = NA, cexRow=1, margins=c(5,10))
title("Custom GS Microglia", adj = 0.5, line = 2)
dev.off()



#----Specific Microglia Gene Sets heatmaps
HM_GeneSet <- geneSetsMat$'Homeostatic Markers'
rownames(genes_sub)=sub("\\|.*", "", rownames(genes_sub))
rownames(genes_sub)=tolower(rownames(genes_sub))
MIG_HM <- genes_sub[which(rownames(genes_sub) %in% HM_GeneSet),]


genes_sub_Z <- MIG_HM %>%
  apply(1,"scale") %>%
  t() %>%
  as_tibble()
colnames(genes_sub_Z) = colnames(MIG_HM)
genes_sub_Z$geneID = rownames(MIG_HM)
genes_sub_Z[is.na(genes_sub_Z)] <- 0   # set any NA values to 0

# Create a matrix of Z-scored subset data
genes_sub_Z_mat <- as.matrix(genes_sub_Z[,1:{ncol(genes_sub_Z)-1}])
rownames(genes_sub_Z_mat) <- str_to_title(genes_sub_Z$geneID)


png(paste0("./Microglia/062524/Microglia Custom GS updated/Heatmap of MIG genes in HM geneSet Male.png"), width=7,height=8,units="in",res=600, pointsize = 14)
#pdf(paste0("./Microglia/062524/Microglia Custom GS updated/Heatmap of MIG genes in HM geneSet Male.pdf"), width=7,height=8,pointsize = 14,useDingbats = FALSE)
heatmap3(genes_sub_Z_mat, ColSideColors =cbind(c(SexColorBar),c(freqColorBar),c(StressColorBar)), ColSideWidth=1, ColSideLabs=NA,
         col=barColors, breaks=breakBarColors, legendfun=function()showLegend(legend=c(NA),col=c(NA)),
         #Rowv=NA, 
         Colv=NA,  scale="none",
         #labRow = NA, 
         labCol = NA, cexRow=1, margins=c(1,5))
title("Homeostatic Markers Microglia Gene Set", adj = 0.5, line = 2)
dev.off()


#--------------------------------------------------------------------DEGs--------------------------------------------------------
#Selecting Condition of interest for volcano Plots
#load your RAW Counts Data 



rawcounts <- read.csv("Data/RawCountMIG.csv", header = TRUE)
geneName=rawcounts[,1]
geneID=paste0(geneName,"|", rawcounts$Gene.ID)
rawcounts=rawcounts[, -c(1,2)]
rownames(rawcounts)=geneID
genes=rawcounts

#filter criteria as recommended on WGCNA FAQ 
# https://horvath.genetics.ucla.edu/html/CoexpressionNetwork/Rpackages/WGCNA/faq.html
#(genes where 50% of samples have count >10)
n=ncol(rawcounts)*0.70 # 50% of the sample
indKeepNew=which(rowSums(rawcounts>50)> n)
datafliter=rawcounts[indKeepNew,]
genes=datafliter
x = rowSums(rawcounts)
numberofFilteredGenes=count(x<50) # run this to see how many genes filterefd

# Reading experimental Group
Expt_data_in <- read_excel("Data/MetaDataMIG.xlsx" , sheet="Sheet1") %>%  # imports the metadata
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



sortby <- which(Expt_data$Condition =="F DarkHz Stress " | Expt_data$Condition  == "F DarkHz Control " )#This  can be changed to any group or condition
Expt_data_sub=Expt_data[sortby,]
genes_sub=genes[,sortby]


# dds <- DESeqDataSetFromMatrix(countData = genes_sub,
#                               colData = Expt_data_sub,
#                               design = ~ Frequency)
# 
# dds$Frequency <- relevel( dds$Frequency, ref = "10Hz")
# 
# 



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
export(outMale,"Microglia/062524/Deseq2/Microglia DEGs_Female Dark Stress vs Female Dark Control DESeq2 ThreshholdPassing FC2 V2.csv")

#####Export Data#####
DesqDataOut=as.matrix(data.frame(DesqData))
DseqExport <-data.frame(rownames(DesqDataOut),DesqDataOut)
export(DseqExport,"Microglia/062524/Deseq2/Microglia DEGs_Female Dark Stress vs Female Dark Control DESeq2 All V2.csv")


# Setting up output folder for the figure
output_folder=output_folder = paste0("Microglia/062524/Deseq2/Volcano/")
dir.create(output_folder,recursive = TRUE, showWarnings=FALSE)


keyvals <- ifelse(
  DesqData$log2FoldChange < -1 & DesqData$pvalue < 5*10e-3, '#0070c0',
  ifelse(DesqData$log2FoldChange > 1 & DesqData$pvalue < 5*10e-3, '#FA0000',
         'gray'))
keyvals[is.na(keyvals)] <- 'gray'
names(keyvals)[keyvals == '#FA0000'] <- 'Up'
names(keyvals)[keyvals == 'gray'] <- 'NS'
names(keyvals)[keyvals == '#0070c0'] <- 'Down'


LabNames <- read_excel("Data/Genes to label for updating volcano plots_12Feb2025.xlsx" , sheet="Female Microglia")
GeneList<-LabNames$`Stress No Stim`
GeneListVolcano<-unique(str_to_title(GeneList))
# LabNames$`Microglia markers`=str_to_title(LabNames$`Microglia markers`)
# LabNames$`Microglia markers`=as.character(LabNames$`Microglia markers`)
# GeneListVolcano=unique(LabNames$`Microglia markers`)

#pdf(paste0(output_folder,"Microglia DEGs Female Dark Stress vs Female Dark Control FC2 v2.pdf"), width=7,height=8,pointsize = 14,useDingbats = FALSE)
png(paste0(output_folder,"Microglia DEGs Female Dark Stress vs Female Dark Control FC2 v2.png"), width=7,height=8,pointsize = 14,units="in",res=600)
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




#-----------------10Hz vs Dark Stress---------------


sortby <- which(Expt_data$Condition =="F 10HzHz Stress " | Expt_data$Condition  == "F DarkHz Stress " )#This  can be changed to any group or condition
Expt_data_sub=Expt_data[sortby,]
genes_sub=genes[,sortby]


 dds <- DESeqDataSetFromMatrix(countData = genes_sub,
                               colData = Expt_data_sub,
                               design = ~ Frequency)
 
 dds$Frequency <- relevel( dds$Frequency, ref = "Dark")
 
# 
# dds <- DESeqDataSetFromMatrix(countData = genes_sub,
#                               colData = Expt_data_sub,
#                               design = ~ Stressed)
# 
# dds$Stressed <- relevel( dds$Stressed, ref = "Control")
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
export(outMale,"Microglia/062524/Deseq2/Microglia DEGs_Female 10Hz Stress vs Female Dark Stress DESeq2 ThreshholdPassing FC2 V2.csv")

#####Export Data#####
DesqDataOut=as.matrix(data.frame(DesqData))
DseqExport <-data.frame(rownames(DesqDataOut),DesqDataOut)
export(DseqExport,"Microglia/062524/Deseq2/Microglia DEGs_Female 10Hz Stress vs Female Dark Stress DESeq2 All V2.csv")


# Setting up output folder for the figure
output_folder=output_folder = paste0("Microglia/062524/Deseq2/Volcano/")
dir.create(output_folder,recursive = TRUE, showWarnings=FALSE)


keyvals <- ifelse(
  DesqData$log2FoldChange < -1 & DesqData$pvalue < 5*10e-3, '#0070c0',
  ifelse(DesqData$log2FoldChange > 1 & DesqData$pvalue < 5*10e-3, '#FA0000',
         'gray'))
keyvals[is.na(keyvals)] <- 'gray'
names(keyvals)[keyvals == '#FA0000'] <- 'Up'
names(keyvals)[keyvals == 'gray'] <- 'NS'
names(keyvals)[keyvals == '#0070c0'] <- 'Down'


LabNames <- read_excel("Data/Genes to label for updating volcano plots_12Feb2025.xlsx" , sheet="Female Microglia")
GeneList<-LabNames$'10Hz'
GeneListVolcano<-unique(str_to_title(GeneList))
# LabNames$`Microglia markers`=str_to_title(LabNames$`Microglia markers`)
# LabNames$`Microglia markers`=as.character(LabNames$`Microglia markers`)
# GeneListVolcano=unique(LabNames$`Microglia markers`)

pdf(paste0(output_folder,"Microglia DEGs Female 10Hz Stress vs Female Dark Control FC2 v2.pdf"), width=7,height=8,pointsize = 14,useDingbats = FALSE)
#png(paste0(output_folder,"Microglia DEGs Female 10Hz Stress vs Female Dark Stress FC2 v2.png"), width=7,height=8,pointsize = 14,units="in",res=600)
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



#-----------------40Hz vs Dark Stress---------------


sortby <- which(Expt_data$Condition =="F 40HzHz Stress " | Expt_data$Condition  == "F DarkHz Stress " )#This  can be changed to any group or condition
Expt_data_sub=Expt_data[sortby,]
genes_sub=genes[,sortby]


dds <- DESeqDataSetFromMatrix(countData = genes_sub,
                              colData = Expt_data_sub,
                              design = ~ Frequency)

dds$Frequency <- relevel( dds$Frequency, ref = "Dark")

# 
# dds <- DESeqDataSetFromMatrix(countData = genes_sub,
#                               colData = Expt_data_sub,
#                               design = ~ Stressed)
# 
# dds$Stressed <- relevel( dds$Stressed, ref = "Control")
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
export(outMale,"Microglia/062524/Deseq2/Microglia DEGs_Female 40Hz Stress vs Female Dark Stress DESeq2 ThreshholdPassing FC2 V2.csv")

#####Export Data#####
DesqDataOut=as.matrix(data.frame(DesqData))
DseqExport <-data.frame(rownames(DesqDataOut),DesqDataOut)
export(DseqExport,"Microglia/062524/Deseq2/Microglia DEGs_Female 40Hz Stress vs Female Dark Stress DESeq2 All V2.csv")


# Setting up output folder for the figure
output_folder=output_folder = paste0("Microglia/062524/Deseq2/Volcano/")
dir.create(output_folder,recursive = TRUE, showWarnings=FALSE)


keyvals <- ifelse(
  DesqData$log2FoldChange < -1 & DesqData$pvalue < 5*10e-3, '#0070c0',
  ifelse(DesqData$log2FoldChange > 1 & DesqData$pvalue < 5*10e-3, '#FA0000',
         'gray'))
keyvals[is.na(keyvals)] <- 'gray'
names(keyvals)[keyvals == '#FA0000'] <- 'Up'
names(keyvals)[keyvals == 'gray'] <- 'NS'
names(keyvals)[keyvals == '#0070c0'] <- 'Down'


LabNames <- read_excel("Data/Genes to label for updating volcano plots_12Feb2025.xlsx" , sheet="Female Microglia")
GeneList<-LabNames$'40Hz'
GeneListVolcano<-unique(str_to_title(GeneList))
# LabNames$`Microglia markers`=str_to_title(LabNames$`Microglia markers`)
# LabNames$`Microglia markers`=as.character(LabNames$`Microglia markers`)
# GeneListVolcano=unique(LabNames$`Microglia markers`)

pdf(paste0(output_folder,"Microglia DEGs Female 40Hz Stress vs Female Dark Control FC2 v2.pdf"), width=7,height=8,pointsize = 14,useDingbats = FALSE)
#png(paste0(output_folder,"Microglia DEGs Female 40Hz Stress vs Female Dark Stress FC2 v2.png"), width=7,height=8,pointsize = 14,units="in",res=600)
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





#-----------------Dark Stress vs Dark Control ---------------

sortby <- which(Expt_data$Condition =="M DarkHz Stress " | Expt_data$Condition  == "M DarkHz Control " )#This  can be changed to any group or condition
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
export(outMale,"Microglia/062524/Deseq2/Microglia DEGs_Male Dark Stress vs Male Dark Control DESeq2 ThreshholdPassing FC2 V2.csv")

#####Export Data#####
DesqDataOut=as.matrix(data.frame(DesqData))
DseqExport <-data.frame(rownames(DesqDataOut),DesqDataOut)
export(DseqExport,"Microglia/062524/Deseq2/Microglia DEGs_Male Dark Stress vs Male Dark Control DESeq2 All V2.csv")


# Setting up output folder for the figure
output_folder=output_folder = paste0("Microglia/062524/Deseq2/Volcano/")
dir.create(output_folder,recursive = TRUE, showWarnings=FALSE)


keyvals <- ifelse(
  DesqData$log2FoldChange < -1 & DesqData$pvalue < 5*10e-3, '#0070c0',
  ifelse(DesqData$log2FoldChange > 1 & DesqData$pvalue < 5*10e-3, '#FA0000',
         'gray'))
keyvals[is.na(keyvals)] <- 'gray'
names(keyvals)[keyvals == '#FA0000'] <- 'Up'
names(keyvals)[keyvals == 'gray'] <- 'NS'
names(keyvals)[keyvals == '#0070c0'] <- 'Down'


LabNames <- read_excel("Data/Genes to label for updating volcano plots_12Feb2025.xlsx" , sheet="Male Microglia")
GeneList<-LabNames$`Stress No Stim`
GeneListVolcano<-unique(str_to_title(GeneList))
# LabNames$`Microglia markers`=str_to_title(LabNames$`Microglia markers`)
# LabNames$`Microglia markers`=as.character(LabNames$`Microglia markers`)
# GeneListVolcano=unique(LabNames$`Microglia markers`)

pdf(paste0(output_folder,"Microglia DEGs Male Dark Stress vs Male Dark Control FC2 v2.pdf"), width=7,height=8,pointsize = 14,useDingbats = FALSE)
#png(paste0(output_folder,"Microglia DEGs Male Dark Stress vs Male Dark Control FC2 v2.png"), width=7,height=8,pointsize = 14,units="in",res=600)
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




#----- Male 40Hz vs Dark


sortby <- which(Expt_data$Condition =="M 40HzHz Stress " | Expt_data$Condition  == "M DarkHz Stress " )#This  can be changed to any group or condition
Expt_data_sub=Expt_data[sortby,]
genes_sub=genes[,sortby]


dds <- DESeqDataSetFromMatrix(countData = genes_sub,
                              colData = Expt_data_sub,
                              design = ~ Frequency)

dds$Frequency <- relevel( dds$Frequency, ref = "Dark")

# 
# dds <- DESeqDataSetFromMatrix(countData = genes_sub,
#                               colData = Expt_data_sub,
#                               design = ~ Stressed)
# 
# dds$Stressed <- relevel( dds$Stressed, ref = "Control")
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
export(outMale,"Microglia/062524/Deseq2/Microglia DEGs_Male 40Hz Stress vs Male Dark Stress DESeq2 ThreshholdPassing FC2 V2.csv")

#####Export Data#####
DesqDataOut=as.matrix(data.frame(DesqData))
DseqExport <-data.frame(rownames(DesqDataOut),DesqDataOut)
export(DseqExport,"Microglia/062524/Deseq2/Microglia DEGs_Male 40Hz Stress vs Male Dark Stress DESeq2 All V2.csv")


# Setting up output folder for the figure
output_folder=output_folder = paste0("Microglia/062524/Deseq2/Volcano/")
dir.create(output_folder,recursive = TRUE, showWarnings=FALSE)


keyvals <- ifelse(
  DesqData$log2FoldChange < -1 & DesqData$pvalue < 5*10e-3, '#0070c0',
  ifelse(DesqData$log2FoldChange > 1 & DesqData$pvalue < 5*10e-3, '#FA0000',
         'gray'))
keyvals[is.na(keyvals)] <- 'gray'
names(keyvals)[keyvals == '#FA0000'] <- 'Up'
names(keyvals)[keyvals == 'gray'] <- 'NS'
names(keyvals)[keyvals == '#0070c0'] <- 'Down'


LabNames <- read_excel("Data/Genes to label for updating volcano plots_12Feb2025.xlsx" , sheet="Male Microglia")
GeneList<-LabNames$'40Hz'
GeneListVolcano<-unique(str_to_title(GeneList))
# LabNames$`Microglia markers`=str_to_title(LabNames$`Microglia markers`)
# LabNames$`Microglia markers`=as.character(LabNames$`Microglia markers`)
# GeneListVolcano=unique(LabNames$`Microglia markers`)

#pdf(paste0(output_folder,"Microglia DEGs Male 40Hz Stress vs Male Dark Control FC2 v2.pdf"), width=7,height=8,pointsize = 14,useDingbats = FALSE)
png(paste0(output_folder,"Microglia DEGs Male 40Hz Stress vs Male Dark Stress FC2 v2.png"), width=7,height=8,pointsize = 14,units="in",res=600)
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



#----- Male 10Hz vs Dark


sortby <- which(Expt_data$Condition =="M 10HzHz Stress " | Expt_data$Condition  == "M DarkHz Stress " )#This  can be changed to any group or condition
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
export(outMale,"Microglia/062524/Deseq2/Microglia DEGs_Male 10Hz Stress vs Male Dark Stress DESeq2 ThreshholdPassing FC2 V2.csv")

#####Export Data#####
DesqDataOut=as.matrix(data.frame(DesqData))
DseqExport <-data.frame(rownames(DesqDataOut),DesqDataOut)
export(DseqExport,"Microglia/062524/Deseq2/Microglia DEGs_Male 10Hz Stress vs Male Dark Stress DESeq2 All V2.csv")


# Setting up output folder for the figure
output_folder=output_folder = paste0("Microglia/062524/Deseq2/Volcano/")
dir.create(output_folder,recursive = TRUE, showWarnings=FALSE)


keyvals <- ifelse(
  DesqData$log2FoldChange < -1 & DesqData$pvalue < 5*10e-3, '#0070c0',
  ifelse(DesqData$log2FoldChange > 1 & DesqData$pvalue < 5*10e-3, '#FA0000',
         'gray'))
keyvals[is.na(keyvals)] <- 'gray'
names(keyvals)[keyvals == '#FA0000'] <- 'Up'
names(keyvals)[keyvals == 'gray'] <- 'NS'
names(keyvals)[keyvals == '#0070c0'] <- 'Down'


LabNames <- read_excel("Data/Genes to label for updating volcano plots_12Feb2025.xlsx" , sheet="Male Microglia")
GeneList<-LabNames$'10Hz'
GeneListVolcano<-unique(str_to_title(GeneList))
# LabNames$`Microglia markers`=str_to_title(LabNames$`Microglia markers`)
# LabNames$`Microglia markers`=as.character(LabNames$`Microglia markers`)
# GeneListVolcano=unique(LabNames$`Microglia markers`)

pdf(paste0(output_folder,"Microglia DEGs Male 10Hz Stress vs Male Dark Control FC2 v2.pdf"), width=7,height=8,pointsize = 14,useDingbats = FALSE)
#png(paste0(output_folder,"Microglia DEGs Male 10Hz Stress vs Male Dark Stress FC2 v2.png"), width=7,height=8,pointsize = 14,units="in",res=600)
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




#---- Generating Volcano Plot---#

#pdf(paste0(output_folder,"Microglia DEGs Male 40Hz vs Male 10Hz  Stress FC2.pdf"), width=7,height=8,pointsize = 14,useDingbats = FALSE)
png(paste0(output_folder," Microglia DEGs Male 40Hz vs Male 10Hz Stress FC2.png"), width=7,height=8,pointsize = 14,units="in",res=600)
EnhancedVolcano(DesqData,
                lab = rownames(DesqData),
                selectLab =rownames(DesqData)[which(names(keyvals) %in% c('Up', 'Down'))],
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






# Count Plots for GO terms


Data <- read_excel("./Microglia/062524/Microglial_Biological Processes_New analysis_Selected GO terms 1.xlsx",
                   sheet= "BP_10Hz_M_all_NC v2")
Data=data.frame(Data)

Name=Data$GO.biological.process.complete
pValFDR=Data$Client.Text.Box.Input..FDR.
FE=Data$Client.Text.Box.Input..fold.Enrichment.
Count=Data$Count
pval=Data$Client.Text.Box.Input..raw.P.value.
FELog=log2(FE)
data2=data.frame(Name,log2(FE),Count,pval)

png(paste0("./Microglia/062524/GO Terms/BP_10Hz_M_all_NC v2.png"), width=10,height=8,pointsize = 20, units="in",res=600)
#pdf("./Microglia/062524/GO Terms/BP_10Hz_M_all_NC v2.pdf", width=10,height=10,pointsize = 6)
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

Data_GeneName= read_excel("Microglia/062524/Stress Mice Portal Transcripts_stress signature comparison_June2024.xlsx",sheet="Microglia DEGs v2") 

DEG_M_10Hz_All=Data_GeneName$`10Hz_M DEG all`
DEG_M_40Hz_All=Data_GeneName$`40Hz_M DEG all`
DEG_M_40Hzvs10Hz_All=Data_GeneName$`40v10Hz_M DEG all`
DEG_F_10Hz_All=Data_GeneName$`10Hz_F DEG all`
DEG_F_40Hz_All=Data_GeneName$`40Hz_F DEG all`
DEG_F_40Hzvs10Hz_All=Data_GeneName$`40v10Hz_F DEG all`

Shared_Two_BioStudies=Data_GeneName$`Supplementary Table 5: Prefrontal Cortex: DEGs shared among 2 BioProjects`
Shared_Three_BioStudies=Data_GeneName$`Table 5: Prefrontal Cortex: DEGs shared among 3 BioProjects`
PublishedGenes=Data_GeneName$`Genes previously associated with stress in either humans, animal models or both.`


indKeep1=match(DEG_M_10Hz_All,Shared_Two_BioStudies)
DEG_M_10Hz_AllSharedWithTwo_BioStudies<- DEG_M_10Hz_All[c(indKeep1)]
export(data.frame(DEG_M_10Hz_AllSharedWithTwo_BioStudies),file="Microglia/062524/Venn Diagrams/Shared Microglial Genes_Shared_Two_BioStudies vs Male 10Hz all .csv")

indKeep2=match(DEG_M_40Hz_All,Shared_Two_BioStudies)
DEG_M_40Hz_AllSharedWithTwo_BioStudies<- DEG_M_40Hz_All[c(indKeep2)]
export(data.frame(DEG_M_40Hz_AllSharedWithTwo_BioStudies),file="Microglia/062524/Venn Diagrams/Shared Microglial Genes_Shared_Two_BioStudies vs Male 40Hz all .csv")


indKeep3=match(DEG_M_40Hz_All,Shared_Three_BioStudies)
DEG_M_40Hz_AllShared_Three_BioStudies<- DEG_M_40Hz_All[c(indKeep3)]
export(data.frame(DEG_M_40Hz_AllShared_Three_BioStudies),file="Microglia/062524/Venn Diagrams/Shared Microglial Genes_Shared_Three_BioStudies vs Male 40Hz all .csv")


indKeep4=match(DEG_M_10Hz_All,Shared_Three_BioStudies)
DEG_M_10Hz_AllShared_Three_BioStudies<- DEG_M_10Hz_All[c(indKeep4)]
export(data.frame(DEG_M_10Hz_AllShared_Three_BioStudies),file="Microglia/062524/Venn Diagrams/Shared Microglial Genes_Shared_Three_BioStudies vs Male 10Hz all .csv")


indKeep5=match(DEG_M_10Hz_All,PublishedGenes)
DEG_M_10Hz_AllSharedWithTwo_BioStudies<- DEG_M_10Hz_All[c(indKeep5)]
export(data.frame(DEG_M_10Hz_AllSharedWithTwo_BioStudies),file="Microglia/062524/Venn Diagrams/Shared Microglial Genes_Shared_PublishedGenes vs Male 10Hz all.csv")

indKeep6=match(DEG_M_40Hz_All,PublishedGenes)
DEG_M_40Hz_All_PublishedGenes<- DEG_M_40Hz_All[c(indKeep6)]
export(data.frame(DEG_M_40Hz_All_PublishedGenes),file="Microglia/062524/Venn Diagrams/Shared Microglial Genes_Shared_PublishedGenes vs Male 40Hz all .csv")


indKeep7=match(DEG_F_10Hz_All,Shared_Two_BioStudies)
DEG_F_10Hz_AllSharedWithTwo_BioStudies<- DEG_F_10Hz_All[c(indKeep7)]
export(data.frame(DEG_F_10Hz_AllSharedWithTwo_BioStudies),file="Microglia/062524/Venn Diagrams/Shared Microglial Genes_Shared_Two_BioStudies vs Female 10Hz all .csv")

indKeep8=match(DEG_F_40Hz_All,Shared_Two_BioStudies)
DEG_F_40Hz_AllSharedWithTwo_BioStudies<- DEG_F_40Hz_All[c(indKeep8)]
export(data.frame(DEG_F_40Hz_AllSharedWithTwo_BioStudies),file="Microglia/062524/Venn Diagrams/Shared Microglial Genes_Shared_Two_BioStudies vs Female 40Hz all .csv")


indKeep9=match(DEG_F_40Hz_All,Shared_Three_BioStudies)
DEG_F_40Hz_AllShared_Three_BioStudies<- DEG_F_40Hz_All[c(indKeep9)]
export(data.frame(DEG_F_40Hz_AllShared_Three_BioStudies),file="Microglia/062524/Venn Diagrams/Shared Microglial Genes_Shared_Three_BioStudies vs Female 40Hz all .csv")


indKeep10=match(DEG_F_10Hz_All,Shared_Three_BioStudies)
DEG_F_10Hz_AllShared_Three_BioStudies<- DEG_F_10Hz_All[c(indKeep10)]
export(data.frame(DEG_F_10Hz_AllShared_Three_BioStudies),file="Microglia/062524/Venn Diagrams/Shared Microglial Genes_Shared_Three_BioStudies vs Female 10Hz all .csv")

indKeep11=match(DEG_F_10Hz_All,PublishedGenes)
DEG_F_10Hz_AllSharedWithTwo_BioStudies<- DEG_F_10Hz_All[c(indKeep11)]
export(data.frame(DEG_F_10Hz_AllSharedWithTwo_BioStudies),file="Microglia/062524/Venn Diagrams/Shared Microglial Genes_Shared_PublishedGenes vs Female 10Hz all.csv")

indKeep12=match(DEG_F_40Hz_All,PublishedGenes)
DEG_F_40Hz_All_PublishedGenes<- DEG_F_40Hz_All[c(indKeep12)]
export(data.frame(DEG_F_40Hz_All_PublishedGenes),file="Microglia/062524/Venn Diagrams/Shared Microglial Genes_Shared_PublishedGenes vs Female 40Hz all .csv")


output_folder = paste0("Microglia/062524/Venn Diagrams/")
#VennFreq1=list(Male10L=DEG_M_10Hz_All, Female40L=DEG_F_40Hz_All)
#VennFreq2=list(Male10L=DEG_M_10Hz_All, Male40L= DEG_M_40Hz_All,Female10L=DEG_F_10Hz_All , Female40L=DEG_F_40Hz_All )


Shared_Two_BioStudies=Data_GeneName$`Supplementary Table 5: Prefrontal Cortex: DEGs shared among 2 BioProjects`
Shared_Three_BioStudies=Data_GeneName$`Table 5: Prefrontal Cortex: DEGs shared among 3 BioProjects`
PublishedGenes=Data_GeneName$`Genes previously associated with stress in either humans, animal models or both.`


VennFreq1=list(Male10L=DEG_M_10Hz_All, Male40HzL=DEG_M_40Hz_All,twoStudy= Shared_Two_BioStudies , ThreeStudy=Shared_Three_BioStudies)
png(paste0(output_folder,"Venn Diagram Male 10Hz vs Male 40Hz vs Two_BioStudies vs three_BioStudies.png"), width=7,height=8,pointsize = 14,units="in",res=600)
#pdf(paste0(output_folder,"Venn Diagram Male 10Hz vs Male 40Hz vs Two_BioStudies vs three_BioStudies.pdf"), width=8,height=8,pointsize = 12)
ggvenn(VennFreq1, 
       #fill_color = c("#872EFE","#D98DFF", ),
       stroke_size = 0.5, set_name_size = 6, text_size = 5)
dev.off()

VennFreq2=list(Male10L=DEG_M_10Hz_All,Male40HzL=DEG_M_40Hz_All, StressStudy=PublishedGenes )
png(paste0(output_folder,"Venn Diagram Male 10Hz vs Male 40Hz vs Published Stress.png"), width=7,height=8,pointsize = 14,units="in",res=600)
#pdf(paste0(output_folder,"Venn Diagram Male 10Hz vs Male 40Hz vs Published Stress.pdf"), width=8,height=8,pointsize = 12)
ggvenn(VennFreq2, 
       #fill_color = c("#872EFE","#D98DFF", ),
       stroke_size = 0.5, set_name_size = 6, text_size = 5)
dev.off()

VennFreq3=list(Female10L=DEG_F_10Hz_All,Female40HzL=DEG_F_40Hz_All, twoStudy= Shared_Two_BioStudies , ThreeStudy=Shared_Three_BioStudies)
png(paste0(output_folder,"Venn Diagram Female 10Hz vs Female 40Hz vs Two_BioStudies vs three_BioStudies.png"), width=7,height=8,pointsize = 14,units="in",res=600)
#pdf(paste0(output_folder,"Venn Diagram Female 10Hz vs Female 40Hz vs Two_BioStudies vs three_BioStudies.pdf"), width=8,height=8,pointsize = 12)
ggvenn(VennFreq3, 
       #fill_color = c("#872EFE","#D98DFF", ),
       stroke_size = 0.5, set_name_size = 6, text_size = 5)
dev.off()


VennFreq4=list(Female10L=DEG_F_10Hz_All, Female40HzL=DEG_F_40Hz_All, StressStudy=PublishedGenes)
png(paste0(output_folder,"Venn Diagram Female 10Hz vs Female 40Hz vs Published Stress.png"), width=7,height=8,pointsize = 14,units="in",res=600)
#pdf(paste0(output_folder,"Venn Diagram Female 10Hz vs Female 40Hz vs Published Stress.pdf"), width=8,height=8,pointsize = 12)
ggvenn(VennFreq4, 
       #fill_color = c("#872EFE","#D98DFF", ),
       stroke_size = 0.5, set_name_size = 6, text_size = 5)
dev.off()


















