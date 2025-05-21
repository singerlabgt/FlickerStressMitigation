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

rawcounts <- read.csv("Data/Astrocye Raw Counts V2.csv", header = TRUE)
geneName=rawcounts[,1]
geneID=paste0(geneName,"|", rawcounts$Gene.ID)
rawcounts=rawcounts[, -c(1,2)]
rownames(rawcounts)=geneID
genes=rawcounts

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
Expt_data_in <- read_excel("Data/MetaDataAst.xlsx" , sheet="Sheet1") %>%  # imports the metadata
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
saveRDS(normcounts,"Astrocyte_NormCount_NewFilter_75PercentAbove50.RDS")



#Reading Norm Count for further analysis

normcounts = readRDS("Astrocyte_NormCount_NewFilter_75PercentAbove50.RDS")
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


#Generating Heatmap 

# Genes Clustered Heatmap
#pdf(paste0("./Astrocyte/062524/All Genes Clustered Heatmap Astrocyte 75percen Above 50.pdf"), width=6,height=8,pointsize = 14, useDingbats = FALSE)
png(paste0("./Astrocyte/062524/All Genes Clustered Heatmap Astrocyte 75percen Above 50.png"), width=7,height=8,pointsize = 14, units="in",res=600)
heatmap3(genes_Z_mat, ColSideColors =cbind(c(SexColorBar),c(freqColorBar),c(StressColorBar)), ColSideWidth=1, ColSideLabs=NA,  
         col=barColors, breaks=breakBarColors, legendfun=function()showLegend(legend=c(NA),col=c(NA)),
         Rowv=as.dendrogram(hrGenes), 
         Colv=NA,  scale="none",
         margins = c(3,10),
         labRow = NA, labCol=NA,
         cexCol=1.3)
title("Genes Cluster", adj = 0.5, line = 2)
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
   pdf(paste0("./Astrocyte/062524/Gene Set ",titles[i],".pdf"),height=5,width=6,useDingbats = FALSE,useKerning = FALSE)
 # png(paste0("./Astrocyte/062524/SanityCheck",titles[i],".png"),height=5,width=6,units="in",res=600)
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
#pdf(paste0("./Astrocyte/062524/Astrocyte Genes Clustered Heatmap  Females Only.pdf"), width=6,height=8,pointsize = 14, useDingbats = FALSE)
png(paste0("./Astrocyte/062524/Astrocyte Genes Clustered Heatmap Females Only.png"), width=7,height=8,pointsize = 14, units="in",res=600)
heatmap3(genes_sub_Z_mat, ColSideColors =cbind(c(SexColorBar),c(StressColorBar),c(freqColorBar)), ColSideWidth=1, ColSideLabs=NA,  
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
saveRDS(GroupsOut2,file="./Astrocyte/062524/Cluster Analysis Genes Female Only Astrocyte.RDS")

export(data.frame(GroupsOut2),file="./Astrocyte/062524/Cluster Analysis Genes Female Only Astrocyte.csv")
Clusters=readRDS("./Astrocyte/062524/Cluster Analysis Genes Female Only Astrocyte.RDS")

png(paste0("Genes Clustered Heatmap  Female Only side color.png"), width=7,height=8,pointsize = 14, units="in",res=600)
heatmap3(genes_sub_Z_mat, RowSideColors=mycolhcCTD2,ColSideColors =cbind(c(SexColorBar),c(freqColorBar),c(StressColorBar)), ColSideWidth=1, ColSideLabs=NA,  
         col=barColors, breaks=breakBarColors, legendfun=function()showLegend(legend=c(NA),col=c(NA)),
         Rowv=as.dendrogram(hrGenes), 
         Colv=NA,  scale="none",
         margins = c(3,10),
         labRow = NA, labCol=NA,
         cexCol=1.3)
title("Genes Cluster", adj = 0.5, line = 2)
dev.off()






# 
# 
# #Cell Type Specific 
# 
# geneSets=read_excel("Reference Files/MyGene-Mouse-SharmaZhangUnion.xlsx") %>%
#   as.list() %>%
#   lapply(function(x) x[!is.na(x)])
# 
# genes_subZ <- genes_sub %>%
#   apply(1,"scale") %>%
#   t() %>%
#   as_tibble()
# colnames(genes_subZ) = colnames(genes_sub)
# genes_subZ$geneID = rownames(genes_sub)
# genes_subZ <- drop_na(genes_subZ)   # remove ANY rows with NA values
# 
# # Create a matrix of Z-scored subset data
# genes_subZ_mat <- as.matrix(genes_subZ[,1:{ncol(genes_subZ)-1}])
# 
# 
# 
# rownames(genes_subZ_mat) <- genes_subZ$geneID
# geneSetEnrich=gsva(genes_subZ_mat, geneSets, mx.diff=TRUE, kcdf="Gaussian")
# geneSetEnrichSortZ=t(apply(geneSetEnrich,1,scale)); #z-score the data
# 
# 
# hrGSVA= hclust(dist((geneSetEnrichSortZ),method = "euclidean"), method="ward.D2")
# 
# 
# 
# png(paste0("./Astrocyte/GSVA_CTE_Female Stressed.png"), width=6,height=8,units="in",res=600, pointsize = 14)
# #pdf(paste0("./Astrocyte/GSVA_CTE_Female Stressed.pdf"), width=6,height=8,pointsize = 14, useDingbats = FALSE)
# heatmap3(geneSetEnrichSortZ, ColSideColors =cbind(c(SexColorBar),c(freqColorBar)), ColSideWidth=1, ColSideLabs=NA,
#          col=barColors, breaks=breakBarColors, legendfun=function()showLegend(legend=c(NA),col=c(NA)),
#          Rowv=as.dendrogram(hrGSVA), Colv=NA,  scale="none",
#          labRow = rownames(geneSetEnrich), labCol = NA, cexRow=1, margins=c(5,10))
# title("GSVA CTE", adj = 0.5, line = 2)
# dev.off()
# 
# 
# 
# # Specific Cell type Gene Sets
# Neuron_gene_set <- geneSets$Neuron 
# Ast_gene_set <- geneSets$Astrocytes 
# Mig_gene_set <- geneSets$Microglia 
# olig_gene_set<-geneSets$Oligodendrocytes
# Endoth_gen_set<-geneSets$Endothelia 
# 
# rownames(genes)=sub("\\|.*", "", rownames(genes))
# 
# #Neuron
# Neuron_gene <- genes_subZ_mat[which(rownames(genes_subZ_mat) %in% Neuron_gene_set),] #1384
# indNeu=match(rownames(Neuron_gene), rownames(genes_subZ_mat))
# NeuGenes=genes_subZ_mat[indNeu,]
# 
# NeuronExpression=rowSums(genes[which(rownames(genes) %in% Neuron_gene_set),])
# export(data.frame(NeuronExpression),file="./Astrocyte/Sum of normliazed Counts in Neuronal CTE_Astrocyte.csv",row.names=T)
# #Astrocyte
# 
# Astro_gene <- genes_subZ_mat[which(rownames(genes_subZ_mat) %in% Ast_gene_set),] #1066
# indAstr=match(rownames(Astro_gene), rownames(genes_subZ_mat))
# astroGenes=genes_subZ_mat[indAstr,]
# AstrocyteExpression=rowSums(genes[which(rownames(genes) %in% Ast_gene_set),])
# export(data.frame(AstrocyteExpression),file="./Astrocyte/Sum of normliazed Counts in Astrocye  CTE_Astrocyte.csv",row.names=T)
# 
# 
# 
# 
# #Microglia
# Mig_gene <- genes_subZ_mat[which(rownames(genes_subZ_mat) %in% Mig_gene_set),] #1125
# indMig=match(rownames(Mig_gene), rownames(genes_subZ_mat))
# MigGenes=genes_subZ_mat[indMig,]
# MicrogliaExpression=rowSums(genes[which(rownames(genes) %in% Mig_gene_set),])
# export(data.frame(MicrogliaExpression),file="./Astrocyte/Sum of normliazed Counts in Microglia  CTE_Astrocyte.csv",row.names=T)
# 
# 
# 
# #Oligodendrocytes
# 
# Olig_gene <- genes_subZ_mat[which(rownames(genes_subZ_mat) %in% olig_gene_set),] #895
# indOlig=match(rownames(Olig_gene), rownames(genes_subZ_mat))
# OligGenes=genes_subZ_mat[indOlig,]
# 
# OligoExpression=rowSums(genes[which(rownames(genes) %in% olig_gene_set),])
# export(data.frame(OligoExpression),file="./Astrocyte/Sum of normliazed Counts in Oligodendrocye CTE_Astrocyte.csv",row.names=T)
# 
# 
# 
# #Endothelial
# 
# Endo_gene <- genes_subZ_mat[which(rownames(genes_subZ_mat) %in% Endoth_gen_set),] #718
# indEndo=match(rownames(Endo_gene), rownames(genes_subZ_mat))
# EndoGenes=genes_subZ_mat[indEndo,]
# 
# EndoExpression=rowSums(genes[which(rownames(genes) %in% Endoth_gen_set),])
# export(data.frame(EndoExpression),file="./Astrocyte/Sum of normliazed Counts in Endothelial CTE_ CTE_Astrocyte.csv",row.names=T)
# 
# 
# 
# 
# 
# 
# AllGenes=rownames(genes_subZ_mat)
# MicrogliaGenes=rownames(MigGenes)
# EnodthelialGenes=rownames(EndoGenes)
# OligodenrocyteGenes=rownames(OligGenes)
# NeuronalGenes=rownames(NeuGenes)
# AstrocyteGenes=rownames(astroGenes)
# 
# 
# #Cell type Distribuation within neurons
# 
# others=length(AllGenes)-sum (length(AstrocyteGenes),length(NeuronalGenes),length(MicrogliaGenes),length(EnodthelialGenes),length(OligodenrocyteGenes))
# value=c(length(AstrocyteGenes),length(NeuronalGenes),length(MicrogliaGenes),length(EnodthelialGenes),length(OligodenrocyteGenes),others)
# group=c("Astrocyte","Neuron","Microglia","Endothelial","Oligodendrocyte","Unclassified")
# labels= c(length(AstrocyteGenes)/length(AllGenes) *100 ,length(NeuronalGenes)/length(AllGenes) *100 ,length(MicrogliaGenes)/length(AllGenes) *100 ,length(EnodthelialGenes)/length(AllGenes) *100 ,length(OligodenrocyteGenes)/length(AllGenes) *100 ,others/length(AllGenes) *100)
# labels=round(labels, digits=2)
# FemaleStressedpie=data.frame(value,group,labels)
# 
# png(paste0("./Astrocyte/Astrocyte Female Stressed Cell Type Distribution.png"), width=6,height=8,units="in",res=600, pointsize = 14)
# #pdf(paste0("./Astrocyte/Female Stressed Cell Type Distribution.pdf"), width=6,height=8,pointsize = 14, useDingbats = FALSE)
# ggplot(FemaleStressedpie, aes(x = "", y = value, fill = group)) +
#   geom_col() +
#   geom_label(aes(label = labels),label.size = 0.25,
#              position = position_stack(vjust = 0.5),
#              show.legend = FALSE) +
#   coord_polar(theta = "y") +
#   scale_fill_brewer(palette="Set2")+
#   theme_void()
# dev.off()
# 
# 
# 
# #Bar Plot
# 
# #png(paste0("./Astrocyte/Astrocyte Female Stressed Cell Type DistributionBar plot all.png"), width=6,height=8,units="in",res=600, pointsize = 14)
# pdf(paste0("./Astrocyte/Astrocyte Female Stressed Cell Type DistributionBar plot all.pdf"), width=6,height=8,pointsize = 14, useDingbats = FALSE)
# ggplot(data=FemaleStressedpie, aes(x=1, y=labels, fill=group)) +
#   geom_bar(stat="identity") + 
#   geom_label(aes(label = labels),label.size = 0.25,
#              position = position_stack(vjust = 0.5),
#              show.legend = FALSE) +
#   scale_fill_manual(values=c("#1a476f","#90353b","#55752f", "#e37e00","#ffd200","#d9e6eb"))+
#   theme_classic()
# dev.off()
# 
# 
# 
# 
# #Heatmaps for each Cell type
# #Clustering Genes 
# 
# # Specific Cell type Gene Sets
# Neuron_gene_set <- geneSets$Neuron 
# Ast_gene_set <- geneSets$Astrocytes 
# Mig_gene_set <- geneSets$Microglia 
# olig_gene_set<-geneSets$Oligodendrocytes
# Endoth_gen_set<-geneSets$Endothelia 
# 
# 
# 
# #Neuron
# Neuron_gene <- genes_sub[which(rownames(genes_sub) %in% Neuron_gene_set),] #1384
# indNeu=match(rownames(Neuron_gene), rownames(genes_sub))
# NeuGenesRaw=genes_sub[indNeu,]
# 
# #Astrocyte
# 
# Astro_gene <- genes_sub[which(rownames(genes_sub) %in% Ast_gene_set),] #1066
# indAstr=match(rownames(Astro_gene), rownames(genes_sub))
# astroGenesRaw=genes_sub[indAstr,]
# 
# #Microglia
# Mig_gene <- genes_sub[which(rownames(genes_sub) %in% Mig_gene_set),] #1125
# indMig=match(rownames(Mig_gene), rownames(genes_sub))
# MigGenesRaw=genes_sub[indMig,]
# 
# 
# #Oligodendrocytes
# 
# Olig_gene <- genes_sub[which(rownames(genes_sub) %in% olig_gene_set),] #895
# indOlig=match(rownames(Olig_gene), rownames(genes_sub))
# OligGenesRaw=genes_sub[indOlig,]
# 
# #Endothelial
# 
# Endo_gene <- genes_sub[which(rownames(genes_sub) %in% Endoth_gen_set),] #718
# indEndo=match(rownames(Endo_gene), rownames(genes_sub))
# EndoGenesRaw=genes_sub[indEndo,]
# 
# 
# 
# genes_subZ <- NeuGenesRaw %>%
#   apply(1,"scale") %>%
#   t() %>%
#   as_tibble()
# colnames(genes_subZ) = colnames(NeuGenesRaw)
# genes_subZ$geneID = rownames(NeuGenesRaw)
# genes_subZ <- drop_na(genes_subZ)   # remove ANY rows with NA values
# 
# # Create a matrix of Z-scored subset data
# genes_subZ_mat <- as.matrix(genes_subZ[,1:{ncol(genes_subZ)-1}])
# 
# 
# 
# rownames(genes_subZ_mat) <- genes_subZ$geneID
# 
# hrGenes= hclust(dist((genes_subZ_mat),method = "euclidean"), method="ward.D2")
# 
# 
# pdf(paste0("./Astrocyte/CTE_Neuron Clustered Heatmap Stressed Subjects Female.pdf"), width=6,height=8,pointsize = 14, useDingbats = FALSE)
# #png(paste0("Astrocyte/CTE_Neuron Clustered Heatmap Stressed Subjects Female.png"), width=7,height=8,pointsize = 14, units="in",res=600)
# heatmap3(genes_subZ_mat, ColSideColors =cbind(c(SexColorBar),c(freqColorBar),c(StressColorBar)), ColSideWidth=1, ColSideLabs=NA,  
#          col=barColors, breaks=breakBarColors, legendfun=function()showLegend(legend=c(NA),col=c(NA)),
#          Rowv=as.dendrogram(hrGenes), 
#          Colv=NA,  scale="none",
#          margins = c(3,10),
#          labRow = NA, labCol=nrow(genes_subZ_mat),
#          cexCol=1.3)
# title("Genes Cluster Neuron", adj = 0.5, line = 2)
# dev.off()
# 
# 
# 
# 
# #Heatmap of All genes with cell type marker
# 
# 
# #Female Data
# 
# sortby <- which(Expt_data$Condition =="F 10HzHz Stress " | Expt_data$Condition  =="F 40HzHz Stress " | Expt_data$Condition  =="F LightHz Stress " )
# #This  can be changed to any group or condition
# Expt_data_sub=Expt_data[sortby,]
# genes_sub=genes[,sortby]
# 
# 
# GeneName=rownames(genes_sub)
# GeneName2=sub("\\|.*", "", GeneName)
# 
# 
# 
# genes_sub_Z <- genes_sub %>%
#   apply(1,"scale") %>%
#   t() %>%
#   as_tibble()
# colnames(genes_sub_Z) = colnames(genes_sub)
# genes_sub_Z$geneID = rownames(genes_sub)
# genes_sub_Z[is.na(genes_sub_Z)] <- 0   # set any NA values to 0
# rownames(genes_sub)=GeneName2
# # Create a matrix of Z-scored subset data
# genes_subZ_mat <- as.matrix(genes_sub_Z[,1:{ncol(genes_sub_Z)-1}])
# rownames(genes_subZ_mat) <- genes_sub_Z$geneID
# 
# 
# 
# 
# 
# 
# # Specific Cell type Gene Sets
# Neuron_gene_set <- geneSets$Neuron 
# Ast_gene_set <- geneSets$Astrocytes 
# Mig_gene_set <- geneSets$Microglia 
# olig_gene_set<-geneSets$Oligodendrocytes
# Endoth_gen_set<-geneSets$Endothelia 
# 
# 
# 
# #Neuron
# Neuron_gene <- genes_subZ_mat[which(rownames(genes_subZ_mat) %in% Neuron_gene_set),] #1384
# indNeu=match(rownames(Neuron_gene), rownames(genes_subZ_mat))
# NeuGenes=genes_subZ_mat[indNeu,]
# NeuGenes=data.frame(NeuGenes)
# NeuGenes <- NeuGenes %>%
#   mutate(Celltype = paste0("Neurons"))
# 
# 
# #Astrocyte
# 
# Astro_gene <- genes_subZ_mat[which(rownames(genes_subZ_mat) %in% Ast_gene_set),] #1066
# indAstr=match(rownames(Astro_gene), rownames(genes_subZ_mat))
# astroGenes=genes_subZ_mat[indAstr,]
# 
# astroGenes=data.frame(astroGenes)
# astroGenes <- astroGenes %>%
#   mutate(Celltype = paste0("Astrocyte"))
# 
# 
# #Microglia
# Mig_gene <- genes_subZ_mat[which(rownames(genes_subZ_mat) %in% Mig_gene_set),] #1125
# indMig=match(rownames(Mig_gene), rownames(genes_subZ_mat))
# MigGenes=genes_subZ_mat[indMig,]
# 
# MigGenes=data.frame(MigGenes)
# MigGenes <- MigGenes %>%
#   mutate(Celltype = paste0("Microglia"))
# 
# #Oligodendrocytes
# 
# Olig_gene <- genes_subZ_mat[which(rownames(genes_subZ_mat) %in% olig_gene_set),] #895
# indOlig=match(rownames(Olig_gene), rownames(genes_subZ_mat))
# OligGenes=genes_subZ_mat[indOlig,]
# OligGenes=data.frame(OligGenes)
# OligGenes <- OligGenes %>%
#   mutate(Celltype = paste0("Oligodendrocyte"))
# #Endothelial
# 
# Endo_gene <- genes_subZ_mat[which(rownames(genes_subZ_mat) %in% Endoth_gen_set),] #718
# indEndo=match(rownames(Endo_gene), rownames(genes_subZ_mat))
# EndoGenes=genes_subZ_mat[indEndo,]
# 
# EndoGenes=data.frame(EndoGenes)
# EndoGenes <- EndoGenes %>%
#   mutate(Celltype = paste0("Endothelial"))
# 
# 
# classified<-genes_subZ_mat[which(rownames(genes_subZ_mat) %in% c(Endoth_gen_set,olig_gene_set,Mig_gene_set,Ast_gene_set,Neuron_gene_set)),]
# 
# indclassified=match(rownames(classified), rownames(genes_subZ_mat))
# unclassified<- genes_subZ_mat[-c(indclassified) ,]
# unclassified=data.frame(unclassified)
# unclassified <- unclassified %>%
#   mutate(Celltype = paste0("Unclassified"))
# 
# 
# 
# Data <- unclassified %>%
#   bind_rows(NeuGenes,astroGenes,MigGenes,OligGenes,EndoGenes) 
# 
# 
# CellTypecolorBar <- case_when(
#   Data$Celltype == "Neurons" ~ "#e37e00",
#   Data$Celltype == "Microglia" ~ "#55752f",
#   Data$Celltype == "Astrocyte" ~ "#1a476f",
#   Data$Celltype == "Endothelial" ~ "#90353b",
#   Data$Celltype == "Oligodendrocyte" ~ "#ffd200",
#   Data$Celltype == "Unclassified" ~ "#d9e6eb"
# )
# 
# export(data.frame(rownames(Data),Data),file="./Astrocyte/Female CTE Data with cell type annotation v2.csv")
# hrGenes= hclust(dist((Data[, c(1:9)]),method = "euclidean"), method="ward.D2")
# 
# pdf(paste0("./Astrocyte/Gene Heatmap Female Stressed CTE Color Bar no cluster.pdf"), width=6,height=8,pointsize = 14, useDingbats = FALSE)
# #png(paste0("./Astrocyte/Gene Heatmap Female Stressed CTE Color Bar no cluster.png"), width=7,height=8,pointsize = 14, units="in",res=600)
# heatmap3(Data[, c(1:9)], RowSideColors=CellTypecolorBar,ColSideColors =cbind(c(SexColorBar),c(freqColorBar),c(StressColorBar)), ColSideWidth=1, ColSideLabs=NA,  
#          col=barColors, breaks=breakBarColors, legendfun=function()showLegend(legend=c(NA),col=c(NA)),
#          Rowv=NA, 
#         # Rowv=as.dendrogram(hrGenes), 
#          Colv=NA,  scale="none",
#          margins = c(3,10),
#          labRow = NA, labCol=NA,
#          cexCol=1.3)
# title("All genes With CTE", adj = 0.5, line = 2)
# dev.off()




###Astrocyte  Specific Gene Set (Custome)


sortby <- which(Expt_data$Condition =="F 10HzHz Stress " | Expt_data$Condition  =="F 40HzHz Stress " | Expt_data$Condition  =="F DarkHz Stress "|Expt_data$Condition  =="F DarkHz Control " )
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
geneSetsMat=read_excel("Reference Files/GeneSets_Astrocytes_FXN1.xlsx")
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



for (ind in 1:20)
{
  
 # png(paste0("Astrocyte/062524/AstrocyteCustom GS/GS Female boxplot",ind,".png"),height=4,width=5,units="in",res=1000)
  pdf(paste0("Astrocyte/062524/AstrocyteCustom GS/GS Female boxplot",ind,".pdf"),height=4,width=5)
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


#png(paste0("./Astrocyte/062524/AstrocyteCustom GS/Custom GS Neuron Female All Heatmap.png"), width=6,height=8,units="in",res=600, pointsize = 14)
pdf(paste0("./Astrocyte/062524/AstrocyteCustom GS/Custom GS Neuron Female All Heatmap.pdf"), width=6,height=8,pointsize = 14, useDingbats = FALSE)
heatmap3(geneSetEnrichSortZ, ColSideColors =cbind(c(SexColorBar),c(freqColorBar),c(StressColorBar)), ColSideWidth=1, ColSideLabs=NA,
         col=barColors, breaks=breakBarColors, legendfun=function()showLegend(legend=c(NA),col=c(NA)),
         Rowv=as.dendrogram(hrGSVA), Colv=NA,  scale="none",
         labRow = rownames(geneSetEnrich), labCol = NA, cexRow=1, margins=c(5,10))
title("Custom GS Astrocyte", adj = 0.5, line = 2)
dev.off()








 
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


#Heatmap for only stressed Data


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
pdf(paste0("./Astrocyte/062524/Astrocyte Genes Clustered Heatmap  Male Only.pdf"), width=6,height=8,pointsize = 14, useDingbats = FALSE)
#png(paste0("./Astrocyte/062524/Astrocyte Genes Clustered Heatmap Male Only.png"), width=7,height=8,pointsize = 14, units="in",res=600)
heatmap3(genes_sub_Z_mat, ColSideColors =cbind(c(SexColorBar),c(StressColorBar),c(freqColorBar)), ColSideWidth=1, ColSideLabs=NA,  
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
saveRDS(GroupsOut2,file="./Astrocyte/062524/Cluster Analysis Genes Male Only Astrocyte.RDS")

export(data.frame(GroupsOut2),file="./Astrocyte/062524/Cluster Analysis Genes Male Only Astrocyte.csv")
Clusters=readRDS("./Astrocyte/062524/Cluster Analysis Genes Male Only Astrocyte.RDS")

png(paste0("Genes Clustered Heatmap  Male Only side color.png"), width=7,height=8,pointsize = 14, units="in",res=600)
heatmap3(genes_sub_Z_mat, RowSideColors=mycolhcCTD2,ColSideColors =cbind(c(SexColorBar),c(freqColorBar),c(StressColorBar)), ColSideWidth=1, ColSideLabs=NA,  
         col=barColors, breaks=breakBarColors, legendfun=function()showLegend(legend=c(NA),col=c(NA)),
         Rowv=as.dendrogram(hrGenes), 
         Colv=NA,  scale="none",
         margins = c(3,10),
         labRow = NA, labCol=NA,
         cexCol=1.3)
title("Genes Cluster", adj = 0.5, line = 2)
dev.off()







###Astrocyte  Specific Gene Set (Custome)


sortby <- which(Expt_data$Condition =="M 10HzHz Stress " | Expt_data$Condition  =="M 40HzHz Stress " | Expt_data$Condition  =="M DarkHz Stress "|Expt_data$Condition  =="M DarkHz Control " )
#This  can be changed to any group or Condition
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
geneSetsMat=read_excel("Reference Files/GeneSets_Astrocytes_FXN1.xlsx")
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



for (ind in 1:20)
{
  
   png(paste0("Astrocyte/062524/AstrocyteCustom GS/GS Male boxplot",ind,".png"),height=4,width=5,units="in",res=1000)
  #pdf(paste0("Astrocyte/062524/AstrocyteCustom GS/GS Male boxplot",ind,".pdf"),height=4,width=5)
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


png(paste0("./Astrocyte/062524/AstrocyteCustom GS/Custom GS Neuron Male All Heatmap.png"), width=6,height=8,units="in",res=600, pointsize = 14)
#pdf(paste0("./Astrocyte/062524/AstrocyteCustom GS/Custom GS Neuron Male All Heatmap.pdf"), width=6,height=8,pointsize = 14, useDingbats = FALSE)
heatmap3(geneSetEnrichSortZ, ColSideColors =cbind(c(SexColorBar),c(freqColorBar),c(StressColorBar)), ColSideWidth=1, ColSideLabs=NA,
         col=barColors, breaks=breakBarColors, legendfun=function()showLegend(legend=c(NA),col=c(NA)),
         Rowv=as.dendrogram(hrGSVA), Colv=NA,  scale="none",
         labRow = rownames(geneSetEnrich), labCol = NA, cexRow=1, margins=c(5,10))
title("Custom GS Astrocyte", adj = 0.5, line = 2)
dev.off()


#------------------------------------062524----------------------------------#



# Count Plots for GO terms


Data <- read_excel("./Astrocyte/062524/Astrocytic_Biological Processes_New analysis_Selected GO terms.xlsx",sheet= "BP_40v10Hz_F_all_NC v2")
Data=data.frame(Data)

Name=Data$GO.biological.process.complete
pValFDR=Data$Client.Text.Box.Input..FDR.
FE=Data$Client.Text.Box.Input..fold.Enrichment.
Count=Data$Count
pval=Data$Client.Text.Box.Input..raw.P.value.
FELog=log2(FE)
data2=data.frame(Name,log2(FE),Count,pval)

png(paste0("./Astrocyte/062524/GOterms/BP_40v10Hz_F_all_NC v2.png"), width=10,height=8,pointsize = 20, units="in",res=600)
#pdf("./Astrocyte/062524/GOterms/BP_40v10Hz_F_all_NC v2.pdf", width=10,height=10,pointsize = 6)
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


#--------------------------------------------------------------------DEGs--------------------------------------------------------
#Selecting Condition of interest for volcano Plots
#load your RAW Counts Data 


rawcounts <- read.csv("Data/Astrocye Raw Counts V2.csv", header = TRUE)
geneName=rawcounts[,1]
geneID=paste0(geneName,"|", rawcounts$Gene.ID)
rawcounts=rawcounts[, -c(1,2)]
rownames(rawcounts)=geneID
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
Expt_data_in <- read_excel("Data/MetaDataAst.xlsx" , sheet="Sheet1") %>%  # imports the metadata
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



sortby <- which(Expt_data$Condition =="F DarkHz Stress " | Expt_data$Condition  =="F DarkHz Control " )#This  can be changed to any group or condition
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
export(outMale,"Astrocyte/062524/Deseq2/Astrocyte DEGs_Stress Female Stress vs Female Control DESeq2 ThreshholdPassing FC2 v2.csv")

#####Export Data#####
DesqDataOut=as.matrix(data.frame(DesqData))
DseqExport <-data.frame(rownames(DesqDataOut),DesqDataOut)
export(DseqExport,"Astrocyte/062524/Deseq2/Astrocyte DEGs_Stress Female Stress vs Female Control DESeq2 All v2.csv")

keyvals <- ifelse(
  DesqData$log2FoldChange < -1 & DesqData$pvalue < 5*10e-3, '#0070c0',
  ifelse(DesqData$log2FoldChange > 1 & DesqData$pvalue < 5*10e-3, '#FA0000',
         'gray'))
keyvals[is.na(keyvals)] <- 'gray'
names(keyvals)[keyvals == '#FA0000'] <- 'Up'
names(keyvals)[keyvals == 'gray'] <- 'NS'
names(keyvals)[keyvals == '#0070c0'] <- 'Down'

# Setting up output folder for the figure
output_folder=output_folder = paste0("Astrocyte/062524/Deseq2/Volcano/")

#---- Generating Volcano Plot---#

LabNames <- read_excel("Data/Genes to label for updating volcano plots_12Feb2025.xlsx" , sheet="Female Astrocytes")
GeneList<-LabNames$`Stress No Stim`
GeneListVolcano<-unique(str_to_title(GeneList))



pdf(paste0(output_folder,"Astrocyte DEGs  Female Stress vs Female Control FC2 v3 040425.pdf"), width=7,height=8,pointsize = 14,useDingbats = FALSE)
#png(paste0(output_folder,"Astrocyte DEGs  Female Stress vs Female Control FC2 v3 040425.png"), width=7,height=8,pointsize = 14,units="in",res=600)
EnhancedVolcano(DesqData,
                lab = rownames(DesqData),
                selectLab=GeneListVolcano,
              #  selectLab =rownames(DesqData)[which(names(keyvals) %in% c('Up', 'Down'))],
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
                xlim=c(-6,6) , ylim=c(0,8),
                # ylim = c(0, max(-log10(DesqData$pvalue), na.rm = TRUE) + 5),
                border = 'full')
dev.off()


#------------------10Hz vs Dark Female----------------#

sortby <- which(Expt_data$Condition =="F 10HzHz Stress " | Expt_data$Condition  =="F DarkHz Stress " )#This  can be changed to any group or condition
Expt_data_sub=Expt_data[sortby,]
genes_sub=genes[,sortby]

# 
# dds <- DESeqDataSetFromMatrix(countData = genes_sub,
#                               colData = Expt_data_sub,
#                               design = ~ Frequency)
# 
# dds$Frequency <- relevel( dds$Frequency, ref = "10Hz")
# 

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
export(outMale,"Astrocyte/062524/Deseq2/Astrocyte DEGs_Stress Female 10Hz vs Stress Female Dark_DESeq2 ThreshholdPassing FC2 v2.csv")

#####Export Data#####
DesqDataOut=as.matrix(data.frame(DesqData))
DseqExport <-data.frame(rownames(DesqDataOut),DesqDataOut)
export(DseqExport,"Astrocyte/062524/Deseq2/Astrocyte DEGs_Stress Female 10Hz vs Stress Female Dark_DESeq2 All v2.csv")

keyvals <- ifelse(
  DesqData$log2FoldChange < -1 & DesqData$pvalue < 5*10e-3, '#0070c0',
  ifelse(DesqData$log2FoldChange > 1 & DesqData$pvalue < 5*10e-3, '#FA0000',
         'gray'))
keyvals[is.na(keyvals)] <- 'gray'
names(keyvals)[keyvals == '#FA0000'] <- 'Up'
names(keyvals)[keyvals == 'gray'] <- 'NS'
names(keyvals)[keyvals == '#0070c0'] <- 'Down'

# Setting up output folder for the figure
output_folder=output_folder = paste0("Astrocyte/062524/Deseq2/Volcano/")

#---- Generating Volcano Plot---#

LabNames <- read_excel("Data/Genes to label for updating volcano plots_12Feb2025.xlsx" , sheet="Female Astrocytes")
GeneList<-LabNames$'10Hz'
GeneListVolcano<-unique(str_to_title(GeneList))



#pdf(paste0(output_folder,"Astrocyte DEGs Stress Female 10Hz vs Dark_FC2 v2.pdf"), width=7,height=8,pointsize = 14,useDingbats = FALSE)
png(paste0(output_folder,"Astrocyte DEGs  Stress Female 10Hz vs Dark FC2 v2.png"), width=7,height=8,pointsize = 14,units="in",res=600)
EnhancedVolcano(DesqData,
                lab = rownames(DesqData),
                selectLab=GeneListVolcano,
                #  selectLab =rownames(DesqData)[which(names(keyvals) %in% c('Up', 'Down'))],
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






#------------------40Hz vs Dark Female----------------#

sortby <- which(Expt_data$Condition =="F 40HzHz Stress " | Expt_data$Condition  =="F DarkHz Stress " )#This  can be changed to any group or condition
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
export(outMale,"Astrocyte/062524/Deseq2/Astrocyte DEGs_Stress Female 40Hz vs Stress Female Dark_DESeq2 ThreshholdPassing FC2 v2.csv")

#####Export Data#####
DesqDataOut=as.matrix(data.frame(DesqData))
DseqExport <-data.frame(rownames(DesqDataOut),DesqDataOut)
export(DseqExport,"Astrocyte/062524/Deseq2/Astrocyte DEGs_Stress Female 40Hz vs Stress Female Dark_DESeq2 All v2.csv")

keyvals <- ifelse(
  DesqData$log2FoldChange < -1 & DesqData$pvalue < 5*10e-3, '#0070c0',
  ifelse(DesqData$log2FoldChange > 1 & DesqData$pvalue < 5*10e-3, '#FA0000',
         'gray'))
keyvals[is.na(keyvals)] <- 'gray'
names(keyvals)[keyvals == '#FA0000'] <- 'Up'
names(keyvals)[keyvals == 'gray'] <- 'NS'
names(keyvals)[keyvals == '#0070c0'] <- 'Down'

# Setting up output folder for the figure
output_folder=output_folder = paste0("Astrocyte/062524/Deseq2/Volcano/")

#---- Generating Volcano Plot---#

LabNames <- read_excel("Data/Genes to label for updating volcano plots_12Feb2025.xlsx" , sheet="Female Astrocytes")
GeneList<-LabNames$'40Hz'
GeneListVolcano<-unique(str_to_title(GeneList))



pdf(paste0(output_folder,"Astrocyte DEGs Stress Female 40Hz vs Dark_FC2 v2.pdf"), width=7,height=8,pointsize = 14,useDingbats = FALSE)
#png(paste0(output_folder,"Astrocyte DEGs  Stress Female 40Hz vs Dark FC2 v2.png"), width=7,height=8,pointsize = 14,units="in",res=600)
EnhancedVolcano(DesqData,
                lab = rownames(DesqData),
                selectLab=GeneListVolcano,
                #  selectLab =rownames(DesqData)[which(names(keyvals) %in% c('Up', 'Down'))],
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


#-----------------Stress vs Ctrl Male------------------#

sortby <- which(Expt_data$Condition =="M DarkHz Stress " | Expt_data$Condition  =="M DarkHz Control " )#This  can be changed to any group or condition
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
export(outMale,"Astrocyte/062524/Deseq2/Astrocyte DEGs_Stress Male Stress vs Male Control DESeq2 ThreshholdPassing FC2 v2.csv")

#####Export Data#####
DesqDataOut=as.matrix(data.frame(DesqData))
DseqExport <-data.frame(rownames(DesqDataOut),DesqDataOut)
export(DseqExport,"Astrocyte/062524/Deseq2/Astrocyte DEGs_Stress Male Stress vs Male Control DESeq2 All v2.csv")

keyvals <- ifelse(
  DesqData$log2FoldChange < -1 & DesqData$pvalue < 5*10e-3, '#0070c0',
  ifelse(DesqData$log2FoldChange > 1 & DesqData$pvalue < 5*10e-3, '#FA0000',
         'gray'))
keyvals[is.na(keyvals)] <- 'gray'
names(keyvals)[keyvals == '#FA0000'] <- 'Up'
names(keyvals)[keyvals == 'gray'] <- 'NS'
names(keyvals)[keyvals == '#0070c0'] <- 'Down'

# Setting up output folder for the figure
output_folder=output_folder = paste0("Astrocyte/062524/Deseq2/Volcano/")

#---- Generating Volcano Plot---#

LabNames <- read_excel("Data/Genes to label for updating volcano plots_12Feb2025.xlsx" , sheet="Male Astrocytes")
GeneList<-LabNames$`Stress No Stim`
GeneListVolcano<-unique(str_to_title(GeneList))



#pdf(paste0(output_folder,"Astrocyte DEGs  Male Stress vs Male Control FC2 v2.pdf"), width=7,height=8,pointsize = 14,useDingbats = FALSE)
png(paste0(output_folder,"Astrocyte DEGs  Male Stress vs Male Control FC2 v2.png"), width=7,height=8,pointsize = 14,units="in",res=600)
EnhancedVolcano(DesqData,
                lab = rownames(DesqData),
                selectLab=GeneListVolcano,
                #  selectLab =rownames(DesqData)[which(names(keyvals) %in% c('Up', 'Down'))],
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



#------------------10Hz vs Dark Male----------------#

sortby <- which(Expt_data$Condition =="M 10HzHz Stress " | Expt_data$Condition  =="M DarkHz Stress " )#This  can be changed to any group or condition
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
export(outMale,"Astrocyte/062524/Deseq2/Astrocyte DEGs_Stress Male 10Hz vs Stress Male Dark_DESeq2 ThreshholdPassing FC2 v2.csv")

#####Export Data#####
DesqDataOut=as.matrix(data.frame(DesqData))
DseqExport <-data.frame(rownames(DesqDataOut),DesqDataOut)
export(DseqExport,"Astrocyte/062524/Deseq2/Astrocyte DEGs_Stress Male 10Hz vs Stress Male Dark_DESeq2 All v2.csv")

keyvals <- ifelse(
  DesqData$log2FoldChange < -1 & DesqData$pvalue < 5*10e-3, '#0070c0',
  ifelse(DesqData$log2FoldChange > 1 & DesqData$pvalue < 5*10e-3, '#FA0000',
         'gray'))
keyvals[is.na(keyvals)] <- 'gray'
names(keyvals)[keyvals == '#FA0000'] <- 'Up'
names(keyvals)[keyvals == 'gray'] <- 'NS'
names(keyvals)[keyvals == '#0070c0'] <- 'Down'

# Setting up output folder for the figure
output_folder=output_folder = paste0("Astrocyte/062524/Deseq2/Volcano/")

#---- Generating Volcano Plot---#

LabNames <- read_excel("Data/Genes to label for updating volcano plots_12Feb2025.xlsx" , sheet="Male Astrocytes")
GeneList<-LabNames$'10Hz'
GeneListVolcano<-unique(str_to_title(GeneList))



pdf(paste0(output_folder,"Astrocyte DEGs Stress Male 10Hz vs Dark_FC2 v2.pdf"), width=7,height=8,pointsize = 14,useDingbats = FALSE)
#png(paste0(output_folder,"Astrocyte DEGs  Stress Male 10Hz vs Dark FC2 v2.png"), width=7,height=8,pointsize = 14,units="in",res=600)
EnhancedVolcano(DesqData,
                lab = rownames(DesqData),
                selectLab=GeneListVolcano,
                #  selectLab =rownames(DesqData)[which(names(keyvals) %in% c('Up', 'Down'))],
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






#------------------40Hz vs Dark Female----------------#

sortby <- which(Expt_data$Condition =="M 40HzHz Stress " | Expt_data$Condition  =="M DarkHz Stress " )#This  can be changed to any group or condition
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
export(outMale,"Astrocyte/062524/Deseq2/Astrocyte DEGs_Stress Male 40Hz vs Stress Male Dark_DESeq2 ThreshholdPassing FC2 v2.csv")

#####Export Data#####
DesqDataOut=as.matrix(data.frame(DesqData))
DseqExport <-data.frame(rownames(DesqDataOut),DesqDataOut)
export(DseqExport,"Astrocyte/062524/Deseq2/Astrocyte DEGs_Stress Male 40Hz vs Stress Male Dark_DESeq2 All v2.csv")

keyvals <- ifelse(
  DesqData$log2FoldChange < -1 & DesqData$pvalue < 5*10e-3, '#0070c0',
  ifelse(DesqData$log2FoldChange > 1 & DesqData$pvalue < 5*10e-3, '#FA0000',
         'gray'))
keyvals[is.na(keyvals)] <- 'gray'
names(keyvals)[keyvals == '#FA0000'] <- 'Up'
names(keyvals)[keyvals == 'gray'] <- 'NS'
names(keyvals)[keyvals == '#0070c0'] <- 'Down'

# Setting up output folder for the figure
output_folder=output_folder = paste0("Astrocyte/062524/Deseq2/Volcano/")

#---- Generating Volcano Plot---#

LabNames <- read_excel("Data/Genes to label for updating volcano plots_12Feb2025.xlsx" , sheet="Male Astrocytes")
GeneList<-LabNames$'40Hz'
GeneListVolcano<-unique(str_to_title(GeneList))



#pdf(paste0(output_folder,"Astrocyte DEGs Stress Male 40Hz vs Dark_FC2 v2.pdf"), width=7,height=8,pointsize = 14,useDingbats = FALSE)
png(paste0(output_folder,"Astrocyte DEGs  Stress Male 40Hz vs Dark FC2 v2.png"), width=7,height=8,pointsize = 14,units="in",res=600)
EnhancedVolcano(DesqData,
                lab = rownames(DesqData),
                selectLab=GeneListVolcano,
                #  selectLab =rownames(DesqData)[which(names(keyvals) %in% c('Up', 'Down'))],
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


















#-----------------Venn Diagram of Neuronal DEGs vs publihsed study markers-------------------#

Data_GeneName= read_excel("Astrocyte/062524/Stress Mice Portal Transcripts_stress signature comparison_June2024.xlsx",sheet="Astrocyte DEGs v2") 

DEG_M_10Hz_All=Data_GeneName$`10Hz_M DEG all`
DEG_M_40Hz_All=Data_GeneName$`40Hz_M DEG all`
DEG_M_40Hzvs10Hz_All=Data_GeneName$`40v10Hz_M DEG all`
DEG_F_10Hz_All=Data_GeneName$`10Hz_F DEG all`
DEG_F_40Hz_All=Data_GeneName$`40Hz_F DEG all`
DEG_F_40Hzvs10Hz_All=Data_GeneName$`40v10Hz_F DEG all`

Shared_Two_BioStudies=Data_GeneName$`Supplementary Table 5: Prefrontal Cortex: DEGs shared among 2 BioProjects`
Shared_Three_BioStudies=Data_GeneName$`Table 5 (see above): Prefrontal Cortex: DEGs shared among 3 BioProjects`
PublishedGenes=Data_GeneName$`Genes previously associated with stress in either humans, animal models or both.`


indKeep1=match(DEG_M_10Hz_All,Shared_Two_BioStudies)
DEG_M_10Hz_AllSharedWithTwo_BioStudies<- DEG_M_10Hz_All[c(indKeep1)]
export(data.frame(DEG_M_10Hz_AllSharedWithTwo_BioStudies),file="Astrocyte/062524/Venn Diagrams/Shared Astrocyte Genes_Shared_Two_BioStudies vs Male 10Hz all .csv")

indKeep2=match(DEG_M_40Hz_All,Shared_Two_BioStudies)
DEG_M_40Hz_AllSharedWithTwo_BioStudies<- DEG_M_40Hz_All[c(indKeep2)]
export(data.frame(DEG_M_40Hz_AllSharedWithTwo_BioStudies),file="Astrocyte/062524/Venn Diagrams/Shared Astrocyte Genes_Shared_Two_BioStudies vs Male 40Hz all .csv")


indKeep3=match(DEG_M_40Hz_All,Shared_Three_BioStudies)
DEG_M_40Hz_AllShared_Three_BioStudies<- DEG_M_40Hz_All[c(indKeep3)]
export(data.frame(DEG_M_40Hz_AllShared_Three_BioStudies),file="Astrocyte/062524/Venn Diagrams/Shared Astrocyte Genes_Shared_Three_BioStudies vs Male 40Hz all .csv")


indKeep4=match(DEG_M_10Hz_All,Shared_Three_BioStudies)
DEG_M_10Hz_AllShared_Three_BioStudies<- DEG_M_10Hz_All[c(indKeep4)]
export(data.frame(DEG_M_10Hz_AllShared_Three_BioStudies),file="Astrocyte/062524/Venn Diagrams/Shared Astrocyte Genes_Shared_Three_BioStudies vs Male 10Hz all .csv")


indKeep5=match(DEG_M_10Hz_All,PublishedGenes)
DEG_M_10Hz_AllSharedWithTwo_BioStudies<- DEG_M_10Hz_All[c(indKeep5)]
export(data.frame(DEG_M_10Hz_AllSharedWithTwo_BioStudies),file="Astrocyte/062524/Venn Diagrams/Shared Astrocyte Genes_Shared_PublishedGenes vs Male 10Hz all.csv")

indKeep6=match(DEG_M_40Hz_All,PublishedGenes)
DEG_M_40Hz_All_PublishedGenes<- DEG_M_40Hz_All[c(indKeep6)]
export(data.frame(DEG_M_40Hz_All_PublishedGenes),file="Astrocyte/062524/Venn Diagrams/Shared Astrocyte Genes_Shared_PublishedGenes vs Male 40Hz all .csv")


indKeep7=match(DEG_F_10Hz_All,Shared_Two_BioStudies)
DEG_F_10Hz_AllSharedWithTwo_BioStudies<- DEG_F_10Hz_All[c(indKeep7)]
export(data.frame(DEG_F_10Hz_AllSharedWithTwo_BioStudies),file="Astrocyte/062524/Venn Diagrams/Shared Astrocyte Genes_Shared_Two_BioStudies vs Female 10Hz all .csv")

indKeep8=match(DEG_F_40Hz_All,Shared_Two_BioStudies)
DEG_F_40Hz_AllSharedWithTwo_BioStudies<- DEG_F_40Hz_All[c(indKeep8)]
export(data.frame(DEG_F_40Hz_AllSharedWithTwo_BioStudies),file="Astrocyte/062524/Venn Diagrams/Shared Astrocyte Genes_Shared_Two_BioStudies vs Female 40Hz all .csv")


indKeep9=match(DEG_F_40Hz_All,Shared_Three_BioStudies)
DEG_F_40Hz_AllShared_Three_BioStudies<- DEG_F_40Hz_All[c(indKeep9)]
export(data.frame(DEG_F_40Hz_AllShared_Three_BioStudies),file="Astrocyte/062524/Venn Diagrams/Shared Astrocyte Genes_Shared_Three_BioStudies vs Female 40Hz all .csv")


indKeep10=match(DEG_F_10Hz_All,Shared_Three_BioStudies)
DEG_F_10Hz_AllShared_Three_BioStudies<- DEG_F_10Hz_All[c(indKeep10)]
export(data.frame(DEG_F_10Hz_AllShared_Three_BioStudies),file="Astrocyte/062524/Venn Diagrams/Shared Astrocyte Genes_Shared_Three_BioStudies vs Female 10Hz all .csv")

indKeep11=match(DEG_F_10Hz_All,PublishedGenes)
DEG_F_10Hz_AllSharedWithTwo_BioStudies<- DEG_F_10Hz_All[c(indKeep11)]
export(data.frame(DEG_F_10Hz_AllSharedWithTwo_BioStudies),file="Astrocyte/062524/Venn Diagrams/Shared Astrocyte Genes_Shared_PublishedGenes vs Female 10Hz all.csv")

indKeep12=match(DEG_F_40Hz_All,PublishedGenes)
DEG_F_40Hz_All_PublishedGenes<- DEG_F_40Hz_All[c(indKeep12)]
export(data.frame(DEG_F_40Hz_All_PublishedGenes),file="Astrocyte/062524/Venn Diagrams/Shared Astrocyte Genes_Shared_PublishedGenes vs Female 40Hz all .csv")


output_folder = paste0("Astrocyte/062524/Venn Diagrams/")
#VennFreq1=list(Male10L=DEG_M_10Hz_All, Female40L=DEG_F_40Hz_All)
#VennFreq2=list(Male10L=DEG_M_10Hz_All, Male40L= DEG_M_40Hz_All,Female10L=DEG_F_10Hz_All , Female40L=DEG_F_40Hz_All )



VennFreq1=list(Male10L=DEG_M_10Hz_All, Male40HzL=DEG_M_40Hz_All,twoStudy= Shared_Two_BioStudies , ThreeStudy=Shared_Three_BioStudies)
#png(paste0(output_folder,"Venn Diagram Male 10Hz vs Male 40Hz vs Two_BioStudies vs three_BioStudies.png"), width=7,height=8,pointsize = 14,units="in",res=600)
pdf(paste0(output_folder,"Venn Diagram Male 10Hz vs Male 40Hz vs Two_BioStudies vs three_BioStudies.pdf"), width=8,height=8,pointsize = 12)
ggvenn(VennFreq1, 
       #fill_color = c("#872EFE","#D98DFF", ),
       stroke_size = 0.5, set_name_size = 6, text_size = 5)
dev.off()

VennFreq2=list(Male10L=DEG_M_10Hz_All,Male40HzL=DEG_M_40Hz_All, StressStudy=PublishedGenes )
#png(paste0(output_folder,"Venn Diagram Male 10Hz vs Male 40Hz vs Published Stress.png"), width=7,height=8,pointsize = 14,units="in",res=600)
pdf(paste0(output_folder,"Venn Diagram Male 10Hz vs Male 40Hz vs Published Stress.pdf"), width=8,height=8,pointsize = 12)
ggvenn(VennFreq2, 
       #fill_color = c("#872EFE","#D98DFF", ),
       stroke_size = 0.5, set_name_size = 6, text_size = 5)
dev.off()

VennFreq3=list(Female10L=DEG_F_10Hz_All,Female40HzL=DEG_F_40Hz_All, twoStudy= Shared_Two_BioStudies , ThreeStudy=Shared_Three_BioStudies)
#png(paste0(output_folder,"Venn Diagram Female 10Hz vs Female 40Hz vs Two_BioStudies vs three_BioStudies.png"), width=7,height=8,pointsize = 14,units="in",res=600)
pdf(paste0(output_folder,"Venn Diagram Female 10Hz vs Female 40Hz vs Two_BioStudies vs three_BioStudies.pdf"), width=8,height=8,pointsize = 12)
ggvenn(VennFreq3, 
       #fill_color = c("#872EFE","#D98DFF", ),
       stroke_size = 0.5, set_name_size = 6, text_size = 5)
dev.off()


VennFreq4=list(Female10L=DEG_F_10Hz_All, Female40HzL=DEG_F_40Hz_All, StressStudy=PublishedGenes)
#png(paste0(output_folder,"Venn Diagram Female 10Hz vs Female 40Hz vs Published Stress.png"), width=7,height=8,pointsize = 14,units="in",res=600)
pdf(paste0(output_folder,"Venn Diagram Female 10Hz vs Female 40Hz vs Published Stress.pdf"), width=8,height=8,pointsize = 12)
ggvenn(VennFreq4, 
       #fill_color = c("#872EFE","#D98DFF", ),
       stroke_size = 0.5, set_name_size = 6, text_size = 5)
dev.off()




#-----------------------------------------------------------------------------------------------------------











































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



png(paste0("./Astrocyte/GSVA_CTE_male Stressed.png"), width=6,height=8,units="in",res=600, pointsize = 14)
#pdf(paste0("./Astrocyte/GSVA_CTE_male Stressed.pdf"), width=6,height=8,pointsize = 14, useDingbats = FALSE)
heatmap3(geneSetEnrichSortZ, ColSideColors =cbind(c(SexColorBar),c(freqColorBar)), ColSideWidth=1, ColSideLabs=NA,
         col=barColors, breaks=breakBarColors, legendfun=function()showLegend(legend=c(NA),col=c(NA)),
         Rowv=as.dendrogram(hrGSVA), Colv=NA,  scale="none",
         labRow = rownames(geneSetEnrich), labCol = NA, cexRow=1, margins=c(5,10))
title("GSVA CTE", adj = 0.5, line = 2)
dev.off()



# Specific Cell type Gene Sets
Neuron_gene_set <- geneSets$Neuron 
Ast_gene_set <- geneSets$Astrocytes 
Mig_gene_set <- geneSets$Microglia 
olig_gene_set<-geneSets$Oligodendrocytes
Endoth_gen_set<-geneSets$Endothelia 

rownames(genes)=sub("\\|.*", "", rownames(genes))

#Neuron
Neuron_gene <- genes_subZ_mat[which(rownames(genes_subZ_mat) %in% Neuron_gene_set),] #1384
indNeu=match(rownames(Neuron_gene), rownames(genes_subZ_mat))
NeuGenes=genes_subZ_mat[indNeu,]

NeuronExpression=rowSums(genes[which(rownames(genes) %in% Neuron_gene_set),])
export(data.frame(NeuronExpression),file="./Astrocyte/Sum of normliazed Counts in Neuronal CTE_Astrocyte Male.csv",row.names=T)
#Astrocyte

Astro_gene <- genes_subZ_mat[which(rownames(genes_subZ_mat) %in% Ast_gene_set),] #1066
indAstr=match(rownames(Astro_gene), rownames(genes_subZ_mat))
astroGenes=genes_subZ_mat[indAstr,]
AstrocyteExpression=rowSums(genes[which(rownames(genes) %in% Ast_gene_set),])
export(data.frame(AstrocyteExpression),file="./Astrocyte/Sum of normliazed Counts in Astrocye  CTE_Astrocyte Male.csv",row.names=T)




#Microglia
Mig_gene <- genes_subZ_mat[which(rownames(genes_subZ_mat) %in% Mig_gene_set),] #1125
indMig=match(rownames(Mig_gene), rownames(genes_subZ_mat))
MigGenes=genes_subZ_mat[indMig,]
MicrogliaExpression=rowSums(genes[which(rownames(genes) %in% Mig_gene_set),])
export(data.frame(MicrogliaExpression),file="./Astrocyte/Sum of normliazed Counts in Microglia  CTE_Astrocyte Male.csv",row.names=T)



#Oligodendrocytes

Olig_gene <- genes_subZ_mat[which(rownames(genes_subZ_mat) %in% olig_gene_set),] #895
indOlig=match(rownames(Olig_gene), rownames(genes_subZ_mat))
OligGenes=genes_subZ_mat[indOlig,]

OligoExpression=rowSums(genes[which(rownames(genes) %in% olig_gene_set),])
export(data.frame(OligoExpression),file="./Astrocyte/Sum of normliazed Counts in Oligodendrocye CTE_Astrocyte Male.csv",row.names=T)



#Endothelial

Endo_gene <- genes_subZ_mat[which(rownames(genes_subZ_mat) %in% Endoth_gen_set),] #718
indEndo=match(rownames(Endo_gene), rownames(genes_subZ_mat))
EndoGenes=genes_subZ_mat[indEndo,]

EndoExpression=rowSums(genes[which(rownames(genes) %in% Endoth_gen_set),])
export(data.frame(EndoExpression),file="./Astrocyte/Sum of normliazed Counts in Endothelial CTE_ CTE_Astrocyte Male.csv",row.names=T)






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
maleStressedpie=data.frame(value,group,labels)

#png(paste0("./Astrocyte/Astrocyte male Stressed Cell Type Distribution.png"), width=6,height=8,units="in",res=600, pointsize = 14)
pdf(paste0("./Astrocyte/Astrocyte male Stressed Cell Type Distribution.pdf"), width=6,height=8,pointsize = 14, useDingbats = FALSE)
ggplot(maleStressedpie, aes(x = "", y = value, fill = group)) +
  geom_col() +
  geom_label(aes(label = labels),label.size = 0.25,
             position = position_stack(vjust = 0.5),
             show.legend = FALSE) +
  coord_polar(theta = "y") +
  scale_fill_brewer(palette="Set2")+
  theme_void()
dev.off()



#Bar Plot

png(paste0("./Astrocyte/Astrocyte male Stressed Cell Type DistributionBar plot all.png"), width=6,height=8,units="in",res=600, pointsize = 14)
#pdf(paste0("./Astrocyte/Astrocyte male Stressed Cell Type DistributionBar plot all.pdf"), width=6,height=8,pointsize = 14, useDingbats = FALSE)
ggplot(data=maleStressedpie, aes(x=1, y=labels, fill=group)) +
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


#pdf(paste0("./Astrocyte/CTE_Endoth Clustered Heatmap Stressed Subjects male.pdf"), width=6,height=8,pointsize = 14, useDingbats = FALSE)
png(paste0("Astrocyte/CTE_Endoth Clustered Heatmap Stressed Subjects male.png"), width=7,height=8,pointsize = 14, units="in",res=600)
heatmap3(genes_subZ_mat, ColSideColors =cbind(c(SexColorBar),c(freqColorBar),c(StressColorBar)), ColSideWidth=1, ColSideLabs=NA,  
         col=barColors, breaks=breakBarColors, legendfun=function()showLegend(legend=c(NA),col=c(NA)),
         Rowv=as.dendrogram(hrGenes), 
         Colv=NA,  scale="none",
         margins = c(3,10),
         labRow = NA, labCol=nrow(genes_subZ_mat),
         cexCol=1.3)
title("Genes Cluster Endothelial", adj = 0.5, line = 2)
dev.off()




#Heatmap of All genes with cell type marker


#male Data

sortby <- which(Expt_data$Condition =="M 10HzHz Stress " | Expt_data$Condition  =="M 40HzHz Stress " | Expt_data$Condition  =="M LightHz Stress " )
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
genes_subZ_mat <- as.matrix(genes_sub_Z[,1:{ncol(genes_sub_Z)-1}])
rownames(genes_subZ_mat) <- genes_sub_Z$geneID






# Specific Cell type Gene Sets
Neuron_gene_set <- geneSets$Neuron 
Ast_gene_set <- geneSets$Astrocytes 
Mig_gene_set <- geneSets$Microglia 
olig_gene_set<-geneSets$Oligodendrocytes
Endoth_gen_set<-geneSets$Endothelia 



#Neuron
Neuron_gene <- genes_subZ_mat[which(rownames(genes_subZ_mat) %in% Neuron_gene_set),] #1384
indNeu=match(rownames(Neuron_gene), rownames(genes_subZ_mat))
NeuGenes=genes_subZ_mat[indNeu,]
NeuGenes=data.frame(NeuGenes)
NeuGenes <- NeuGenes %>%
  mutate(Celltype = paste0("Neurons"))


#Astrocyte

Astro_gene <- genes_subZ_mat[which(rownames(genes_subZ_mat) %in% Ast_gene_set),] #1066
indAstr=match(rownames(Astro_gene), rownames(genes_subZ_mat))
astroGenes=genes_subZ_mat[indAstr,]

astroGenes=data.frame(astroGenes)
astroGenes <- astroGenes %>%
  mutate(Celltype = paste0("Astrocyte"))


#Microglia
Mig_gene <- genes_subZ_mat[which(rownames(genes_subZ_mat) %in% Mig_gene_set),] #1125
indMig=match(rownames(Mig_gene), rownames(genes_subZ_mat))
MigGenes=genes_subZ_mat[indMig,]

MigGenes=data.frame(MigGenes)
MigGenes <- MigGenes %>%
  mutate(Celltype = paste0("Microglia"))

#Oligodendrocytes

Olig_gene <- genes_subZ_mat[which(rownames(genes_subZ_mat) %in% olig_gene_set),] #895
indOlig=match(rownames(Olig_gene), rownames(genes_subZ_mat))
OligGenes=genes_subZ_mat[indOlig,]
OligGenes=data.frame(OligGenes)
OligGenes <- OligGenes %>%
  mutate(Celltype = paste0("Oligodendrocyte"))
#Endothelial

Endo_gene <- genes_subZ_mat[which(rownames(genes_subZ_mat) %in% Endoth_gen_set),] #718
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

export(data.frame(rownames(Data),Data),file="./Astrocyte/Male CTE Data with cell type annotation v2.csv")
hrGenes= hclust(dist((Data[, c(1:9)]),method = "euclidean"), method="ward.D2")

pdf(paste0("./Astrocyte/Gene Heatmap male Stressed CTE Color Bar with cluster.pdf"), width=6,height=8,pointsize = 14, useDingbats = FALSE)
#png(paste0("./Astrocyte/Gene Heatmap male Stressed CTE Color Bar with cluster.png"), width=7,height=8,pointsize = 14, units="in",res=600)
heatmap3(Data[, c(1:9)], RowSideColors=CellTypecolorBar,ColSideColors =cbind(c(SexColorBar),c(freqColorBar),c(StressColorBar)), ColSideWidth=1, ColSideLabs=NA,  
         col=barColors, breaks=breakBarColors, legendfun=function()showLegend(legend=c(NA),col=c(NA)),
        # Rowv=NA, 
         Rowv=as.dendrogram(hrGenes), 
         Colv=NA,  scale="none",
         margins = c(3,10),
         labRow = NA, labCol=NA,
         cexCol=1.3)
title("All genes With CTE", adj = 0.5, line = 2)
dev.off()




###Astrocyte  Specific Gene Set (Custome)

sortby <- which(Expt_data$Condition =="M 10HzHz Stress " | Expt_data$Condition  =="M 40HzHz Stress " | Expt_data$Condition  =="M LightHz Stress " )
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
geneSetsMat=read_excel("Reference Files/GeneSets_Astrocytes_FXN1.xlsx")
geneSets=as.list(geneSetsMat)
geneSets=lapply(geneSets, function(x) x[!is.na(x)])


geneSetEnrich=gsva(genes_sub_Z_mat, geneSets, mx.diff=TRUE, kcdf="Gaussian")
geneSetEnrichSortZ=t(apply(geneSetEnrich,1,scale)); #z-score the data
geneSetEnrichSortZ[is.na(geneSetEnrichSortZ)]=0
colnames(geneSetEnrichSortZ)=colnames(genes_sub_Z_mat)





Freq=Expt_data_sub$Frequency
datagroupVars=data.frame(t(geneSetEnrichSortZ), Freq)


for (ind in 1:20)
{
  
  #png(paste0("Astrocyte/Astrocyte Custom GS/GS male boxplot",ind,".png"),height=4,width=5,units="in",res=1000)
  pdf(paste0("Astrocyte/Astrocyte Custom GS/GS male boxplot",ind,".pdf"),height=4,width=5)
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



#--------------------------------------------------------------------DEGs--------------------------------------------------------
#Selecting Condition of interest for volcano Plots
#load your RAW Counts Data 


rawcounts <- read.csv("Data/Astrocye Raw Counts V2.csv", header = TRUE)
geneName=rawcounts[,1]
geneID=paste0(geneName,"|", rawcounts$Gene.ID)
rawcounts=rawcounts[, -c(1,2)]
rownames(rawcounts)=geneID
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
Expt_data_in <- read_excel("Data/MetaDataAst.xlsx" , sheet="Sheet1") %>%  # imports the metadata
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






sortby <- which(Expt_data$Condition =="M 10HzHz Stress " | Expt_data$Condition  =="M LightHz Stress " )#This  can be changed to any group or condition
Expt_data_sub=Expt_data[sortby,]
genes_sub=genes[,sortby]


dds <- DESeqDataSetFromMatrix(countData = genes_sub,
                              colData = Expt_data_sub,
                              design = ~ Frequency)

dds$Frequency <- relevel( dds$Frequency, ref = "Light")


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
export(outMale,"Astrocyte/Deseq2/Astrocyte DEGs_Stress male 10Hz vs male Light DESeq2 ThreshholdPassing FC2.csv")

#####Export Data#####
DesqDataOut=as.matrix(data.frame(DesqData))
DseqExport <-data.frame(rownames(DesqDataOut),DesqDataOut)
export(DseqExport,"Astrocyte/Deseq2/Astrocyte DEGs_Stress male 10Hz vs male Light  DESeq2 All.csv")


keyvals <- ifelse(
  DesqData$log2FoldChange < -1 & DesqData$pvalue < 5*10e-3, 'royalblue',
  ifelse(DesqData$log2FoldChange > 1 & DesqData$pvalue < 5*10e-3, 'red',
         'gray'))
keyvals[is.na(keyvals)] <- 'gray'
names(keyvals)[keyvals == 'red'] <- 'Up'
names(keyvals)[keyvals == 'gray'] <- 'NS'
names(keyvals)[keyvals == 'royalblue'] <- 'Down'

# Setting up output folder for the figure
output_folder=output_folder = paste0("Astrocyte/Deseq2/")
dir.create(output_folder,recursive = TRUE, showWarnings=FALSE)


#---- Generating Volcano Plot---#

pdf(paste0(output_folder,"Astrocyte DEGs male 10Hz vs male Light Stress FC2.pdf"), width=7,height=8,pointsize = 14,useDingbats = FALSE)
#png(paste0(output_folder," Astrocyte DEGs male 10Hz vs male Light Stress FC2.png"), width=7,height=8,pointsize = 14,units="in",res=600)
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


#-------------------------DEGs Venn Diagrams---------------------------#

Expt_data_in <- read_excel("Data/MetaDataAst.xlsx" , sheet="Sheet1") %>%  # imports the metadata
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



#Female

sortby1 <- which(Expt_data$Condition =="F 10HzHz Stress " | Expt_data$Condition  =="F 40HzHz Stress " | Expt_data$Condition  =="F LightHz Stress " )
#This  can be changed to any group or condition
Expt_data_sub_F=Expt_data[sortby1,]
genes_sub_F=genes[,sortby1]

rownames(genes_sub_F)=sub("\\|.*", "", rownames(genes_sub_F))

genes_sub_Z_F <- genes_sub_F %>%
  apply(1,"scale") %>%
  t() %>%
  as_tibble()
colnames(genes_sub_Z_F) = colnames(genes_sub_F)
genes_sub_Z_F$geneID = rownames(genes_sub_F)
genes_sub_Z_F[is.na(genes_sub_Z_F)] <- 0   # set any NA values to 0

# Create a matrix of Z-scored subset data
genes_subZ_mat_F <- as.matrix(genes_sub_Z_F[,1:{ncol(genes_sub_Z_F)-1}])
rownames(genes_subZ_mat_F) <- genes_sub_Z_F$geneID


#Male

sortby2 <- which(Expt_data$Condition =="M 10HzHz Stress " | Expt_data$Condition  =="M 40HzHz Stress " | Expt_data$Condition  =="M LightHz Stress ")
#This  can be changed to any group or condition
Expt_data_sub_M=Expt_data[sortby2,]
genes_sub_M=genes[,sortby2]

GeneName_M=rownames(genes_sub_M)

rownames(genes_sub_M)=sub("\\|.*", "", GeneName_M)
genes_sub_Z_M <- genes_sub_M %>%
  apply(1,"scale") %>%
  t() %>%
  as_tibble()
colnames(genes_sub_Z_M) = colnames(genes_sub_M)
genes_sub_Z_M$geneID = rownames(genes_sub_M)
genes_sub_Z_M[is.na(genes_sub_Z_M)] <- 0   # set any NA values to 0

# Create a matrix of Z-scored subset data

genes_subZ_mat_M <- as.matrix(genes_sub_Z_M[,1:{ncol(genes_sub_Z_M)-1}])
rownames(genes_subZ_mat_M) <- genes_sub_Z_M$geneID

#Venn Diagram of all DEGs
#FC =2 , p<0.05


#Female
#"Astrocyte/Deseq2/Astrocyte DEGs_Stress male 10Hz vs male Light DESeq2 ThreshholdPassing FC2.csv")
Female40 <- read_excel("Astrocyte/Deseq2/Astrocyte DEGs_Stress Female 40Hz vs Female Light DESeq2 ThreshholdPassing FC2.xlsx",sheet="Up") 
FemalesigOE40 <- data.frame(subset( Female40, threshold==TRUE))

Female10 <- read_excel("Astrocyte/Deseq2/Astrocyte DEGs_Stress Female 10Hz vs Female Light DESeq2 ThreshholdPassing FC2.xlsx",sheet="Up") 
FemalesigOE10 <- data.frame(subset( Female10, threshold==TRUE))
# 
# Female1040 <- read_excel("Astrocyte/Deseq2/Tina Data Stress 40Hz vs 10 Female DESeq2 All FC2 111323SB.xlsx",sheet="Up") 
# FemalesigOE1040 <- data.frame(subset( Female1040, threshold==TRUE))
# 

#Male


Male40 <- read_excel("Astrocyte/Deseq2/Astrocyte DEGs_Stress male 40Hz vs male Light DESeq2 ThreshholdPassing FC2.xlsx",sheet="Up") 
MalesigOE40 <- data.frame(subset(Male40, threshold==TRUE))


Male10 <- read_excel("Astrocyte/Deseq2/Astrocyte DEGs_Stress male 10Hz vs male Light DESeq2 ThreshholdPassing FC2.xlsx",sheet="Up") 
MalesigOE10 <- data.frame(subset( Male10, threshold==TRUE))

#Cutoff Threshhold
pvalue.cutoff <- 5*10e-3
lfc.cutoff <- 1



Male40L=MalesigOE40$rownames.sigOE_ordered.
Male10L=MalesigOE10$rownames.sigOE_ordered.
Female10L=FemalesigOE10$rownames.sigOE_ordered.
Female40L=FemalesigOE40$rownames.sigOE_ordered.

#Male down 40Hz specific 
MaleData40=genes_subZ_mat_M[match(MalesigOE40$rownames.sigOE_ordered. , sub("\\|.*", "",rownames(genes_subZ_mat_M))) ,]
MaleDataOthers<-genes_subZ_mat_M[which(Male40L %in% c(Male10L,Female10L,Female40L)),]

indRemove=match(rownames(MaleDataOthers), rownames(genes_subZ_mat_M))

MaleData40only<- MaleData40[-c(indRemove) ,]

export(data.frame(rownames(MaleData40only),MaleData40only),file="Astrocyte/Deseq2/Male Astrocyte 40Hz specific DEGs on venn diagram Up.csv")

#Male down 10Hz specific 

MaleData10=genes_subZ_mat_M[match(MalesigOE10$rownames.sigOE_ordered. , sub("\\|.*", "",rownames(genes_subZ_mat_M))) ,]
MaleDataOthers<-genes_subZ_mat_M[which(Male10L %in% c(Male40L,Female10L,Female40L)),]


indRemove=match(rownames(MaleDataOthers), rownames(genes_subZ_mat_M))

MaleData10only<- MaleData10[-c(indRemove) ,]

export(data.frame(rownames(MaleData10only),MaleData10only),file="Astrocyte/Deseq2/Male Astrocyte 10Hz specific DEGs on venn diagram Up.csv")



#Female down 10Hz specific 

FemaleData10=genes_subZ_mat_F[match(FemalesigOE10$rownames.sigOE_ordered. , sub("\\|.*", "",rownames(genes_subZ_mat_F))) ,]
FemaleDataOthers<-genes_subZ_mat_F[which(Female10L %in% c(Male40L,Male10L,Female40L)),]


indRemove=match(rownames(FemaleDataOthers), rownames(genes_subZ_mat_F))

FemaleData10only<- FemaleData10[-c(indRemove) ,]

export(data.frame(rownames(FemaleData10only),FemaleData10only),file="Astrocyte/Deseq2/Female Astrocyte 10Hz specific DEGs on venn diagram Up.csv")




#Female down 40Hz specific 

FemaleData40=genes_subZ_mat_F[match(FemalesigOE40$rownames.sigOE_ordered. , sub("\\|.*", "",rownames(genes_subZ_mat_F))) ,]
FemaleDataOthers<-genes_subZ_mat_F[which(Female40L %in% c(Male40L,Male10L,Female10L)),]


indRemove=match(rownames(FemaleDataOthers), rownames(genes_subZ_mat_F))

FemaleData40only<- FemaleData40[-c(indRemove) ,]

export(data.frame(rownames(FemaleData40only),FemaleData40only),file="Astrocyte/Deseq2/Female Astrocyte 40Hz specific DEGs on venn diagram Up.csv")


output_folder = paste0("Astrocyte/Deseq2/")
VennFreq=list(Male40L=MalesigOE40$rownames.sigOE_ordered.,Male10L=MalesigOE10$rownames.sigOE_ordered.,Female10L=FemalesigOE10$rownames.sigOE_ordered., Female40L=FemalesigOE40$rownames.sigOE_ordered.)

png(paste0(output_folder,"Venn Diagram DEGs Astrocyte Up Only.png"), width=7,height=8,pointsize = 14,units="in",res=600)
#pdf(paste0(output_folder,"Venn Diagram DEGs Astrocyte Up Only.pdf"), width=8,height=8,pointsize = 12)
ggvenn(VennFreq, 
       #fill_color = c("#E41A1C","#377EB8"),
       stroke_size = 0.5, set_name_size = 6, text_size = 5)
dev.off()








#Female
#"Astrocyte/Deseq2/Astrocyte DEGs_Stress male 10Hz vs male Light DESeq2 ThreshholdPassing FC2.csv")
Female40 <- read_excel("Astrocyte/Deseq2/Astrocyte DEGs_Stress Female 40Hz vs Female Light DESeq2 ThreshholdPassing FC2.xlsx",sheet="Down") 
FemalesigOE40 <- data.frame(subset( Female40, threshold==TRUE))

Female10 <- read_excel("Astrocyte/Deseq2/Astrocyte DEGs_Stress Female 10Hz vs Female Light DESeq2 ThreshholdPassing FC2.xlsx",sheet="Down") 
FemalesigOE10 <- data.frame(subset( Female10, threshold==TRUE))
# 
# Female1040 <- read_excel("Astrocyte/Deseq2/Tina Data Stress 40Hz vs 10 Female DESeq2 All FC2 111323SB.xlsx",sheet="Down") 
# FemalesigOE1040 <- data.frame(subset( Female1040, threshold==TRUE))
# 

#Male


Male40 <- read_excel("Astrocyte/Deseq2/Astrocyte DEGs_Stress male 40Hz vs male Light DESeq2 ThreshholdPassing FC2.xlsx",sheet="Down") 
MalesigOE40 <- data.frame(subset(Male40, threshold==TRUE))


Male10 <- read_excel("Astrocyte/Deseq2/Astrocyte DEGs_Stress male 10Hz vs male Light DESeq2 ThreshholdPassing FC2.xlsx",sheet="Down") 
MalesigOE10 <- data.frame(subset( Male10, threshold==TRUE))

#Cutoff Threshhold
pvalue.cutoff <- 5*10e-3
lfc.cutoff <- 1


















#Cutoff Threshhold
pvalue.cutoff <- 5*10e-3
lfc.cutoff <- 1
output_folder = paste0("Astrocyte/Deseq2/")
VennFreq=list(Male10L=MalesigOE10$rownames.sigOE_ordered., Female40L=FemalesigOE40$rownames.sigOE_ordered.)



png(paste0(output_folder,"Venn Diagram Male10 Female 40 Down.png"), width=7,height=8,pointsize = 14,units="in",res=600)
#pdf(paste0(output_folder,"Venn Diagram Male10 Female 40 Down.pdf"), width=8,height=8,pointsize = 12)
ggvenn(VennFreq, 
       fill_color = c("dodgerblue2","hotpink"),
       stroke_size = 0.5, set_name_size = 6, text_size = 5)
dev.off()





