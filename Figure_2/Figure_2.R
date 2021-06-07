###################################################
# Matt Regner
# Franco Lab 
# 2020-2021
# Description: plot P2G histograms, plot P2G heatmap,
# re-run overlap analysis, plot browser track
###################################################
source("P2G_Heatmap_Distal.R")
source("Archr_Peak_RawPval.R")
source("Archr_Peak_Null_Permute.R")
library(ggplot2)
library(Seurat)
library(scales)
library(forcats)
library(RColorBrewer)
library(ArchR)
library(stringr)
library(stringi)
library(dplyr)
library(tidyr)
addArchRThreads(32)

####################################################################
# PART 1: summary distribution histograms (correlation and p-values)
# 1) Observed data
# 2) Null data example
####################################################################
# Read in P2G data
p2g <- readRDS("./All_P2G_Observed.rds")
p2g$RawPVal <- as.numeric(p2g$RawPVal)
p2g$Correlation <- as.numeric(p2g$Correlation)

# Plot correlation histogram
p1 <- ggplot(p2g,aes(x=Correlation))+geom_histogram(color="black",fill="gray90",size=0.3,boundary = 0,bins=15)+
  theme_bw()+
  scale_y_continuous(expand = c(0,0),limits=c(0,1500000)) +
  scale_x_continuous(expand= c(0,0),limits = c(-1,1))

# Plot p-value histogram
p2 <- ggplot(p2g,aes(x=RawPVal))+geom_histogram(color="black",fill="gray90",size=0.3,boundary = 0,bins=15)+
  theme_bw()+
  scale_y_continuous(expand = c(0,0),
                     limits = c(0,1.75e+06)) +
  scale_x_continuous(expand = c(0,0),
                     limits = c(0,1))


# Read in P2G data (null example)
p2g <- readRDS("./All_P2G_Null_example.rds")
p2g$RawPVal <- as.numeric(p2g$RawPVal)
p2g$Correlation <- as.numeric(p2g$Correlation)

# Plot correlation histogram
p3 <- ggplot(p2g,aes(x=Correlation))+geom_histogram(color="black",fill="gray90",size=0.3,boundary = 0,bins=15)+
  theme_bw()+
  scale_y_continuous(limits=c(0,1500000),expand=c(0,0)) +
  scale_x_continuous(expand = c(0,0),limits = c(-1,1))

# Plot p-value histogram
p4 <- ggplot(p2g,aes(x=RawPVal))+geom_histogram(color="black",fill="gray90",size=0.3,boundary = 0,bins=15)+
  theme_bw()+
  scale_y_continuous(expand = c(0,0),
                     limits = c(0,1.75e+06)) +
  scale_x_continuous(expand = c(0,0),limits = c(0,1))
# Plot side-by-side
CombinePlots(list(p1,p2,p3,p4),ncol=2)+ggsave("Histograms_obs_null.pdf",height =3,width = 4)


####################################################################
# PART 2: Plot histograms for number of genes regulated by peaks
####################################################################
# Number of peaks per number of genes
p2g <- readRDS("./Cancer_specific_P2G_table.rds")
p2g$RawPVal <- as.numeric(p2g$RawPVal)
p2g$Correlation <- as.numeric(p2g$Correlation)

df <- data.frame(num.genes = table(p2g$peakName))

df$cat <- ifelse(df$num.genes.Freq < 3,"1to2",df$num.genes.Freq)
df$cat <- ifelse(df$num.genes.Freq < 6 &df$num.genes.Freq >2,"3to5",df$cat)
df$cat <- ifelse(df$num.genes.Freq < 50 &df$num.genes.Freq >5,"6plus",df$cat)

df <- data.frame(table(df$cat))

p1 <- ggplot(df,aes(x=Var1,y=Freq))+geom_bar(stat = "identity")+theme_classic()+
  geom_text(aes(label=Freq), vjust=-0.3, size=3.5)



# Number of genes per number of peaks
df <- data.frame(num.peaks = table(p2g$geneName))

df$cat <- ifelse(df$num.peaks.Freq < 3,"1to2",df$num.peaks.Freq)
df$cat <- ifelse(df$num.peaks.Freq < 6 &df$num.peaks.Freq >2,"3to5",df$cat)
df$cat <- ifelse(df$num.peaks.Freq < 50 &df$num.peaks.Freq >5,"6plus",df$cat)

df <- data.frame(table(df$cat))

p2 <- ggplot(df,aes(x=Var1,y=Freq))+geom_bar(stat = "identity")+theme_classic()+
  geom_text(aes(label=Freq), vjust=-0.3, size=3.5)

CombinePlots(list(p1,p2))+ggsave("GenesPer_PeaksPer_Histograms.pdf",width = 10,height = 4)

####################################################################
# PART 3: browser track for RHEB enhancers
# 1) Plot browser track
# 2) Verify cancer-specific enhancers have differential accessibility
####################################################################

# Read in other annotation features:
encode.all <- read.delim("./GRCh38-ccREs.bed",header =F)
colnames(encode.all)[1:3] <- c("seqnames","start","end")
encode.all <- makeGRangesFromDataFrame(encode.all)

ft.peaks <- readRDS("./Fallopian_Tube_Cell_line_Peaks.rds")
ft.peaks <- ft.peaks[,3:5]
colnames(ft.peaks)[1:3] <- c("seqnames","start","end")
ft.peaks <- makeGRangesFromDataFrame(ft.peaks)

ov.peaks <- readRDS("./Ovarian_Epithelial_Cell_line_Peaks.rds")
ov.peaks <- ov.peaks[,3:5]
colnames(ov.peaks)[1:3] <- c("seqnames","start","end")
ov.peaks <- makeGRangesFromDataFrame(ov.peaks)

atac <- readRDS("final_archr_proj_archrGS-P2Gs.rds")
# ATAC
levels(factor(atac$predictedGroup_ArchR))
my_levels <- as.character(c(11,20,21,22,31,
               19,34,
               3,
               9,10,
               16,17,
               0,27,
               6,8,12,14,15,18,24,25,26,29,
               7,23,
               1,33,
               2,4,30,
               5,13,
               32,
               28,35))

for ( i in levels(factor(atac$predictedGroup_ArchR))){
  num <-  gsub("-.*","",i)
  idx <- match(num,my_levels)
  atac$predictedGroup_ArchR <- str_replace(atac$predictedGroup_ArchR,pattern = i,replacement = paste0(idx,"_",atac$predictedGroup_ArchR))
  print("iter complete")
}
test <- as.data.frame(atac@cellColData)
test <- test[test$Sample == "38FE7L",]
levels(factor(test$predictedGroup_ArchR))

new.idents <- setdiff(levels(factor(atac$predictedGroup_ArchR)),
        c("16_8-Stromal fibroblasts",
          "17_12-Stromal fibroblasts",
          "18_14-Stromal fibroblasts",
          "19_15-Stromal fibroblasts",
          "20_18-Fibroblast",
          "21_24-Fibroblast",
          "23_215_23_26-Fibroblast",
          "24_29-Fibroblast",
          "26_23-Stromal fibroblasts",
          "31_30-T cell",
          "30_4-Lymphocytes",
          "34_32-Mast cell",
          "33_13-Macrophages",
          "36_35-B cell",
          "32_5-Macrophage",
          "35_28-B cell" ))

idxSample <- BiocGenerics::which(atac$predictedGroup_ArchR %in% new.idents)
cellsSample <- atac$cellNames[idxSample]
atac.sub <- atac[cellsSample, ]

#########################################################################
# Make modified getP2G function:
getPeak2GeneLinks.mod <- function(
  ArchRProj = NULL, 
  corCutOff = 0.45, 
  PValCutOff = 0.0001,
  varCutOffATAC = 0.25,
  varCutOffRNA = 0.25,
  resolution = 1, 
  returnLoops = TRUE
){
  
  .validInput(input = ArchRProj, name = "ArchRProj", valid = "ArchRProject")
  .validInput(input = corCutOff, name = "corCutOff", valid = "numeric")
  .validInput(input = PValCutOff, name = "PValCutOff", valid = "numeric")
  .validInput(input = varCutOffATAC, name = "varCutOffATAC", valid = "numeric")
  .validInput(input = varCutOffRNA, name = "varCutOffRNA", valid = "numeric")
  .validInput(input = resolution, name = "resolution", valid = c("integer", "null"))
  .validInput(input = returnLoops, name = "returnLoops", valid = "boolean")
  
  if(is.null(ArchRProj@peakSet)){
    return(NULL)
  }
  
  if(is.null(metadata(ArchRProj@peakSet)$Peak2GeneLinks)){
    
    return(NULL)
    
  }else{
    
    p2g <- metadata(ArchRProj@peakSet)$Peak2GeneLinks
    p2g <- p2g[which(p2g$Correlation >= corCutOff & p2g$RawPVal <= PValCutOff), ,drop=FALSE]
    
    if(!is.null(varCutOffATAC)){
      p2g <- p2g[which(p2g$VarQATAC > varCutOffATAC),]
    }
    
    if(!is.null(varCutOffRNA)){
      p2g <- p2g[which(p2g$VarQRNA > varCutOffRNA),]
    }
    
    if(returnLoops){
      
      peakSummits <- resize(metadata(p2g)$peakSet, 1, "center")
      geneStarts <- resize(metadata(p2g)$geneSet, 1, "start")
      
      if(!is.null(resolution)){
        summitTiles <- floor(start(peakSummits) / resolution) * resolution + floor(resolution / 2)
        geneTiles <- floor(start(geneStarts) / resolution) * resolution + floor(resolution / 2)
      }else{
        summitTiles <- start(peakSummits)
        geneTiles <- start(geneTiles)
      }
      
      loops <- .constructGR(
        seqnames = seqnames(peakSummits[p2g$idxATAC]),
        start = summitTiles[p2g$idxATAC],
        end = geneTiles[p2g$idxRNA]
      )
      mcols(loops)$value <- p2g$Correlation
      mcols(loops)$FDR <- p2g$FDR
      
      loops <- loops[order(mcols(loops)$value, decreasing=TRUE)]
      loops <- unique(loops)
      loops <- loops[width(loops) > 0]
      loops <- sort(sortSeqlevels(loops))
      
      loops <- SimpleList(Peak2GeneLinks = loops)
      
      return(loops)
      
    }else{
      
      return(p2g)
      
    }
    
  }
  
}



#########################################################################

plot <- plotBrowserTrack(atac.sub,geneSymbol ="RHEB", groupBy = "predictedGroup_ArchR",
                         features = GRangesList(TrackA = encode.all,TrackB = ft.peaks,TrackC = ov.peaks), 
                         loops = getPeak2GeneLinks.mod(atac,corCutOff = 0.45,
                                                   PValCutOff = 1e-12,varCutOffATAC = 0,
                                                   varCutOffRNA = 0),upstream = 7000,downstream = 35000)

pdf("RHEB_final.pdf",width = 6,height = 8)
grid::grid.draw(plot[[1]])
dev.off()

# 
# # Make violin plots for RHEB expression and mTOR pathway expression
# # Relevel object@ident
# rna@active.ident <- factor(x =rna$RNA_snn_res.0.7, levels = rev(my_levels))
# p1 <- VlnPlot(rna,features = "RHEB",pt.size = 0)+coord_flip()+NoLegend()+ggsave("RHEB_vln_new.pdf",width = 4,height = 8)
# 
# p1 <- p1$data
# colnames(p1)
# kruskal.test(RHEB~ident,data = p1)
# 
# 
# gset <- read.delim("./BIOCARTA_MTOR_PATHWAY.txt",header = T)
# gset <- as.character(gset$BIOCARTA_MTOR_PATHWAY[2:length(gset$BIOCARTA_MTOR_PATHWAY)])
# rna <- AddModuleScore(rna,features = list(gset),name = "mTOR_members",search = T)
# 
# 
# p1 <- VlnPlot(rna,features = "mTOR_members1",pt.size = 0)+coord_flip()+NoLegend()+
#   ggsave("mTOR_pathway_vln_new.pdf",width = 4,height = 8)
# 
# p1 <- p1$data
# colnames(p1)
# kruskal.test(mTOR_members1~ident,data = p1)
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# ###########################################################
# #############################################################
# # Read in P2G data
# p2g <- readRDS("./All_P2G_Observed.rds")
# p2g$RawPVal <- as.numeric(p2g$RawPVal)
# p2g$Correlation <- as.numeric(p2g$Correlation)
# 
# # Plot histogram of Correlations and Pvalues
# sub <- p2g[p2g$RawPVal <=1e-12,]
# sub.neg <- sub[sub$Correlation < 0,]
# lower <- last(sort(sub.neg$Correlation))
# sub.pos <- sub[sub$Correlation > 0,]
# upper <- first(sort(sub.pos$Correlation))
# 
# min.obs.cor <- min(p2g$Correlation)
# max.obs.cor <- max(p2g$Correlation)
# 
# p1 <- ggplot(sub,aes(x=Correlation))+geom_histogram(color="black",fill="gray90",size=0.3,boundary = 0,bins = 12)+
#   theme_bw()+geom_vline(xintercept = c(lower,upper),linetype="dashed",size=0.3,color="red")+
#   scale_y_continuous(expand = c(0,0),limits=c(0,1.5e+06)) +
#   scale_x_continuous(expand= c(0,0),limits = c(min.obs.cor,max.obs.cor))
# 
# min.obs.p <- min(sub$RawPVal)
# max.obs.p <- max(sub$RawPVal)
# 
# p2 <- ggplot(sub,aes(x=RawPVal))+geom_histogram(color="black",fill="gray90",size=0.3,boundary = 0,bins=10)+
#   theme_bw()+geom_vline(xintercept = 1e-12,linetype="dashed",color="red")+
#   scale_y_continuous(expand = c(0,0),
#                      limits = c(0,3.5e+05)) +
#   scale_x_continuous(expand = c(0,0),
#                      limits = c(min.obs.p,max.obs.p))
# 
# CombinePlots(list(p1,p2),ncol=2)+ggsave("Histrograms.pdf",height =2,width = 6)
# 
# # Read in P2G data (null example)
# # Read in P2G data
# p2g <- readRDS("./All_P2G_Null_example.rds")
# p2g$RawPVal <- as.numeric(p2g$RawPVal)
# p2g$Correlation <- as.numeric(p2g$Correlation)
# 
# # Plot histogram of Correlations and Pvalues
# sub <- p2g[p2g$RawPVal <=1e-12,]
# sub.neg <- sub[sub$Correlation < 0,]
# lower <- last(sort(sub.neg$Correlation))
# sub.pos <- sub[sub$Correlation > 0,]
# upper <- first(sort(sub.pos$Correlation))
# 
# p1 <- ggplot(sub,aes(x=Correlation))+geom_histogram(color="black",fill="gray90",size=0.3,boundary = 0,bins = 12)+
#   theme_bw()+geom_vline(xintercept = c(lower,upper),linetype="dashed",size=0.3,color="red")+
#   scale_y_continuous(expand = c(0,0),limits=c(0,100)) +
#   scale_x_continuous(expand= c(0,0),limits = c(min.obs.cor,max.obs.cor))
# 
# p2 <- ggplot(sub,aes(x=RawPVal))+geom_histogram(color="black",fill="gray90",size=0.3,boundary = 0,bins=10)+
#   theme_bw()+geom_vline(xintercept = 1e-12,linetype="dashed",color="red")+
#   scale_y_continuous(expand = c(0,0),
#                      limits = c(0,100)) +
#   scale_x_continuous(expand = c(0,0),
#                      limits = c(0,1e-12))
# 
# CombinePlots(list(p1,p2),ncol=2)+ggsave("Histrograms-null.pdf",height =2,width = 6)
# 


