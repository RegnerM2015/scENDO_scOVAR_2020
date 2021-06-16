###################################################
# Matt Regner
# Franco Lab 
# 2020-2021
# Description: plot figures for Fig. 2
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
h5disableFileLocking()
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




############################################################################
# PART 1.5.0.5: How many cancer specific peaks are deferentially accessible?
############################################################################

# Marker Peaks for malignant clusters with 100% patient specificity 
atac <- readRDS("final_archr_proj_archrGS-P2Gs.rds")
markersPeaks <- getMarkerFeatures(
  ArchRProj = atac, 
  useMatrix = "PeakMatrix", 
  groupBy = "predictedGroup_ArchR",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon"
)
levels(factor(atac$predictedGroup_ArchR))
labels <- c("0-Fibroblast",
            "10-Epithelial cell" ,
            "11-Unciliated epithelia 1",
            "16-Fibroblast",
            "17-Epithelial cell" ,
            "19-Epithelial cell" ,
            "20-Ciliated",
            "21-Unciliated epithelia 1",
            "22-Unciliated epithelia 2",
            "3-Epithelial cell",
            "31-Unciliated epithelia 1",
            "34-Epithelial cell",
            "9-Epithelial cell",
            "27-Fibroblast" 
            
            )

markerList <- getMarkers(markersPeaks, cutOff = "FDR <= 0.1 & Log2FC >= 0.5")
markerList.sub <- markerList[names(markerList) %in% labels]
markers.peaks <- unlist(markerList.sub)
markers.peaks <- paste0(markers.peaks$seqnames,":",markers.peaks$start,"-",markers.peaks$end)
markers.peaks <- unique(markers.peaks)

p2g.cancer <- readRDS("./Cancer_specific_P2G_table.rds")
p2g.cancer.diff <- p2g.cancer[p2g.cancer$peakName %in% markers.peaks,]
print(paste0(round(nrow(p2g.cancer.diff)/nrow(p2g.cancer)*100,3),"% of cancer-specifc P2Gs have differentially accessible peaks in the malignant cell populations! (FDR <= 0.1 & Log2FC >= 0.5)"))
print(paste0(round(length(unique(p2g.cancer.diff$peakName))/length(unique(p2g.cancer$peakName))*100,3),"% of cancer-specific peaks are differentially accessible in the malignant cell populations! (FDR <= 0.1 & Log2FC >= 0.5)"))

p2g.cancer <- readRDS("./Cancer_enriched_P2G_table.rds")
p2g.cancer.diff <- p2g.cancer[p2g.cancer$peakName %in% markers.peaks,]
print(paste0(round(nrow(p2g.cancer.diff)/nrow(p2g.cancer)*100,3),"% of cancer-enriched P2Gs have differentially accessible peaks in the malignant cell populations! (FDR <= 0.1 & Log2FC >= 0.5)"))
print(paste0(round(length(unique(p2g.cancer.diff$peakName))/length(unique(p2g.cancer$peakName))*100,3),"% of cancer-enriched peaks are differentially accessible in the malignant cell populations! (FDR <= 0.1 & Log2FC >= 0.5)"))


####################################################################
# PART 1.5.1: Investigate the near neighboring gene rule in P2Gs:
####################################################################
p2g <- readRDS("./All_P2G_Observed.rds")
p2g <- p2g[p2g$Correlation >= 0.45,]
p2g <- p2g[p2g$RawPVal <= 1e-12,]
p2g <- p2g[p2g$peakType == "Distal",]# Subset to distal P2Gs
p2g$idx <- paste0(p2g$idxATAC,"-",p2g$idxRNA)
p2g.cancer <- readRDS("./Cancer_specific_P2G_table.rds")
p2g$Cancer.Specific <- ifelse(p2g$idx %in% p2g.cancer$idx,"Cancer","Normal")

p2g$p2g.pair <- paste0(p2g$peakName,"::",p2g$geneName)

atac <- readRDS("final_archr_proj_archrGS-P2Gs.rds")
peak.info <- getPeakSet(atac)
p2g.nearest.pair <- paste0(peak.info@seqnames,":",peak.info@ranges,"::",peak.info$nearestGene)

p2g$Nearest <- ifelse(p2g$p2g.pair %in% p2g.nearest.pair,"Yes","No")


p2g.norm <- p2g[p2g$Cancer.Specific == "Normal",]
p2g.cancer <- p2g[p2g$Cancer.Specific == "Cancer",]


p2g.yes <- p2g[p2g$Nearest == "Yes",]
print(length(unique(p2g.yes$peakName))/length(unique(p2g$peakName)))
print(nrow(p2g.yes)/nrow(p2g))

ggplot(p2g) +
  aes(x = Cancer.Specific, fill = factor(Nearest),boundary=0) +
  geom_bar(position = "fill")+
  scale_y_continuous(expand = c(0,0))+
  scale_x_discrete(expand=c(0,0))+
  theme_classic()+ggsave("ProportionofP2Gs_predicted_by_NearestRule.pdf",width = 4,height = 3)

# Proportions differ significantly by fisher exact test:

cancer.nearest <- p2g.cancer[p2g.cancer$Nearest == "Yes",]
norm.nearest <- p2g.norm[p2g.norm$Nearest == "Yes",]

res <- fisher.test(matrix(c(nrow(cancer.nearest), nrow(norm.nearest), 
                            nrow(p2g.cancer)-nrow(cancer.nearest), nrow(p2g.norm)-nrow(norm.nearest)),
                            ncol = 2))
print(res)
print(paste0("The proportion of cancer-specific P2Gs predictable by nearest neighboring rule is significantly less relative to the normal P2Gs (p=",res$p.value,")"))
####################################################################
# PART 1.5.2: Investigate the distances between peaks and linked genes
####################################################################
# p2g <- readRDS("./All_P2G_Observed.rds")
# p2g <- p2g[p2g$Correlation >= 0.45,]
# p2g <- p2g[p2g$RawPVal <= 1e-12,]
# p2g <- p2g[p2g$peakType == "Distal",]# Subset to distal P2Gs
# p2g$idx <- paste0(p2g$idxATAC,"-",p2g$idxRNA)
# p2g.cancer <- readRDS("./Cancer_specific_P2G_table.rds")
# p2g$Cancer.Specific <- ifelse(p2g$idx %in% p2g.cancer$idx,"Cancer","Normal")
# 
# p2g$p2g.pair <- paste0(p2g$peakName,"::",p2g$geneName)
# 
# atac <- readRDS("final_archr_proj_archrGS-P2Gs.rds")
# gene.info <- getGenes(atac)
# gene.info <- data.frame(geneName=gene.info$symbol,
#                         geneCoord=paste0(gene.info@seqnames,":",gene.info@ranges))
# 
# p2g <- merge(p2g,gene.info,by="geneName")
# p2g.check.trans <- tidyr::separate(p2g, col="peakName",into=c("chrom.peak","ranges.peak"),sep=":")
# p2g.check.trans <- tidyr::separate(p2g.check.trans, col="geneCoord",into=c("chrom.gene","ranges.gene"),sep=":")
# length(which(p2g.check.trans$chrom.gene == p2g.check.trans$chrom.peak))# No trans interactions are in the data
# 
# peak.gr <- GRanges(p2g$peakName)
# gene.gr <- GRanges(p2g$geneCoord)
# distance.gr <- GenomicRanges::distance(peak.gr,gene.gr)
# 
# p2g$p2g.distance <- distance.gr
# 
# 
# p2g.norm <- p2g[p2g$Cancer.Specific == "Normal",]
# p2g.cancer <- p2g[p2g$Cancer.Specific == "Cancer",]
# 
# 
# p <- ggplot(p2g,aes(x=p2g.distance))
# p+geom_histogram(aes(fill=Cancer.Specific),position = "identity",alpha=0.5)+
#   scale_fill_manual(values = c("orange", "gray"))+
#   geom_vline( aes(xintercept=mean(p2g.cancer$p2g.distance)),
#               linetype="dashed",col="orange")+
#   geom_vline( aes(xintercept=mean(p2g.norm$p2g.distance)),
#               linetype="dashed",col="gray")+
#   theme_classic()
# 
# test <- wilcox.test(p2g.norm$p2g.distance,p2g.cancer$p2g.distance)
# test <- wilcox.test(p2g.cancer$p2g.distance,p2g.norm$p2g.distance)

##################################################################
# PART 2: Plot proportion of peaks per number of peaks and 
# average number of target genes per cancer and normal groups
####################################################################
# Number of peaks per number of genes
p2g.cancer <- readRDS("./Cancer_specific_P2G_table.rds")
p2g.normal <- readRDS("./All_P2G_Observed.rds")
p2g.normal <- p2g.normal[p2g.normal$Correlation >= 0.45,]
p2g.normal <- p2g.normal[p2g.normal$RawPVal <= 1e-12,]
p2g.normal <- p2g.normal[p2g.normal$peakType == "Distal",]
p2g.normal <- p2g.normal[p2g.normal$peakName %ni% p2g.cancer$peakName,]

df.cancer <- data.frame(num.genes = table(p2g.cancer$peakName))
df.cancer$cat <- ifelse(df.cancer$num.genes.Freq < 3,"1-2",df.cancer$num.genes.Freq)
df.cancer$cat <- ifelse(df.cancer$num.genes.Freq >2,"3+",df.cancer$cat)

df.normal <- data.frame(num.genes = table(p2g.normal$peakName))
df.normal$cat <- ifelse(df.normal$num.genes.Freq < 3,"1-2",df.normal$num.genes.Freq)
df.normal$cat <- ifelse(df.normal$num.genes.Freq >2,"3+",df.normal$cat)

head(df.normal)
head(df.cancer)

df.normal$type <- "normal"
df.cancer$type <- "cancer"

df<- rbind(df.normal,df.cancer)


p1 <- ggplot(df, aes(x=type, y=num.genes.Freq, fill = type)) +
  stat_summary(geom = "bar", fun.y = mean, position = "dodge",width=0.3) +
  stat_summary(geom = "errorbar", fun.data = mean_se, position = "dodge",width=0.15)+
  theme_classic()+  coord_cartesian(ylim=c(1.44,1.65))+NoLegend()

df.cancer %>% 
  count(cat) %>% 
  mutate(perc = (n / nrow(df.cancer)*100)) -> cancer
df.normal %>% 
  count(cat) %>% 
  mutate(perc = (n / nrow(df.normal)*100)) -> normal

cancer$type <- "cancer"
normal$type <- "normal"

comb <- rbind(cancer,normal)
print(comb)
p2 <- ggplot(comb,aes(x=cat,y=perc,fill=type))+
  geom_bar(stat="identity",position = position_dodge(width=0.75),width=0.6)+
  theme_classic()+
  geom_text(aes(label=round(perc,3)), vjust=-0.4, size=3.5,position = position_dodge(width=0.75))+
  scale_y_continuous(expand = c(0,0),limits=c(0,100))+NoLegend()


CombinePlots(list(p2,p1),ncol=2)+ggsave("Barcharts_P2G_plots.pdf",width=6,height=3)

# Cancer specific peaks link to more genes on average with statistical significance:
test <- wilcox.test(df.cancer$num.genes.Freq,df.normal$num.genes.Freq,correct = F)
print(test)
print(test$p.value)
print(paste0("Cancer-specific peaks link to more genes on average (1.57 v. 1.44 genes, p=",test$p.value,")"))

# Proportion of 3+ peaks is greater in cancer-specific peaks relative to normal peaks:
res <- fisher.test(matrix(c(414, 1979, 
                            3274, 20187),
                          ncol = 2))
print(res)
print(res$p.value)
print(paste0("The proportion of cancer-specific peaks linking to 3 or more genes is significantly greater relative to the normal peaks (p=",res$p.value,")"))

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
my_levels <- as.character(c(3,9,10,16,17,
                            11,20,21,22,31,19,34,
                            0,27,
                            6,8,12,14,15,18,24,25,26,29,7,23,
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
# Color rows:

epithelial.cols <- colorRampPalette(c("#A0E989", "#245719"))
epithelial.cols <- epithelial.cols(14)

fibro.cols <-colorRampPalette(c("#FABFD2", "#B1339E"))
fibro.cols <- fibro.cols(10)

smooth.cols <- c("#b47fe5","#d8b7e8")

endo.cols <- c("#93CEFF","#4A99FF")

t.cols <- c("gray60","gray20","gray40")

macro.cols <- c("#ff6600","#ff9d5c")

mast.cols <- "gold3"

b.cols <- c("#B22222","#CD5C5C")


cols <- c(epithelial.cols,fibro.cols,smooth.cols,endo.cols,t.cols,macro.cols,mast.cols,b.cols)
cols <- cols[-c(16:21,23:24,26,30:36)]
cols <- cols[c(8:12,1:7,13:20)]
plot <- plotBrowserTrack(atac.sub,geneSymbol ="RHEB", groupBy = "predictedGroup_ArchR",
                         features = GRangesList(TrackA = encode.all,TrackB = ft.peaks,TrackC = ov.peaks), 
                         loops = getPeak2GeneLinks.mod(atac,corCutOff = 0.45,
                                                   PValCutOff = 1e-12,varCutOffATAC = 0,
                                                   varCutOffRNA = 0),upstream = 6000,downstream = 35000,
                         pal=cols)

pdf("RHEB_final.pdf",width = 6,height = 8)
grid::grid.draw(plot[[1]])
dev.off()


plot <- plotBrowserTrack(atac.sub,geneSymbol ="MUC16", groupBy = "predictedGroup_ArchR",
                         features = GRangesList(TrackA = encode.all,TrackB = ft.peaks,TrackC = ov.peaks), 
                         loops = getPeak2GeneLinks.mod(atac,corCutOff = 0.45,
                                                       PValCutOff = 1e-12,varCutOffATAC = 0,
                                                       varCutOffRNA = 0),upstream = 300000,downstream = 80000,
                         pal=cols)

pdf("MUC16_final.pdf",width = 6,height = 8)
grid::grid.draw(plot[[1]])
dev.off()

names <- gsub(".*_","",atac.sub$predictedGroup_ArchR)
saveRDS(names,"names.rds")
rm(atac.sub)
rm(atac)
####################################################################
# PART 4: plot matching violin plots for RHEB expression and mTOR
####################################################################
names <- readRDS("./names.rds")
rna <- readRDS("./endo_ovar_All_scRNA_processed.rds")
rna.sub <- rna[,rna$cell.type %in% intersect(levels(factor(rna$cell.type)),levels(factor(names)))]

my_levels <- as.character(c(3,9,10,16,17,
                            11,20,21,22,31,19,34,
                            0,27,
                            6,25,7,
                            1,33,
                            2))

# Make violin plots for RHEB expression and mTOR pathway expression
# Relevel object@ident
rna.sub@active.ident <- factor(x =rna.sub$RNA_snn_res.0.7, levels = rev(my_levels))
p1 <- VlnPlot(rna.sub,features = "RHEB",pt.size = 0)+coord_flip()+NoLegend()
p1 <- ggplot(p1$data,aes(y=ident,x=RHEB))+geom_boxplot(aes(fill=ident),lwd=0.45,outlier.size = 0.95,fatten = 0.95)+NoLegend()+
  theme_classic()+NoLegend()+ylab("Cluster #")+xlab("RHEB expression")


# Compute mTOR pathway member expression
gset <- read.delim("./BIOCARTA_MTOR_PATHWAY.txt",header = T)
gset <- as.character(gset$BIOCARTA_MTOR_PATHWAY[2:length(gset$BIOCARTA_MTOR_PATHWAY)])
rna.sub <- AddModuleScore(rna.sub,features = list(gset),name = "mTOR_members",search = T)


p2 <- VlnPlot(rna.sub,features = "mTOR_members1",pt.size = 0)+coord_flip()+NoLegend()
p2 <- ggplot(p2$data,aes(y=ident,x=mTOR_members1))+geom_boxplot(aes(fill=ident),lwd=0.45,outlier.size = 0.95,fatten = 0.95)+NoLegend()+
  theme_classic()+NoLegend()+ylab("Cluster #")+xlab("mTOR member expression")
CombinePlots(list(p1,p2),ncol=2)+ggsave("VlnPlots.pdf",width = 6,height = 8)

# Differential enrichment of mTOR pathway members
res <- kruskal.test(data=p2$data,mTOR_members1 ~ ident)
print(res)
print(res$p.value)
# Differential expression of RHEB in cluster 3
grouped.markers <- FindMarkers(rna.sub,ident.1 = "3",
                               ident.2 = as.character(c(9,10,16,17,
                                           11,20,21,22,31,19,34,
                                           0,27,
                                           6,25,7,
                                           1,33,
                                           2)),only.pos = T)
grouped.markers$gene <- rownames(grouped.markers)
grouped.markers.RHEB <- grouped.markers[grouped.markers$gene == "RHEB",]
print(grouped.markers.RHEB)# Significant up-regulation in cluster 3

writeLines(capture.output(sessionInfo()), "sessionInfo.txt")
