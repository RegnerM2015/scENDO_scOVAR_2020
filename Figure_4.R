###########################################################
# Matt Regner
# Franco Lab
# Date: April-December 2020
#
# Description: Figure 4, HGSOC enhancer discovery 
# RNA and ATAC data
###########################################################
source("./P2G_Heatmap_Distal.R")
library(ggplot2)
library(Seurat)
library(scales)
library(forcats)
library(dplyr)
library(tidyverse)
library(tidyr)

# Make patient sample metadata and color assignments 

sampleColors <- RColorBrewer::brewer.pal(11,"Paired")
sampleColors[11] <- "#8c8b8b"
pie(rep(1,11), col=sampleColors) 


# Color patient tumors to resemble the cancer ribbon color 
sampleColors <- c(sampleColors[5],sampleColors[7],sampleColors[6],sampleColors[8],sampleColors[10],sampleColors[9],sampleColors[4],sampleColors[3],sampleColors[2],sampleColors[11],sampleColors[1])



sampleAnnot <- data.frame(Sample = c("3533EL","3571DL","36186L","36639L",
                                     "366C5L","37EACL","38FE7L","3BAE2L","3CCF1L","3E4D1L","3E5CFL"),
                          Color = sampleColors,
                          Cancer = c("endometrial","endometrial","endometrial","endometrial","endometrial","endometrial",
                                     "ovarian","ovarian","ovarian","ovarian","ovarian"),
                          Histology = c("endometrioid","endometrioid","endometrioid","endometrioid","endometrioid",
                                        "serous","endometrioid","serous","carcinosarcoma","GIST","serous"),
                          BMI = c(39.89,30.5,38.55,55.29,49.44,29.94,34.8,22.13,23.72,33.96,22.37),
                          Age = c(70,70,70,49,62,74,76,61,69,59,59),
                          Race = c("AA","CAU","CAU","CAU","CAU","CAU","CAU","CAU","CAU","CAU","AS"),
                          Stage = c("IA","IA","IA","IA","IA","IIIA","IA","IIB","IVB","IV","IIIC"),
                          Site = c("Endometrium","Endometrium","Endometrium","Endometrium","Endometrium","Ovary","Ovary","Ovary","Ovary","Ovary","Ovary"),
                          Type = c("Endometrial","Endometrial","Endometrial","Endometrial","Endometrial","Endometrial","Ovarian","Ovarian","Ovarian","Gastric","Ovarian"))

sampleColors[11] <- "#8fb2c4"# Make blue slightly darker for Figure 4

###########################################################
# Read in RNA data for full cohort
rna <- readRDS("ovar_HGSOC_scRNA_processed.rds")


# Read in matching ATAC data for full cohort (ArchR project after label transfer)
atac <- readRDS("proj_LSI_GeneScores_Annotations_Int.rds")
###########################################################

# Plot Patient UMAPs for RNA/ATAC
###########################################################
# RNA UMAP first
rna.df <- as.data.frame(rna@reductions$umap@cell.embeddings)
length(which(rownames(rna.df)==rownames(rna@meta.data)))
rna.df$Sample <- rna$Sample

rna.sample.plot <-ggplot(rna.df,aes(x = UMAP_1,y=UMAP_2,color = Sample))+
  geom_point(size = .1)+
  theme_classic()+
  theme(plot.title = element_text(face = "bold"))+
  xlab("UMAP_1")+
  ylab("UMAP_2")+ 
  theme(legend.key.size = unit(0.2, "cm"))+
  scale_color_manual(values = sampleColors[c(8,11)])+
  guides(colour = guide_legend(override.aes = list(size=6)))
rna.sample.plot +ggsave("Sample_RNA.pdf",width = 8,height = 6)

# ATAC UMAP second:
atac.df <- plotEmbedding(atac,colorBy = "cellColData",name = "Sample",embedding = "UMAP")
atac.df <- as.data.frame(atac.df$data)
atac.df$Sample <- gsub(".*-","",atac.df$color)

ggplot(atac.df,aes_string(x = "x",y="y",color = "Sample"))+
  geom_point(size = .1)+
  theme_classic()+
  theme(plot.title = element_text(face = "bold"))+
  xlab("UMAP_1")+
  ylab("UMAP_2")+
  theme(legend.key.size = unit(0.2, "cm"))+
  scale_color_manual(values = sampleColors[c(8,11)])+
  guides(colour = guide_legend(override.aes = list(size=6)))+
  ggsave(paste0("Sample_ATAC.pdf"),width = 8,height = 6)
###########################################################




# Plot cell type UMAPs for RNA/ATAC
###########################################################
# RNA UMAP first
rna.df <- as.data.frame(rna@reductions$umap@cell.embeddings)
length(which(rownames(rna.df)==rownames(rna@meta.data)))
rna.df$cell.type <- rna$cell.type

rna.df$cell.type <- factor(rna.df$cell.type,levels = c("0-Epithelial cell",
                                                       "2-Epithelial cell",
                                                       "3-Epithelial cell",
                                                       "7-Epithelial cell",
                                                       "11-Epithelial cell",
                                                       "15-Empty/Epithelial cell",
                                                       "16-Epithelial cell",
                                                       "19-Epithelial cell",
                                                       "1-Fibroblast",
                                                       "4-Fibroblast",
                                                       "9-Fibroblast",
                                                       "10-Fibroblast",
                                                       "22-Fibroblast",
                                                       "14-Endothelial cell",
                                                       "5-T cell",
                                                       "12-NK cell",
                                                       "18-T cell",
                                                       "21-T cell",
                                                       "6-Macrophage",
                                                       "8-Macrophage",
                                                       "13-Macrophage" ,
                                                       "20-Macrophage",
                                                       "23-Macrophage","17-B cell" 
                                                       ))


epithelial.cols <- colorRampPalette(c("#A0E989", "#337B24"))
epithelial.cols <- epithelial.cols(8)

fibro.cols <- c("#FB8CAB","#E65C9C","#CF268A","#Af1281","#6B0772")

endo.cols <- "#377EB8"

t.cols <- c("#CCCCCC","#333333","#666666","#999999")

macro.cols <- colorRampPalette(c("#FECC51", "#FD6A02"))
macro.cols <- macro.cols(5)


b.cols <- "#E41A1C"


cols <- c(epithelial.cols,fibro.cols,endo.cols,t.cols,macro.cols,b.cols)

rna.sample.plot <-ggplot(rna.df,aes(x = UMAP_1,y=UMAP_2,color =cell.type))+
  geom_point(size = .1)+
  theme_classic()+
  theme(plot.title = element_text(face = "bold"))+
  xlab("UMAP_1")+
  ylab("UMAP_2")+ 
  theme(legend.key.size = unit(0.2, "cm"))+
  scale_color_manual(values = cols)+
  guides(colour = guide_legend(override.aes = list(size=3)))
rna.sample.plot +ggsave("Cell_type_RNA.pdf",width = 8,height = 6)

# ATAC UMAP second:
atac.df <- plotEmbedding(atac,colorBy = "cellColData",name = "predictedGroup_ArchR",embedding = "UMAP")
atac.df <- as.data.frame(atac.df$data)
atac.df$cell.type <- sub(".*?-","",atac.df$color)


atac.df$cell.type <- factor(atac.df$cell.type,levels = c("0-Epithelial cell",
                                                       "2-Epithelial cell",
                                                       "3-Epithelial cell",
                                                       "7-Epithelial cell",
                                                       "11-Epithelial cell",
                                                       "15-Empty_Epithelial cell",
                                                       "16-Epithelial cell",
                                                       "19-Epithelial cell",
                                                       "1-Fibroblast",
                                                       "4-Fibroblast",
                                                       "9-Fibroblast",
                                                       "10-Fibroblast",
                                                       "22-Fibroblast",
                                                       "14-Endothelial cell",
                                                       "5-T cell",
                                                       "12-NK cell",
                                                       "18-T cell",
                                                       "21-T cell",
                                                       "6-Macrophage",
                                                       "8-Macrophage",
                                                       "13-Macrophage" ,
                                                       "20-Macrophage",
                                                       "23-Macrophage","17-B cell" 
))



epithelial.cols <- colorRampPalette(c("#A0E989", "#337B24"))
epithelial.cols <- epithelial.cols(8)

fibro.cols <- c("#FB8CAB","#E65C9C","#CF268A","#Af1281","#6B0772")

endo.cols <- "#377EB8"

t.cols <- c("#CCCCCC","#333333","#666666","#999999")

macro.cols <- colorRampPalette(c("#FECC51", "#FD6A02"))
macro.cols <- macro.cols(5)


b.cols <- "#E41A1C"


cols <- c(epithelial.cols,fibro.cols,endo.cols,t.cols,macro.cols,b.cols)


ggplot(atac.df,aes_string(x = "x",y="y",color = "cell.type"))+
  geom_point(size = .1)+
  theme_classic()+
  theme(plot.title = element_text(face = "bold"))+
  xlab("UMAP_1")+
  ylab("UMAP_2")+
  theme(legend.key.size = unit(0.2, "cm"))+
  scale_color_manual(values = cols)+
  guides(colour = guide_legend(override.aes = list(size=3)))+
  ggsave(paste0("Cell_type_ATAC.pdf"),width = 8,height = 6)
###########################################################


# Patient proportion per subcluster in RNA:
###########################################################
meta <- rna@meta.data

df <- meta %>% dplyr::group_by(RNA_snn_res.0.7) %>% dplyr::count(Sample)
colnames(df) <- c("Cluster","Sample","Cells")

# Reorder cluster factor levels to group by cell type 
levels(factor(rna$cell.type))
df %>% 
  dplyr::mutate(cell.type = factor(Cluster,levels = c("0","2","3","7","11","15","16","19",
                                                      "1","4","9","10","22",
                                                      "14",
                                                      "5","12","18","21",
                                                      "6","8","13","20","23",
                                                      "17"))) %>% 
  ggplot(aes(fill=Sample, y=Cells, x= fct_rev(cell.type))) + 
  scale_fill_manual(values = sampleColors[c(8,11)])+
  geom_bar(position="fill", stat="identity")+
  coord_flip()+theme_classic()+xlab("Clusters")+ylab("# of cells")+
  ggsave("Cell_Type_Prop_RNA.pdf",width = 4,height = 11)


# Patient proportion per subcluster in ATAC:
meta <- as.data.frame(atac@cellColData)
meta$predictedGroup_ArchR <- gsub("-.*", "", meta$predictedGroup_ArchR)

df <- meta %>% dplyr::group_by(predictedGroup_ArchR) %>% dplyr::count(Sample)
colnames(df) <- c("Cluster","Sample","Cells")

# Reorder cluster factor levels to group by cell type 
levels(factor(atac$predictedGroup_ArchR))
df %>% 
  dplyr::mutate(cell.type = factor(Cluster,levels = c("0","2","3","7","11","15","16","19",
                                                      "1","4","9","10","22",
                                                      "14",
                                                      "5","12","18","21",
                                                      "6","8","13","20","23",
                                                      "17"))) %>% 
  ggplot(aes(fill=Sample, y=Cells, x= fct_rev(cell.type))) + 
  scale_fill_manual(values = sampleColors[c(8,11)])+
  geom_bar(position="fill", stat="identity")+
  coord_flip()+theme_classic()+xlab("Clusters")+ylab("# of cells")+
  ggsave("Cell_Type_Prop_ATAC.pdf",width = 4,height = 11)
###########################################################




# Immune cell proportion barcharts
###########################################################
# ATAC 
meta.atac <- as.data.frame(atac@cellColData)

meta.atac$predictedGroup_ArchR <- gsub(".*-","",meta.atac$predictedGroup_ArchR)


meta.atac$immune <- ifelse(meta.atac$predictedGroup_ArchR == "T cell" |
                             meta.atac$predictedGroup_ArchR == "B cell"|
                             meta.atac$predictedGroup_ArchR == "NK cell" |
                             meta.atac$predictedGroup_ArchR == "Macrophage","Immune","Non-immune")


BAE2L.prop <- nrow(meta.atac[meta.atac$immune == "Immune" & meta.atac$Sample == '3BAE2L',])/nrow(meta.atac[meta.atac$Sample == '3BAE2L',])
E5CFL.prop <- nrow(meta.atac[meta.atac$immune == "Immune" & meta.atac$Sample == '3E5CFL',])/nrow(meta.atac[meta.atac$Sample == '3E5CFL',])

df <- data.frame(Sample = c("3BAE2L","3E5CFL"),
                 Proportion = c(BAE2L.prop,E5CFL.prop))

ggplot(df,aes(x= Sample,y = Proportion,fill = Sample))+
  geom_bar(stat = "identity")+
  scale_y_continuous(labels = scales::percent,limits = c(0,0.4))+
  scale_fill_manual(values = sampleColors[c(8,11)])+
  ylab("% immune cells")+
  theme_classic()+ggsave("Immune_Prop_ATAC.pdf",width = 3,height = 4)
  
# RNA
meta.rna <- as.data.frame(rna@meta.data)

meta.rna$cell.type <- gsub(".*-","",meta.rna$cell.type)


meta.rna$immune <- ifelse(meta.rna$cell.type == "T cell" |
                            meta.rna$cell.type == "B cell"|
                            meta.rna$cell.type == "NK cell" |
                            meta.rna$cell.type == "Macrophage","Immune","Non-immune")


BAE2L.prop <- nrow(meta.rna[meta.rna$immune == "Immune" & meta.rna$Sample == '3BAE2L',])/nrow(meta.rna[meta.rna$Sample == '3BAE2L',])
E5CFL.prop <- nrow(meta.rna[meta.rna$immune == "Immune" & meta.rna$Sample == '3E5CFL',])/nrow(meta.rna[meta.rna$Sample == '3E5CFL',])

df <- data.frame(Sample = c("3BAE2L","3E5CFL"),
                 Proportion = c(BAE2L.prop,E5CFL.prop))

ggplot(df,aes(x= Sample,y = Proportion,fill = Sample))+
  geom_bar(stat = "identity")+
  scale_y_continuous(labels = scales::percent,limits = c(0,0.4))+
  scale_fill_manual(values = sampleColors[c(8,11)])+
  ylab("% immune cells")+
  theme_classic()+ggsave("Immune_Prop_RNA.pdf",width = 3,height = 4)
###########################################################


# Write cluster total # of cells to output files
###########################################################
total.atac <- as.data.frame(table(atac$predictedGroup_ArchR))
colnames(total.atac) <- c("Cluster","ATAC cells")

total.rna <- as.data.frame(table(rna$cell.type))
colnames(total.rna) <- c("Cluster","RNA cells")
total.rna$Cluster <- gsub("/","_",total.rna$Cluster)

total <- merge(total.rna,total.atac,by = "Cluster")
WriteXLS::WriteXLS(total,"./Total_Cell_Numbers_per_Cluster.xlsx")
###########################################################



# Plot peak2gene heatmap
###########################################################
atac.peaks <- readRDS("./final_archr_proj_archrGS.rds")
p2g.HGSOC <- getPeak2GeneLinks(ArchRProj = atac.peaks,returnLoops = F,corCutOff = 0,FDRCutOff = 1)


hist(p2g.HGSOC$EmpFDR)
summary(p2g.HGSOC$EmpFDR)
hist(p2g.HGSOC$FDR)
summary(p2g.HGSOC$FDR)
hist(p2g.HGSOC$Correlation)
summary(p2g.HGSOC$Correlation)


#####################################################################################
# Check to see if Peak to Gene links overlap with statistically significant P2G links 
# in full cohort by EmpFDR:
##################################################################################
p2g.coords <- metadata(p2g.HGSOC)[[1]]
ranges <- as.data.frame(p2g.coords@ranges)

HGSOC <- data.frame(chrom = p2g.coords@seqnames , start = ranges$start,end = ranges$end)

sig.distal <- readRDS("./Peak2GeneLinks_EmpFDR_Ranked_Distal_ArchR.rds")
sig.distal <- dplyr::filter(sig.distal,EmpFDR <= 0.05)

full <- data.frame(chrom = sig.distal$seqnames.x,start = sig.distal$start.x,end = sig.distal$end.x)
full.gr <- makeGRangesFromDataFrame(full)

# Use bedtoolsr to intersect Peak2Gene links with full cohort Peak2Gene links:
intersect <- bedtoolsr::bt.intersect(HGSOC,full.gr,wa = T)
colnames(intersect) <- c("seqnames","start","end")

empirical.sig.HGSOC.peaks.from.full.cohort <- makeGRangesFromDataFrame(intersect)
###################################################################################

# Subset colors to new cell cluster labels
names(cols) <- levels(atac.df$cell.type)
cols <- cols[names(cols) %in% levels(factor(atac.peaks$predictedGroup_ArchR))]

p2g <- plotPeak2GeneHeatmap.distal(atac.peaks,
                            corCutOff = 0.45,
                            FDRCutOff = 0.0001,
                            groupBy = "predictedGroup_ArchR",
                            k=length(cols),palGroup = cols,returnMatrices = F,nPlot = 25000,
                            palRNA =  paletteContinuous("solarExtra"),
                            palATAC =  paletteContinuous("solarExtra"))

pdf("Peak2Gene_Heatmap_Legend.pdf",width = 14,height = 12)
draw(p2g, heatmap_legend_side = "right", annotation_legend_side = "right")
dev.off()

pdf("Peak2Gene_Heatmap.pdf",width = 8,height = 10)
draw(p2g, heatmap_legend_side = "bot", annotation_legend_side = "bot")
dev.off()



# For every gene in user defined list, plot P2G linkage browser track with 
# 1) All ENCODE cCREs
# 2) Distal ENCODE cCREs
# 3) HGSOC statistically significant peaks that overlap with EmpFDR significant peaks from full cohort
# 4) HGSOC statistically signiicant peaks (EmpFDR is calculated within HGSOC cohort only)
# Next, plot matching scRNA-seq expression in violin plot
########################################################################################################

# Read in ENCODE cCRE information:
encode.all <- read.delim("./GRCh38-ccREs.bed",header =F)
colnames(encode.all)[1:3] <- c("seqnames","start","end")
encode.all <- makeGRangesFromDataFrame(encode.all)

encode.distal <- read.delim("./GRCh38-ccREs.dELS.bed",header =F)
colnames(encode.distal)[1:3] <- c("seqnames","start","end")
encode.distal <- makeGRangesFromDataFrame(encode.distal)

# Read in ENCODE DN-ase epithelium 
encode.epithelium <- read.delim("./ENCODE_Epithelium_DNase.bed",header = T,sep = "\t")
colnames(encode.epithelium)[1:3] <- c("seqnames","start","end")
encode.epithelium <- makeGRangesFromDataFrame(encode.epithelium)


# Read in all P2G peaks 
# Find HGSOC statistically signiicant peaks (EmpFDR is calculated within HGSOC cohort only)
p2g <- plotPeak2GeneHeatmap.distal(atac.peaks,
                            corCutOff = 0.45,
                            FDRCutOff = 0.0001,
                            groupBy = "predictedGroup_ArchR",
                            k=length(cols),palGroup = cols,returnMatrices = T,nPlot = 500000)# Make Heatmap object with ALL P2Gs!

mat <- p2g$RNA$matrix
colnames(mat) <- make.unique(p2g$RNA$colData$groupBy)
rownames(mat) <- make.unique(p2g$Peak2GeneLinks$gene)

kmeans <- p2g$RNA$kmeansId

# list genes of interest:

genes.of.interest <- c("LAPTM4B","KLK6","RHOD","FOLR1","GAL")

for (i in genes.of.interest){
  
  # Save P2G peaknames and Kmeans cluster for the gene being iterated
  idx <- grep(i,rownames(mat))
  kmeans.idx <- kmeans[idx]
  
  peaks.genes <- p2g$Peak2GeneLinks
  peaks.genes <- peaks.genes[idx,]
  
  peaks.genes <- as.data.frame(peaks.genes)
  peaks.genes$kmeans.group <- kmeans.idx
  
  saveRDS(peaks.genes,paste0("P2G_Hits_",i,"_data.rds"))
  #peaks.genes <- dplyr::filter(peaks.genes,EmpFDR <= 0.05)
  
  sig.peaks <- as.data.frame(peaks.genes$peak)
  sig.peaks$start <- sig.peaks$`peaks.genes$peak`
  sig.peaks$end <- sig.peaks$`peaks.genes$peak`
  colnames(sig.peaks) <- c("seqnames","start","end")
  
  sig.peaks$seqnames <- gsub(":.*","",sig.peaks$seqnames)
  sig.peaks$start <- gsub(".*:","",sig.peaks$start)
  sig.peaks <- sig.peaks[,-3]
  
  sig.peaks <- tidyr::separate(data = sig.peaks, col = start, into = c("start", "end"), sep = "\\-")
  
  sig.peaks <- makeGRangesFromDataFrame(sig.peaks)
  ###########################################################
  
  
  # Plot Peak2Gene links in browser track with the following annotations:
  # 1) All ENCODE cCREs
  # 2) Distal ENCODE cCREs
  # 3) HGSOC statistically significant peaks that overlap with EmpFDR significant peaks from full cohort
  # 4) HGSOC statistically signiicant peaks (EmpFDR is calculated within HGSOC cohort only)
  ###########################################################
  
  # Annotation 1:
  plot <- plotBrowserTrack(atac.peaks,geneSymbol =i, groupBy = "predictedGroup_ArchR",
                           features = GRangesList(TrackA = encode.all,TrackB = encode.distal), 
                           pal = cols,
                           loops = getPeak2GeneLinks(atac.peaks,corCutOff = 0.45,
                                                     FDRCutOff = 0.0001))
  
  pdf(paste0("Peak2Gene_",i,"_cCREs_plus_Distal.pdf"),width = 6,height = 8)
  grid::grid.draw(plot[[1]])
  dev.off()
  
  
  # Annotation 2:
  plot <- plotBrowserTrack(atac.peaks,geneSymbol =i, groupBy = "predictedGroup_ArchR",
                           features = GRangesList(TrackA = encode.all,TrackB = encode.distal,TrackC = getPeakSet(atac.peaks),TrackD = empirical.sig.HGSOC.peaks.from.full.cohort,TrackE=sig.peaks), 
                           pal = cols,
                           loops = getPeak2GeneLinks(atac.peaks,corCutOff = 0.45,
                                                     FDRCutOff = 0.0001))
  
  pdf(paste0("Peak2Gene_",i,"_cCREs_plus_Distal_plus_HGSOC_sig_p2gs.pdf"),width = 6,height = 8)
  grid::grid.draw(plot[[1]])
  dev.off()
  
  
  # Finally plot matching scRNA-seq expression in violin plots:
  
  # Plot violin plots for matching scRNA-seq expression:
  ###########################################################
  Idents(rna) <- "RNA_snn_res.0.7"
  
  my_levels <- rev(c(0,2,3,7,11,16,
                     1,4,9,10,
                     14,
                     5,12,18,21,
                     6,8,13,
                     17,15,19,22,20,23
  ))
  
  # Relevel object@ident
  rna$cluster.new <- factor(x = rna$RNA_snn_res.0.7, levels = my_levels)
  Idents(rna) <- "cluster.new"
  
  VlnPlot(rna,features = i,group.by = "cluster.new",pt.size = 0,
          idents = c(0,2,3,7,11,16,
                     1,4,9,10,
                     14,
                     5,12,18,21,
                     6,8,13,
                     17
          ))+coord_flip()+NoLegend()+ggsave(paste0(i,"_vln.pdf"),width = 4,height = 8)
  
  
  
  
}






# Follow up and customize genomic coordinate view for each:
#######################################################

peaks.genes <- readRDS("./P2G_Hits_LAPTM4B_data.rds")
sig.peaks <- as.data.frame(peaks.genes$peak)
sig.peaks$start <- sig.peaks$`peaks.genes$peak`
sig.peaks$end <- sig.peaks$`peaks.genes$peak`
colnames(sig.peaks) <- c("seqnames","start","end")

sig.peaks$seqnames <- gsub(":.*","",sig.peaks$seqnames)
sig.peaks$start <- gsub(".*:","",sig.peaks$start)
sig.peaks <- sig.peaks[,-3]

sig.peaks <- tidyr::separate(data = sig.peaks, col = start, into = c("start", "end"), sep = "\\-")

sig.peaks <- makeGRangesFromDataFrame(sig.peaks)
plot <- plotBrowserTrack(atac.peaks,geneSymbol ="LAPTM4B", groupBy = "predictedGroup_ArchR",
                         features = GRangesList(TrackA = encode.all,TrackB = encode.distal,TrackC=encode.epithelium), 
                         pal = cols,upstream = 45000,downstream = 50000,
                         loops = getPeak2GeneLinks(atac.peaks,corCutOff = 0.45,
                                                   FDRCutOff = 0.0001))

pdf("LAPTM4B_final.pdf",width = 6,height = 8)
grid::grid.draw(plot[[1]])
dev.off()


# Annotation 2:
plot <- plotBrowserTrack(atac.peaks,geneSymbol ="LAPTM4B", groupBy = "predictedGroup_ArchR",
                         features = GRangesList(TrackA = encode.all,TrackB = encode.distal,TrackC = getPeakSet(atac.peaks),TrackD = empirical.sig.HGSOC.peaks.from.full.cohort,TrackE=sig.peaks), 
                         pal = cols,upstream = 45000,downstream = 50000,
                         loops = getPeak2GeneLinks(atac.peaks,corCutOff = 0.45,
                                                   FDRCutOff = 0.0001))

pdf("LAPTM4B_final_sig_peaks_anno.pdf",width = 6,height = 8)
grid::grid.draw(plot[[1]])
dev.off()




peaks.genes <- readRDS("./P2G_Hits_KLK11_data.rds")
sig.peaks <- as.data.frame(peaks.genes$peak)
sig.peaks$start <- sig.peaks$`peaks.genes$peak`
sig.peaks$end <- sig.peaks$`peaks.genes$peak`
colnames(sig.peaks) <- c("seqnames","start","end")

sig.peaks$seqnames <- gsub(":.*","",sig.peaks$seqnames)
sig.peaks$start <- gsub(".*:","",sig.peaks$start)
sig.peaks <- sig.peaks[,-3]

sig.peaks <- tidyr::separate(data = sig.peaks, col = start, into = c("start", "end"), sep = "\\-")

sig.peaks <- makeGRangesFromDataFrame(sig.peaks)
plot <- plotBrowserTrack(atac.peaks,geneSymbol ="KLK11", groupBy = "predictedGroup_ArchR",
                         features = GRangesList(TrackA = encode.all,TrackB = encode.distal,TrackC=encode.epithelium), 
                         pal = cols,upstream = 112000,downstream = 30000,
                         loops = getPeak2GeneLinks(atac.peaks,corCutOff = 0.45,
                                                   FDRCutOff = 0.0001,resolution = 6))

pdf("KLK_locus_final.pdf",width = 6,height = 8)
grid::grid.draw(plot[[1]])
dev.off()


# Annotation 2:
plot <- plotBrowserTrack(atac.peaks,geneSymbol ="KLK11", groupBy = "predictedGroup_ArchR",
                         features = GRangesList(TrackA = encode.all,TrackB = encode.distal,TrackC = getPeakSet(atac.peaks),TrackD = empirical.sig.HGSOC.peaks.from.full.cohort,TrackE=sig.peaks), 
                         pal = cols,upstream = 112000,downstream =30000,
                         loops = getPeak2GeneLinks(atac.peaks,corCutOff = 0.45,
                                                   FDRCutOff = 0.0001,resolution = 6))

pdf("KLK_locus_final_sig_peaks_anno.pdf",width = 6,height = 8)
grid::grid.draw(plot[[1]])
dev.off()



peaks.genes <- readRDS("./P2G_Hits_RHOD_data.rds")
sig.peaks <- as.data.frame(peaks.genes$peak)
sig.peaks$start <- sig.peaks$`peaks.genes$peak`
sig.peaks$end <- sig.peaks$`peaks.genes$peak`
colnames(sig.peaks) <- c("seqnames","start","end")

sig.peaks$seqnames <- gsub(":.*","",sig.peaks$seqnames)
sig.peaks$start <- gsub(".*:","",sig.peaks$start)
sig.peaks <- sig.peaks[,-3]

sig.peaks <- tidyr::separate(data = sig.peaks, col = start, into = c("start", "end"), sep = "\\-")

sig.peaks <- makeGRangesFromDataFrame(sig.peaks)
plot <- plotBrowserTrack(atac.peaks,geneSymbol ="RHOD", groupBy = "predictedGroup_ArchR",
                         features = GRangesList(TrackA = encode.all,TrackB = encode.distal,TrackC=encode.epithelium),  
                         pal = cols,upstream = 10000,downstream = 45000,
                         loops = getPeak2GeneLinks(atac.peaks,corCutOff = 0.45,
                                                   FDRCutOff = 0.0001))

pdf("RHOD_final.pdf",width = 6,height = 8)
grid::grid.draw(plot[[1]])
dev.off()


# Annotation 2:
plot <- plotBrowserTrack(atac.peaks,geneSymbol ="RHOD", groupBy = "predictedGroup_ArchR",
                         features = GRangesList(TrackA = encode.all,TrackB = encode.distal,TrackC = getPeakSet(atac.peaks),TrackD = empirical.sig.HGSOC.peaks.from.full.cohort,TrackE=sig.peaks), 
                         pal = cols,upstream = 10000,downstream = 45000,
                         loops = getPeak2GeneLinks(atac.peaks,corCutOff = 0.45,
                                                   FDRCutOff = 0.0001))

pdf("RHOD_final_sig_peaks_anno.pdf",width = 6,height = 8)
grid::grid.draw(plot[[1]])
dev.off()




peaks.genes <- readRDS("./P2G_Hits_FOLR1_data.rds")
sig.peaks <- as.data.frame(peaks.genes$peak)
sig.peaks$start <- sig.peaks$`peaks.genes$peak`
sig.peaks$end <- sig.peaks$`peaks.genes$peak`
colnames(sig.peaks) <- c("seqnames","start","end")

sig.peaks$seqnames <- gsub(":.*","",sig.peaks$seqnames)
sig.peaks$start <- gsub(".*:","",sig.peaks$start)
sig.peaks <- sig.peaks[,-3]

sig.peaks <- tidyr::separate(data = sig.peaks, col = start, into = c("start", "end"), sep = "\\-")

sig.peaks <- makeGRangesFromDataFrame(sig.peaks)
plot <- plotBrowserTrack(atac.peaks,geneSymbol ="FOLR1", groupBy = "predictedGroup_ArchR",
                         features = GRangesList(TrackA = encode.all,TrackB = encode.distal,TrackC=encode.epithelium), 
                         pal = cols,upstream = 60000,downstream = 30000,
                         loops = getPeak2GeneLinks(atac.peaks,corCutOff = 0.45,
                                                   FDRCutOff = 0.0001))

pdf("FOLR1_final.pdf",width = 6,height = 8)
grid::grid.draw(plot[[1]])
dev.off()


# Annotation 2:
plot <- plotBrowserTrack(atac.peaks,geneSymbol ="FOLR1", groupBy = "predictedGroup_ArchR",
                         features = GRangesList(TrackA = encode.all,TrackB = encode.distal,TrackC = getPeakSet(atac.peaks),TrackD = empirical.sig.HGSOC.peaks.from.full.cohort,TrackE=sig.peaks), 
                         pal = cols,upstream = 60000,downstream = 30000,
                         loops = getPeak2GeneLinks(atac.peaks,corCutOff = 0.45,
                                                   FDRCutOff = 0.0001))

pdf("FOLR1_final_sig_peaks_anno.pdf",width = 6,height = 8)
grid::grid.draw(plot[[1]])
dev.off()



peaks.genes <- readRDS("./P2G_Hits_GAL_data.rds")
sig.peaks <- as.data.frame(peaks.genes$peak)
sig.peaks$start <- sig.peaks$`peaks.genes$peak`
sig.peaks$end <- sig.peaks$`peaks.genes$peak`
colnames(sig.peaks) <- c("seqnames","start","end")

sig.peaks$seqnames <- gsub(":.*","",sig.peaks$seqnames)
sig.peaks$start <- gsub(".*:","",sig.peaks$start)
sig.peaks <- sig.peaks[,-3]

sig.peaks <- tidyr::separate(data = sig.peaks, col = start, into = c("start", "end"), sep = "\\-")

sig.peaks <- makeGRangesFromDataFrame(sig.peaks)
plot <- plotBrowserTrack(atac.peaks,geneSymbol ="GAL", groupBy = "predictedGroup_ArchR",
                         features = GRangesList(TrackA = encode.all,TrackB = encode.distal,TrackC=encode.epithelium), 
                         pal = cols,upstream = 60000,downstream = 25000,
                         loops = getPeak2GeneLinks(atac.peaks,corCutOff = 0.45,
                                                   FDRCutOff = 0.0001))

pdf("GAL_final.pdf",width = 6,height = 8)
grid::grid.draw(plot[[1]])
dev.off()


# Annotation 2:
plot <- plotBrowserTrack(atac.peaks,geneSymbol ="GAL", groupBy = "predictedGroup_ArchR",
                         features = GRangesList(TrackA = encode.all,TrackB = encode.distal,TrackC = getPeakSet(atac.peaks),TrackD = empirical.sig.HGSOC.peaks.from.full.cohort,TrackE=sig.peaks), 
                         pal = cols,upstream = 60000,downstream = 25000,
                         loops = getPeak2GeneLinks(atac.peaks,corCutOff = 0.45,
                                                   FDRCutOff = 0.0001))

pdf("GAL_final_sig_peaks_anno.pdf",width = 6,height = 8)
grid::grid.draw(plot[[1]])
dev.off()



rna@active.ident <- factor(rna$RNA_snn_res.0.7,levels = rev(c("0","2","3","7","11","15","16","19",
                          "1","4","9","10","22",
                          "14",
                          "5","12","18","21",
                          "6","8","13","20","23",
                          "17")))

VlnPlot(rna,features = "Total_CNVs",pt.size = 0)+coord_flip()+NoLegend()+ggsave("CNV_vln_new.pdf",width = 4,height = 8)



# 
# Idents(rna) <- "RNA_snn_res.0.7"
# 
# my_levels <- rev(c(0,2,3,7,11,16,
#                    1,4,9,10,
#                    14,
#                    5,12,18,21,
#                    6,8,13,
#                    17,15,19,22,20,23
# ))
# 
# # Relevel object@ident
# rna$cluster.new <- factor(x = rna$RNA_snn_res.0.7, levels = my_levels)
# Idents(rna) <- "cluster.new"
# 
# VlnPlot(rna,features = "KLK5",group.by = "cluster.new",pt.size = 0,
#         idents = c(0,2,3,7,11,16,
#                    1,4,9,10,
#                    14,
#                    5,12,18,21,
#                    6,8,13,
#                    17))+coord_flip()+NoLegend()+ggsave("KLK5.pdf",width = 4,height = 8)
# 
# VlnPlot(rna,features = "KLK6",group.by = "cluster.new",pt.size = 0,
#         idents = c(0,2,3,7,11,16,
#                    1,4,9,10,
#                    14,
#                    5,12,18,21,
#                    6,8,13,
#                    17))+coord_flip()+NoLegend()+ggsave("KLK6.pdf",width = 4,height = 8)
# 
#                            