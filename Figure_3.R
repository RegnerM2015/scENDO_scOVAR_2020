###########################################################
# Matt Regner
# Franco Lab
# Date: April-December 2020
#
# Description: Figure 3, EEC enhancer discovery 
# RNA and ATAC data
###########################################################
source("./P2G_Heatmap_Distal.R")
library(ggplot2)
library(Seurat)
library(scales)
library(forcats)
library(bedtoolsr)
library(dplyr)

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


###########################################################
# Read in RNA data for full cohort
rna <- readRDS("endo_EEC_scRNA_processed.rds")


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
  scale_color_manual(values = sampleColors[1:5])+
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
  scale_color_manual(values = sampleColors[1:5])+
  guides(colour = guide_legend(override.aes = list(size=6)))+
  ggsave(paste0("Sample_ATAC.pdf"),width = 8,height = 6)
###########################################################




# Plot cell type UMAPs for RNA/ATAC
###########################################################
# RNA UMAP first
rna.df <- as.data.frame(rna@reductions$umap@cell.embeddings)
length(which(rownames(rna.df)==rownames(rna@meta.data)))
rna.df$cell.type <- rna$cell.type

rna.df$cell.type <- factor(rna.df$cell.type,levels = c("6-Unciliated epithelia 1",
                                                       "8-Unciliated epithelia 1",
                                                       "10-Unciliated epithelia 2",
                                                       "14-Ciliated",
                                                       "15-Unciliated epithelia 1",
                                                       "19-Ciliated",
                                                       "21-Unciliated epithelia 2",
                                                       "22-Unciliated epithelia 1",
                                                       "1-Stromal fibroblasts",
                                                       "2-Stromal fibroblasts",
                                                       "9-Stromal fibroblasts",
                                                       "12-Stromal fibroblasts",
                                                       "13-Stromal fibroblasts",
                                                       "16-Stromal fibroblasts",
                                                       "17-Stromal fibroblasts",
                                                       "20-Stromal fibroblasts","23-Stromal fibroblasts",
                                                       "4-Smooth muscle cells" ,"5-Smooth muscle cells",
                                                       "3-Endothelia","11-Endothelia","26-Endothelia",
                                                       "0-Lymphocytes","18-Lymphocytes",
                                                       "24-Macrophages","7-Macrophages",
                                                       "25-Mast cell",
                                                       "27-B cell","28-B cell"
))


epithelial.cols <- colorRampPalette(c("#A0E989", "#337B24"))
epithelial.cols <- epithelial.cols(8)

fibro.cols <- c("#f593f1","#f272ed","#ef52e9","#ec31e5","#e415dd","#c412bd","#a30f9d","#820c7e","#62095e")
smooth.cols <- c("#b47fe5","#9900cc")

endo.cols <- c("#93CEFF","#4A99FF","#01579B")

t.cols <- c("gray40","gray60")

macro.cols <- c("#ff6600","#ff9d5c")

mast.cols <- "gold3"

b.cols <- c("#B22222","#CD5C5C")


cols <- c(epithelial.cols,fibro.cols,smooth.cols,endo.cols,t.cols,macro.cols,mast.cols,b.cols)

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


atac.df$cell.type <- factor(atac.df$cell.type,levels = c("6-Unciliated epithelia 1",
                                                       "8-Unciliated epithelia 1",
                                                       "10-Unciliated epithelia 2",
                                                       "14-Ciliated",
                                                       "15-Unciliated epithelia 1",
                                                       "19-Ciliated",
                                                       "21-Unciliated epithelia 2",
                                                       "22-Unciliated epithelia 1",
                                                       "1-Stromal fibroblasts",
                                                       "2-Stromal fibroblasts",
                                                       "9-Stromal fibroblasts",
                                                       "12-Stromal fibroblasts",
                                                       "13-Stromal fibroblasts",
                                                       "16-Stromal fibroblasts",
                                                       "17-Stromal fibroblasts",
                                                       "20-Stromal fibroblasts","23-Stromal fibroblasts",
                                                       "4-Smooth muscle cells" ,"5-Smooth muscle cells",
                                                       "3-Endothelia","11-Endothelia","26-Endothelia",
                                                       "0-Lymphocytes","18-Lymphocytes",
                                                       "24-Macrophages","7-Macrophages",
                                                       "25-Mast cell",
                                                       "27-B cell","28-B cell"
))



epithelial.cols <- colorRampPalette(c("#A0E989", "#337B24"))
epithelial.cols <- epithelial.cols(8)

fibro.cols <- c("#f593f1","#f272ed","#ef52e9","#ec31e5","#e415dd","#c412bd","#a30f9d","#820c7e","#62095e")
smooth.cols <- c("#b47fe5","#9900cc")

endo.cols <- c("#93CEFF","#4A99FF","#01579B")

t.cols <- c("gray40","gray60")

macro.cols <- c("#ff6600","#ff9d5c")

mast.cols <- "gold3"

b.cols <- c("#B22222","#CD5C5C")


cols <- c(epithelial.cols,fibro.cols,smooth.cols,endo.cols,t.cols,macro.cols,mast.cols,b.cols)

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

df <- meta %>% group_by(RNA_snn_res.0.7) %>% count(Sample)
colnames(df) <- c("Cluster","Sample","Cells")

# Reorder cluster factor levels to group by cell type 
levels(factor(rna$cell.type))
df %>% 
  dplyr::mutate(cell.type = factor(Cluster,levels = c("6",
                                                      "8",
                                                      "10",
                                                      "14",
                                                      "15",
                                                      "19",
                                                      "21",
                                                      "22",
                                                      "1",
                                                      "2",
                                                      "9",
                                                      "12",
                                                      "13",
                                                      "16",
                                                      "17",
                                                      "20","23",
                                                      "4" ,"5",
                                                      "3","11","26",
                                                      "0","18",
                                                      "7","24",
                                                      "25",
                                                      "27","28"))) %>% 
  ggplot(aes(fill=Sample, y=Cells, x= fct_rev(cell.type))) + 
  scale_fill_manual(values = sampleColors[1:5])+
  geom_bar(position="fill", stat="identity")+
  coord_flip()+theme_classic()+xlab("Clusters")+ylab("# of cells")+
  ggsave("Cell_Type_Prop_RNA.pdf",width =4,height = 11)


# Patient proportion per subcluster in ATAC:
meta <- as.data.frame(atac@cellColData)
meta$predictedGroup_ArchR <- gsub("-.*", "", meta$predictedGroup_ArchR)

df <- meta %>% group_by(predictedGroup_ArchR) %>% count(Sample)
colnames(df) <- c("Cluster","Sample","Cells")

# Reorder cluster factor levels to group by cell type 
levels(factor(atac$predictedGroup_ArchR))
df %>% 
  dplyr::mutate(cell.type = factor(Cluster,levels = c("6",
                                                      "8",
                                                      "10",
                                                      "14",
                                                      "15",
                                                      "19",
                                                      "21",
                                                      "22",
                                                      "1",
                                                      "2",
                                                      "9",
                                                      "12",
                                                      "13",
                                                      "16",
                                                      "17",
                                                      "20","23",
                                                      "4" ,"5",
                                                      "3","11","26",
                                                      "0","18",
                                                      "7","24",
                                                      "25",
                                                      "27","28"))) %>% 
  ggplot(aes(fill=Sample, y=Cells, x= fct_rev(cell.type))) + 
  scale_fill_manual(values = sampleColors[1:5])+
  geom_bar(position="fill", stat="identity")+
  coord_flip()+theme_classic()+xlab("Clusters")+ylab("# of cells")+
  ggsave("Cell_Type_Prop_ATAC.pdf",width = 4,height =11)
###########################################################




# Immune cell proportion barcharts
###########################################################
# ATAC 
meta.atac <- as.data.frame(atac@cellColData)

meta.atac$predictedGroup_ArchR <- gsub(".*-","",meta.atac$predictedGroup_ArchR)


meta.atac$immune <- ifelse(meta.atac$predictedGroup_ArchR == "Lymphocytes" |
                             meta.atac$predictedGroup_ArchR == "B cell"|
                             meta.atac$predictedGroup_ArchR == "Mast cell" |
                             meta.atac$predictedGroup_ArchR == "Macrophages","Immune","Non-immune")


endo.1.prop <- nrow(meta.atac[meta.atac$immune == "Immune" & meta.atac$Sample == '3533EL',])/nrow(meta.atac[meta.atac$Sample == '3533EL',])
endo.2.prop <- nrow(meta.atac[meta.atac$immune == "Immune" & meta.atac$Sample == '3571DL',])/nrow(meta.atac[meta.atac$Sample == '3571DL',])
endo.3.prop <- nrow(meta.atac[meta.atac$immune == "Immune" & meta.atac$Sample == '36186L',])/nrow(meta.atac[meta.atac$Sample == '36186L',])
endo.4.prop <- nrow(meta.atac[meta.atac$immune == "Immune" & meta.atac$Sample == '36639L',])/nrow(meta.atac[meta.atac$Sample == '36639L',])
endo.5.prop <- nrow(meta.atac[meta.atac$immune == "Immune" & meta.atac$Sample == '366C5L',])/nrow(meta.atac[meta.atac$Sample == '366C5L',])

df <- data.frame(Sample = c("3533EL","3571DL","36186L","36639L","366C5L"),
                 Proportion = c(endo.1.prop,endo.2.prop,endo.3.prop,endo.4.prop,endo.5.prop))

ggplot(df,aes(x= Sample,y = Proportion,fill = Sample))+
  geom_bar(stat = "identity")+
  scale_y_continuous(labels = scales::percent,limits = c(0,0.4))+
  scale_fill_manual(values = sampleColors[1:5])+
  ylab("% immune cells")+
  theme_classic()+ggsave("Immune_Prop_ATAC.pdf",width = 4,height = 6)

# RNA
meta.rna <- as.data.frame(rna@meta.data)

meta.rna$cell.type <- gsub(".*-","",meta.rna$cell.type)


meta.rna$immune <- ifelse(meta.rna$cell.type == "Lymphocytes" |
                            meta.rna$cell.type == "B cell"|
                            meta.rna$cell.type == "Mast cell" |
                            meta.rna$cell.type == "Macrophages","Immune","Non-immune")


endo.1.prop <- nrow(meta.rna[meta.rna$immune == "Immune" & meta.rna$Sample == '3533EL',])/nrow(meta.rna[meta.rna$Sample == '3533EL',])
endo.2.prop <- nrow(meta.rna[meta.rna$immune == "Immune" & meta.rna$Sample == '3571DL',])/nrow(meta.rna[meta.rna$Sample == '3571DL',])
endo.3.prop <- nrow(meta.rna[meta.rna$immune == "Immune" & meta.rna$Sample == '36186L',])/nrow(meta.rna[meta.rna$Sample == '36186L',])
endo.4.prop <- nrow(meta.rna[meta.rna$immune == "Immune" & meta.rna$Sample == '36639L',])/nrow(meta.rna[meta.rna$Sample == '36639L',])
endo.5.prop <- nrow(meta.rna[meta.rna$immune == "Immune" & meta.rna$Sample == '366C5L',])/nrow(meta.rna[meta.rna$Sample == '366C5L',])

df <- data.frame(Sample = c("3533EL","3571DL","36186L","36639L","366C5L"),
                 Proportion = c(endo.1.prop,endo.2.prop,endo.3.prop,endo.4.prop,endo.5.prop))

ggplot(df,aes(x= Sample,y = Proportion,fill = Sample))+
  geom_bar(stat = "identity")+
  scale_y_continuous(labels = scales::percent,limits = c(0,0.4))+
  scale_fill_manual(values = sampleColors[1:5])+
  ylab("% immune cells")+
  theme_classic()+ggsave("Immune_Prop_RNA.pdf",width = 4,height = 6)
###########################################################


# Write cluster total # of cells to output files
###########################################################
total.atac <- as.data.frame(table(atac$predictedGroup_ArchR))
colnames(total.atac) <- c("Cluster","ATAC cells")

total.rna <- as.data.frame(table(rna$cell.type))
colnames(total.rna) <- c("Cluster","RNA cells")


total <- merge(total.rna,total.atac,by = "Cluster")
WriteXLS::WriteXLS(total,"./Total_Cell_Numbers_per_Cluster.xlsx")
###########################################################


# Plot peak2gene heatmap
###########################################################
atac.peaks <- readRDS("./final_archr_proj_archrGS.rds")
p2g.EEC <- getPeak2GeneLinks(ArchRProj = atac.peaks,returnLoops = F,corCutOff = 0,FDRCutOff = 1)

hist(p2g.EEC$EmpFDR)
summary(p2g.EEC$EmpFDR)
hist(p2g.EEC$FDR)
summary(p2g.EEC$FDR)
hist(p2g.EEC$Correlation)
summary(p2g.EEC$Correlation)


#####################################################################################
# Check to see if Peak to Gene links overlap with statistically significant P2G links 
# in full cohort by EmpFDR:
##################################################################################
p2g.coords <- metadata(p2g.EEC)[[1]]
ranges <- as.data.frame(p2g.coords@ranges)

EEC <- data.frame(chrom = p2g.coords@seqnames , start = ranges$start,end = ranges$end)

sig.distal <- readRDS("./Peak2GeneLinks_EmpFDR_Ranked_Distal_ArchR.rds")
sig.distal <- dplyr::filter(sig.distal,EmpFDR <= 0.05)

full <- data.frame(chrom = sig.distal$seqnames.x,start = sig.distal$start.x,end = sig.distal$end.x)
full.gr <- makeGRangesFromDataFrame(full)

# Use bedtoolsr to intersect Peak2Gene links with full cohort Peak2Gene links:
intersect <- bedtoolsr::bt.intersect(EEC,full.gr,wa = T)
colnames(intersect) <- c("seqnames","start","end")

empirical.sig.EEC.peaks.from.full.cohort <- makeGRangesFromDataFrame(intersect)
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
# 3) EEC statistically significant peaks that overlap with EmpFDR significant peaks from full cohort
# 4) EEC statistically signiicant peaks (EmpFDR is calculated within EEC cohort only)
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
# Find EEC statistically signiicant peaks (EmpFDR is calculated within EEC cohort only)
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

genes.of.interest <- c("ELF3","LAMC2","TACSTD2","CLDN3","CLDN4")

for (i in genes.of.interest){
  
  # Save P2G peaknames and Kmeans cluster for the gene being iterated
  idx <- grep(i,rownames(mat))
  kmeans.idx <- kmeans[idx]
  
  peaks.genes <- p2g$Peak2GeneLinks
  peaks.genes <- peaks.genes[idx,]
  
  peaks.genes <- as.data.frame(peaks.genes)
  peaks.genes$kmeans.group <- kmeans.idx
  
  saveRDS(peaks.genes,paste0("P2G_Hits_",i,"_data.rds"))
  peaks.genes <- dplyr::filter(peaks.genes,EmpFDR <= 0.05)

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
  # 3) EEC statistically significant peaks that overlap with EmpFDR significant peaks from full cohort
  # 4) EEC statistically signiicant peaks (EmpFDR is calculated within EEC cohort only)
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
                           features = GRangesList(TrackA = encode.all,TrackB = encode.distal,TrackC = getPeakSet(atac.peaks),TrackD = empirical.sig.EEC.peaks.from.full.cohort,TrackE=sig.peaks), 
                           pal = cols,
                           loops = getPeak2GeneLinks(atac.peaks,corCutOff = 0.45,
                                                     FDRCutOff = 0.0001))
  
  pdf(paste0("Peak2Gene_",i,"_cCREs_plus_Distal_plus_EEC_sig_p2gs.pdf"),width = 6,height = 8)
  grid::grid.draw(plot[[1]])
  dev.off()
  
  
  # Finally plot matching scRNA-seq expression in violin plots:
  
  # Plot violin plots for matching scRNA-seq expression:
  ###########################################################
  Idents(rna) <- "RNA_snn_res.0.7"
  
  my_levels <- rev(c(6,
                     8,
                     10,
                     14,
                     15,
                     19,
                     21,
                     22,
                     1,
                     2,
                     9,
                     12,
                     13,
                     16,
                     17,
                     23,
                     4 ,5,
                     3,11,26,
                     0,18,
                     7,24,
                     25,
                     27,28,20))
  
  # Relevel object@ident
  rna$cluster.new <- factor(x = rna$RNA_snn_res.0.7, levels = my_levels)
  Idents(rna) <- "cluster.new"
  
  VlnPlot(rna,features = i,group.by = "cluster.new",pt.size = 0,
          idents = c("6",
                     "8",
                     "10","14",
                     "15",
                     "19",
                     "21",
                     "22",
                     "1",
                     "2",
                     "9",
                     "12",
                     "13",
                     "16",
                     "17",
                     "23",
                     "4" ,"5",
                     "3","11","26",
                     "0","18",
                     "24","7",
                     "25",
                     "27","28"))+coord_flip()+NoLegend()+ggsave(paste0(i,"_vln.pdf"),width = 4,height = 8)
  
  
}



# Follow up and customize genomic coordinate view for each:
#######################################################

peaks.genes <- readRDS("./P2G_Hits_ELF3_data.rds")
sig.peaks <- as.data.frame(peaks.genes$peak)
sig.peaks$start <- sig.peaks$`peaks.genes$peak`
sig.peaks$end <- sig.peaks$`peaks.genes$peak`
colnames(sig.peaks) <- c("seqnames","start","end")

sig.peaks$seqnames <- gsub(":.*","",sig.peaks$seqnames)
sig.peaks$start <- gsub(".*:","",sig.peaks$start)
sig.peaks <- sig.peaks[,-3]

sig.peaks <- tidyr::separate(data = sig.peaks, col = start, into = c("start", "end"), sep = "\\-")

sig.peaks <- makeGRangesFromDataFrame(sig.peaks)
plot <- plotBrowserTrack(atac.peaks,geneSymbol ="ELF3", groupBy = "predictedGroup_ArchR",
                         features = GRangesList(TrackA = encode.all,TrackB = encode.distal,TrackC = encode.epithelium), 
                         pal = cols,upstream = 12000,downstream = 110000,
                         loops = getPeak2GeneLinks(atac.peaks,corCutOff = 0.45,
                                                   FDRCutOff = 0.0001))

pdf("ELF3_final.pdf",width = 6,height =8)
grid::grid.draw(plot[[1]])
dev.off()


# Annotation 2:
plot <- plotBrowserTrack(atac.peaks,geneSymbol ="ELF3", groupBy = "predictedGroup_ArchR",
                         features = GRangesList(TrackA = encode.all,TrackB = encode.distal,TrackC = getPeakSet(atac.peaks),TrackD = empirical.sig.EEC.peaks.from.full.cohort,TrackE=sig.peaks), 
                         pal = cols,upstream = 12000,downstream = 110000,
                         loops = getPeak2GeneLinks(atac.peaks,corCutOff = 0.45,
                                                   FDRCutOff = 0.0001))

pdf("ELF3_final_sig_peaks_anno.pdf",width = 6,height =8)
grid::grid.draw(plot[[1]])
dev.off()


peaks.genes <- readRDS("./P2G_Hits_LAMC2_data.rds")
sig.peaks <- as.data.frame(peaks.genes$peak)
sig.peaks$start <- sig.peaks$`peaks.genes$peak`
sig.peaks$end <- sig.peaks$`peaks.genes$peak`
colnames(sig.peaks) <- c("seqnames","start","end")

sig.peaks$seqnames <- gsub(":.*","",sig.peaks$seqnames)
sig.peaks$start <- gsub(".*:","",sig.peaks$start)
sig.peaks <- sig.peaks[,-3]

sig.peaks <- tidyr::separate(data = sig.peaks, col = start, into = c("start", "end"), sep = "\\-")

sig.peaks <- makeGRangesFromDataFrame(sig.peaks)
plot <- plotBrowserTrack(atac.peaks,geneSymbol ="LAMC2", groupBy = "predictedGroup_ArchR",
                         features = GRangesList(TrackA = encode.all,TrackB = encode.distal,TrackC = encode.epithelium),
                         pal = cols,upstream = 30000,downstream = 10000,
                         loops = getPeak2GeneLinks(atac.peaks,corCutOff = 0.45,
                                                   FDRCutOff = 0.0001))

pdf("LAMC2_final.pdf",width = 6,height =8)
grid::grid.draw(plot[[1]])
dev.off()


# Annotation 2:
plot <- plotBrowserTrack(atac.peaks,geneSymbol ="LAMC2", groupBy = "predictedGroup_ArchR",
                         features = GRangesList(TrackA = encode.all,TrackB = encode.distal,TrackC = getPeakSet(atac.peaks),TrackD = empirical.sig.EEC.peaks.from.full.cohort,TrackE=sig.peaks), 
                         pal = cols,upstream = 30000,downstream = 10000,
                         loops = getPeak2GeneLinks(atac.peaks,corCutOff = 0.45,
                                                   FDRCutOff = 0.0001))

pdf("LAMC2_final_sig_peaks_anno.pdf",width = 6,height =8)
grid::grid.draw(plot[[1]])
dev.off()




peaks.genes <- readRDS("./P2G_Hits_TACSTD2_data.rds")
sig.peaks <- as.data.frame(peaks.genes$peak)
sig.peaks$start <- sig.peaks$`peaks.genes$peak`
sig.peaks$end <- sig.peaks$`peaks.genes$peak`
colnames(sig.peaks) <- c("seqnames","start","end")

sig.peaks$seqnames <- gsub(":.*","",sig.peaks$seqnames)
sig.peaks$start <- gsub(".*:","",sig.peaks$start)
sig.peaks <- sig.peaks[,-3]

sig.peaks <- tidyr::separate(data = sig.peaks, col = start, into = c("start", "end"), sep = "\\-")

sig.peaks <- makeGRangesFromDataFrame(sig.peaks)
plot <- plotBrowserTrack(atac.peaks,geneSymbol ="TACSTD2", groupBy = "predictedGroup_ArchR",
                         features = GRangesList(TrackA = encode.all,TrackB = encode.distal,TrackC = encode.epithelium), 
                         pal = cols,upstream = 10000,downstream =45000,
                         loops = getPeak2GeneLinks(atac.peaks,corCutOff = 0.40,
                                                   FDRCutOff = 0.0001))

pdf("TACSTD2_final.pdf",width = 6,height =8)
grid::grid.draw(plot[[1]])
dev.off()


# Annotation 2:
plot <- plotBrowserTrack(atac.peaks,geneSymbol ="TACSTD2", groupBy = "predictedGroup_ArchR",
                         features = GRangesList(TrackA = encode.all,TrackB = encode.distal,TrackC = getPeakSet(atac.peaks),TrackD = empirical.sig.EEC.peaks.from.full.cohort,TrackE=sig.peaks), 
                         pal = cols,upstream = 10000,downstream = 45000,
                         loops = getPeak2GeneLinks(atac.peaks,corCutOff = 0.40,
                                                   FDRCutOff = 0.0001))

pdf("TACSTD2_final_sig_peaks_anno.pdf",width = 6,height =8)
grid::grid.draw(plot[[1]])
dev.off()





peaks.genes <- readRDS("./P2G_Hits_CLDN3_data.rds")
sig.peaks <- as.data.frame(peaks.genes$peak)
sig.peaks$start <- sig.peaks$`peaks.genes$peak`
sig.peaks$end <- sig.peaks$`peaks.genes$peak`
colnames(sig.peaks) <- c("seqnames","start","end")

sig.peaks$seqnames <- gsub(":.*","",sig.peaks$seqnames)
sig.peaks$start <- gsub(".*:","",sig.peaks$start)
sig.peaks <- sig.peaks[,-3]

sig.peaks <- tidyr::separate(data = sig.peaks, col = start, into = c("start", "end"), sep = "\\-")

sig.peaks <- makeGRangesFromDataFrame(sig.peaks)
plot <- plotBrowserTrack(atac.peaks,geneSymbol ="CLDN3", groupBy = "predictedGroup_ArchR",
                         features = GRangesList(TrackA = encode.all,TrackB = encode.distal,TrackC = encode.epithelium), 
                         pal = cols,upstream =10000,downstream =70000,
                         loops = getPeak2GeneLinks(atac.peaks,corCutOff = 0.45,
                                                   FDRCutOff = 0.0001))

pdf("CLDN3_final.pdf",width = 6,height =8)
grid::grid.draw(plot[[1]])
dev.off()


# Annotation 2:
plot <- plotBrowserTrack(atac.peaks,geneSymbol ="CLDN3", groupBy = "predictedGroup_ArchR",
                         features = GRangesList(TrackA = encode.all,TrackB = encode.distal,TrackC = getPeakSet(atac.peaks),TrackD = empirical.sig.EEC.peaks.from.full.cohort,TrackE=sig.peaks), 
                         pal = cols,upstream = 10000,downstream = 70000,
                         loops = getPeak2GeneLinks(atac.peaks,corCutOff = 0.45,
                                                   FDRCutOff = 0.0001))

pdf("CLDN3_final_sig_peaks_anno.pdf",width = 6,height =8)
grid::grid.draw(plot[[1]])
dev.off()



peaks.genes <- readRDS("./P2G_Hits_CLDN4_data.rds")
sig.peaks <- as.data.frame(peaks.genes$peak)
sig.peaks$start <- sig.peaks$`peaks.genes$peak`
sig.peaks$end <- sig.peaks$`peaks.genes$peak`
colnames(sig.peaks) <- c("seqnames","start","end")

sig.peaks$seqnames <- gsub(":.*","",sig.peaks$seqnames)
sig.peaks$start <- gsub(".*:","",sig.peaks$start)
sig.peaks <- sig.peaks[,-3]

sig.peaks <- tidyr::separate(data = sig.peaks, col = start, into = c("start", "end"), sep = "\\-")

sig.peaks <- makeGRangesFromDataFrame(sig.peaks)
plot <- plotBrowserTrack(atac.peaks,geneSymbol ="CLDN4", groupBy = "predictedGroup_ArchR",
                         features = GRangesList(TrackA = encode.all,TrackB = encode.distal,TrackC = encode.epithelium), 
                         pal = cols,upstream =70000,downstream =30000,
                         loops = getPeak2GeneLinks(atac.peaks,corCutOff = 0.7,
                                                   FDRCutOff = 0.0001))

pdf("CLDN4_final.pdf",width = 6,height =8)
grid::grid.draw(plot[[1]])
dev.off()


# Annotation 2:
plot <- plotBrowserTrack(atac.peaks,geneSymbol ="CLDN4", groupBy = "predictedGroup_ArchR",
                         features = GRangesList(TrackA = encode.all,TrackB = encode.distal,TrackC = getPeakSet(atac.peaks),TrackD = empirical.sig.EEC.peaks.from.full.cohort,TrackE=sig.peaks), 
                         pal = cols,upstream = 70000,downstream = 30000,
                         loops = getPeak2GeneLinks(atac.peaks,corCutOff = 0.7,
                                                   FDRCutOff = 0.0001))

pdf("CLDN4_final_sig_peaks_anno.pdf",width = 6,height =8)
grid::grid.draw(plot[[1]])
dev.off()



# Plot CNV violin: 
rna@active.ident <- factor(rna$RNA_snn_res.0.7,levels = rev(c("6",
                          "8",
                          "10",
                          "14",
                          "15",
                          "19",
                          "21",
                          "22",
                          "1",
                          "2",
                          "9",
                          "12",
                          "13",
                          "16",
                          "17",
                          "20","23",
                          "4" ,"5",
                          "3","11","26",
                          "0","18",
                          "7","24",
                          "25",
                          "27","28")))

VlnPlot(rna,features = "Total_CNVs",pt.size = 0)+coord_flip()+NoLegend()+ggsave("CNV_vln_new.pdf",width = 4,height = 8)



