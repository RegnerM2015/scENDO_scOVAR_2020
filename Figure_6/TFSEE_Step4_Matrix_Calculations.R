###########################################################
# Matt Regner
# Franco Lab
# Date: May-December 2020
# 
# Sample: All 
# Description: This script performs the following tasks  
#     1) Read in TF motif prediction matrix 
#     2) Pseudobulk TF expression (DESeq rlog transformation)
#     3) Pseudobulk enhancer activity (DESeq rlog transformation)
#     4) Scaling of each matrix from 0-1
#     5) TFSEE matrix operations 
#################################################################
library(ArchR)
library(ComplexHeatmap)
library(DESeq2)
library(Seurat)
library(ggplot2)
library(dplyr)

###############################################
# Part 1: Read in TF motif prediction matrix 
###############################################
# Read in TF motif prediction matrix
enhancer.motifs<- read.delim("./motif_enhancers_scaled.txt",sep = "\t")
rownames(enhancer.motifs) <- enhancer.motifs[,1]
enhancer.motifs <- enhancer.motifs[,-1]
# Rename syntax of peaks
colnames(enhancer.motifs) <- sub("\\.",":",colnames(enhancer.motifs))
colnames(enhancer.motifs) <- sub("\\.","-",colnames(enhancer.motifs))


# Set labels of cell types of interest

# Only use clusters with 100% patient specificity 
labels <- c("31-Unciliated epithelia 1",
            "21-Unciliated epithelia 1",
            "19-Epithelial cell",
            "34-Epithelial cell",
            "3-Epithelial cell",
            "10-Epithelial cell",
            "16-Fibroblast",
            "17-Epithelial cell",
            "0-Fibroblast",
            "27-Fibroblast",
            "9-Epithelial cell")# Screen for "pure" patient clusters


###############################################
# Part 2: Pseudobulk TF expression 
###############################################

# Pseudobulk RNA TF expression
rna <- readRDS("./endo_ovar_All_scRNA_processed.rds")
rna.counts <- rna@assays$RNA@counts

rna.pseudobulk <- data.frame(rownames(rna))
for (i in labels){
  cells <- rownames(dplyr::filter(rna@meta.data,cell.type == i))
  
  rna.counts.sub <- rna.counts[,colnames(rna.counts) %in% cells]
  
  rna.counts.bulk <- rowSums(as.matrix(rna.counts.sub))
  
  rna.pseudobulk$i <- rna.counts.bulk
  colnames(rna.pseudobulk)[dim(rna.pseudobulk)[2]] <- i
  
}
rownames(rna.pseudobulk) <- rna.pseudobulk[,1]
rna.pseudobulk <- rna.pseudobulk[,-1]
dim(rna.pseudobulk)

# Remove any features/rows that have zero counts in ANY sample
dim(rna.pseudobulk)
rna.pseudobulk[rna.pseudobulk == 0] <- NA
rna.pseudobulk <- rna.pseudobulk[complete.cases(rna.pseudobulk),]
dim(rna.pseudobulk)


# rlog normalize data
library(DESeq2)
dds.rna <- DESeqDataSetFromMatrix(countData = rna.pseudobulk,
                              colData = data.frame(meta=colnames(rna.pseudobulk)),design = ~ 1)

dds.rna <- estimateSizeFactors(dds.rna)
sizeFactors(dds.rna)

rlog.rna <- rlog(dds.rna,blind = T)
saveRDS(dds.rna,"dds_rna.rds")
saveRDS(rlog.rna,"rlog_rna.rds")


###############################################
# Part 3: Pseudobulk enhancer activity
###############################################

# Pseudobulk ATAC enhancer matrix

atac <- readRDS("./final_archr_proj_archrGS.rds")

source("./Archr_Peak_Null_Permute.R")
source("./Archr_Peak_RawPval.R")
source("./getMatrixFromProject.R")
peaks.mat <- getMatrixFromProject.mod(atac,useMatrix = "PeakMatrix")
cell.names <- colnames(assay(peaks.mat))
peak.names <- paste0(seqnames(rowRanges(peaks.mat)),":", start(ranges(rowRanges(peaks.mat))),"-", end(ranges(rowRanges(peaks.mat))))


peaks.mat <- assay(peaks.mat)

dim(peaks.mat)
length(cell.names)
length(peak.names)

colnames(peaks.mat) <- cell.names
rownames(peaks.mat) <- peak.names
head(peaks.mat[,1:4])


# Construct pseudobulk peak/enhancer matrix
peaks.pseudobulk <- data.frame(rownames(peaks.mat))
for (i in labels){
  cells <- rownames(dplyr::filter(as.data.frame(atac@cellColData),predictedGroup_ArchR == i))
  
  peaks.sub <- peaks.mat[,colnames(peaks.mat) %in% cells]
  peaks.sub <- as(peaks.sub, "dgCMatrix")
  peaks.bulk <- Matrix::rowSums(peaks.sub)
  
  peaks.pseudobulk$i <- peaks.bulk
  colnames(peaks.pseudobulk)[dim(peaks.pseudobulk)[2]] <- i
  
  print("Iteration complete")
}
rownames(peaks.pseudobulk) <- peaks.pseudobulk$rownames.peaks.
peaks.pseudobulk <- peaks.pseudobulk[,-1]

# Remove features/peaks that have zero counts in ANY sample
dim(peaks.pseudobulk)
peaks.pseudobulk[peaks.pseudobulk == 0] <- NA
peaks.pseudobulk <- peaks.pseudobulk[complete.cases(peaks.pseudobulk),]
dim(peaks.pseudobulk)

#DESeq ATAC

library(DESeq2)
dds.atac <- DESeqDataSetFromMatrix(countData = peaks.pseudobulk,
                              colData = data.frame(meta=colnames(peaks.pseudobulk)),design = ~ 1)

dds.atac <- estimateSizeFactors(dds.atac)
sizeFactors(dds.atac)

rlog.atac <- rlog(dds.atac,blind = T)

saveRDS(rlog.atac,"rlog_atac.rds")
saveRDS(dds.atac,"dds_atac.rds")
rlog.rna <- readRDS("rlog_rna.rds")
rlog.atac <- readRDS("rlog_atac.rds")

###############################################
# Part 4: Subsetting and scaling from 0 to 1
###############################################

# Subset and reorder motif matrix, peak matrix and expression matrix
enhancer.mat <- assay(rlog.atac)
enhancer.mat <- as.data.frame(enhancer.mat)
enhancer.mat <- enhancer.mat[rownames(enhancer.mat) %in% colnames(enhancer.motifs),]

# Check order
motif.mat <- enhancer.motifs
motif.mat <- motif.mat[,colnames(motif.mat) %in% rownames(enhancer.mat)]
motif.mat <- motif.mat[ ,order(match(colnames(motif.mat), rownames(enhancer.mat))) ]
length(which(colnames(motif.mat) == rownames(enhancer.mat)))

expr.mat <- assay(rlog.rna)
expr.mat <- as.data.frame(expr.mat)

# Check order
tfs <- intersect(rownames(motif.mat),rownames(expr.mat))
expr.mat <- expr.mat[rownames(expr.mat) %in% tfs,]
motif.mat <- motif.mat[rownames(motif.mat) %in% tfs,]
expr.mat <- expr.mat[order(match(rownames(expr.mat), rownames(motif.mat))), ]

length(which(rownames(expr.mat)==rownames(motif.mat)))

# Rescale peak matrix and rna expr matrix
library(scales)
for (i in 1:ncol(expr.mat)){
  expr.mat[,i] <- rescale(expr.mat[,i],to = c(0,1))
}
library(scales)# Scale cell types
for (i in 1:ncol(enhancer.mat)){
  enhancer.mat[,i] <- rescale(enhancer.mat[,i],to = c(0,1))
}


dim(enhancer.mat)
dim(motif.mat)
dim(expr.mat)



length(which(rownames(motif.mat) == rownames(expr.mat)))

length(which(rownames(enhancer.mat) == colnames(motif.mat)))



length(which(rownames(motif.mat) == rownames(expr.mat)))

length(which(rownames(enhancer.mat) == colnames(motif.mat)))

dim(enhancer.mat)
dim(motif.mat)
dim(expr.mat)


###############################################
# Part 5: TFSEE matrix operations
###############################################

# Peform matrix multiplication of enhancer activity by motif prediction
int <- t(enhancer.mat) %*% t(motif.mat)

dim(int)
dim(expr.mat)

length(which(colnames(int)==rownames(expr.mat)))



# Multiply intermediate matrix by TF expression matrix 
tfsee <- int*t(expr.mat)
dim(tfsee)
head(tfsee)

saveRDS(tfsee,"tfsee_matrix.rds")
saveRDS(int,"int_matrix.rds")
saveRDS(expr.mat,"expr_matrix.rds")
saveRDS(motif.mat,"motif_matrix.rds")

tfsee <- readRDS("tfsee_matrix.rds")
# Write table of TFs for canSAR search:
write.csv(data.frame(colnames(tfsee)),"TFs_for_canSAR.csv")

##################################################################
# END OF TFSEE COMPUTATION
##################################################################



##################################################################
# Plotting of TFSEE heatmap 
##################################################################
# Enrich for highest signal TFs:
tfsee.t <- t(tfsee)
tfsee.t <- as.data.frame(tfsee.t)
tfsee.t$rowMeans <- rowMeans(as.matrix(tfsee.t))
tfsee.t <- dplyr::filter(tfsee.t,rowMeans > summary(tfsee.t$rowMeans)[2])
tfsee.t <- tfsee.t[,-12]# Remove rowMeans column
tfsee <- t(tfsee.t)


# Output TF names to csv file for canSAR search:
write.csv(data.frame(TF=colnames(tfsee)),quote = F,sep = ",","TFs_for_canSAR.csv")

# Make heatmap annotation
ha = HeatmapAnnotation(
  Site = factor(c(rep("1",5),rep("2",6))),
  Type = factor(c(rep("1",5),rep("2",5),rep("3",1))),
  Histology = factor(c(rep("1",5),rep("2",4),rep("3",1),rep("4",1))),
  Stage = factor(c(rep("1",1),rep("2",2),rep("3",1),rep("4",1),rep("5",4),rep("6",2))),
  Patient = factor(as.character(seq(1:11))),
  size.rna = DESeq2::sizeFactors(dds.rna),
  size.atac = DESeq2::sizeFactors(dds.atac))


library(ComplexHeatmap)
# Z-score the TFSEEs 
pdf(paste0("TFSEE_TF_scaled.pdf"), width = 5, height = 9)
Heatmap(t(scale(tfsee)),show_row_names =T,clustering_method_columns = "complete",top_annotation = ha,
        clustering_method_rows="complete",
        clustering_distance_columns = "euclidean",
        clustering_distance_rows = "euclidean")
dev.off()
# Z-score the cell types
pdf(paste0("TFSEE_celltype_scaled.pdf"), width = 5, height =9)
Heatmap(scale(t(tfsee)),top_annotation = ha, 
        show_row_names =T,clustering_method_rows = "complete",clustering_method_columns = "complete",clustering_distance_rows = "euclidean",
        clustering_distance_columns = "euclidean")
dev.off()
# No z-score normalization
pdf(paste0("TFSEE_celltype_uncscaled.pdf"), width = 5, height = 9)
Heatmap(t(tfsee),top_annotation = ha,
        show_row_names =T,clustering_method_rows = "complete",clustering_method_columns = "complete",
        clustering_distance_columns="euclidean",
        clustering_distance_row="euclidean")
dev.off()


writeLines(capture.output(sessionInfo()), "sessionInfo_TFSEE_Step4.txt")
