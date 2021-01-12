library(ComplexHeatmap)
library(tidyverse)
library(DESeq2)
###############################################
# Part 4: Scaling of each matrix from 0-1
###############################################
rlog.atac <- readRDS("rlog_atac.rds")
rlog.rna <- readRDS("rlog_rna.rds")
enhancer.motifs <- readRDS("enhancer_motifs.rds")

# Scale Enhancer matrix :
enhancer.mat <- assay(rlog.atac)
enhancer.mat <- as.data.frame(enhancer.mat)
enhancer.mat$mean <- rowMeans(enhancer.mat)
enhancer.mat <- dplyr::filter(enhancer.mat,mean > summary(enhancer.mat$mean)[2])
enhancer.mat <- enhancer.mat[,-12]
enhancer.mat <- as.matrix(enhancer.mat)
library(scales)
for (i in 1:ncol(enhancer.mat)){
  enhancer.mat[,i] <- rescale(enhancer.mat[,i],to = c(0,1))
}

#saveRDS(enhancer.motifs,"enhancer_motifs.rds")
# Scale TF motif matrix (invert the pvalue):
motif.mat <- enhancer.motifs
motif.mat <- motif.mat[,colnames(motif.mat) %in% rownames(enhancer.mat)]
length(which(colnames(motif.mat) == rownames(enhancer.mat)))
library(scales)
for (i in 1:ncol(motif.mat)){
  motif.mat[,i] <- rescale(-log10(motif.mat[,i]+0.5),to = c(0,1))
}
# library(scales)
# for (i in 1:ncol(motif.mat)){
#   motif.mat[,i] <- rescale(motif.mat[,i],to = c(0,1))
# }


# Scale TF expression matrix:
expr.mat <- assay(rlog.rna)
expr.mat <- as.data.frame(expr.mat)
expr.mat <- expr.mat[rownames(expr.mat) %in% TFs.to.use,]
expr.mat$mean <- rowMeans(expr.mat)
expr.mat <- dplyr::filter(expr.mat,mean > summary(expr.mat$mean)[2])
expr.mat <- expr.mat[,-12]
expr.mat <- as.matrix(expr.mat)
library(scales)
for (i in 1:ncol(expr.mat)){
  expr.mat[,i] <- rescale(expr.mat[,i],to = c(0,1))
}

TFs <- intersect(rownames(motif.mat),rownames(expr.mat))
expr.mat <- expr.mat[TFs,]
motif.mat <- motif.mat[TFs,]
length(which(rownames(motif.mat) == rownames(expr.mat)))

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
saveRDS(expr.mat,"expr_matrix.rds")
saveRDS(motif.mat,"motif_matrix.rds")
##################################################################
# END OF TFSEE COMPUTATION
##################################################################



##################################################################
# Plotting of TFSEE heatmap and log2FC scatter plot 
##################################################################

# Make heatmap annotation
ha = HeatmapAnnotation(
  Site = factor(c(rep("1",5),rep("2",6))),
  Type = factor(c(rep("1",5),rep("2",5),rep("3",1))),
  Histology = factor(c(rep("1",5),rep("2",4),rep("3",1),rep("4",1))),
  Stage = factor(c(rep("1",1),rep("2",2),rep("3",1),rep("4",1),rep("5",4),rep("6",2))),
  Patient = factor(as.character(seq(1:11)))
)


# Z-score the TFSEEs to reveal patterns/groups of cell types 

pdf(paste0("TFSEE_TF_scaled_Zscore_v2.pdf"), width = 5, height = 6)

Heatmap(t(scale(tfsee)),show_row_names =T)

dev.off()



# Z-score the cell types to reveal patterns/groups of TF/enhancers

pdf(paste0("TFSEE_celltype_scaled_Zscore_v2.pdf"), width = 5, height = 6)

Heatmap(scale(t(tfsee)),show_row_names =T)

dev.off()
