###########################################################
# Matt Regner
# Franco Lab
# Date: May 2020-May 2021
# 
# Sample: All 
# Description: This script performs the following tasks  
#         1) Find motifs in peaks
############################################################
source("TFSEE_script.R")
library(data.table)
library(ggplot2)
library(Seurat)
library(Signac)
library(JASPAR2018)
library(motifmatchr)
library(GenomicRanges)
library(tidyr)
library(JASPAR2020)
library(TFBSTools)
library(BSgenome.Hsapiens.UCSC.hg38)
library(TFBSTools)
library(JASPAR2018)
library(GenomicRanges)
library(stringi)
library(stringr)
library(scales)
library(Seurat)
library(Signac)
library(Rtsne)
library(BSgenome.Mmusculus.UCSC.mm10)
library(patchwork)
library(ggplot2)
library(tidyr)
library(dplyr)
library(Rfast)
library(readxl)
library(motifmatchr)
# Get cluster labels from peak matrix
peak.mat <- readRDS("./Peak_matrix.rds")
labels <- colnames(peak.mat)



###############################################
# Part 1: Pseudobulk TF expression 
###############################################
# Read in RNA data
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
rna.pseudobulk


library(DESeq2)
dds <- DESeqDataSetFromMatrix(countData = as.matrix(rna.pseudobulk),
                              colData = data.frame(meta=colnames(rna.pseudobulk)),design = ~ 1)

dds <- estimateSizeFactors(dds)
sizeFactors(dds)

rlog.rna <- rlog(dds,blind = T)
saveRDS(rlog.rna,"rlog_rna.rds")



###############################################
# Part 2: Pseudobulk peak activity
###############################################
peaks.mat <- readRDS("./Peak_matrix.rds")

library(DESeq2)
dds <- DESeqDataSetFromMatrix(countData = peaks.mat,
                              colData = data.frame(meta=colnames(peaks.mat)),design = ~ 1)

dds <- estimateSizeFactors(dds)
sizeFactors(dds)

rlog.atac <- rlog(dds,blind = T)
saveRDS(rlog.atac,"rlog_atac.rds")


###############################################
# Part 3: Run TFSEE
###############################################
# Read in JASPAR2020 motif list:
opts <- list()
opts[["species"]] <- 9606
opts[["all_versions"]] <- TRUE
PFMatrixList <- getMatrixSet(JASPAR2020, opts)
PFMatrixList


res <- TFSEE(peak.mat=assay(rlog.atac),
             expr.mat=assay(rlog.rna),
             activity.mat=assay(rlog.rna),
             motif.list = PFMatrixList,
             motif.p.val = 1e-05,
             genome="hg38",
             return.names=T
             )


##################################################################
# Plotting of TFSEE heatmap and log2FC scatter plot 
##################################################################
library(ComplexHeatmap)
# Make heatmap annotation
ha = HeatmapAnnotation(
  Site = factor(c(rep("1",5),rep("2",6))),
  Type = factor(c(rep("1",5),rep("2",5),rep("3",1))),
  Histology = factor(c(rep("1",5),rep("2",4),rep("3",1),rep("4",1))),
  Stage = factor(c(rep("1",1),rep("2",2),rep("3",1),rep("4",1),rep("5",4),rep("6",2))),
  Patient = factor(as.character(seq(1:11)))
)


# Z-score the columns

pdf(paste0("TFSEE_Column_Scaled.pdf"), width = 5, height = 6)

Heatmap(scale(res),show_row_names =T,clustering_method_columns = "ward.D")

dev.off()

# Z-score the rows

pdf(paste0("TFSEE_Row_Scaled.pdf"), width = 5, height = 6)

Heatmap(t(scale(t(res))),show_row_names =T,clustering_method_columns = "ward.D")

dev.off()



# Plot scatter plot
tfsee <- readRDS("./TFSEE_Outputs/TFSEE_final.rds")
expr.mat <- readRDS("./TFSEE_Outputs/expression_matrix-postscaling.rds")
# # Fold in druggability information for TFs
# drug.scores <- read.csv("./TFs_CanSAR.csv")
# drug.scores <- as.data.frame(cbind(drug.scores$gene_name,drug.scores$ligand_druggability_score))
# colnames(drug.scores) <- c("TF","Score")
# dim(drug.scores)
# 
# drug.scores <- drug.scores[duplicated(drug.scores) == FALSE,]
# dim(drug.scores)
# drug.scores <- dplyr::arrange(drug.scores,TF)
# drug.scores <- drug.scores[-18,]# Remove duplicate DDIT3
# drug.scores <- drug.scores[drug.scores$TF %in% colnames(tfsee),]
# 
# drug.scores.sub <- dplyr::filter(drug.scores,Score > 0)
# 
# length(which(drug.scores$TF== colnames(tfsee)))

# # Scatter plot of log2FC TFSEE score versus log2FC TF expression (Sarcoma v. Epithelial of carsarc)

# Subset data to the most robustly expressed TFs across pseudobulk clusters


length(which(rownames(tfsee)==rownames(expr.mat)))
tfsee <- tfsee[ order(match(rownames(tfsee),rownames(expr.mat))), ]
expr.mat <- expr.mat[ order(match(rownames(expr.mat),rownames(tfsee))), ]
length(which(rownames(tfsee)==rownames(expr.mat)))
length(which(colnames(tfsee)==colnames(expr.mat)))
# tfsee <- tfsee[summary(rowMeans(tfsee)) >2,]
# expr.mat <- expr.mat[rownames(expr.mat) %in% rownames(tfsee),]

# Build TFSEE fold change vector
tfsee.sarcoma <- as.data.frame(tfsee[,3])
rownames(tfsee.sarcoma) <- make.unique(rownames(tfsee))
colnames(tfsee.sarcoma) <- "Sarcoma"
tfsee.sarcoma$Sarcoma <- tfsee.sarcoma$Sarcoma+0.1

tfsee.carcinoma <- as.data.frame(tfsee[,4])
rownames(tfsee.carcinoma) <- make.unique(rownames(tfsee))
colnames(tfsee.carcinoma) <- "carcinoma"
tfsee.carcinoma$carcinoma <- tfsee.carcinoma$carcinoma+0.1


FC.tfsee <- as.data.frame(log2(tfsee.sarcoma$Sarcoma/tfsee.carcinoma$carcinoma))
rownames(FC.tfsee) <- make.unique(rownames(tfsee))
colnames(FC.tfsee) <- "TFSEE"
hist(FC.tfsee)

# Build Expression fold change vector
expr.mat.sarcoma <- as.data.frame(expr.mat[,3])
rownames(expr.mat.sarcoma) <- make.unique(rownames(expr.mat))
colnames(expr.mat.sarcoma) <- "Sarcoma"
expr.mat.sarcoma$Sarcoma <- expr.mat.sarcoma$Sarcoma+0.1

expr.mat.carcinoma <- as.data.frame(expr.mat[,4])
rownames(expr.mat.carcinoma) <- make.unique(rownames(expr.mat))
colnames(expr.mat.carcinoma) <- "carcinoma"
expr.mat.carcinoma$carcinoma <- expr.mat.carcinoma$carcinoma+0.1


FC.expr.mat <- as.data.frame(log2(expr.mat.sarcoma$Sarcoma/expr.mat.carcinoma$carcinoma))
rownames(FC.expr.mat) <- make.unique(rownames(expr.mat))
colnames(FC.expr.mat) <- "Expr"

hist(FC.expr.mat)

plot(FC.tfsee$TFSEE,FC.expr.mat$Expr)

total <- as.data.frame(cbind(FC.expr.mat,FC.tfsee))
rownames(total) <- make.unique(rownames(tfsee))
total$sum <- rowSums(total)

total.sarcoma <- dplyr::filter(total,sum > 0)
total.carcinoma <- dplyr::filter(total,sum < 0)

drug.serous <- intersect(rownames(total.serous),drug.scores.sub$TF)
drug.endo <- intersect(rownames(total.endo),drug.scores.sub$TF)

total.serous <- total.serous[rownames(total.serous) %in% drug.serous,]
total.endo <- total.endo[rownames(total.endo) %in% drug.endo,]



# Make Scatter plot for FC TFSEE versus FC Expression

# Add druggability to DF:

drug.scores$Druggable <- ifelse(drug.scores$Score > 0, "Yes", "No")

drug.scores.new <- dplyr::filter(drug.scores,TF %in% rownames(total))

length(which(drug.scores.new$TF == rownames(total)))

total$drug.score <- drug.scores.new$Druggable
ggplot(total,aes(x=FC.tfsee,y = FC.expr))+
  geom_point(aes(color = drug.score),size = 2)+
  ylab("log2FC(TF Expression)")+theme_bw()+
  xlab("log2FC(Raw TFSEE Score)")+
  scale_color_manual(values = c("gray60","black"))+
  geom_vline(aes(xintercept = 0.5),linetype = "dashed") + 
  geom_vline(aes(xintercept = -0.5),linetype = "dashed") + 
  geom_hline(aes(yintercept = 0.5),linetype = "dashed") + 
  geom_hline(aes(yintercept = -0.5),linetype = "dashed") + ggsave("TFSEE_FC_Scar_Epi_Ov.pdf",width =5,height = 4)















