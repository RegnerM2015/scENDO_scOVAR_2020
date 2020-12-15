###########################################################
# Matt Regner
# Franco Lab
# Date: May-December 2020
# 

# Description: This script performs the following tasks  
#         1) scRNA-seq processing before doublet removal
#         2) Doublet detection and removal
#         3) scRNA-seq processing after doublet removal 
#         4) SingleR cell typing 
###########################################################

# Change to your working directory 
dir <- "."
setwd(dir) 

#Change to your R library path
########################################################################
.libPaths('/home/regnerm/anaconda3/envs/r-environment/lib/R/library')
source("./rowr.R")
########################################################################

library(scater)
library(dplyr)
library(Seurat)
library(patchwork)
library(SingleCellExperiment)
library(scater)
library(ComplexHeatmap)
library(ConsensusClusterPlus)
library(msigdbr)
library(fgsea)
library(dplyr)
library(tibble)
library(DoubletFinder)
library(Signac)
library(ggplot2)
library(stringr)
library(SingleR)

# Define filepaths/variables/signature sets:
###############################################

# PanglaoDB 
tsv=gzfile("./PanglaoDB_markers_27_Mar_2020.tsv.gz")  
panglaodb <- read.csv(tsv,header=T,sep = "\t") 
panglaodb <- dplyr::filter(panglaodb,species == "Hs" | species == "Mm Hs")# Human subset 
panglaodb <- dplyr::filter(panglaodb,organ == "Connective tissue" |
                             organ == "Epithelium" |
                             organ == "Immune system" |
                             organ == "Reproductive"|
                             organ == "Vasculature" |
                             organ == "Smooth muscle"
)
panglaodb <- split(as.character(panglaodb$official.gene.symbol), panglaodb$cell.type)

# ESTIMATE signatures 
ESTIMATE.signatures <- "./ESTIMATE_signatures.csv"

GRCH38.annotations <- "./Homo_sapiens.GRCh38.86.txt"

doublet.rate = 0.0760

SAMPLE.ID = "ovar_3E4D1L"

###########################################################
# Part 1: scRNA-seq processing before doublet removal
###########################################################

# Load the RNA dataset
counts <- Read10X_h5(filename= "./filtered_feature_bc_matrix.h5")
# Initialize the Seurat object with the raw (non-normalized data).
rna <- CreateSeuratObject(counts = counts, min.cells = 3)# genes not present in at least 3 cells are removed 
rna

PreQCNumCells <- length(colnames(rna))

# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
rna[["percent.mt"]] <- PercentageFeatureSet(rna, pattern = "^MT-")


# QC metrics: nCount_RNA, nFeature_RNA, and percent mitochrondial counts
# Outliers are >2 MADs
rna@meta.data$nCount_RNA_outlier_2mad <- isOutlier(log(rna@meta.data$nCount_RNA),log = F,type = "lower",nmads = 2)
rna@meta.data$nFeature_RNA_outlier_2mad <- isOutlier(log(rna@meta.data$nFeature_RNA),log = F,type = "lower",nmads = 2)
rna@meta.data$percent_mt_outlier_2mad <- isOutlier(log1p(rna@meta.data$percent.mt),log = F,type = "higher",nmads = 2)

rna <- subset(rna, subset = nCount_RNA_outlier_2mad == "FALSE" &
                nFeature_RNA_outlier_2mad == 'FALSE' & 
                percent_mt_outlier_2mad == "FALSE")

PostQCNumCells <- length(colnames(rna))

# Default Seurat processing
rna <- NormalizeData(rna, normalization.method = "LogNormalize", scale.factor = 10000)

#Feature Selection
rna <- FindVariableFeatures(rna, selection.method = "vst", nfeatures = 2000)

#Scaling & PCA
all.genes <- rownames(rna)
rna <- ScaleData(rna, features = all.genes)
rna <- RunPCA(rna)


# Score cells for immune/stromal/fibroblast/endothelial signatures
############################################################################################

immune.stromal <- read.csv(ESTIMATE.signatures,header = F)

stromal <- immune.stromal$V1[1:141]
immune <- immune.stromal$V1[142:282]
fibroblast <- panglaodb$Fibroblasts
endothelial <- panglaodb$`Endothelial cells`
epithelial <- panglaodb$`Epithelial cells`
smooth <- panglaodb$`Smooth muscle cells`
plasma <- panglaodb$`Plasma cells`

feats <- list(stromal,immune,fibroblast,endothelial,epithelial,smooth,plasma)

rna <- AddModuleScore(rna,features = feats,name = c("stromal.","immune.","fibroblast.","endothelial.",
                                                    "epithelial.","smooth.","plasma."),search = T)
#########################################################################

#######################################################################

# Add PC1 to metadata
rna@meta.data$PC1 <- rna@reductions$pca@cell.embeddings[,1]

count_cor_PC1 <- cor(rna$PC1,rna$nCount_RNA,method = "spearman")

stromal.cor <- cor(rna$PC1,rna$stromal.1,method = "spearman")
immune.cor <- cor(rna$PC1,rna$immune.2,method = "spearman")
fibroblast.cor <- cor(rna$PC1,rna$fibroblast.3,method = "spearman")
endothelial.cor <- cor(rna$PC1,rna$endothelial.4,method = "spearman")
epithelial.cor <- cor(rna$PC1,rna$epithelial.5,method = "spearman")
smooth.cor <- cor(rna$PC1,rna$smooth.6,method = "spearman")
plasma.cor <- cor(rna$PC1,rna$plasma.7,method = "spearman")


# Make JackStraw Plot
rna <- JackStraw(rna, num.replicate = 100,dims = 50)
rna <- ScoreJackStraw(rna,dims = 1:50)
JackStrawPlot(rna, dims = 1:50)+ggsave("JackStraw_predoublet.png")
# If PC1 is correalted with read depth, check to see if biological variation is corralted to PC1
if (round(abs(count_cor_PC1),2) > 0.5){
  
  if( round(abs(stromal.cor),2) >= 0.5 |
      round(abs(immune.cor),2) >= 0.5 |
      round(abs(fibroblast.cor),2) >= 0.5 |
      round(abs(endothelial.cor),2) >= 0.5 |
      round(abs(epithelial.cor),2) >= 0.5 |
      round(abs(smooth.cor),2) >= 0.5 |
      round(abs(plasma.cor),2) >= 0.5){
    rna <- FindNeighbors(rna,dims = 1:50)
    rna <- FindClusters(rna,resolution = 0.7)
    rna <- RunUMAP(rna,dims = 1:50)
    Idents(rna) <- "RNA_snn_res.0.7"
    
    # Verify with inferCNV: is PC1 correlated with CNV events/Malignancy?
    #########################################################################
    # inferCNV: does PC1 also correlated with CNV/malignancy status?
    library(infercnv)
    library(stringr)
    library(Seurat)
    counts_matrix = GetAssayData(rna, slot="counts")
    
    # Identify immune clusters
    #######################################################
    # Find immune cells by relative enrichment of ESTIMATE immune signature
    library(psych)
    test <- VlnPlot(rna,features = "immune.2")
    data <- describeBy(test$data$immune.2, test$data$ident, mat = TRUE)
    data.immune <- dplyr::filter(data,median > 0.1)
    
    test <- VlnPlot(rna,features = "plasma.7")
    data <- describeBy(test$data$plasma.7, test$data$ident, mat = TRUE)
    data.plasma <- dplyr::filter(data,median > 0.1)
    
    immune.clusters <- intersect(data.immune$group1,levels(Idents(rna)))
    plasma.clusters <- intersect(data.plasma$group1,levels(Idents(rna)))
    
    immune.clusters <- unique(append(immune.clusters,plasma.clusters))
    
    for (i in 1:length(immune.clusters)){
      j <- which(levels(Idents(rna)) == immune.clusters[i])
      levels(Idents(rna))[j] <- paste0("immune.",immune.clusters[i])
    }
    rna@meta.data$predoublet.idents <- Idents(rna)
    idents <- data.frame(rownames(rna@meta.data),rna@meta.data$predoublet.idents)
    
    colnames(idents) <- c("V1","V2")
    
    saveRDS(rna,"./rna_predoublet_preinferCNV.rds")
    # Make inferCNV input files
    rownames(idents) <- NULL
    colnames(idents) <- NULL
    write.table(idents,"./sample_annotation_file_inferCNV.txt",sep = "\t",row.names = FALSE)
    idents <- read.delim("./sample_annotation_file_inferCNV.txt",header = F)
    
    
    gtf <- read.delim(GRCH38.annotations,header = F)
    library(EnsDb.Hsapiens.v86)
    convert.symbol = function(Data){
      ensembls <- Data$V1
      ensembls <- gsub("\\.[0-9]*$", "", ensembls)
      geneIDs1 <- ensembldb::select(EnsDb.Hsapiens.v86, keys= ensembls, keytype = "GENEID", columns = "SYMBOL")
      Data <- cbind.fill(Data, geneIDs1, fill = NA)
      Data <- na.omit(Data)
      Data$feature <- Data$SYMBOL
      Data.new <- data.frame(Data$SYMBOL,Data$V2,Data$V3,Data$V4)
      Data.new$Data.V2 <- paste("chr",Data.new$Data.V2,sep = "")
      Data.new$Data.SYMBOL <- make.unique(Data.new$Data.SYMBOL)
      return(Data.new)
    }
    gtf <- convert.symbol(gtf)
    head(gtf)
    
    write.table(gtf,"./Homo_sapiens.GRCh38.86.symbol.txt",sep = "\t",row.names = FALSE,col.names = FALSE)
    
    
    num.immune.clusters = length(immune.clusters)
    # create the infercnv object
    if ( num.immune.clusters == 1) {
      infercnv_obj = CreateInfercnvObject(raw_counts_matrix=as.matrix(counts_matrix),
                                          annotations_file="./sample_annotation_file_inferCNV.txt",
                                          gene_order_file="./Homo_sapiens.GRCh38.86.symbol.txt",
                                          ref_group_names=paste0("immune.",immune.clusters[1]))
      
    } else if (num.immune.clusters == 2) {
      infercnv_obj = CreateInfercnvObject(raw_counts_matrix=as.matrix(counts_matrix),
                                          annotations_file="./sample_annotation_file_inferCNV.txt",
                                          gene_order_file="./Homo_sapiens.GRCh38.86.symbol.txt",
                                          ref_group_names=c(paste0("immune.",immune.clusters[1]),
                                                            paste0("immune.",immune.clusters[2])))
      
    } else if ( num.immune.clusters == 3 ) {
      infercnv_obj = CreateInfercnvObject(raw_counts_matrix=as.matrix(counts_matrix),
                                          annotations_file="./sample_annotation_file_inferCNV.txt",
                                          gene_order_file="./Homo_sapiens.GRCh38.86.symbol.txt",
                                          ref_group_names=c(paste0("immune.",immune.clusters[1]),
                                                            paste0("immune.",immune.clusters[2]),
                                                            paste0("immune.",immune.clusters[3])))
    } else if (num.immune.clusters == 4) {
      infercnv_obj = CreateInfercnvObject(raw_counts_matrix=as.matrix(counts_matrix),
                                          annotations_file="./sample_annotation_file_inferCNV.txt",
                                          gene_order_file="./Homo_sapiens.GRCh38.86.symbol.txt",
                                          ref_group_names=c(paste0("immune.",immune.clusters[1]),
                                                            paste0("immune.",immune.clusters[2]),
                                                            paste0("immune.",immune.clusters[3]),
                                                            paste0("immune.",immune.clusters[4])))
    } else if (num.immune.clusters == 5) {
      infercnv_obj = CreateInfercnvObject(raw_counts_matrix=as.matrix(counts_matrix),
                                          annotations_file="./sample_annotation_file_inferCNV.txt",
                                          gene_order_file="./Homo_sapiens.GRCh38.86.symbol.txt",
                                          ref_group_names=c(paste0("immune.",immune.clusters[1]),
                                                            paste0("immune.",immune.clusters[2]),
                                                            paste0("immune.",immune.clusters[3]),
                                                            paste0("immune.",immune.clusters[4]),
                                                            paste0("immune.",immune.clusters[5])))
    } else if (num.immune.clusters == 6) {
      infercnv_obj = CreateInfercnvObject(raw_counts_matrix=as.matrix(counts_matrix),
                                          annotations_file="./sample_annotation_file_inferCNV.txt",
                                          gene_order_file="./Homo_sapiens.GRCh38.86.symbol.txt",
                                          ref_group_names=c(paste0("immune.",immune.clusters[1]),
                                                            paste0("immune.",immune.clusters[2]),
                                                            paste0("immune.",immune.clusters[3]),
                                                            paste0("immune.",immune.clusters[4]),
                                                            paste0("immune.",immune.clusters[5]),
                                                            paste0("immune.",immune.clusters[6])))
    }else if (num.immune.clusters == 7) {
      infercnv_obj = CreateInfercnvObject(raw_counts_matrix=as.matrix(counts_matrix),
                                          annotations_file="./sample_annotation_file_inferCNV.txt",
                                          gene_order_file="./Homo_sapiens.GRCh38.86.symbol.txt",
                                          ref_group_names=c(paste0("immune.",immune.clusters[1]),
                                                            paste0("immune.",immune.clusters[2]),
                                                            paste0("immune.",immune.clusters[3]),
                                                            paste0("immune.",immune.clusters[4]),
                                                            paste0("immune.",immune.clusters[5]),
                                                            paste0("immune.",immune.clusters[6]),
                                                            paste0("immune.",immune.clusters[7])))
    }else if (num.immune.clusters == 8) {
      infercnv_obj = CreateInfercnvObject(raw_counts_matrix=as.matrix(counts_matrix),
                                          annotations_file="./sample_annotation_file_inferCNV.txt",
                                          gene_order_file="./Homo_sapiens.GRCh38.86.symbol.txt",
                                          ref_group_names=c(paste0("immune.",immune.clusters[1]),
                                                            paste0("immune.",immune.clusters[2]),
                                                            paste0("immune.",immune.clusters[3]),
                                                            paste0("immune.",immune.clusters[4]),
                                                            paste0("immune.",immune.clusters[5]),
                                                            paste0("immune.",immune.clusters[6]),
                                                            paste0("immune.",immune.clusters[7]),
                                                            paste0("immune.",immune.clusters[8])))
    }else if (num.immune.clusters == 9) {
      infercnv_obj = CreateInfercnvObject(raw_counts_matrix=as.matrix(counts_matrix),
                                          annotations_file="./sample_annotation_file_inferCNV.txt",
                                          gene_order_file="./Homo_sapiens.GRCh38.86.symbol.txt",
                                          ref_group_names=c(paste0("immune.",immune.clusters[1]),
                                                            paste0("immune.",immune.clusters[2]),
                                                            paste0("immune.",immune.clusters[3]),
                                                            paste0("immune.",immune.clusters[4]),
                                                            paste0("immune.",immune.clusters[5]),
                                                            paste0("immune.",immune.clusters[6]),
                                                            paste0("immune.",immune.clusters[7]),
                                                            paste0("immune.",immune.clusters[8]),
                                                            paste0("immune.",immune.clusters[9])))
    }else if (num.immune.clusters == 10) {
      infercnv_obj = CreateInfercnvObject(raw_counts_matrix=as.matrix(counts_matrix),
                                          annotations_file="./sample_annotation_file_inferCNV.txt",
                                          gene_order_file="./Homo_sapiens.GRCh38.86.symbol.txt",
                                          ref_group_names=c(paste0("immune.",immune.clusters[1]),
                                                            paste0("immune.",immune.clusters[2]),
                                                            paste0("immune.",immune.clusters[3]),
                                                            paste0("immune.",immune.clusters[4]),
                                                            paste0("immune.",immune.clusters[5]),
                                                            paste0("immune.",immune.clusters[6]),
                                                            paste0("immune.",immune.clusters[7]),
                                                            paste0("immune.",immune.clusters[8]),
                                                            paste0("immune.",immune.clusters[9]),
                                                            paste0("immune.",immune.clusters[10])))
    }
    
    # perform infercnv operations to reveal cnv signal
    infercnv_obj = infercnv::run(infercnv_obj,
                                 cutoff=0.1,  # use 1 for smart-seq, 0.1 for 10x-genomics
                                 out_dir="./output_dir_CNV_predoublet",  # dir is auto-created for storing outputs
                                 cluster_by_groups=T,   # cluster
                                 denoise=T,scale_data = T,
                                 HMM=T,HMM_type = "i6",analysis_mode = "samples",min_cells_per_gene = 10,
                                 BayesMaxPNormal = 0.4, num_threads = 8
                                 
    )
    
    
    
    regions <- read.delim("./output_dir_CNV_predoublet/HMM_CNV_predictions.HMMi6.hmm_mode-samples.Pnorm_0.4.pred_cnv_regions.dat")
    probs <- read.delim("./output_dir_CNV_predoublet/BayesNetOutput.HMMi6.hmm_mode-samples/CNV_State_Probabilities.dat")
    
    probs <- as.data.frame(t(probs[3,]))
    colnames(probs) <- "Prob.Normal"
    probs <- dplyr::filter(probs,Prob.Normal < 0.05)
    cnvs <- rownames(probs)
    cnvs <- gsub("\\.","-",cnvs)
    
    regions <- regions[regions$cnv_name %in% cnvs, ]
    
    
    cnv.groups <- sub("\\..*", "", regions$cell_group_name)
    
    
    length(which(rownames(rna@reductions$pca@cell.embeddings) == rownames(rna@meta.data)))
    rna$PC1.loading <- rna@reductions$pca@cell.embeddings[,1]
    rna$cell.barcode <- rownames(rna@meta.data)
    rna$CNV.Pos <- ifelse(as.character(rna$predoublet.idents) %in% cnv.groups,1,0)
    
    
    cnv.freq <- data.frame(table(regions$cell_group_name))
    cnv.freq$Var1 <- sub("\\..*", "", cnv.freq$Var1)
    
    rna$Total_CNVs <- ifelse(as.character(rna$predoublet.idents) %in% cnv.freq$Var1,cnv.freq$Freq,0)

    boxplot.cnv <- ggplot(rna@meta.data,aes(x= predoublet.idents,y=PC1.loading,color = as.factor(CNV.Pos)))+geom_boxplot()
    boxplot.cnv+ggsave("Predoublet_CNV_PC1_boxplot.png")
    
    data <- describeBy(boxplot.cnv$data$PC1.loading, boxplot.cnv$data$predoublet.idents, mat = TRUE)
    data$CNV <- ifelse(data$group1 %in% cnv.groups,1,0)
    
    wilcox <- wilcox.test(data = rna@meta.data,PC1.loading~CNV.Pos)
    
    if (wilcox$p.value < 0.05){
      rna <- rna
      library(stringr)
      levels(Idents(rna)) <- str_remove(levels(Idents(rna)),"immune.")
      saveRDS(rna,"./rna_predoublet_PassedPC1Checks.rds")
    }else{
      all.genes <- rownames(rna)
      rna <- ScaleData(rna, features = all.genes,vars.to.regress = "nCount_RNA")
      rna <- FindNeighbors(rna,dims = 1:50)
      rna <- FindClusters(rna,resolution = 0.7)
      rna <- RunUMAP(rna,dims = 1:50)
      Idents(rna) <- "RNA_snn_res.0.7"
      
      saveRDS(rna,"./rna_predoublet_FailedCNVTest.rds")
    }

  }else{
    all.genes <- rownames(rna)
    rna <- ScaleData(rna, features = all.genes,vars.to.regress = "nCount_RNA")
    rna <- FindNeighbors(rna,dims = 1:50)
    rna <- FindClusters(rna,resolution = 0.7)
    rna <- RunUMAP(rna,dims = 1:50)
    Idents(rna) <- "RNA_snn_res.0.7"
    saveRDS(rna,"./rna_predoublet_FailedCorTest.rds")
  }
}else{
  rna <- FindNeighbors(rna,dims = 1:50)
  rna <- FindClusters(rna,resolution = 0.7)
  rna <- RunUMAP(rna,dims = 1:50)
  Idents(rna) <- "RNA_snn_res.0.7"
  saveRDS(rna,"./rna_predoublet_SkipChecks.rds")
}

###########################################################
# Part 2: Doublet detection and removal
###########################################################

#' # Doublet detection: 1) DoubletDecon 2) DoubletFinder 3) Take intersection of calls
######################################################################################

# # Doublet Decon
# # Add if else statement to regress out nCount RNA if needed before running DoubletDecon
library(DoubletDecon)
idents.length <- length(levels(Idents(rna)))
for (i in levels(Idents(rna))){
  i <- as.numeric(i)
  levels(Idents(rna))[i] <- i -1
}

#Improved_Seurat_Pre_Process()
#Idents(rna) <- as.factor(Idents(rna))
seuratObject=rna
newFiles=Improved_Seurat_Pre_Process(seuratObject, num_genes=50, write_files=F)

#dev.off()
results=Main_Doublet_Decon(rawDataFile=newFiles$newExpressionFile,
                           groupsFile=newFiles$newGroupsFile,
                           filename=SAMPLE.ID,
                           location="./DoubletDecon",
                           fullDataFile=NULL,
                           removeCC=FALSE,
                           species="hsa",
                           rhop=0.9,
                           write=F,
                           PMF=TRUE,
                           useFull=FALSE,
                           heatmap=FALSE,
                           centroids=TRUE,
                           num_doubs=100,
                           only50=FALSE,
                           min_uniq=4,
                           nCores=4)

decon.doublets <- rownames(results$Final_doublets_groups)
decon.doublets <- gsub("\\.","-",decon.doublets)

# DoubletFinder
# Add modfieid Parameter sweep function to regress out nCOunt RNA if needed
# Add if else statement to regress out nCount RNA if needed before running DoubletFinder
#########################################################################################
library(DoubletFinder)
# Estimate Doublets
#################################################################################################
## pK Identification (no ground-truth) ---------------------------------------------------------------------------------------
sweep.res.list <- paramSweep_v3(rna, PCs = 1:50, sct = FALSE)

sweep.stats <- summarizeSweep(sweep.res.list , GT = FALSE)

bcmvn<- find.pK(sweep.stats)

pK.1 <- as.numeric(unfactor(dplyr::arrange(bcmvn,desc(BCmetric))$pK[1]))


nExp_poi <- round(doublet.rate*length(colnames(counts)))  ## Assuming 4.6% doublet formation rate - tailor for your dataset

# Run doubletFinder with pk.1
rna <- doubletFinder_v3(rna, PCs = 1:50, pN = 0.25, pK = pK.1, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)
meta <- rna@meta.data


doublet.column <- paste0("DF.classifications_0.25_",pK.1,"_",nExp_poi)


doublet.calls <- rna[[doublet.column]]
colnames(doublet.calls) <- "Call"

rna.dub <- dplyr::filter(doublet.calls, Call == "Doublet")
rna.singlet <- dplyr::filter(doublet.calls, Call == "Singlet")

DF.doublets <- rownames(rna.dub)


# Intersect doublet calls
###########################################################################

head(DF.doublets)
head(decon.doublets)

doublets <- intersect(decon.doublets,DF.doublets)
rna@meta.data$Doublet.Call <- ifelse(rownames(rna@meta.data) %in% doublets,"TRUE","FALSE")
FeatureScatter(rna,feature1 = "nCount_RNA",feature2 = "nFeature_RNA",group.by = "Doublet.Call")+ggsave("Doublet_calls_FeatureScatter.png")
DimPlot(rna,group.by = "Doublet.Call")+ggsave("Doublet_calls_UMAP.png")
saveRDS(doublets,"Intersection_DoubletDecon_DoubletFinder_doublets.rds")
cells <- colnames(rna)
##################################################################################################################

###########################################################
# Part 3: scRNA-seq processing after doublet removal
###########################################################
setwd(dir)

# Load the RNA dataset
counts.init <- Read10X_h5(filename= "./filtered_feature_bc_matrix.h5")
# Initialize the Seurat object with the raw (non-normalized data).
rna <- CreateSeuratObject(counts = counts.init, min.cells = 3)# genes not present in at least 3 cells are removed
rna

rna <- rna[,cells]
rna <- rna[,!(colnames(rna) %in% doublets)]
dim(rna)


#Normalize
rna <- NormalizeData(rna, normalization.method = "LogNormalize", scale.factor = 10000)

#Feature Selection
rna <- FindVariableFeatures(rna, selection.method = "vst", nfeatures = 2000)

#Scaling
all.genes <- rownames(rna)
rna <- ScaleData(rna, features = all.genes)
rna <- RunPCA(rna)


feats <- list(stromal,immune,fibroblast,endothelial,epithelial,smooth,plasma)

rna <- AddModuleScore(rna,features = feats,name = c("stromal.","immune.","fibroblast.","endothelial.",
                                                    "epithelial.","smooth.","plasma."),search = T)


# Add PC1 to metadata
rna@meta.data$PC1 <- rna@reductions$pca@cell.embeddings[,1]

count_cor_PC1 <- cor(rna$PC1,rna$nCount_RNA,method = "spearman")

stromal.cor <- cor(rna$PC1,rna$stromal.1,method = "spearman")
immune.cor <- cor(rna$PC1,rna$immune.2,method = "spearman")
fibroblast.cor <- cor(rna$PC1,rna$fibroblast.3,method = "spearman")
endothelial.cor <- cor(rna$PC1,rna$endothelial.4,method = "spearman")
epithelial.cor <- cor(rna$PC1,rna$epithelial.5,method = "spearman")
smooth.cor <- cor(rna$PC1,rna$smooth.6,method = "spearman")
plasma.cor <- cor(rna$PC1,rna$plasma.7,method = "spearman")


# Make JackStraw Plot
rna <- JackStraw(rna, num.replicate = 100,dims = 50)
rna <- ScoreJackStraw(rna,dims = 1:50)
JackStrawPlot(rna, dims = 1:50)+ggsave("JackStraw_postdoublet.png")


# If PC1 is correalted with read depth, check to see if biological variation is corralted to PC1
if (round(abs(count_cor_PC1),2) > 0.5){
  
  if( round(abs(stromal.cor),2) >= 0.5 |
      round(abs(immune.cor),2) >= 0.5 |
      round(abs(fibroblast.cor),2) >= 0.5 |
      round(abs(endothelial.cor),2) >= 0.5 |
      round(abs(epithelial.cor),2) >= 0.5 |
      round(abs(smooth.cor),2) >= 0.5 |
      round(abs(plasma.cor),2) >= 0.5){
    
    rna <- FindNeighbors(rna,dims = 1:50)
    rna <- FindClusters(rna,resolution = 0.7)
    rna <- RunUMAP(rna,dims = 1:50)
    Idents(rna) <- "RNA_snn_res.0.7"
    
    
    # Verify with inferCNV: is PC1 correlated with CNV events/Malignancy?
    #########################################################################
    # inferCNV: does PC1 also correlated with CNV/malignancy status?
    library(infercnv)
    library(stringr)
    library(Seurat)
    counts_matrix = GetAssayData(rna, slot="counts")
    
    # Identify immune clusters
    #######################################################
    # Find immune cells by relative enrichment of ESTIMATE immune signature
    library(psych)
    test <- VlnPlot(rna,features = "immune.2")
    data <- describeBy(test$data$immune.2, test$data$ident, mat = TRUE)
    data.immune <- dplyr::filter(data,median > 0.1)
    
    test <- VlnPlot(rna,features = "plasma.7")
    data <- describeBy(test$data$plasma.7, test$data$ident, mat = TRUE)
    data.plasma <- dplyr::filter(data,median > 0.1)
    
    immune.clusters <- intersect(data.immune$group1,levels(Idents(rna)))
    plasma.clusters <- intersect(data.plasma$group1,levels(Idents(rna)))
    
    immune.clusters <- unique(append(immune.clusters,plasma.clusters))
    
    for (i in 1:length(immune.clusters)){
      j <- which(levels(Idents(rna)) == immune.clusters[i])
      levels(Idents(rna))[j] <- paste0("immune.",immune.clusters[i])
    }
    rna@meta.data$postdoublet.idents <- Idents(rna)
    idents <- data.frame(rownames(rna@meta.data),rna@meta.data$postdoublet.idents)
    
    
    colnames(idents) <- c("V1","V2")
    saveRDS(rna,"./rna_postdoublet_preinferCNV.rds")
    # Make inferCNV inputs
    
    rownames(idents) <- NULL
    colnames(idents) <- NULL
    write.table(idents,"./sample_annotation_file_inferCNV.txt",sep = "\t",row.names = FALSE)
    idents <- read.delim("./sample_annotation_file_inferCNV.txt",header = F)
    
    
    
    gtf <- read.delim(GRCH38.annotations,header = F)
    library(EnsDb.Hsapiens.v86)
    convert.symbol = function(Data){
      ensembls <- Data$V1
      ensembls <- gsub("\\.[0-9]*$", "", ensembls)
      geneIDs1 <- ensembldb::select(EnsDb.Hsapiens.v86, keys= ensembls, keytype = "GENEID", columns = "SYMBOL")
      Data <- cbind.fill(Data, geneIDs1, fill = NA)
      Data <- na.omit(Data)
      Data$feature <- Data$SYMBOL
      Data.new <- data.frame(Data$SYMBOL,Data$V2,Data$V3,Data$V4)
      Data.new$Data.V2 <- paste("chr",Data.new$Data.V2,sep = "")
      Data.new$Data.SYMBOL <- make.unique(Data.new$Data.SYMBOL)
      return(Data.new)
    }
    gtf <- convert.symbol(gtf)
    head(gtf)
    
    write.table(gtf,"./Homo_sapiens.GRCh38.86.symbol.txt",sep = "\t",row.names = FALSE,col.names = FALSE)
    
    
    num.immune.clusters = length(immune.clusters)
    # create the infercnv object
    if ( num.immune.clusters == 1) {
      infercnv_obj = CreateInfercnvObject(raw_counts_matrix=as.matrix(counts_matrix),
                                          annotations_file="./sample_annotation_file_inferCNV.txt",
                                          gene_order_file="./Homo_sapiens.GRCh38.86.symbol.txt",
                                          ref_group_names=paste0("immune.",immune.clusters[1]))
      
    } else if (num.immune.clusters == 2) {
      infercnv_obj = CreateInfercnvObject(raw_counts_matrix=as.matrix(counts_matrix),
                                          annotations_file="./sample_annotation_file_inferCNV.txt",
                                          gene_order_file="./Homo_sapiens.GRCh38.86.symbol.txt",
                                          ref_group_names=c(paste0("immune.",immune.clusters[1]),
                                                            paste0("immune.",immune.clusters[2])))
      
    } else if ( num.immune.clusters == 3 ) {
      infercnv_obj = CreateInfercnvObject(raw_counts_matrix=as.matrix(counts_matrix),
                                          annotations_file="./sample_annotation_file_inferCNV.txt",
                                          gene_order_file="./Homo_sapiens.GRCh38.86.symbol.txt",
                                          ref_group_names=c(paste0("immune.",immune.clusters[1]),
                                                            paste0("immune.",immune.clusters[2]),
                                                            paste0("immune.",immune.clusters[3])))
    } else if (num.immune.clusters == 4) {
      infercnv_obj = CreateInfercnvObject(raw_counts_matrix=as.matrix(counts_matrix),
                                          annotations_file="./sample_annotation_file_inferCNV.txt",
                                          gene_order_file="./Homo_sapiens.GRCh38.86.symbol.txt",
                                          ref_group_names=c(paste0("immune.",immune.clusters[1]),
                                                            paste0("immune.",immune.clusters[2]),
                                                            paste0("immune.",immune.clusters[3]),
                                                            paste0("immune.",immune.clusters[4])))
    } else if (num.immune.clusters == 5) {
      infercnv_obj = CreateInfercnvObject(raw_counts_matrix=as.matrix(counts_matrix),
                                          annotations_file="./sample_annotation_file_inferCNV.txt",
                                          gene_order_file="./Homo_sapiens.GRCh38.86.symbol.txt",
                                          ref_group_names=c(paste0("immune.",immune.clusters[1]),
                                                            paste0("immune.",immune.clusters[2]),
                                                            paste0("immune.",immune.clusters[3]),
                                                            paste0("immune.",immune.clusters[4]),
                                                            paste0("immune.",immune.clusters[5])))
    } else if (num.immune.clusters == 6) {
      infercnv_obj = CreateInfercnvObject(raw_counts_matrix=as.matrix(counts_matrix),
                                          annotations_file="./sample_annotation_file_inferCNV.txt",
                                          gene_order_file="./Homo_sapiens.GRCh38.86.symbol.txt",
                                          ref_group_names=c(paste0("immune.",immune.clusters[1]),
                                                            paste0("immune.",immune.clusters[2]),
                                                            paste0("immune.",immune.clusters[3]),
                                                            paste0("immune.",immune.clusters[4]),
                                                            paste0("immune.",immune.clusters[5]),
                                                            paste0("immune.",immune.clusters[6])))
    }else if (num.immune.clusters == 7) {
      infercnv_obj = CreateInfercnvObject(raw_counts_matrix=as.matrix(counts_matrix),
                                          annotations_file="./sample_annotation_file_inferCNV.txt",
                                          gene_order_file="./Homo_sapiens.GRCh38.86.symbol.txt",
                                          ref_group_names=c(paste0("immune.",immune.clusters[1]),
                                                            paste0("immune.",immune.clusters[2]),
                                                            paste0("immune.",immune.clusters[3]),
                                                            paste0("immune.",immune.clusters[4]),
                                                            paste0("immune.",immune.clusters[5]),
                                                            paste0("immune.",immune.clusters[6]),
                                                            paste0("immune.",immune.clusters[7])))
    }else if (num.immune.clusters == 8) {
      infercnv_obj = CreateInfercnvObject(raw_counts_matrix=as.matrix(counts_matrix),
                                          annotations_file="./sample_annotation_file_inferCNV.txt",
                                          gene_order_file="./Homo_sapiens.GRCh38.86.symbol.txt",
                                          ref_group_names=c(paste0("immune.",immune.clusters[1]),
                                                            paste0("immune.",immune.clusters[2]),
                                                            paste0("immune.",immune.clusters[3]),
                                                            paste0("immune.",immune.clusters[4]),
                                                            paste0("immune.",immune.clusters[5]),
                                                            paste0("immune.",immune.clusters[6]),
                                                            paste0("immune.",immune.clusters[7]),
                                                            paste0("immune.",immune.clusters[8])))
    }else if (num.immune.clusters == 9) {
      infercnv_obj = CreateInfercnvObject(raw_counts_matrix=as.matrix(counts_matrix),
                                          annotations_file="./sample_annotation_file_inferCNV.txt",
                                          gene_order_file="./Homo_sapiens.GRCh38.86.symbol.txt",
                                          ref_group_names=c(paste0("immune.",immune.clusters[1]),
                                                            paste0("immune.",immune.clusters[2]),
                                                            paste0("immune.",immune.clusters[3]),
                                                            paste0("immune.",immune.clusters[4]),
                                                            paste0("immune.",immune.clusters[5]),
                                                            paste0("immune.",immune.clusters[6]),
                                                            paste0("immune.",immune.clusters[7]),
                                                            paste0("immune.",immune.clusters[8]),
                                                            paste0("immune.",immune.clusters[9])))
    }else if (num.immune.clusters == 10) {
      infercnv_obj = CreateInfercnvObject(raw_counts_matrix=as.matrix(counts_matrix),
                                          annotations_file="./sample_annotation_file_inferCNV.txt",
                                          gene_order_file="./Homo_sapiens.GRCh38.86.symbol.txt",
                                          ref_group_names=c(paste0("immune.",immune.clusters[1]),
                                                            paste0("immune.",immune.clusters[2]),
                                                            paste0("immune.",immune.clusters[3]),
                                                            paste0("immune.",immune.clusters[4]),
                                                            paste0("immune.",immune.clusters[5]),
                                                            paste0("immune.",immune.clusters[6]),
                                                            paste0("immune.",immune.clusters[7]),
                                                            paste0("immune.",immune.clusters[8]),
                                                            paste0("immune.",immune.clusters[9]),
                                                            paste0("immune.",immune.clusters[10])))
    }
    # perform infercnv operations to reveal cnv signal
    infercnv_obj = infercnv::run(infercnv_obj,
                                 cutoff=0.1,  # use 1 for smart-seq, 0.1 for 10x-genomics
                                 out_dir="./output_dir_CNV_postdoublet_PassedPC1Checks",  # dir is auto-created for storing outputs
                                 cluster_by_groups=T,   # cluster
                                 denoise=T,scale_data = T,
                                 HMM=T,HMM_type = "i6",analysis_mode = "samples",min_cells_per_gene = 10,
                                 BayesMaxPNormal = 0.4, num_threads = 8
                                 
    )
    
    
    
    regions <- read.delim("./output_dir_CNV_postdoublet_PassedPC1Checks/HMM_CNV_predictions.HMMi6.hmm_mode-samples.Pnorm_0.4.pred_cnv_regions.dat")
    probs <- read.delim("./output_dir_CNV_postdoublet_PassedPC1Checks/BayesNetOutput.HMMi6.hmm_mode-samples/CNV_State_Probabilities.dat")
    
    probs <- as.data.frame(t(probs[3,]))
    colnames(probs) <- "Prob.Normal"
    probs <- dplyr::filter(probs,Prob.Normal < 0.05)
    cnvs <- rownames(probs)
    cnvs <- gsub("\\.","-",cnvs)
    
    regions <- regions[regions$cnv_name %in% cnvs, ]
    
    
    cnv.groups <- sub("\\..*", "", regions$cell_group_name)
    
    
    length(which(rownames(rna@reductions$pca@cell.embeddings) == rownames(rna@meta.data)))
    rna$PC1.loading <- rna@reductions$pca@cell.embeddings[,1]
    rna$cell.barcode <- rownames(rna@meta.data)
    rna$CNV.Pos <- ifelse(as.character(rna$postdoublet.idents) %in% cnv.groups,1,0)
    
    
    cnv.freq <- data.frame(table(regions$cell_group_name))
    cnv.freq$Var1 <- sub("\\..*", "", cnv.freq$Var1)
    
    rna$Total_CNVs <- ifelse(as.character(rna$postdoublet.idents) %in% cnv.freq$Var1,cnv.freq$Freq,0)
    
    boxplot.cnv <- ggplot(rna@meta.data,aes(x= postdoublet.idents,y=PC1.loading,color = as.factor(CNV.Pos)))+geom_boxplot()
    boxplot.cnv+ggsave("postdoublet_CNV_PC1_boxplot.png")
    
    data <- describeBy(boxplot.cnv$data$PC1.loading, boxplot.cnv$data$postdoublet.idents, mat = TRUE)
    data$CNV <- ifelse(data$group1 %in% cnv.groups,1,0)
    
    wilcox <- wilcox.test(data = rna@meta.data,PC1.loading~CNV.Pos)
    
    if (wilcox$p.value < 0.05){
      rna <- rna
      library(stringr)
      levels(Idents(rna)) <- str_remove(levels(Idents(rna)),"immune.")
      saveRDS(rna,"./rna_postdoublet_PassedPC1Checks.rds")
    }else{
      all.genes <- rownames(rna)
      rna <- ScaleData(rna, features = all.genes,vars.to.regress = "nCount_RNA")
      rna <- FindNeighbors(rna,dims = 1:50)
      rna <- FindClusters(rna,resolution = 0.7)
      rna <- RunUMAP(rna,dims = 1:50)
      Idents(rna) <- "RNA_snn_res.0.7"
      
      library(infercnv)
      library(stringr)
      library(Seurat)
      counts_matrix = GetAssayData(rna, slot="counts")
      
      
      # Identify immune clusters
      #######################################################
      # Find immune cells by relative enrichment of ESTIMATE immune signature
      library(psych)
      test <- VlnPlot(rna,features = "immune.2")
      data <- describeBy(test$data$immune.2, test$data$ident, mat = TRUE)
      data.immune <- dplyr::filter(data,median > 0.1)
      
      test <- VlnPlot(rna,features = "plasma.7")
      data <- describeBy(test$data$plasma.7, test$data$ident, mat = TRUE)
      data.plasma <- dplyr::filter(data,median > 0.1)
      
      immune.clusters <- intersect(data.immune$group1,levels(Idents(rna)))
      plasma.clusters <- intersect(data.plasma$group1,levels(Idents(rna)))
      
      immune.clusters <- unique(append(immune.clusters,plasma.clusters))
      
      for (i in 1:length(immune.clusters)){
        j <- which(levels(Idents(rna)) == immune.clusters[i])
        levels(Idents(rna))[j] <- paste0("immune.",immune.clusters[i])
      }
      rna@meta.data$postdoublet.idents <- Idents(rna)
      idents <- data.frame(rownames(rna@meta.data),rna@meta.data$postdoublet.idents)
      saveRDS(rna,"./rna_postdoublet_preinferCNV.rds")
      #MAke inferCNV input
      
      colnames(idents) <- c("V1","V2")
      
      
      rownames(idents) <- NULL
      colnames(idents) <- NULL
      write.table(idents,"./sample_annotation_file_inferCNV.txt",sep = "\t",row.names = FALSE)
      idents <- read.delim("./sample_annotation_file_inferCNV.txt",header = F)
      
      
      
      gtf <- read.delim(GRCH38.annotations,header = F)
      library(EnsDb.Hsapiens.v86)
      convert.symbol = function(Data){
        ensembls <- Data$V1
        ensembls <- gsub("\\.[0-9]*$", "", ensembls)
        geneIDs1 <- ensembldb::select(EnsDb.Hsapiens.v86, keys= ensembls, keytype = "GENEID", columns = "SYMBOL")
        Data <- cbind.fill(Data, geneIDs1, fill = NA)
        Data <- na.omit(Data)
        Data$feature <- Data$SYMBOL
        Data.new <- data.frame(Data$SYMBOL,Data$V2,Data$V3,Data$V4)
        Data.new$Data.V2 <- paste("chr",Data.new$Data.V2,sep = "")
        Data.new$Data.SYMBOL <- make.unique(Data.new$Data.SYMBOL)
        return(Data.new)
      }
      gtf <- convert.symbol(gtf)
      head(gtf)
      
      write.table(gtf,"./Homo_sapiens.GRCh38.86.symbol.txt",sep = "\t",row.names = FALSE,col.names = FALSE)
      
      
      num.immune.clusters = length(immune.clusters)
      # create the infercnv object
      if ( num.immune.clusters == 1) {
        infercnv_obj = CreateInfercnvObject(raw_counts_matrix=as.matrix(counts_matrix),
                                            annotations_file="./sample_annotation_file_inferCNV.txt",
                                            gene_order_file="./Homo_sapiens.GRCh38.86.symbol.txt",
                                            ref_group_names=paste0("immune.",immune.clusters[1]))
        
      } else if (num.immune.clusters == 2) {
        infercnv_obj = CreateInfercnvObject(raw_counts_matrix=as.matrix(counts_matrix),
                                            annotations_file="./sample_annotation_file_inferCNV.txt",
                                            gene_order_file="./Homo_sapiens.GRCh38.86.symbol.txt",
                                            ref_group_names=c(paste0("immune.",immune.clusters[1]),
                                                              paste0("immune.",immune.clusters[2])))
        
      } else if ( num.immune.clusters == 3 ) {
        infercnv_obj = CreateInfercnvObject(raw_counts_matrix=as.matrix(counts_matrix),
                                            annotations_file="./sample_annotation_file_inferCNV.txt",
                                            gene_order_file="./Homo_sapiens.GRCh38.86.symbol.txt",
                                            ref_group_names=c(paste0("immune.",immune.clusters[1]),
                                                              paste0("immune.",immune.clusters[2]),
                                                              paste0("immune.",immune.clusters[3])))
      } else if (num.immune.clusters == 4) {
        infercnv_obj = CreateInfercnvObject(raw_counts_matrix=as.matrix(counts_matrix),
                                            annotations_file="./sample_annotation_file_inferCNV.txt",
                                            gene_order_file="./Homo_sapiens.GRCh38.86.symbol.txt",
                                            ref_group_names=c(paste0("immune.",immune.clusters[1]),
                                                              paste0("immune.",immune.clusters[2]),
                                                              paste0("immune.",immune.clusters[3]),
                                                              paste0("immune.",immune.clusters[4])))
      } else if (num.immune.clusters == 5) {
        infercnv_obj = CreateInfercnvObject(raw_counts_matrix=as.matrix(counts_matrix),
                                            annotations_file="./sample_annotation_file_inferCNV.txt",
                                            gene_order_file="./Homo_sapiens.GRCh38.86.symbol.txt",
                                            ref_group_names=c(paste0("immune.",immune.clusters[1]),
                                                              paste0("immune.",immune.clusters[2]),
                                                              paste0("immune.",immune.clusters[3]),
                                                              paste0("immune.",immune.clusters[4]),
                                                              paste0("immune.",immune.clusters[5])))
      } else if (num.immune.clusters == 6) {
        infercnv_obj = CreateInfercnvObject(raw_counts_matrix=as.matrix(counts_matrix),
                                            annotations_file="./sample_annotation_file_inferCNV.txt",
                                            gene_order_file="./Homo_sapiens.GRCh38.86.symbol.txt",
                                            ref_group_names=c(paste0("immune.",immune.clusters[1]),
                                                              paste0("immune.",immune.clusters[2]),
                                                              paste0("immune.",immune.clusters[3]),
                                                              paste0("immune.",immune.clusters[4]),
                                                              paste0("immune.",immune.clusters[5]),
                                                              paste0("immune.",immune.clusters[6])))
      }else if (num.immune.clusters == 7) {
        infercnv_obj = CreateInfercnvObject(raw_counts_matrix=as.matrix(counts_matrix),
                                            annotations_file="./sample_annotation_file_inferCNV.txt",
                                            gene_order_file="./Homo_sapiens.GRCh38.86.symbol.txt",
                                            ref_group_names=c(paste0("immune.",immune.clusters[1]),
                                                              paste0("immune.",immune.clusters[2]),
                                                              paste0("immune.",immune.clusters[3]),
                                                              paste0("immune.",immune.clusters[4]),
                                                              paste0("immune.",immune.clusters[5]),
                                                              paste0("immune.",immune.clusters[6]),
                                                              paste0("immune.",immune.clusters[7])))
      }else if (num.immune.clusters == 8) {
        infercnv_obj = CreateInfercnvObject(raw_counts_matrix=as.matrix(counts_matrix),
                                            annotations_file="./sample_annotation_file_inferCNV.txt",
                                            gene_order_file="./Homo_sapiens.GRCh38.86.symbol.txt",
                                            ref_group_names=c(paste0("immune.",immune.clusters[1]),
                                                              paste0("immune.",immune.clusters[2]),
                                                              paste0("immune.",immune.clusters[3]),
                                                              paste0("immune.",immune.clusters[4]),
                                                              paste0("immune.",immune.clusters[5]),
                                                              paste0("immune.",immune.clusters[6]),
                                                              paste0("immune.",immune.clusters[7]),
                                                              paste0("immune.",immune.clusters[8])))
      }else if (num.immune.clusters == 9) {
        infercnv_obj = CreateInfercnvObject(raw_counts_matrix=as.matrix(counts_matrix),
                                            annotations_file="./sample_annotation_file_inferCNV.txt",
                                            gene_order_file="./Homo_sapiens.GRCh38.86.symbol.txt",
                                            ref_group_names=c(paste0("immune.",immune.clusters[1]),
                                                              paste0("immune.",immune.clusters[2]),
                                                              paste0("immune.",immune.clusters[3]),
                                                              paste0("immune.",immune.clusters[4]),
                                                              paste0("immune.",immune.clusters[5]),
                                                              paste0("immune.",immune.clusters[6]),
                                                              paste0("immune.",immune.clusters[7]),
                                                              paste0("immune.",immune.clusters[8]),
                                                              paste0("immune.",immune.clusters[9])))
      }else if (num.immune.clusters == 10) {
        infercnv_obj = CreateInfercnvObject(raw_counts_matrix=as.matrix(counts_matrix),
                                            annotations_file="./sample_annotation_file_inferCNV.txt",
                                            gene_order_file="./Homo_sapiens.GRCh38.86.symbol.txt",
                                            ref_group_names=c(paste0("immune.",immune.clusters[1]),
                                                              paste0("immune.",immune.clusters[2]),
                                                              paste0("immune.",immune.clusters[3]),
                                                              paste0("immune.",immune.clusters[4]),
                                                              paste0("immune.",immune.clusters[5]),
                                                              paste0("immune.",immune.clusters[6]),
                                                              paste0("immune.",immune.clusters[7]),
                                                              paste0("immune.",immune.clusters[8]),
                                                              paste0("immune.",immune.clusters[9]),
                                                              paste0("immune.",immune.clusters[10])))
      }
      # perform infercnv operations to reveal cnv signal
      infercnv_obj = infercnv::run(infercnv_obj,
                                   cutoff=0.1,  # use 1 for smart-seq, 0.1 for 10x-genomics
                                   out_dir="./output_dir_CNV_postdoublet_FailedCNVTest",  # dir is auto-created for storing outputs
                                   cluster_by_groups=T,   # cluster
                                   denoise=T,scale_data = T,
                                   HMM=T,HMM_type = "i6",analysis_mode = "samples",min_cells_per_gene = 10,
                                   BayesMaxPNormal = 0.4, num_threads = 8
                                   
      )
      
      
      
      regions <- read.delim("./output_dir_CNV_postdoublet_FailedCNVTest/HMM_CNV_predictions.HMMi6.hmm_mode-samples.Pnorm_0.4.pred_cnv_regions.dat")
      probs <- read.delim("./output_dir_CNV_postdoublet_FailedCNVTest/BayesNetOutput.HMMi6.hmm_mode-samples/CNV_State_Probabilities.dat")
      
      probs <- as.data.frame(t(probs[3,]))
      colnames(probs) <- "Prob.Normal"
      probs <- dplyr::filter(probs,Prob.Normal < 0.05)
      cnvs <- rownames(probs)
      cnvs <- gsub("\\.","-",cnvs)
      
      regions <- regions[regions$cnv_name %in% cnvs, ]
      
      
      cnv.groups <- sub("\\..*", "", regions$cell_group_name)
      
      
      length(which(rownames(rna@reductions$pca@cell.embeddings) == rownames(rna@meta.data)))
      rna$PC1.loading <- rna@reductions$pca@cell.embeddings[,1]
      rna$cell.barcode <- rownames(rna@meta.data)
      rna$CNV.Pos <- ifelse(as.character(rna$postdoublet.idents) %in% cnv.groups,1,0)
      
      
      cnv.freq <- data.frame(table(regions$cell_group_name))
      cnv.freq$Var1 <- sub("\\..*", "", cnv.freq$Var1)
      
      rna$Total_CNVs <- ifelse(as.character(rna$postdoublet.idents) %in% cnv.freq$Var1,cnv.freq$Freq,0)
      
      saveRDS(rna,"./rna_postdoublet_FailedCNVTest.rds")
    }

  }else{
    all.genes <- rownames(rna)
    rna <- ScaleData(rna, features = all.genes,vars.to.regress = "nCount_RNA")
    rna <- FindNeighbors(rna,dims = 1:50)
    rna <- FindClusters(rna,resolution = 0.7)
    rna <- RunUMAP(rna,dims = 1:50)
    Idents(rna) <- "RNA_snn_res.0.7"
    
    # Proceed with CNV
    library(infercnv)
    library(stringr)
    library(Seurat)
    counts_matrix = GetAssayData(rna, slot="counts")
    
    # Identify immune clusters
    #######################################################
    # Find immune cells by relative enrichment of ESTIMATE immune signature
    library(psych)
    test <- VlnPlot(rna,features = "immune.2")
    data <- describeBy(test$data$immune.2, test$data$ident, mat = TRUE)
    data.immune <- dplyr::filter(data,median > 0.1)
    
    test <- VlnPlot(rna,features = "plasma.7")
    data <- describeBy(test$data$plasma.7, test$data$ident, mat = TRUE)
    data.plasma <- dplyr::filter(data,median > 0.1)
    
    immune.clusters <- intersect(data.immune$group1,levels(Idents(rna)))
    plasma.clusters <- intersect(data.plasma$group1,levels(Idents(rna)))
    
    immune.clusters <- unique(append(immune.clusters,plasma.clusters))
    
    for (i in 1:length(immune.clusters)){
      j <- which(levels(Idents(rna)) == immune.clusters[i])
      levels(Idents(rna))[j] <- paste0("immune.",immune.clusters[i])
    }
    rna@meta.data$postdoublet.idents <- Idents(rna)
    idents <- data.frame(rownames(rna@meta.data),rna@meta.data$postdoublet.idents)
    
    
    
    colnames(idents) <- c("V1","V2")
    
    saveRDS(rna,"./rna_postdoublet_preinferCNV.rds")
    # Make inferCNV inputs
    
    rownames(idents) <- NULL
    colnames(idents) <- NULL
    write.table(idents,"./sample_annotation_file_inferCNV.txt",sep = "\t",row.names = FALSE)
    idents <- read.delim("./sample_annotation_file_inferCNV.txt",header = F)
    
    
    
    gtf <- read.delim(GRCH38.annotations,header = F)
    library(EnsDb.Hsapiens.v86)
    convert.symbol = function(Data){
      ensembls <- Data$V1
      ensembls <- gsub("\\.[0-9]*$", "", ensembls)
      geneIDs1 <- ensembldb::select(EnsDb.Hsapiens.v86, keys= ensembls, keytype = "GENEID", columns = "SYMBOL")
      Data <- cbind.fill(Data, geneIDs1, fill = NA)
      Data <- na.omit(Data)
      Data$feature <- Data$SYMBOL
      Data.new <- data.frame(Data$SYMBOL,Data$V2,Data$V3,Data$V4)
      Data.new$Data.V2 <- paste("chr",Data.new$Data.V2,sep = "")
      Data.new$Data.SYMBOL <- make.unique(Data.new$Data.SYMBOL)
      return(Data.new)
    }
    gtf <- convert.symbol(gtf)
    head(gtf)
    
    write.table(gtf,"./Homo_sapiens.GRCh38.86.symbol.txt",sep = "\t",row.names = FALSE,col.names = FALSE)
    
    
    num.immune.clusters = length(immune.clusters)
    # create the infercnv object
    if ( num.immune.clusters == 1) {
      infercnv_obj = CreateInfercnvObject(raw_counts_matrix=as.matrix(counts_matrix),
                                          annotations_file="./sample_annotation_file_inferCNV.txt",
                                          gene_order_file="./Homo_sapiens.GRCh38.86.symbol.txt",
                                          ref_group_names=paste0("immune.",immune.clusters[1]))
      
    } else if (num.immune.clusters == 2) {
      infercnv_obj = CreateInfercnvObject(raw_counts_matrix=as.matrix(counts_matrix),
                                          annotations_file="./sample_annotation_file_inferCNV.txt",
                                          gene_order_file="./Homo_sapiens.GRCh38.86.symbol.txt",
                                          ref_group_names=c(paste0("immune.",immune.clusters[1]),
                                                            paste0("immune.",immune.clusters[2])))
      
    } else if ( num.immune.clusters == 3 ) {
      infercnv_obj = CreateInfercnvObject(raw_counts_matrix=as.matrix(counts_matrix),
                                          annotations_file="./sample_annotation_file_inferCNV.txt",
                                          gene_order_file="./Homo_sapiens.GRCh38.86.symbol.txt",
                                          ref_group_names=c(paste0("immune.",immune.clusters[1]),
                                                            paste0("immune.",immune.clusters[2]),
                                                            paste0("immune.",immune.clusters[3])))
    } else if (num.immune.clusters == 4) {
      infercnv_obj = CreateInfercnvObject(raw_counts_matrix=as.matrix(counts_matrix),
                                          annotations_file="./sample_annotation_file_inferCNV.txt",
                                          gene_order_file="./Homo_sapiens.GRCh38.86.symbol.txt",
                                          ref_group_names=c(paste0("immune.",immune.clusters[1]),
                                                            paste0("immune.",immune.clusters[2]),
                                                            paste0("immune.",immune.clusters[3]),
                                                            paste0("immune.",immune.clusters[4])))
    } else if (num.immune.clusters == 5) {
      infercnv_obj = CreateInfercnvObject(raw_counts_matrix=as.matrix(counts_matrix),
                                          annotations_file="./sample_annotation_file_inferCNV.txt",
                                          gene_order_file="./Homo_sapiens.GRCh38.86.symbol.txt",
                                          ref_group_names=c(paste0("immune.",immune.clusters[1]),
                                                            paste0("immune.",immune.clusters[2]),
                                                            paste0("immune.",immune.clusters[3]),
                                                            paste0("immune.",immune.clusters[4]),
                                                            paste0("immune.",immune.clusters[5])))
    } else if (num.immune.clusters == 6) {
      infercnv_obj = CreateInfercnvObject(raw_counts_matrix=as.matrix(counts_matrix),
                                          annotations_file="./sample_annotation_file_inferCNV.txt",
                                          gene_order_file="./Homo_sapiens.GRCh38.86.symbol.txt",
                                          ref_group_names=c(paste0("immune.",immune.clusters[1]),
                                                            paste0("immune.",immune.clusters[2]),
                                                            paste0("immune.",immune.clusters[3]),
                                                            paste0("immune.",immune.clusters[4]),
                                                            paste0("immune.",immune.clusters[5]),
                                                            paste0("immune.",immune.clusters[6])))
    }else if (num.immune.clusters == 7) {
      infercnv_obj = CreateInfercnvObject(raw_counts_matrix=as.matrix(counts_matrix),
                                          annotations_file="./sample_annotation_file_inferCNV.txt",
                                          gene_order_file="./Homo_sapiens.GRCh38.86.symbol.txt",
                                          ref_group_names=c(paste0("immune.",immune.clusters[1]),
                                                            paste0("immune.",immune.clusters[2]),
                                                            paste0("immune.",immune.clusters[3]),
                                                            paste0("immune.",immune.clusters[4]),
                                                            paste0("immune.",immune.clusters[5]),
                                                            paste0("immune.",immune.clusters[6]),
                                                            paste0("immune.",immune.clusters[7])))
    }else if (num.immune.clusters == 8) {
      infercnv_obj = CreateInfercnvObject(raw_counts_matrix=as.matrix(counts_matrix),
                                          annotations_file="./sample_annotation_file_inferCNV.txt",
                                          gene_order_file="./Homo_sapiens.GRCh38.86.symbol.txt",
                                          ref_group_names=c(paste0("immune.",immune.clusters[1]),
                                                            paste0("immune.",immune.clusters[2]),
                                                            paste0("immune.",immune.clusters[3]),
                                                            paste0("immune.",immune.clusters[4]),
                                                            paste0("immune.",immune.clusters[5]),
                                                            paste0("immune.",immune.clusters[6]),
                                                            paste0("immune.",immune.clusters[7]),
                                                            paste0("immune.",immune.clusters[8])))
    }else if (num.immune.clusters == 9) {
      infercnv_obj = CreateInfercnvObject(raw_counts_matrix=as.matrix(counts_matrix),
                                          annotations_file="./sample_annotation_file_inferCNV.txt",
                                          gene_order_file="./Homo_sapiens.GRCh38.86.symbol.txt",
                                          ref_group_names=c(paste0("immune.",immune.clusters[1]),
                                                            paste0("immune.",immune.clusters[2]),
                                                            paste0("immune.",immune.clusters[3]),
                                                            paste0("immune.",immune.clusters[4]),
                                                            paste0("immune.",immune.clusters[5]),
                                                            paste0("immune.",immune.clusters[6]),
                                                            paste0("immune.",immune.clusters[7]),
                                                            paste0("immune.",immune.clusters[8]),
                                                            paste0("immune.",immune.clusters[9])))
    }else if (num.immune.clusters == 10) {
      infercnv_obj = CreateInfercnvObject(raw_counts_matrix=as.matrix(counts_matrix),
                                          annotations_file="./sample_annotation_file_inferCNV.txt",
                                          gene_order_file="./Homo_sapiens.GRCh38.86.symbol.txt",
                                          ref_group_names=c(paste0("immune.",immune.clusters[1]),
                                                            paste0("immune.",immune.clusters[2]),
                                                            paste0("immune.",immune.clusters[3]),
                                                            paste0("immune.",immune.clusters[4]),
                                                            paste0("immune.",immune.clusters[5]),
                                                            paste0("immune.",immune.clusters[6]),
                                                            paste0("immune.",immune.clusters[7]),
                                                            paste0("immune.",immune.clusters[8]),
                                                            paste0("immune.",immune.clusters[9]),
                                                            paste0("immune.",immune.clusters[10])))
    }
    # perform infercnv operations to reveal cnv signal
    infercnv_obj = infercnv::run(infercnv_obj,
                                 cutoff=0.1,  # use 1 for smart-seq, 0.1 for 10x-genomics
                                 out_dir="./output_dir_CNV_postdoublet_FailedCorTest",  # dir is auto-created for storing outputs
                                 cluster_by_groups=T,   # cluster
                                 denoise=T,scale_data = T,
                                 HMM=T,HMM_type = "i6",analysis_mode = "samples",min_cells_per_gene = 10,
                                 BayesMaxPNormal = 0.4, num_threads = 8
                                 
    )
    
    
    
    regions <- read.delim("./output_dir_CNV_postdoublet_FailedCorTest/HMM_CNV_predictions.HMMi6.hmm_mode-samples.Pnorm_0.4.pred_cnv_regions.dat")
    probs <- read.delim("./output_dir_CNV_postdoublet_FailedCorTest/BayesNetOutput.HMMi6.hmm_mode-samples/CNV_State_Probabilities.dat")
    
    probs <- as.data.frame(t(probs[3,]))
    colnames(probs) <- "Prob.Normal"
    probs <- dplyr::filter(probs,Prob.Normal < 0.05)
    cnvs <- rownames(probs)
    cnvs <- gsub("\\.","-",cnvs)
    
    regions <- regions[regions$cnv_name %in% cnvs, ]
    
    
    cnv.groups <- sub("\\..*", "", regions$cell_group_name)
    
    
    length(which(rownames(rna@reductions$pca@cell.embeddings) == rownames(rna@meta.data)))
    rna$PC1.loading <- rna@reductions$pca@cell.embeddings[,1]
    rna$cell.barcode <- rownames(rna@meta.data)
    rna$CNV.Pos <- ifelse(as.character(rna$postdoublet.idents) %in% cnv.groups,1,0)
    
    
    cnv.freq <- data.frame(table(regions$cell_group_name))
    cnv.freq$Var1 <- sub("\\..*", "", cnv.freq$Var1)
    
    rna$Total_CNVs <- ifelse(as.character(rna$postdoublet.idents) %in% cnv.freq$Var1,cnv.freq$Freq,0)
    saveRDS(rna,"./rna_postdoublet_FailedCorTest.rds")
  }
}else{
  
  rna <- FindNeighbors(rna,dims = 1:50)
  rna <- FindClusters(rna,resolution = 0.7)
  rna <- RunUMAP(rna,dims = 1:50)
  Idents(rna) <- "RNA_snn_res.0.7"
  
  # Proceed with CNV
  library(infercnv)
  library(stringr)
  library(Seurat)
  counts_matrix = GetAssayData(rna, slot="counts")
  
  
  # Identify immune clusters
  #######################################################
  # Find immune cells by relative enrichment of ESTIMATE immune signature
  library(psych)
  test <- VlnPlot(rna,features = "immune.2")
  data <- describeBy(test$data$immune.2, test$data$ident, mat = TRUE)
  data.immune <- dplyr::filter(data,median > 0.1)
  
  test <- VlnPlot(rna,features = "plasma.7")
  data <- describeBy(test$data$plasma.7, test$data$ident, mat = TRUE)
  data.plasma <- dplyr::filter(data,median > 0.1)
  
  immune.clusters <- intersect(data.immune$group1,levels(Idents(rna)))
  plasma.clusters <- intersect(data.plasma$group1,levels(Idents(rna)))
  
  immune.clusters <- unique(append(immune.clusters,plasma.clusters))
  
  for (i in 1:length(immune.clusters)){
    j <- which(levels(Idents(rna)) == immune.clusters[i])
    levels(Idents(rna))[j] <- paste0("immune.",immune.clusters[i])
  }
  rna@meta.data$postdoublet.idents <- Idents(rna)
  idents <- data.frame(rownames(rna@meta.data),rna@meta.data$postdoublet.idents)
  
  
  colnames(idents) <- c("V1","V2")
  saveRDS(rna,"./rna_postdoublet_preinferCNV.rds")
  # Make inferCNV inputs
  rownames(idents) <- NULL
  colnames(idents) <- NULL
  write.table(idents,"./sample_annotation_file_inferCNV.txt",sep = "\t",row.names = FALSE)
  idents <- read.delim("./sample_annotation_file_inferCNV.txt",header = F)
  
  
  
  gtf <- read.delim(GRCH38.annotations,header = F)
  library(EnsDb.Hsapiens.v86)
  convert.symbol = function(Data){
    ensembls <- Data$V1
    ensembls <- gsub("\\.[0-9]*$", "", ensembls)
    geneIDs1 <- ensembldb::select(EnsDb.Hsapiens.v86, keys= ensembls, keytype = "GENEID", columns = "SYMBOL")
    Data <- cbind.fill(Data, geneIDs1, fill = NA)
    Data <- na.omit(Data)
    Data$feature <- Data$SYMBOL
    Data.new <- data.frame(Data$SYMBOL,Data$V2,Data$V3,Data$V4)
    Data.new$Data.V2 <- paste("chr",Data.new$Data.V2,sep = "")
    Data.new$Data.SYMBOL <- make.unique(Data.new$Data.SYMBOL)
    return(Data.new)
  }
  gtf <- convert.symbol(gtf)
  head(gtf)
  
  write.table(gtf,"./Homo_sapiens.GRCh38.86.symbol.txt",sep = "\t",row.names = FALSE,col.names = FALSE)
  
  
  num.immune.clusters = length(immune.clusters)
  # create the infercnv object
  # create the infercnv object
  if ( num.immune.clusters == 1) {
    infercnv_obj = CreateInfercnvObject(raw_counts_matrix=as.matrix(counts_matrix),
                                        annotations_file="./sample_annotation_file_inferCNV.txt",
                                        gene_order_file="./Homo_sapiens.GRCh38.86.symbol.txt",
                                        ref_group_names=paste0("immune.",immune.clusters[1]))
    
  } else if (num.immune.clusters == 2) {
    infercnv_obj = CreateInfercnvObject(raw_counts_matrix=as.matrix(counts_matrix),
                                        annotations_file="./sample_annotation_file_inferCNV.txt",
                                        gene_order_file="./Homo_sapiens.GRCh38.86.symbol.txt",
                                        ref_group_names=c(paste0("immune.",immune.clusters[1]),
                                                          paste0("immune.",immune.clusters[2])))
    
  } else if ( num.immune.clusters == 3 ) {
    infercnv_obj = CreateInfercnvObject(raw_counts_matrix=as.matrix(counts_matrix),
                                        annotations_file="./sample_annotation_file_inferCNV.txt",
                                        gene_order_file="./Homo_sapiens.GRCh38.86.symbol.txt",
                                        ref_group_names=c(paste0("immune.",immune.clusters[1]),
                                                          paste0("immune.",immune.clusters[2]),
                                                          paste0("immune.",immune.clusters[3])))
  } else if (num.immune.clusters == 4) {
    infercnv_obj = CreateInfercnvObject(raw_counts_matrix=as.matrix(counts_matrix),
                                        annotations_file="./sample_annotation_file_inferCNV.txt",
                                        gene_order_file="./Homo_sapiens.GRCh38.86.symbol.txt",
                                        ref_group_names=c(paste0("immune.",immune.clusters[1]),
                                                          paste0("immune.",immune.clusters[2]),
                                                          paste0("immune.",immune.clusters[3]),
                                                          paste0("immune.",immune.clusters[4])))
  } else if (num.immune.clusters == 5) {
    infercnv_obj = CreateInfercnvObject(raw_counts_matrix=as.matrix(counts_matrix),
                                        annotations_file="./sample_annotation_file_inferCNV.txt",
                                        gene_order_file="./Homo_sapiens.GRCh38.86.symbol.txt",
                                        ref_group_names=c(paste0("immune.",immune.clusters[1]),
                                                          paste0("immune.",immune.clusters[2]),
                                                          paste0("immune.",immune.clusters[3]),
                                                          paste0("immune.",immune.clusters[4]),
                                                          paste0("immune.",immune.clusters[5])))
  } else if (num.immune.clusters == 6) {
    infercnv_obj = CreateInfercnvObject(raw_counts_matrix=as.matrix(counts_matrix),
                                        annotations_file="./sample_annotation_file_inferCNV.txt",
                                        gene_order_file="./Homo_sapiens.GRCh38.86.symbol.txt",
                                        ref_group_names=c(paste0("immune.",immune.clusters[1]),
                                                          paste0("immune.",immune.clusters[2]),
                                                          paste0("immune.",immune.clusters[3]),
                                                          paste0("immune.",immune.clusters[4]),
                                                          paste0("immune.",immune.clusters[5]),
                                                          paste0("immune.",immune.clusters[6])))
  }else if (num.immune.clusters == 7) {
    infercnv_obj = CreateInfercnvObject(raw_counts_matrix=as.matrix(counts_matrix),
                                        annotations_file="./sample_annotation_file_inferCNV.txt",
                                        gene_order_file="./Homo_sapiens.GRCh38.86.symbol.txt",
                                        ref_group_names=c(paste0("immune.",immune.clusters[1]),
                                                          paste0("immune.",immune.clusters[2]),
                                                          paste0("immune.",immune.clusters[3]),
                                                          paste0("immune.",immune.clusters[4]),
                                                          paste0("immune.",immune.clusters[5]),
                                                          paste0("immune.",immune.clusters[6]),
                                                          paste0("immune.",immune.clusters[7])))
  }else if (num.immune.clusters == 8) {
    infercnv_obj = CreateInfercnvObject(raw_counts_matrix=as.matrix(counts_matrix),
                                        annotations_file="./sample_annotation_file_inferCNV.txt",
                                        gene_order_file="./Homo_sapiens.GRCh38.86.symbol.txt",
                                        ref_group_names=c(paste0("immune.",immune.clusters[1]),
                                                          paste0("immune.",immune.clusters[2]),
                                                          paste0("immune.",immune.clusters[3]),
                                                          paste0("immune.",immune.clusters[4]),
                                                          paste0("immune.",immune.clusters[5]),
                                                          paste0("immune.",immune.clusters[6]),
                                                          paste0("immune.",immune.clusters[7]),
                                                          paste0("immune.",immune.clusters[8])))
  }else if (num.immune.clusters == 9) {
    infercnv_obj = CreateInfercnvObject(raw_counts_matrix=as.matrix(counts_matrix),
                                        annotations_file="./sample_annotation_file_inferCNV.txt",
                                        gene_order_file="./Homo_sapiens.GRCh38.86.symbol.txt",
                                        ref_group_names=c(paste0("immune.",immune.clusters[1]),
                                                          paste0("immune.",immune.clusters[2]),
                                                          paste0("immune.",immune.clusters[3]),
                                                          paste0("immune.",immune.clusters[4]),
                                                          paste0("immune.",immune.clusters[5]),
                                                          paste0("immune.",immune.clusters[6]),
                                                          paste0("immune.",immune.clusters[7]),
                                                          paste0("immune.",immune.clusters[8]),
                                                          paste0("immune.",immune.clusters[9])))
  }else if (num.immune.clusters == 10) {
    infercnv_obj = CreateInfercnvObject(raw_counts_matrix=as.matrix(counts_matrix),
                                        annotations_file="./sample_annotation_file_inferCNV.txt",
                                        gene_order_file="./Homo_sapiens.GRCh38.86.symbol.txt",
                                        ref_group_names=c(paste0("immune.",immune.clusters[1]),
                                                          paste0("immune.",immune.clusters[2]),
                                                          paste0("immune.",immune.clusters[3]),
                                                          paste0("immune.",immune.clusters[4]),
                                                          paste0("immune.",immune.clusters[5]),
                                                          paste0("immune.",immune.clusters[6]),
                                                          paste0("immune.",immune.clusters[7]),
                                                          paste0("immune.",immune.clusters[8]),
                                                          paste0("immune.",immune.clusters[9]),
                                                          paste0("immune.",immune.clusters[10])))
  }
  # perform infercnv operations to reveal cnv signal
  infercnv_obj = infercnv::run(infercnv_obj,
                               cutoff=0.1,  # use 1 for smart-seq, 0.1 for 10x-genomics
                               out_dir="./output_dir_CNV_postdoublet_SkipChecks",  # dir is auto-created for storing outputs
                               cluster_by_groups=T,   # cluster
                               denoise=T,scale_data = T,
                               HMM=T,HMM_type = "i6",analysis_mode = "samples",min_cells_per_gene = 10,
                               BayesMaxPNormal = 0.4, num_threads = 8
                               
  )
  
  
  
  regions <- read.delim("./output_dir_CNV_postdoublet_SkipChecks/HMM_CNV_predictions.HMMi6.hmm_mode-samples.Pnorm_0.4.pred_cnv_regions.dat")
  probs <- read.delim("./output_dir_CNV_postdoublet_SkipChecks/BayesNetOutput.HMMi6.hmm_mode-samples/CNV_State_Probabilities.dat")
  
  probs <- as.data.frame(t(probs[3,]))
  colnames(probs) <- "Prob.Normal"
  probs <- dplyr::filter(probs,Prob.Normal < 0.05)
  cnvs <- rownames(probs)
  cnvs <- gsub("\\.","-",cnvs)
  
  regions <- regions[regions$cnv_name %in% cnvs, ]
  
  
  cnv.groups <- sub("\\..*", "", regions$cell_group_name)
  
  
  length(which(rownames(rna@reductions$pca@cell.embeddings) == rownames(rna@meta.data)))
  rna$PC1.loading <- rna@reductions$pca@cell.embeddings[,1]
  rna$cell.barcode <- rownames(rna@meta.data)
  rna$CNV.Pos <- ifelse(as.character(rna$postdoublet.idents) %in% cnv.groups,1,0)
  
  
  cnv.freq <- data.frame(table(regions$cell_group_name))
  cnv.freq$Var1 <- sub("\\..*", "", cnv.freq$Var1)
  
  rna$Total_CNVs <- ifelse(as.character(rna$postdoublet.idents) %in% cnv.freq$Var1,cnv.freq$Freq,0)
  
  
  boxplot.cnv <- ggplot(rna@meta.data,aes(x= postdoublet.idents,y=PC1.loading,color = as.factor(CNV.Pos)))+geom_boxplot()
  boxplot.cnv+ggsave("Postdoublet_CNV_PC1_boxplot.png")
  saveRDS(rna,"./rna_postdoublet_SkipChecks.rds")
}

Idents(rna)<- "RNA_snn_res.0.7"

# DEG analysis with Wilcox
##########################################################################

# Wilcox
Wilcox.markers <- FindAllMarkers(object =rna, min.pct = 0.25,only.pos = F,
                                 test.use = "wilcox")
saveRDS(Wilcox.markers,"./wilcox_DEGs.rds")


###########################################################
# Part 4: SingleR cell typing 
###########################################################

# SingleR labeling of celltypes
##########################################################################
rna.sce <- as.SingleCellExperiment(rna)


# Set reference paths:
# Slyper et al. GSM4186985 (sample ID: HTAPP-624-SMP-3212)
# Processed using default Seurat parameters 
ref.data.ovar <- readRDS("./HTAPP_ovarian_reference.rds")
# SingleR annotation
#####################################################

# Read in reference datasets for SingleR annotation 

# 1) Slyper et al. Nat. Medicine 2020 scRNA-seq ovarian tumor 
ref.data.ovar <- as.SingleCellExperiment(ref.data.ovar)

# 2) Human Primary Cell Atlas Data (microarray)
ref.data.HPCA <- readRDS("./HPCA_celldex.rds")
# 
# # 3) BluePrint Encode (bulk RNA-seq)
ref.data.BED <- readRDS("./BluePrintEncode_celldex.rds")



# 2) Single-cell level label transfer: 
predictions.HPCA.sc <- SingleR(test=rna.sce, assay.type.test="logcounts", assay.type.ref="logcounts",
                               ref=ref.data.HPCA, labels=ref.data.HPCA$label.main)
predictions.BED.sc <- SingleR(test=rna.sce, assay.type.test="logcounts", assay.type.ref="logcounts",
                              ref=ref.data.BED, labels=ref.data.BED$label.main)
# Use de.method wilcox with scRNA-seq reference b/c the reference data is more sparse
predictions.ovar.sc <- SingleR(test=rna.sce, assay.type.test="logcounts", assay.type.ref="logcounts",
                               ref=ref.data.ovar, labels=ref.data.ovar$annotate,de.method = "wilcox")


rna$SingleR.HPCA <- predictions.HPCA.sc$pruned.labels
rna$SingleR.BED <- predictions.BED.sc$pruned.labels
rna$SingleR.ovar <- predictions.ovar.sc$pruned.labels


# Save Seurat object 
#date <- Sys.Date()
saveRDS(rna,paste0(SAMPLE.ID,"_scRNA_processed.rds"))
# ###########################################################################################################
# # Starting cells, PostQC cells, doublets, Post doublet/QC cells, Cluster #
output.meta <- data.frame(StartingNumCells=length(colnames(counts.init)),
                          nMADLogCounts =2,
                          nMADLogFeatures = 2,
                          nMADLog1pMito =2,
                          PostQCNumCells=PostQCNumCells,
                          ExpectedDoubletFraction=doublet.rate,
                          ObservedDoubletFraction=length(doublets)/length(colnames(counts.init)),
                          PostDoubletNumCells=length(colnames(rna)),
                          NumClusters=length(levels(Idents(rna))),
                          DoubletFinderpK = pK.1,
                          MinNumCounts=min(rna$nCount_RNA),
                          MaxNumCounts= max(rna$nCount_RNA),
                          MedianNumbCounts = median(rna$nCount_RNA),
                          MinNumFeats=min(rna$nFeature_RNA),
                          MaxNumFeats= max(rna$nFeature_RNA),
                          MedianNumbFeats = median(rna$nFeature_RNA),
                          stringsAsFactors=FALSE)
output <- as.data.frame(t(output.meta))
colnames(output) <- SAMPLE.ID
xlsx::write.xlsx(output, "scRNA_pipeline_summary.xlsx",
                 row.names = T, col.names = TRUE)




#########################################################################################################
# END OF SCRIPT
#########################################################################################################
