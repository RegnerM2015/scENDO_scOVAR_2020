###########################################################
# Matt Regner
# Franco Lab
# Date: May-December 2020
# 
# Sample: EEC
# Description: This script performs the following tasks  
#         1) scATAC-seq processing
#         2) scRNA-seq/scATAC-seq integration
#         3) Peak2Gene linkage analysis/Co-accessiblity 
#         4) Make input files for TFSEE analysis 
###########################################################
.libPaths('/home/regnerm/anaconda3/envs/r-environment/lib/R/library')

source("./Archr_Peak_Null.R")
source("./filterDoublets_modified.R")
###############################################################
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
library(Signac)
library(ggplot2)
library(stringr)
library(EnsDb.Hsapiens.v86)
library(Seurat)
library(Signac)
library(ggplot2)
library(ensembldb)
library(EnsDb.Hsapiens.v86)
library(ArchR)
library(SingleR)
library(viridis)

set.seed(2)
addArchRThreads(threads = 16) 
addArchRGenome("hg38")

# Set up directories and file variables:
##################################################################################################
SAMPLE.ID <- "EEC"

output.dir <- "."
####################
setwd(output.dir)
####################

inputFiles <- list.files(pattern = "\\.gz$")

sampleNames <- c("3533EL","3571DL","36186L","36639L","366C5L")


encode.all <- read.delim("./GRCh38-ccREs.bed",header =F)
colnames(encode.all)[1:3] <- c("seqnames","start","end")


# Read in matching scRNA-seq
###################################################################################################
rna <- readRDS("./endo_EEC_scRNA_processed.rds")

# Name cancer cell populations of interest (mostly epithelial with nuclear CNV)
labels <- c("22-Unciliated epithelia 1","10-Unciliated epithelia 2","19-Ciliated","21-Unciliated epithelia 2",
            "8-Unciliated epithelia 1", "14-Ciliated", "6-Unciliated epithelia 1","15-Unciliated epithelia 1")

# Redo differential expression with new cell type markers 
Wilcox.markers <- readRDS("./wilcox_DEGs.rds")


# Create Arrow and ArchR project
##########################################################################
ArrowFiles <- createArrowFiles(
  inputFiles = inputFiles,
  sampleNames = sampleNames,
  filterTSS = 0, #Dont set this too high because you can always increase later
  filterFrags = 0,
  addTileMat = T,
  addGeneScoreMat = F
)

doubScores <- addDoubletScores(
  input = ArrowFiles,
  k = 10, #Refers to how many cells near a "pseudo-doublet" to count.
  knnMethod = "UMAP",useMatrix = "TileMatrix",nTrials=5,LSIMethod = 1,scaleDims = F,
  corCutOff = 0.75,UMAPParams = list(n_neighbors =30, min_dist = 0.3, metric = "cosine", verbose =FALSE),
  dimsToUse = 1:50
)

proj <- ArchRProject(
  ArrowFiles = ArrowFiles,
  outputDirectory =SAMPLE.ID,
  copyArrows = T #This is recommened so that you maintain an unaltered copy for later usage.
)

# Filter out outlier low quality cells and doublets
###############################################################################
# GMM for fragments per cell
library(mclust)

for (i in sampleNames){
  proj.i <- proj[proj$Sample == i]
  
  # GMM for fragments per cell
  depth.clust <- Mclust(log10(proj.i$nFrags),G = 2)
  proj.i$depth.cluster <- depth.clust$classification
  proj.i$depth.cluster.uncertainty <- depth.clust$uncertainty
  
  ggPoint(
    x = log10(proj.i$nFrags),
    y = log10(proj.i$TSSEnrichment+1),
    color = as.character(proj.i$depth.cluster),
    xlabel = "log10(unique fragments)",
    ylabel = "log10(TSS Enrichment+1)"
  ) + ggtitle(paste0("GMM classification:\n",i," log10(fragments)"))+
    ggsave(paste0(i,"_depth.pdf"),width = 4,height = 4)
  
  # GMM for TSS per cell
  TSS.clust <- Mclust(log10(proj.i$TSSEnrichment+1),G = 2)
  proj.i$TSS.cluster <- TSS.clust$classification
  proj.i$TSS.cluster.uncertainty <- TSS.clust$uncertainty
  
  ggPoint(
    x = log10(proj.i$nFrags),
    y = log10(proj.i$TSSEnrichment+1),
    color = as.character(proj.i$TSS.cluster),
    discrete = T,
    xlabel = "log10(unique fragments)",
    ylabel = "log10(TSS Enrichment+1)"
  ) + ggtitle(paste0("GMM classification:\n",i," TSS Enrichment"))+
    ggsave(paste0(i,"_TSS.pdf"),width = 4,height = 4)
  
  
  df.TSS <- data.frame(proj.i$cellNames,proj.i$TSS.cluster,proj.i$TSS.cluster.uncertainty,proj.i$TSSEnrichment)
  df.TSS <- dplyr::filter(df.TSS,proj.i.TSS.cluster == "2")
  df.TSS <- dplyr::filter(df.TSS,proj.i.TSS.cluster.uncertainty <= 0.05)
  saveRDS(df.TSS,paste0("df_TSS_",i,".rds"))
  
  df.depth <- data.frame(proj.i$cellNames,proj.i$depth.cluster,proj.i$depth.cluster.uncertainty,proj.i$nFrags)
  df.depth <- dplyr::filter(df.depth,proj.i.depth.cluster == "2")
  df.depth <- dplyr::filter(df.depth,proj.i.depth.cluster.uncertainty <= 0.05)
  saveRDS(df.depth,paste0("df_depth_",i,".rds"))
  
  ggPoint(
    x = log10(proj.i$nFrags),
    y = log10(proj.i$TSSEnrichment+1),
    colorDensity = T,
    continuousSet = "sambaNight",
    xlabel = "log10(unique fragments)",
    ylabel = "log10(TSS Enrichment+1)"
  ) +geom_hline(yintercept = log10(min(df.TSS$proj.i.TSSEnrichment)+1),linetype = "dashed")+
    geom_vline(xintercept = min(log10(df.depth$proj.i.nFrags)),linetype = "dashed")+
    ggtitle(paste0("QC thresholds:\n",i))+
    ggsave(paste0(i,"_QC.pdf"),width = 4,height = 4)
  
  ggPoint(
    x = log10(proj.i$nFrags),
    y = log10(proj.i$TSSEnrichment+1),
    color = proj.i$DoubletEnrichment,
    discrete = F,
    continuousSet = "sambaNight",
    xlabel = "log10(unique fragments)",
    ylabel = "log10(TSS Enrichment+1)"
  ) +geom_hline(yintercept = min(log10(df.TSS$proj.i.TSSEnrichment+1)),linetype = "dashed")+
    geom_vline(xintercept = min(log10(df.depth$proj.i.nFrags)),linetype = "dashed")+
    ggtitle(paste0("Doublet Enrichment:\n",i))+
    ggsave(paste0(i,"_doublets.pdf"),width = 4,height = 4)
  
}


for (i in sampleNames[2]){
  proj.i <- proj[proj$Sample == i]
  
  # GMM for fragments per cell
  depth.clust <- Mclust(log10(proj.i$nFrags),G = 2)
  proj.i$depth.cluster <- depth.clust$classification
  proj.i$depth.cluster.uncertainty <- depth.clust$uncertainty
  
  ggPoint(
    x = log10(proj.i$nFrags),
    y = log10(proj.i$TSSEnrichment+1),
    color = as.character(proj.i$depth.cluster),
    xlabel = "log10(unique fragments)",
    ylabel = "log10(TSS Enrichment+1)"
  ) + ggtitle(paste0("GMM classification:\n",i," log10(fragments)"))+
    ggsave(paste0(i,"_depth.pdf"),width = 4,height = 4)
  
  # Manually set TSS threshold 
  #TSS.clust <- Mclust(log10(proj.i$TSSEnrichment+1),G = 2)
  proj.i$TSS.cluster <- ifelse(log10(proj.i$TSSEnrichment+1) >= 0.80,"2","1")
  proj.i$TSS.cluster.uncertainty <- rep(NA,nrow(proj.i@cellColData))

  ggPoint(
    x = log10(proj.i$nFrags),
    y = log10(proj.i$TSSEnrichment+1),
    color = as.character(proj.i$TSS.cluster),
    discrete = T,
    xlabel = "log10(unique fragments)",
    ylabel = "log10(TSS Enrichment+1)"
  ) + ggtitle(paste0("GMM classification:\n",i," TSS Enrichment"))+
    ggsave(paste0(i,"_TSS.pdf"),width = 4,height = 4)
  
  
  df.TSS <- data.frame(proj.i$cellNames,proj.i$TSS.cluster,proj.i$TSS.cluster.uncertainty,proj.i$TSSEnrichment)
  df.TSS <- dplyr::filter(df.TSS,proj.i.TSS.cluster == "2")
  #df.TSS <- dplyr::filter(df.TSS,proj.i.TSS.cluster.uncertainty <= 0.05)
  saveRDS(df.TSS,paste0("df_TSS_",i,".rds"))
  
  df.depth <- data.frame(proj.i$cellNames,proj.i$depth.cluster,proj.i$depth.cluster.uncertainty,proj.i$nFrags)
  df.depth <- dplyr::filter(df.depth,proj.i.depth.cluster == "2")
  df.depth <- dplyr::filter(df.depth,proj.i.depth.cluster.uncertainty <= 0.05)
  saveRDS(df.depth,paste0("df_depth_",i,".rds"))
  
  ggPoint(
    x = log10(proj.i$nFrags),
    y = log10(proj.i$TSSEnrichment+1),
    colorDensity = T,
    continuousSet = "sambaNight",
    xlabel = "log10(unique fragments)",
    ylabel = "log10(TSS Enrichment+1)"
  ) +geom_hline(yintercept = log10(min(df.TSS$proj.i.TSSEnrichment)+1),linetype = "dashed")+
    geom_vline(xintercept = min(log10(df.depth$proj.i.nFrags)),linetype = "dashed")+
    ggtitle(paste0("QC thresholds:\n",i))+
    ggsave(paste0(i,"_QC.pdf"),width = 4,height = 4)
  
  ggPoint(
    x = log10(proj.i$nFrags),
    y = log10(proj.i$TSSEnrichment+1),
    color = proj.i$DoubletEnrichment,
    discrete = F,
    continuousSet = "sambaNight",
    xlabel = "log10(unique fragments)",
    ylabel = "log10(TSS Enrichment+1)"
  ) +geom_hline(yintercept = min(log10(df.TSS$proj.i.TSSEnrichment+1)),linetype = "dashed")+
    geom_vline(xintercept = min(log10(df.depth$proj.i.nFrags)),linetype = "dashed")+
    ggtitle(paste0("Doublet Enrichment:\n",i))+
    ggsave(paste0(i,"_doublets.pdf"),width = 4,height = 4)
  
}


###############################################################################
dev.off()

# Filter out low quality cells, and remove doublets
##############################################################################

list.depth <- list.files(pattern = "^df_depth")

df.depth <-  data.frame(cellNames=character(),
                        cluster=character(),
                        cluster.uncertainty=character(),
                        nFrags = character())
for (i in list.depth){
  df <- readRDS(i)
  colnames(df) <- c("cellNames","cluster","cluster.uncertainty","nFrags")
  df.depth <- rbind(df.depth,df)
}

list.TSS <- list.files(pattern = "^df_TSS")

df.TSS <-  data.frame(cellNames=character(),
                      cluster=character(),
                      cluster.uncertainty=character(),
                      TSSEnrichment = character())
for (i in list.TSS){
  df <- readRDS(i)
  colnames(df) <- c("cellNames","cluster","cluster.uncertainty","TSSEnrichment")
  df.TSS <- rbind(df.TSS,df)
}


colnames(df.TSS) <- c("cellNames","TSS.cluster","TSS.cluster.uncertainty","TSSEnrichment")
colnames(df.depth) <- c("cellNames","depth.cluster","depth.cluster.uncertainty","nFrags")

cellsPass <- intersect(df.TSS$cellNames,df.depth$cellNames)

cellsFail <-  proj$cellNames[!(proj$cellNames %in% cellsPass)]

# Screen for high quality barcodes (remove non cellular barcodes)
proj.filter <- proj[proj$cellNames %in% cellsPass]


proj <- filterDoublets(proj.filter,filterRatio = 1,cutEnrich = 1,cutScore = -Inf)


plotFragmentSizes(proj)+ggtitle("Fragment Size Histogram")+ggsave("Frags_hist.pdf",width = 6,height = 4)
plotTSSEnrichment(proj)+ggtitle("TSS Enrichment")+ggsave("TSS.pdf",width = 6,height = 4)
###############################################################################################################


# Perform LSI reduction and clustering with ATAC data only
#######################################################################
# Add LSI dimreduc
proj <- addIterativeLSI(
  ArchRProj = proj,
  useMatrix = "TileMatrix",
  name = "IterativeLSI",
  iterations = 2,
  LSIMethod = 2,
  scaleDims = T,
  clusterParams = list( #See Seurat::FindClusters
    resolution = c(0.2),
    sampleCells = 10000,
    n.start = 10
  ),
  UMAPParams = list(n_neighbors =30, 
                    min_dist = 0.3, 
                    metric = "cosine", 
                    verbose =FALSE),
  varFeatures = 25000,
  dimsToUse = 1:50,
  binarize = T,
  corCutOff = 0.75,
  force = T
)

proj <- addClusters(
  input = proj,
  reducedDims = "IterativeLSI",
  method = "Seurat",
  name = "ATAC_clusters",
  resolution = 0.4,
  dimsToUse = 1:50,force = T
)

# Add UMAP based on LSI dims
proj <- addUMAP(proj,nNeighbors = 30,minDist = 0.3,dimsToUse = 1:50,metric = "cosine",force = T)

###################################################################################################


# Estimate gene activity in ATAC data and perform cell type annotation: 

# Add Gene activity matrix using ArchR model 
proj <- addGeneScoreMatrix(proj,matrixName = "ArchRGeneScore",force = T)

# Add Gene activity matrix using Signac model 
# gene.ranges <- ArchR::geneAnnoHg38
# gene.ranges <- gene.ranges$genes
# genebodyandpromoter.coords <- Signac::Extend(x = gene.ranges, upstream = 2000, downstream = 0)
# genebodyandpromoter.coords$name <- genebodyandpromoter.coords$symbol
# proj <- addFeatureMatrix(proj,matrixName = "SignacGeneScore",features = genebodyandpromoter.coords,force = T)


getAvailableMatrices(proj)
saveRDS(proj,"proj_LSI_AND_GeneScores.rds")



#  Constrained Integration to only align cells from the same patient tumor

groupList <- SimpleList()
for (i in levels(factor(proj$Sample))){
  
  rna.sub <- rna[,rna$Sample == i]
  RNA.cells <- colnames(rna.sub)

  idxSample <- BiocGenerics::which(proj$Sample == i)
  cellsSample <- proj$cellNames[idxSample]
  proj.filter <- proj[cellsSample, ]
  ATAC.cells <- proj.filter$cellNames

  groupList[[i]] <- SimpleList(
    ATAC = ATAC.cells,
    RNA = RNA.cells
  )
}



# ###########################################################################################
# # Perform Seurat v3 label transfer between RNA/ATAC
# proj <- addGeneIntegrationMatrix(
#   ArchRProj = proj,
#   useMatrix = "SignacGeneScore",
#   matrixName = "GeneIntegrationMatrix_Signac",
#   reducedDims = "IterativeLSI",
#   seRNA = rna,
#   groupList = groupList,
#   addToArrow = T,
#   force= TRUE,
#   groupRNA = "cell.type",
#   nameCell = "predictedCell_Signac",
#   nameGroup = "predictedGroup_Signac",
#   nameScore = "predictedScore_Signac",
#   plotUMAP = F,
#   useImputation = F,
#   transferParams = list(dims = 1:50)
# )

proj <- addGeneIntegrationMatrix(
  ArchRProj = proj,
  useMatrix = "ArchRGeneScore",
  matrixName = "GeneIntegrationMatrix_ArchR",
  reducedDims = "IterativeLSI",
  seRNA = rna,
  groupList = groupList,
  addToArrow = T,
  force= TRUE,
  groupRNA = "cell.type",
  nameCell = "predictedCell_ArchR",
  nameGroup = "predictedGroup_ArchR",
  nameScore = "predictedScore_ArchR",
  plotUMAP = F,
  useImputation = F,
  transferParams = list(dims = 1:50)
)
getAvailableMatrices(proj)
saveRDS(proj,"proj_LSI_GeneScores_Annotations_Int.rds")



############################################################################################
# Begin Downstream Analysis
#      1)  Plotting RNA/ATAC by sample, by cluster, by predicted label 
#      2)  Marker Gene (RNA/ATAC) intersection
#      3)  Peak2GeneLinks/Coaccessiblity
######################################################

# PART 1: Plotting 
######################################################################################

# Make embedding highlighting by 1) Predicted group ArchR 2) Predicted group Signac
# 3) Sample 4) ATAC-only clusters 
atac.archr <- plotEmbedding(proj,colorBy = "cellColData",name = "predictedGroup_ArchR")
atac.archr.emb <- as.data.frame(atac.archr$data)
atac.archr.emb$cell.type.archr <- atac.archr.emb$color
atac.archr.emb$cell.type.archr <- sub("-", ":", atac.archr.emb$cell.type.archr)
atac.archr.emb$cell.type.archr <- gsub(".*:", "", atac.archr.emb$cell.type.archr)
head(atac.archr.emb)
atac.archr.emb$cell.type.archr <- factor(atac.archr.emb$cell.type.archr, levels = levels(as.factor(rna$cell.type)))
head(atac.archr.emb)

# atac.signac <- plotEmbedding(proj,colorBy = "cellColData",name = "predictedGroup_Signac")
# atac.signac.emb <- as.data.frame(atac.signac$data)
# atac.signac.emb$cell.type.signac <- atac.signac.emb$color
# atac.signac.emb$cell.type.signac <- sub("-", ":", atac.signac.emb$cell.type.signac)
# atac.signac.emb$cell.type.signac <- gsub(".*:", "", atac.signac.emb$cell.type.signac)
# head(atac.signac.emb)
# atac.signac.emb$cell.type.signac <- factor(atac.signac.emb$cell.type.signac, levels = levels(as.factor(rna$cell.type)))
# head(atac.signac.emb)

atac <- plotEmbedding(proj,colorBy = "cellColData",name = "Sample")
atac.emb.sample <- as.data.frame(atac$data)
atac.emb.sample$sample <- atac.emb.sample$color
atac.emb.sample$sample <- sub("-", ":", atac.emb.sample$sample )
atac.emb.sample$sample  <- gsub(".*:", "", atac.emb.sample$sample )
head(atac.emb.sample)
head(atac.emb.sample)


atac <- plotEmbedding(proj,colorBy = "cellColData",name = "ATAC_clusters")
atac.emb.cluster <- as.data.frame(atac$data)
atac.emb.cluster$sample <- atac.emb.cluster$color
atac.emb.cluster$sample <- sub("-", ":", atac.emb.cluster$sample )
atac.emb.cluster$sample  <- gsub(".*:", "", atac.emb.cluster$sample )
head(atac.emb.cluster)

atac.emb.all <- cbind(atac.archr.emb[,c(1:2,4)],
                      atac.emb.sample[,4],
                      atac.emb.cluster[,4])

atac.emb.all$plain <- "Plain"

colnames(atac.emb.all) <- c("UMAP1","UMAP2","Predicted.Group.ArchR",
                            "Sample","ATAC_clusters","Blank")

head(atac.emb.all)

var.list <- colnames(atac.emb.all)[3:4]

for (i in 1:length(var.list)){

  ggplot(atac.emb.all,aes_string(x = "UMAP1",y="UMAP2",color = var.list[i]))+
    geom_point(size = .1)+
    theme_classic()+
    ggtitle(paste0("scATAC-seq: ",var.list[i]))+
    theme(plot.title = element_text(face = "bold"))+
    xlab("UMAP_1")+
    ylab("UMAP_2")+
    theme(legend.key.size = unit(0.2, "cm"))+
    guides(colour = guide_legend(override.aes = list(size=3)))+
    ggsave(paste0(var.list[i],"_ATAC.pdf"),width = 8,height = 6)
  
  
}


prediction.scores <- data.frame(ArchR= proj$predictedScore_ArchR)
var.list <- colnames(prediction.scores)

for (i in 1:length(var.list)){
  
  ggplot(prediction.scores,aes_string(x = var.list[i]))+
    geom_histogram(binwidth = 0.025,fill="#000000", color="#e9ecef", alpha=0.9)+
    theme_classic()+
    ggsave(paste0(var.list[i],"_ATAC.pdf"),width = 8,height = 6)

}



# Plot matching scRNA-seq plots:
#####################################################################################
rna.emb <- as.data.frame(rna@reductions$umap@cell.embeddings)
rna.emb$cell.type <- as.factor(rna$cell.type)
rna.emb$sample <- rna$Sample

rna.emb$cell.type <- factor(rna.emb$cell.type,levels = levels(atac.emb.all$Predicted.Group.ArchR))
rna.cell.plot <- ggplot(rna.emb,aes(x = UMAP_1,y=UMAP_2,color = cell.type))+
  geom_point(size = .1)+
  theme_classic()+
  ggtitle("scRNAseq")+
  theme(plot.title = element_text(face = "bold"))+
  xlab("UMAP_1")+
  ylab("UMAP_2")+ 
  theme(legend.key.size = unit(0.2, "cm"))+
  guides(colour = guide_legend(override.aes = list(size=3)))
rna.cell.plot +ggsave("RNA_All_labels.pdf",width = 8,height = 6)
  


rna.emb$sample <- factor(rna.emb$sample,levels = levels(factor(atac.emb.all$Sample)))
rna.sample.plot <-ggplot(rna.emb,aes(x = UMAP_1,y=UMAP_2,color = sample))+
  geom_point(size = .1)+
  theme_classic()+
  ggtitle("scRNAseq")+
  theme(plot.title = element_text(face = "bold"))+
  xlab("UMAP_1")+
  ylab("UMAP_2")+ 
  theme(legend.key.size = unit(0.2, "cm"))+
  guides(colour = guide_legend(override.aes = list(size=3)))+ggsave("RNA_All_sample.pdf")
rna.sample.plot +ggsave("RNA_All_samples.pdf",width = 8,height = 6)

#########################################################################################################


# PART 2: Marker gene ovarlap RNA/ATAC


# Differential pseudo-gene activity analysis: 
############################

######################################################################################
#ArchR
########################################################################################
# DEGs using ATAC labels
markersGS.archr <- getMarkerFeatures(
  ArchRProj = proj,
  useMatrix = "ArchRGeneScore",
  groupBy = "ATAC_clusters",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon"
)


heatmapGS.archr <- markerHeatmap(
  seMarker = markersGS.archr,
  cutOff = "FDR <= 0.01 & Log2FC >= 1.25",
  labelMarkers =NULL,
  transpose = F,
  pal =  viridis(n=256),
  limits = c(-2,2)
)
ComplexHeatmap::draw(heatmapGS.archr, heatmap_legend_side = "bot", annotation_legend_side = "bot")
plotPDF(heatmapGS.archr , name = "GeneScores-Marker-Heatmap_ArchR", width = 8, height = 6, ArchRProj = proj, addDOC = FALSE)




# DEGs using predicted labels (removing small groups)
idxSample <- BiocGenerics::which(proj$predictedScore_ArchR > 0.5)
cellsSample <- proj$cellNames[idxSample]
proj.filter <- proj[cellsSample, ]

popular.groups <- summary(factor(proj.filter$predictedGroup_ArchR))
popular.groups <- popular.groups[popular.groups > 10]
proj.filter$Mode.Label <- ifelse(proj.filter$predictedGroup_ArchR %in% names(popular.groups),TRUE,FALSE)

idxSample <- BiocGenerics::which(proj.filter$Mode.Label == TRUE)
cellsSample <- proj.filter$cellNames[idxSample]
proj.filter <- proj.filter[cellsSample, ]

# DEGs using predicted labels
markersGS.archr.pred <- getMarkerFeatures(
  ArchRProj = proj.filter,
  useMatrix = "ArchRGeneScore",
  groupBy = "predictedGroup_ArchR",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon"
)
heatmapGS.archr.pred <- markerHeatmap(
  seMarker = markersGS.archr.pred,
  cutOff = "FDR <= 0.01 & Log2FC >= 1.25",
  labelMarkers =NULL,
  transpose = F,
  pal =  viridis(n=256),
  limits = c(-2,2)
)
ComplexHeatmap::draw(heatmapGS.archr.pred, heatmap_legend_side = "bot", annotation_legend_side = "bot")
plotPDF(heatmapGS.archr.pred, name = "GeneScores-Marker-Heatmap_ArchR_pred", width = 8, height = 6, ArchRProj = proj.filter, addDOC = FALSE)

########################################################################################################################
# 
# 
# # Signac
# ########################################################################################################################
# # DEGs using ATAC labels
# markersGS.signac <- getMarkerFeatures(
#   ArchRProj = proj,
#   useMatrix = "SignacGeneScore",
#   groupBy = "ATAC_clusters",
#   bias = c("TSSEnrichment", "log10(nFrags)"),
#   testMethod = "wilcoxon"
# )
# 
# 
# heatmapGS.signac <- markerHeatmap(
#   seMarker = markersGS.signac,
#   cutOff = "FDR <= 0.01 & Log2FC >= 1.25",
#   labelMarkers =NULL,
#   transpose = F,
#   pal =  viridis(n=256),
#   limits = c(-2,2)
# )
# ComplexHeatmap::draw(heatmapGS.signac, heatmap_legend_side = "bot", annotation_legend_side = "bot")
# plotPDF(heatmapGS.signac , name = "GeneScores-Marker-Heatmap_signac", width = 8, height = 6, ArchRProj = proj, addDOC = FALSE)
# 
# 
# 
# 
# # DEGs using predicted labels (removing small groups)
# idxSample <- BiocGenerics::which(proj$predictedScore_Signac > 0.5)
# cellsSample <- proj$cellNames[idxSample]
# proj.filter <- proj[cellsSample, ]
# 
# popular.groups <- summary(factor(proj.filter$predictedGroup_Signac))
# popular.groups <- popular.groups[popular.groups > 10]
# proj.filter$Mode.Label <- ifelse(proj.filter$predictedGroup_Signac %in% names(popular.groups),TRUE,FALSE)
# 
# idxSample <- BiocGenerics::which(proj.filter$Mode.Label == TRUE)
# cellsSample <- proj.filter$cellNames[idxSample]
# proj.filter <- proj.filter[cellsSample, ]
# 
# # DEGs using predicted labels
# markersGS.signac.pred <- getMarkerFeatures(
#   ArchRProj = proj.filter,
#   useMatrix = "SignacGeneScore",
#   groupBy = "predictedGroup_Signac",
#   bias = c("TSSEnrichment", "log10(nFrags)"),
#   testMethod = "wilcoxon"
# )
# heatmapGS.signac.pred <- markerHeatmap(
#   seMarker = markersGS.signac.pred,
#   cutOff = "FDR <= 0.01 & Log2FC >= 1.25",
#   labelMarkers =NULL,
#   transpose = F,
#   pal =  viridis(n=256),
#   limits = c(-2,2)
# )
# ComplexHeatmap::draw(heatmapGS.signac.pred, heatmap_legend_side = "bot", annotation_legend_side = "bot")
# plotPDF(heatmapGS.signac.pred, name = "GeneScores-Marker-Heatmap_Signac_pred", width = 8, height = 6, ArchRProj = proj.filter, addDOC = FALSE)

# Differential peak analysis:
############################
# ATAC clusters
proj <- addGroupCoverages(ArchRProj = proj,groupBy = "ATAC_clusters",force = T)

pathToMacs2 <- findMacs2()

proj <- addReproduciblePeakSet(
  ArchRProj = proj,
  groupBy = "ATAC_clusters",
  pathToMacs2 = pathToMacs2,force = T
)
proj <- addPeakMatrix(proj,force = T)
proj <- addBgdPeaks(proj)

markersPeaks <- getMarkerFeatures(
  ArchRProj = proj,
  useMatrix = "PeakMatrix",
  groupBy = "ATAC_clusters",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon"
)

heatmapPeaks<- markerHeatmap(
  seMarker = markersPeaks,
  cutOff = "FDR <= 0.01 & Log2FC >= 1.25",
  labelMarkers =NULL,
  transpose = F,
  pal =  viridis(n=256),
  limits = c(-2,2)
)
ComplexHeatmap::draw(heatmapPeaks, heatmap_legend_side = "bot", annotation_legend_side = "bot")
plotPDF(heatmapPeaks, name = "Markers_peaks_ATAC_clusters", width = 8, height = 6, ArchRProj = proj, addDOC = FALSE)


# ArchR predicted labels


# DEGs using predicted labels (removing small groups)
idxSample <- BiocGenerics::which(proj$predictedScore_ArchR >= 0.5)
cellsSample <- proj$cellNames[idxSample]
proj.filter <- proj[cellsSample, ]

popular.groups <- summary(factor(proj.filter$predictedGroup_ArchR))
popular.groups <- popular.groups[popular.groups > 10]
proj.filter$Mode.Label <- ifelse(proj.filter$predictedGroup_ArchR %in% names(popular.groups),TRUE,FALSE)

idxSample <- BiocGenerics::which(proj.filter$Mode.Label == TRUE)
cellsSample <- proj.filter$cellNames[idxSample]
proj.filter <- proj.filter[cellsSample, ]

proj.archr <- addGroupCoverages(ArchRProj = proj.filter,groupBy = "predictedGroup_ArchR",force = T)

pathToMacs2 <- findMacs2()

proj.archr <- addReproduciblePeakSet(
  ArchRProj = proj.archr,
  groupBy = "predictedGroup_ArchR",
  pathToMacs2 = pathToMacs2,force = T
)
proj.archr <- addPeakMatrix(proj.archr,force = T)
proj.archr <- addBgdPeaks(proj.archr)
markersPeaks.archr <- getMarkerFeatures(
  ArchRProj = proj.archr,
  useMatrix = "PeakMatrix",
  groupBy = "predictedGroup_ArchR",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon"
)

heatmapPeaks.archr<- markerHeatmap(
  seMarker = markersPeaks.archr,
  cutOff = "FDR <= 0.01 & Log2FC >= 1.25",
  labelMarkers =NULL,
  transpose = F,
  pal =  viridis(n=256),
  limits = c(-2,2)
)
ComplexHeatmap::draw(heatmapPeaks.archr, heatmap_legend_side = "bot", annotation_legend_side = "bot")
plotPDF(heatmapPeaks.archr , name = "Markers_peaks_Archr_Predicted_labels", width = 8, height = 6, ArchRProj =  proj.archr, addDOC = FALSE)

# 
# 
# # Signac predicted labels
# # DEGs using predicted labels (removing small groups)
# idxSample <- BiocGenerics::which(proj$predictedScore_Signac > 0.5)
# cellsSample <- proj$cellNames[idxSample]
# proj.filter <- proj[cellsSample, ]
# 
# popular.groups <- summary(factor(proj.filter$predictedGroup_Signac))
# popular.groups <- popular.groups[popular.groups > 10]
# proj.filter$Mode.Label <- ifelse(proj.filter$predictedGroup_Signac %in% names(popular.groups),TRUE,FALSE)
# 
# idxSample <- BiocGenerics::which(proj.filter$Mode.Label == TRUE)
# cellsSample <- proj.filter$cellNames[idxSample]
# proj.filter <- proj.filter[cellsSample, ]
# 
# proj.signac <- addGroupCoverages(ArchRProj = proj.filter,groupBy = "predictedGroup_Signac",force = T)
# 
# pathToMacs2 <- findMacs2()
# 
# proj.signac <- addReproduciblePeakSet(
#   ArchRProj = proj.signac,
#   groupBy = "predictedGroup_Signac",
#   pathToMacs2 = pathToMacs2,force = T
# )
# proj.signac <- addPeakMatrix( proj.signac,force = T)
# proj.signac <- addBgdPeaks(proj.signac)
# markersPeaks.signac <- getMarkerFeatures(
#   ArchRProj = proj.signac,
#   useMatrix = "PeakMatrix",
#   groupBy = "predictedGroup_Signac",
#   bias = c("TSSEnrichment", "log10(nFrags)"),
#   testMethod = "wilcoxon"
# )
# 
# heatmapPeaks.signac <- markerHeatmap(
#   seMarker = markersPeaks.signac,
#   cutOff = "FDR <= 0.01 & Log2FC >= 1.25",
#   labelMarkers =NULL,
#   transpose = F,
#   pal =  viridis(n=256),
#   limits = c(-2,2)
# )
# ComplexHeatmap::draw(heatmapPeaks.signac, heatmap_legend_side = "bot", annotation_legend_side = "bot")
# plotPDF(heatmapPeaks.signac, name = "Markers_peaks_Signac_Predicted_labels", width = 8, height = 6, ArchRProj =  proj.signac, addDOC = FALSE)
# 


# RNA heatmap
####################################################################################
topN <-Wilcox.markers %>% group_by(cluster) %>% dplyr::filter(p_val_adj <= 0.01) %>% top_n(30, desc(avg_logFC))

# Seurat Heatmap (top N DEGs per cluster)
#rna.sub <- subset(rna,downsample = 300)
#DoHeatmap(rna.sub,features = top30$gene,draw.lines = F,size = 2,disp.min = -2,disp.max = 2)+scale_fill_viridis()

## All 2000 variable features
#DoHeatmap(rna.sub,features = VariableFeatures(rna),draw.lines = F,size = 2,disp.min = -2,disp.max = 2)+scale_fill_viridis()


# Downsample cells from each cluster
rna.sub <- subset(rna,downsample =300)
rna.sub <- NormalizeData(rna.sub)
rna.sub <- rna.sub[rownames(rna.sub) %in% topN$gene,]
rna.sub <- ScaleData(rna.sub,features = rownames(rna.sub))

mat <- rna.sub@assays$RNA@scale.data

cluster_anno<- rna.sub@meta.data$cell.type

# gg <- ggplot_build(rna.cell.plot)
#
# data <- as.data.frame(gg$data)
# names <- data$group
# names <- names -1
# data$group <- names
# #data$group <- data$colour
#
#
# cols <- unique(data$colour)
# names <- as.numeric(unique(data$group))
# names(cols) <- names
#
# cols <- cols[order(factor(names(cols), levels=levels(cluster_anno)))]

col_fun = circlize::colorRamp2(c(-2, 0, 2),viridis(n = 3))

heatmapRNA <- Heatmap(mat, name = "Expression",
        column_split = factor(cluster_anno),
        cluster_columns =T,
        show_column_dend = F,
        cluster_column_slices = T,
        column_title_gp = gpar(fontsize = 8),
        column_gap = unit(0.1, "mm"),
        cluster_rows = T,
        show_row_dend = FALSE,
        col = col_fun,
        column_title_rot = 90,
        show_column_names = F)
plotPDF(heatmapRNA, name = "Heatmap_RNA", width = 8, height = 6)



##################################################################################################

source("./Archr_Peak_Null.R")
# Overlap RNA/ATAC DEG hits 1) ArchR 2) Signac
markerList.atac <- getMarkers(markersGS.archr.pred, cutOff = "FDR <= 0.01 & Log2FC >= 1.25")

Idents(rna) <- "cell.type"
Wilcox.markers <- readRDS("./wilcox_DEGs.rds")
# Make marker RNA list
for ( i in levels(factor(rna$cell.type))){
  rna.cluster <- Wilcox.markers[Wilcox.markers$cluster == i &
                                  Wilcox.markers$p_val_adj <= 0.01 &
                                  Wilcox.markers$avg_logFC >= 1.0,]
  saveRDS(rna.cluster,paste0("cluster_markers_",i,".rds"))
}

my.files <- list.files(pattern = "cluster_markers")
my.files


# now importing them all at once as entries in a list
# for my files, I have to provide some options to read.csv
markerList.rna <- lapply(my.files,
               readRDS)

names <- str_remove(my.files,"cluster_markers_")
names(markerList.rna) <- str_remove(names,".rds")


int.names <- intersect(names(markerList.rna),names(markerList.atac))


gene.hits <- list()
for (i in int.names){

  atac.clust <- markerList.atac[[i]]
  rna.clust <- markerList.rna[[i]]
  cluster.hits <- intersect(atac.clust$name,rna.clust$gene)

  gene.hits[[i]] <- cluster.hits
}
########################################################################################

# Add Coaccessiblity and Peak2Gene links:
proj.archr <- addCoAccessibility(
  ArchRProj = proj.archr,
  reducedDims = "IterativeLSI",
  dimsToUse = 1:50,
  scaleDims = NULL,
  corCutOff = 0.75,
  k = 100,
  knnIteration = 500,
  overlapCutoff = 0.8,
  maxDist = 1e+05,
  scaleTo = 10^4,
  log2Norm = TRUE,
  seed = 1,
  threads = getArchRThreads(),
  verbose = TRUE,
  logFile = createLogFile("addCoAccessibility")
)

proj.archr <- addPeak2GeneLinks(
  ArchRProj = proj.archr ,
  reducedDims = "IterativeLSI",
  useMatrix = "GeneIntegrationMatrix_ArchR",
  dimsToUse = 1:50,
  scaleDims = NULL,
  corCutOff = 0.75,
  k = 100,
  knnIteration = 500,
  overlapCutoff = 0.8,
  maxDist = 250000,
  scaleTo = 10^4,
  log2Norm = TRUE,
  predictionCutoff = 0.5,
  seed = 1,
  threads = max(floor(getArchRThreads()/2), 1),
  verbose = TRUE,
  logFile = createLogFile("addPeak2GeneLinks")
)


# For every cluster, visualize the peaks coaccessible with marker genes

gg <- ggplot_build(rna.cell.plot)

data <- as.data.frame(gg$data)
data$cell.type <- rna.emb$cell.type


cols <- unique(data$colour)
names(cols) <- unique(data$cell.type)

for (i in names(gene.hits)){
  genes <- gene.hits[[i]]

  p <- plotBrowserTrack(
    ArchRProj = proj.archr,
    groupBy = "predictedGroup_ArchR",
    geneSymbol = genes,
    upstream = 50000,
    downstream = 50000,
    loops = getCoAccessibility(
      ArchRProj = proj.archr,
      corCutOff = 0.70,
      resolution = 1,
      returnLoops = TRUE
    ),
    features = makeGRangesFromDataFrame(encode.all),
    pal = cols
  )
  plotPDF(plotList = p,
          name = paste0("ArchR_CoAscblty_Hits",i,".pdf"),
          ArchRProj = proj.archr,
          addDOC = T, width = 5, height = 5)


  p <- plotBrowserTrack(
    ArchRProj = proj.archr,
    groupBy = "predictedGroup_ArchR",
    geneSymbol = genes,
    upstream = 50000,
    downstream = 50000,
    loops = getPeak2GeneLinks(
      ArchRProj = proj.archr,
      corCutOff = 0.70,
      FDRCutOff = 1e-04,
      resolution = 1,
      returnLoops = TRUE
    ),
    features = makeGRangesFromDataFrame(encode.all),
    pal = cols

  )
  loops <- getPeak2GeneLinks(
    ArchRProj = proj.archr,
    corCutOff = 0.70,
    FDRCutOff = 1e-04,
    resolution = 1,
    returnLoops = F
  )
  saveRDS(loops,paste0(i,"_empirical_pval_archr.rds"))

  plotPDF(plotList = p,
          name = paste0("ArchR_Peak2Gene_Hits_",i,".pdf"),
          ArchRProj = proj.archr,
          addDOC = T, width = 6, height = 8)

}

# source("./Archr_Peak_Null.R")
# # Overlap RNA/ATAC DEG hits 
# markerList.atac <- getMarkers(markersGS.signac.pred, cutOff = "FDR <= 0.01 & Log2FC >= 1.25")
# 
# # Make marker RNA list
# for ( i in levels(factor(rna$cell.type))){
#   rna.cluster <- Wilcox.markers[Wilcox.markers$cluster == i &
#                                   Wilcox.markers$p_val_adj <= 0.01 &
#                                   Wilcox.markers$avg_logFC >= 1.0,]
#   saveRDS(rna.cluster,paste0("cluster_markers_signac_",i,".rds"))
# }
# 
# my.files <- list.files(pattern = "cluster_markers_signac_")
# my.files
# 
# 
# # now importing them all at once as entries in a list
# # for my files, I have to provide some options to read.csv
# markerList.rna <- lapply(my.files,
#                          readRDS)
# 
# names <- str_remove(my.files,"cluster_markers_signac_")
# names(markerList.rna) <- str_remove(names,".rds")
# 
# 
# int.names <- intersect(names(markerList.rna),names(markerList.atac))
# 
# 
# gene.hits <- list()
# for (i in int.names){
# 
#   atac.clust <- markerList.atac[[i]]
#   rna.clust <- markerList.rna[[i]]
#   cluster.hits <- intersect(atac.clust$name,rna.clust$gene)
# 
#   gene.hits[[i]] <- cluster.hits
# }
# ########################################################################################
# 
# # Add Coaccessiblity and Peak2Gene links:
# proj.signac <- addCoAccessibility(
#   ArchRProj = proj.signac,
#   reducedDims = "IterativeLSI",
#   dimsToUse = 1:50,
#   scaleDims = NULL,
#   corCutOff = 0.75,
#   k = 100,
#   knnIteration = 500,
#   overlapCutoff = 0.8,
#   maxDist = 1e+05,
#   scaleTo = 10^4,
#   log2Norm = TRUE,
#   seed = 1,
#   threads = getArchRThreads(),
#   verbose = TRUE,
#   logFile = createLogFile("addCoAccessibility")
# )
# 
# proj.signac <- addPeak2GeneLinks(
#   ArchRProj = proj.signac ,
#   reducedDims = "IterativeLSI",
#   useMatrix = "GeneIntegrationMatrix_Signac",
#   dimsToUse = 1:50,
#   scaleDims = NULL,
#   corCutOff = 0.75,
#   k = 100,
#   knnIteration = 500,
#   overlapCutoff = 0.8,
#   maxDist = 250000,
#   scaleTo = 10^4,
#   log2Norm = TRUE,
#   predictionCutoff = 0.5,
#   seed = 1,
#   threads = max(floor(getArchRThreads()/2), 1),
#   verbose = TRUE,
#   logFile = createLogFile("addPeak2GeneLinks")
# )
# 
# 
# # For every cluster, visualize the peaks coaccessible with marker genes
# 
# 
# gg <- ggplot_build(rna.cell.plot)
# 
# data <- as.data.frame(gg$data)
# data$cell.type <- rna.emb$cell.type
# 
# 
# cols <- unique(data$colour)
# names(cols) <- unique(data$cell.type)
# 
# 
# for (i in names(gene.hits)){
#   genes <- gene.hits[[i]]
# 
#   p <- plotBrowserTrack(
#     ArchRProj = proj.signac,
#     groupBy = "predictedGroup_Signac",
#     geneSymbol = genes,
#     upstream = 50000,
#     downstream = 50000,
#     loops = getCoAccessibility(
#       ArchRProj = proj.signac,
#       corCutOff = 0.70,
#       resolution = 1,
#       returnLoops = TRUE
#     ),
#     features = makeGRangesFromDataFrame(encode.all),
#     pal = cols
#   )
#   plotPDF(plotList = p,
#           name = paste0("signac_CoAscblty_Hits",i,".pdf"),
#           ArchRProj = proj.signac,
#           addDOC = T, width = 5, height = 5)
# 
# 
#   p <- plotBrowserTrack(
#     ArchRProj = proj.signac,
#     groupBy = "predictedGroup_Signac",
#     geneSymbol = genes,
#     upstream = 50000,
#     downstream = 50000,
#     loops = getPeak2GeneLinks(
#       ArchRProj = proj.signac,
#       corCutOff = 0.70,
#       FDRCutOff = 1e-04,
#       resolution = 1,
#       returnLoops = TRUE
#     ),
#     features = makeGRangesFromDataFrame(encode.all),
#     pal = cols
# 
#   )
# 
#   loops <- getPeak2GeneLinks(
#     ArchRProj = proj.signac,
#     corCutOff = 0.70,
#     FDRCutOff = 1e-04,
#     resolution = 1,
#     returnLoops = F
#   )
#   saveRDS(loops,paste0(i,"_empirical_pval_signac.rds"))
# 
# 
#   plotPDF(plotList = p,
#           name = paste0("signac_Peak2Gene_Hits",i,".pdf"),
#           ArchRProj = proj.signac,
#           addDOC = T, width = 6, height = 8)
# 
# }


saveRDS(proj.archr,"final_archr_proj_archrGS.rds")
#saveRDS(proj.signac,"final_archr_proj_signacGS.rds")


################
# Archr Version
################
dir.create("./BAM_Subset_ArchR")
setwd("./BAM_Subset_ArchR")

for (i in levels(factor(proj.archr$Sample))){
  barcode.df.archr <- dplyr::filter(as.data.frame(proj.archr@cellColData),Sample == i)
  
  print(paste0("Starting Sample:",i))
  
  for (k in levels(factor(barcode.df.archr$predictedGroup_ArchR))){
    barcode.df.archr.sub <- dplyr::filter(barcode.df.archr,predictedGroup_ArchR == k)
    csv <- data.frame(barcode=rep(NA,1000000))
    
    print(paste0("Starting Sample:",k))
    for (j in 1:nrow(barcode.df.archr.sub)){
      barcode.df.archr.sub.sub <- barcode.df.archr.sub[j,]
      
      if (barcode.df.archr.sub.sub$predictedGroup_ArchR %in% labels){
        csv$barcode[j] <-rownames(barcode.df.archr.sub.sub)
      }else{
        print("Cell not of interest")
      }
      
    }
    
    csv <- na.omit(csv)
    csv$barcode <- gsub(".*#","",csv$barcode)
    write.table(csv$barcode,paste0(i,"_",k,"_barcodes.csv"),row.names = F,col.names = FALSE,sep = ",",quote = F)
  }
  
}
ff <- dir(getwd(), recursive=TRUE, full.names=TRUE)
## Extract vector of empty files' names
eff <- ff[file.info(ff)[["size"]]==0]
## Remove empty files
unlink(eff, recursive=TRUE, force=FALSE)


# Note that GIST 3E4D1L causes error in for loop, still outputs necessary files
for (i in levels(factor(proj.archr$Sample))){
  my.files <- list.files(pattern = "_barcodes.csv")
  idx <- grep(i,my.files)
  my.files <- my.files[idx]
  
  csv_config <- data.frame(library_id = str_remove(my.files,"_barcodes.csv"),
                           barcodes_csv = paste0(getwd(),"/",my.files))
  #csv_config$library_id <- sub("\\ .*", "", my.files)
  write.csv(csv_config,paste0(SAMPLE.ID,"_",i,"_config.csv"),row.names = F,quote = F)
}



peak2genes <- getPeak2GeneLinks(proj.archr,returnLoops = F,corCutOff = 0,FDRCutOff = 1)
p2g <- as.data.frame(peak2genes)

coords <- as.data.frame(metadata(peak2genes)[[1]])
coords$idxATAC <- 1:nrow(coords)

coords$coord <- paste0(coords$seqnames,":",coords$start,"-",coords$end)

peaks.atac <- readRDS(paste0("../",SAMPLE.ID,"/Peak2GeneLinks/seATAC-Group-KNN.rds"))
peaks <- as.data.frame(rowRanges(peaks.atac))

peaks <- peaks[,c(1:3,13)]
peaks$coord <- paste0(peaks$seqnames,":",peaks$start,"-",peaks$end)


head(coords)
head(peaks)
dim(coords)
dim(peaks)


peaks.coords <- merge(coords,peaks,by ="coord")


p2g.coords <- merge(peaks.coords,p2g,by="idxATAC")


p2g.coords <- dplyr::filter(p2g.coords,peakType == "Distal")
p2g.coords <- dplyr::arrange(p2g.coords,EmpFDR)


saveRDS(p2g.coords,"Peak2GeneLinks_EmpFDR_Ranked_Distal_ArchR.rds")


# Write ALL distal peaks to BED file:

#p2g.coords.sub <- dplyr::filter(p2g.coords,EmpFDR <= 0.05)
p2g.coords.sub <- dplyr::arrange(p2g.coords,EmpFDR)
bed <- p2g.coords.sub[,3:5]
write.table(bed,"Peak2GeneLinks_EmpFDR_0.05_Ranked_Distal_ArchR.bed",row.names = F,col.names = F,quote = F,sep = "\t")



##################################################################

# Write marker distal peaks to BED file(s):
markerList <- getMarkers(markersPeaks.archr, cutOff = "FDR <= 0.10 & Log2FC >= 1.25")

names <- intersect(names(markerList),labels)
markerList <- as.list(markerList[names])

marker.index <- numeric()
for (i in 1:length(markerList)){
  df <- as.data.frame(markerList[[i]])
  if(nrow(df) > 0 ){
    marker.index <- append(marker.index,i)
  }else{
    print("No marker peaks are significant")
  }
}
markerList <- markerList[marker.index]# 

# Subset marker list to cell types of interest
for (i in 1:length(markerList)){
  df <- as.data.frame(markerList[[i]])
  df$coord <- paste0(df$seqnames,":",df$start,"-",df$end)
  markerList[[i]] <- df
}

# Write marker enhancers to BED files 
for (i in names(markerList)){
  marker.peaks <- as.data.frame(markerList[[i]])
  peaks <- as.data.frame(intersect(p2g.coords.sub$coord,marker.peaks$coord))
  colnames(peaks) <- "coord"
  
  peaks.signif <- merge(peaks,marker.peaks,by = "coord")
  peaks.signif <- dplyr::arrange(peaks.signif,FDR)
  
  peaks <- as.data.frame(peaks.signif[,1])
  colnames(peaks) <- "Peak"
  
  peaks <- as.data.frame(str_split_fixed(peaks$Peak,":",2))
  colnames(peaks) <- c("chr","range")
  peaks.old <- peaks
  peaks.new <- as.data.frame(str_split_fixed(peaks$range,"-",2))
  
  peaks <- cbind(peaks.old,peaks.new)
  
  colnames(peaks) <- c("chr","range","start","end")
  
  peaks <- peaks[,-2]
  
  
  bed <- peaks
  if (nrow(bed) < 1){
    print("No marker peaks overlapped with peak2gene links")
  }else{
    write.table(bed,paste0("Marker_Enhancers_ArchR_",i,".bed"),row.names = F,col.names = F,quote = F,sep = "\t")
  }
  
}


# ################
# # Signac Version
# ################
# dir.create("./BAM_Subset_Signac")
# setwd("./BAM_Subset_Signac")
# 
# for (i in levels(factor(proj.signac$Sample))){
#   barcode.df.signac <- dplyr::filter(as.data.frame(proj.signac@cellColData),Sample == i)
#   
#   print(paste0("Starting Sample:",i))
#   
#   for (k in levels(factor(barcode.df.signac$predictedGroup_Signac))){
#     barcode.df.signac.sub <- dplyr::filter(barcode.df.signac,predictedGroup_Signac == k)
#     csv <- data.frame(barcode=rep(NA,1000000))
#     
#     print(paste0("Starting Sample:",k))
#     for (j in 1:nrow(barcode.df.signac.sub)){
#       barcode.df.signac.sub.sub <- barcode.df.signac.sub[j,]
#       
#       if (barcode.df.signac.sub.sub$predictedGroup_Signac %in% labels){
#         csv$barcode[j] <-rownames(barcode.df.signac.sub.sub)
#       }else{
#         print("Cell not of interest")
#       }
#       
#     }
#     
#     csv <- na.omit(csv)
#     csv$barcode <- gsub(".*#","",csv$barcode)
#     write.table(csv$barcode,paste0(i,"_",k,"_barcodes.csv"),row.names = F,col.names = FALSE,sep = ",",quote = F)
#   }
#   
# }
# ff <- dir(getwd(), recursive=TRUE, full.names=TRUE)
# ## Extract vector of empty files' names
# eff <- ff[file.info(ff)[["size"]]==0]
# ## Remove empty files
# unlink(eff, recursive=TRUE, force=FALSE)
# 
# 
# # Note that GIST 3E4D1L causes error in for loop, still outputs necessary files
# for (i in levels(factor(proj.signac$Sample))){
#   my.files <- list.files(pattern = "_barcodes.csv")
#   idx <- grep(i,my.files)
#   my.files <- my.files[idx]
#   
#   csv_config <- data.frame(library_id = str_remove(my.files,"_barcodes.csv"),
#                            barcodes_csv = paste0(getwd(),"/",my.files))
#   #csv_config$library_id <- sub("\\ .*", "", my.files)
#   write.csv(csv_config,paste0(SAMPLE.ID,"_",i,"_config.csv"),row.names = F,quote = F)
# }
# 
# 
# 
# peak2genes <- getPeak2GeneLinks(proj.signac,returnLoops = F,corCutOff = 0,FDRCutOff = 1)
# p2g <- as.data.frame(peak2genes)
# 
# coords <- as.data.frame(metadata(peak2genes)[[1]])
# coords$idxATAC <- 1:nrow(coords)
# 
# coords$coord <- paste0(coords$seqnames,":",coords$start,"-",coords$end)
# 
# peaks.atac <- readRDS(paste0("../",SAMPLE.ID,"/Peak2GeneLinks/seATAC-Group-KNN.rds"))
# peaks <- as.data.frame(rowRanges(peaks.atac))
# 
# peaks <- peaks[,c(1:3,13)]
# peaks$coord <- paste0(peaks$seqnames,":",peaks$start,"-",peaks$end)
# 
# 
# head(coords)
# head(peaks)
# dim(coords)
# dim(peaks)
# 
# 
# peaks.coords <- merge(coords,peaks,by ="coord")
# 
# 
# p2g.coords <- merge(peaks.coords,p2g,by="idxATAC")
# 
# 
# p2g.coords <- dplyr::filter(p2g.coords,peakType == "Distal")
# p2g.coords <- dplyr::arrange(p2g.coords,EmpFDR)
# 
# 
# saveRDS(p2g.coords,"Peak2GeneLinks_EmpFDR_Ranked_Distal_Signac.rds")
# 
# 
# # Write ALL distal peaks to BED file:
# 
# p2g.coords.sub <- dplyr::filter(p2g.coords,EmpFDR <= 0.05)
# p2g.coords.sub <- dplyr::arrange(p2g.coords.sub,EmpFDR)
# bed <- p2g.coords.sub[,3:5]
# write.table(bed,"Peak2GeneLinks_EmpFDR_0.05_Ranked_Distal_Signac.bed",row.names = F,col.names = F,quote = F,sep = "\t")
# 
# 
# 
# ##################################################################
# 
# # Write marker distal peaks to BED file(s):
# 
# markerList <- getMarkers(markersPeaks.signac, cutOff = "FDR <= 0.01 & Log2FC >= 1.25")
# 
# names <- intersect(names(markerList),labels)
# markerList <- as.list(markerList[labels])
# 
# 
# # Subset marker list to cell types of interest
# for (i in 1:length(markerList)){
#   df <- as.data.frame(markerList[[i]])
#   df$coord <- paste0(df$seqnames,":",df$start,"-",df$end)
#   markerList[[i]] <- df
# }
# 
# # Write marker enhancers to BED files 
# for (i in labels){
#   marker.peaks <- as.data.frame(markerList[[i]])
#   peaks <- as.data.frame(intersect(p2g.coords.sub$coord,marker.peaks$coord))
#   colnames(peaks) <- "Peak"
#   peaks <- as.data.frame(str_split_fixed(peaks$Peak,":",2))
#   colnames(peaks) <- c("chr","range")
#   peaks.old <- peaks
#   peaks.new <- as.data.frame(str_split_fixed(peaks$range,"-",2))
#   
#   peaks <- cbind(peaks.old,peaks.new)
#   
#   colnames(peaks) <- c("chr","range","start","end")
#   
#   peaks <- peaks[,-2]
#   
#   
#   bed <- peaks
#   if (nrow(bed) < 1){
#     print("No marker peaks overlapped with peak2gene links")
#   }else{
#     write.table(bed,paste0("Marker_Enhancers_Signac_",i,".bed"),row.names = F,col.names = F,quote = F,sep = "\t")
#   }
#   
# }
# 
# 
# 


######################################################################################################
##########################################################################
# END OF SCRIPT
##########################################################################

