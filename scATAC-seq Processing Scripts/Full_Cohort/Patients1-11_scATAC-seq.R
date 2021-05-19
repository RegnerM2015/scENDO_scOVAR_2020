###########################################################
# Matt Regner
# Franco Lab
# Date: May-December 2020
# 
# Sample: All 
# Description: This script performs the following tasks  
#         1) scATAC-seq processing
#         2) scRNA-seq/scATAC-seq integration
#         3) Peak2Gene linkage analysis/Co-accessiblity 
#         4) Make input files for TFSEE analysis 
###########################################################
.libPaths('/home/regnerm/anaconda3/envs/scENDO_scOVAR/lib/R/library')


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

set.seed(123)
addArchRThreads(threads = 24) 
addArchRGenome("hg38")

# Set up directories and file variables:
##################################################################################################
SAMPLE.ID <- "All"

output.dir <- "."
####################
setwd(output.dir)
####################

inputFiles <- list.files(pattern = "\\.gz$")

sampleNames <- c("3533EL","3571DL","36186L","36639L","366C5L","37EACL","38FE7L","3BAE2L","3CCF1L","3E4D1L","3E5CFL")


encode.all <- read.delim("./GRCh38-ccREs.bed",header =F)
colnames(encode.all)[1:3] <- c("seqnames","start","end")


# Store patient metadata and colors:

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

# 
# # Read in matching scRNA-seq
# ###################################################################################################
rna <- readRDS("./endo_ovar_All_scRNA_processed.rds")
rna$cell.type <- str_replace(rna$cell.type,"/","_")
Idents(rna) <- "cell.type"

# Redo differential expression with new cell type markers 
Wilcox.markers <- readRDS("./wilcox_DEGs.rds")
Wilcox.markers$cluster <- str_replace(Wilcox.markers$cluster,"/","_")

# 
# Create Arrow and ArchR project
# ##########################################################################
# ArrowFiles <- createArrowFiles(
#  inputFiles = inputFiles,
#  sampleNames = sampleNames,
#  filterTSS = 0, #Dont set this too high because you can always increase later
#  filterFrags = 0,
#  addTileMat = T,
#  addGeneScoreMat = F
# )
# ArrowFiles <- list.files(pattern=".arrow")
# doubScores <- addDoubletScores(
#   input = ArrowFiles,
#   k = 10, #Refers to how many cells near a "pseudo-doublet" to count.
#   knnMethod = "UMAP",useMatrix = "TileMatrix",nTrials=5,LSIMethod = 1,scaleDims = F,
#   corCutOff = 0.75,UMAPParams = list(n_neighbors =30, min_dist = 0.3, metric = "cosine", verbose =FALSE),
#   dimsToUse = 1:50
# )
# 
# proj <- ArchRProject(
#   ArrowFiles = ArrowFiles,
#   outputDirectory = "All",
#   copyArrows = T #This is recommened so that you maintain an unaltered copy for later usage.
# )
# 
# # Filter out outlier low quality cells and doublets
# ###############################################################################
# # GMM for fragments per cell
# library(mclust)
# 
# for (i in sampleNames){
#   proj.i <- proj[proj$Sample == i]
# 
#   # GMM for fragments per cell
#   depth.clust <- Mclust(log10(proj.i$nFrags),G = 2)
#   proj.i$depth.cluster <- depth.clust$classification
#   proj.i$depth.cluster.uncertainty <- depth.clust$uncertainty
# 
#   ggPoint(
#     x = log10(proj.i$nFrags),
#     y = log10(proj.i$TSSEnrichment+1),
#     color = as.character(proj.i$depth.cluster),
#     xlabel = "log10(unique fragments)",
#     ylabel = "log10(TSS Enrichment+1)"
#   ) + ggtitle(paste0("GMM classification:\n",i," log10(fragments)"))+
#     ggsave(paste0(i,"_depth.pdf"),width = 4,height = 4)
# 
#   # GMM for TSS per cell
#   TSS.clust <- Mclust(log10(proj.i$TSSEnrichment+1),G = 2)
#   proj.i$TSS.cluster <- TSS.clust$classification
#   proj.i$TSS.cluster.uncertainty <- TSS.clust$uncertainty
# 
#   ggPoint(
#     x = log10(proj.i$nFrags),
#     y = log10(proj.i$TSSEnrichment+1),
#     color = as.character(proj.i$TSS.cluster),
#     discrete = T,
#     xlabel = "log10(unique fragments)",
#     ylabel = "log10(TSS Enrichment+1)"
#   ) + ggtitle(paste0("GMM classification:\n",i," TSS Enrichment"))+
#     ggsave(paste0(i,"_TSS.pdf"),width = 4,height = 4)
# 
# 
#   df.TSS <- data.frame(proj.i$cellNames,proj.i$TSS.cluster,proj.i$TSS.cluster.uncertainty,proj.i$TSSEnrichment)
#   df.TSS <- dplyr::filter(df.TSS,proj.i.TSS.cluster == "2")
#   df.TSS <- dplyr::filter(df.TSS,proj.i.TSS.cluster.uncertainty <= 0.05)
#   saveRDS(df.TSS,paste0("df_TSS_",i,".rds"))
# 
#   df.depth <- data.frame(proj.i$cellNames,proj.i$depth.cluster,proj.i$depth.cluster.uncertainty,proj.i$nFrags)
#   df.depth <- dplyr::filter(df.depth,proj.i.depth.cluster == "2")
#   df.depth <- dplyr::filter(df.depth,proj.i.depth.cluster.uncertainty <= 0.05)
#   saveRDS(df.depth,paste0("df_depth_",i,".rds"))
# 
#   ggPoint(
#     x = log10(proj.i$nFrags),
#     y = log10(proj.i$TSSEnrichment+1),
#     colorDensity = T,
#     continuousSet = "sambaNight",
#     xlabel = "log10(unique fragments)",
#     ylabel = "log10(TSS Enrichment+1)"
#   ) +geom_hline(yintercept = log10(min(df.TSS$proj.i.TSSEnrichment)+1),linetype = "dashed")+
#     geom_vline(xintercept = min(log10(df.depth$proj.i.nFrags)),linetype = "dashed")+
#     ggtitle(paste0("QC thresholds:\n",i))+
#     ggsave(paste0(i,"_QC.pdf"),width = 4,height = 4)
# 
#   ggPoint(
#     x = log10(proj.i$nFrags),
#     y = log10(proj.i$TSSEnrichment+1),
#     color = proj.i$DoubletEnrichment,
#     discrete = F,
#     continuousSet = "sambaNight",
#     xlabel = "log10(unique fragments)",
#     ylabel = "log10(TSS Enrichment+1)"
#   ) +geom_hline(yintercept = min(log10(df.TSS$proj.i.TSSEnrichment+1)),linetype = "dashed")+
#     geom_vline(xintercept = min(log10(df.depth$proj.i.nFrags)),linetype = "dashed")+
#     ggtitle(paste0("Doublet Enrichment:\n",i))+
#     ggsave(paste0(i,"_doublets.pdf"),width = 4,height = 4)
# 
# }
# 
# 
# for (i in sampleNames[2]){
#   proj.i <- proj[proj$Sample == i]
# 
#   # GMM for fragments per cell
#   depth.clust <- Mclust(log10(proj.i$nFrags),G = 2)
#   proj.i$depth.cluster <- depth.clust$classification
#   proj.i$depth.cluster.uncertainty <- depth.clust$uncertainty
# 
#   ggPoint(
#     x = log10(proj.i$nFrags),
#     y = log10(proj.i$TSSEnrichment+1),
#     color = as.character(proj.i$depth.cluster),
#     xlabel = "log10(unique fragments)",
#     ylabel = "log10(TSS Enrichment+1)"
#   ) + ggtitle(paste0("GMM classification:\n",i," log10(fragments)"))+
#     ggsave(paste0(i,"_depth.pdf"),width = 4,height = 4)
# 
#   # Manually set TSS threshold
#   #TSS.clust <- Mclust(log10(proj.i$TSSEnrichment+1),G = 2)
#   proj.i$TSS.cluster <- ifelse(log10(proj.i$TSSEnrichment+1) >= 0.80,"2","1")
#   proj.i$TSS.cluster.uncertainty <- rep(NA,nrow(proj.i@cellColData))
# 
#   ggPoint(
#     x = log10(proj.i$nFrags),
#     y = log10(proj.i$TSSEnrichment+1),
#     color = as.character(proj.i$TSS.cluster),
#     discrete = T,
#     xlabel = "log10(unique fragments)",
#     ylabel = "log10(TSS Enrichment+1)"
#   ) + ggtitle(paste0("GMM classification:\n",i," TSS Enrichment"))+
#     ggsave(paste0(i,"_TSS.pdf"),width = 4,height = 4)
# 
# 
#   df.TSS <- data.frame(proj.i$cellNames,proj.i$TSS.cluster,proj.i$TSS.cluster.uncertainty,proj.i$TSSEnrichment)
#   df.TSS <- dplyr::filter(df.TSS,proj.i.TSS.cluster == "2")
#   #df.TSS <- dplyr::filter(df.TSS,proj.i.TSS.cluster.uncertainty <= 0.05)
#   saveRDS(df.TSS,paste0("df_TSS_",i,".rds"))
# 
#   df.depth <- data.frame(proj.i$cellNames,proj.i$depth.cluster,proj.i$depth.cluster.uncertainty,proj.i$nFrags)
#   df.depth <- dplyr::filter(df.depth,proj.i.depth.cluster == "2")
#   df.depth <- dplyr::filter(df.depth,proj.i.depth.cluster.uncertainty <= 0.05)
#   saveRDS(df.depth,paste0("df_depth_",i,".rds"))
# 
#   ggPoint(
#     x = log10(proj.i$nFrags),
#     y = log10(proj.i$TSSEnrichment+1),
#     colorDensity = T,
#     continuousSet = "sambaNight",
#     xlabel = "log10(unique fragments)",
#     ylabel = "log10(TSS Enrichment+1)"
#   ) +geom_hline(yintercept = log10(min(df.TSS$proj.i.TSSEnrichment)+1),linetype = "dashed")+
#     geom_vline(xintercept = min(log10(df.depth$proj.i.nFrags)),linetype = "dashed")+
#     ggtitle(paste0("QC thresholds:\n",i))+
#     ggsave(paste0(i,"_QC.pdf"),width = 4,height = 4)
# 
#   ggPoint(
#     x = log10(proj.i$nFrags),
#     y = log10(proj.i$TSSEnrichment+1),
#     color = proj.i$DoubletEnrichment,
#     discrete = F,
#     continuousSet = "sambaNight",
#     xlabel = "log10(unique fragments)",
#     ylabel = "log10(TSS Enrichment+1)"
#   ) +geom_hline(yintercept = min(log10(df.TSS$proj.i.TSSEnrichment+1)),linetype = "dashed")+
#     geom_vline(xintercept = min(log10(df.depth$proj.i.nFrags)),linetype = "dashed")+
#     ggtitle(paste0("Doublet Enrichment:\n",i))+
#     ggsave(paste0(i,"_doublets.pdf"),width = 4,height = 4)
# 
# }
# 
# for (i in sampleNames[7]){
#   proj.i <- proj[proj$Sample == i]
# 
#   # GMM for fragments per cell
#   depth.clust <- Mclust(log10(proj.i$nFrags),G = 2)
#   proj.i$depth.cluster <- depth.clust$classification
#   proj.i$depth.cluster.uncertainty <- depth.clust$uncertainty
# 
#   ggPoint(
#     x = log10(proj.i$nFrags),
#     y = log10(proj.i$TSSEnrichment+1),
#     color = as.character(proj.i$depth.cluster),
#     xlabel = "log10(unique fragments)",
#     ylabel = "log10(TSS Enrichment+1)"
#   ) + ggtitle(paste0("GMM classification:\n",i," log10(fragments)"))+
#     ggsave(paste0(i,"_depth.pdf"),width = 4,height = 4)
# 
#   # Manually set TSS threshold
#   #TSS.clust <- Mclust(log10(proj.i$TSSEnrichment+1),G = 2)
#   proj.i$TSS.cluster <- ifelse(log10(proj.i$TSSEnrichment+1) >= 0.80,"2","1")
#   proj.i$TSS.cluster.uncertainty <- rep(NA,nrow(proj.i@cellColData))
# 
#   ggPoint(
#     x = log10(proj.i$nFrags),
#     y = log10(proj.i$TSSEnrichment+1),
#     color = as.character(proj.i$TSS.cluster),
#     discrete = T,
#     xlabel = "log10(unique fragments)",
#     ylabel = "log10(TSS Enrichment+1)"
#   ) + ggtitle(paste0("GMM classification:\n",i," TSS Enrichment"))+
#     ggsave(paste0(i,"_TSS.pdf"),width = 4,height = 4)
# 
# 
#   df.TSS <- data.frame(proj.i$cellNames,proj.i$TSS.cluster,proj.i$TSS.cluster.uncertainty,proj.i$TSSEnrichment)
#   df.TSS <- dplyr::filter(df.TSS,proj.i.TSS.cluster == "2")
#   #df.TSS <- dplyr::filter(df.TSS,proj.i.TSS.cluster.uncertainty <= 0.05)
#   saveRDS(df.TSS,paste0("df_TSS_",i,".rds"))
# 
#   df.depth <- data.frame(proj.i$cellNames,proj.i$depth.cluster,proj.i$depth.cluster.uncertainty,proj.i$nFrags)
#   df.depth <- dplyr::filter(df.depth,proj.i.depth.cluster == "2")
#   df.depth <- dplyr::filter(df.depth,proj.i.depth.cluster.uncertainty <= 0.05)
#   saveRDS(df.depth,paste0("df_depth_",i,".rds"))
# 
#   ggPoint(
#     x = log10(proj.i$nFrags),
#     y = log10(proj.i$TSSEnrichment+1),
#     colorDensity = T,
#     continuousSet = "sambaNight",
#     xlabel = "log10(unique fragments)",
#     ylabel = "log10(TSS Enrichment+1)"
#   ) +geom_hline(yintercept = log10(min(df.TSS$proj.i.TSSEnrichment)+1),linetype = "dashed")+
#     geom_vline(xintercept = min(log10(df.depth$proj.i.nFrags)),linetype = "dashed")+
#     ggtitle(paste0("QC thresholds:\n",i))+
#     ggsave(paste0(i,"_QC.pdf"),width = 4,height = 4)
# 
#   ggPoint(
#     x = log10(proj.i$nFrags),
#     y = log10(proj.i$TSSEnrichment+1),
#     color = proj.i$DoubletEnrichment,
#     discrete = F,
#     continuousSet = "sambaNight",
#     xlabel = "log10(unique fragments)",
#     ylabel = "log10(TSS Enrichment+1)"
#   ) +geom_hline(yintercept = min(log10(df.TSS$proj.i.TSSEnrichment+1)),linetype = "dashed")+
#     geom_vline(xintercept = min(log10(df.depth$proj.i.nFrags)),linetype = "dashed")+
#     ggtitle(paste0("Doublet Enrichment:\n",i))+
#     ggsave(paste0(i,"_doublets.pdf"),width = 4,height = 4)
# 
# }
# 
# 
# 
# ###############################################################################
# dev.off()
# 
# # Filter out low quality cells, and remove doublets
# ##############################################################################
# 
# list.depth <- list.files(pattern = "^df_depth")
# 
# df.depth <-  data.frame(cellNames=character(),
#                         cluster=character(),
#                         cluster.uncertainty=character(),
#                         nFrags = character())
# for (i in list.depth){
#   df <- readRDS(i)
#   colnames(df) <- c("cellNames","cluster","cluster.uncertainty","nFrags")
#   df.depth <- rbind(df.depth,df)
# }
# 
# list.TSS <- list.files(pattern = "^df_TSS")
# 
# df.TSS <-  data.frame(cellNames=character(),
#                       cluster=character(),
#                       cluster.uncertainty=character(),
#                       TSSEnrichment = character())
# for (i in list.TSS){
#   df <- readRDS(i)
#   colnames(df) <- c("cellNames","cluster","cluster.uncertainty","TSSEnrichment")
#   df.TSS <- rbind(df.TSS,df)
# }
# 
# 
# colnames(df.TSS) <- c("cellNames","TSS.cluster","TSS.cluster.uncertainty","TSSEnrichment")
# colnames(df.depth) <- c("cellNames","depth.cluster","depth.cluster.uncertainty","nFrags")
# 
# cellsPass <- intersect(df.TSS$cellNames,df.depth$cellNames)
# 
# cellsFail <-  proj$cellNames[!(proj$cellNames %in% cellsPass)]
# 
# # Screen for high quality barcodes (remove non cellular barcodes)
# proj.filter <- proj[proj$cellNames %in% cellsPass]
# 
# 
# proj <- filterDoublets(proj.filter,filterRatio = 1,cutEnrich = 1,cutScore = -Inf)
# 
# 
# plotFragmentSizes(proj)+ggtitle("Fragment Size Histogram")+ggsave("Frags_hist.pdf",width = 6,height = 4)
# plotTSSEnrichment(proj)+ggtitle("TSS Enrichment")+ggsave("TSS.pdf",width = 6,height = 4)
# ###############################################################################################################

proj <- readRDS("./proj_LSI_AND_UMAP.rds")
# Perform LSI reduction and clustering with ATAC data only
#######################################################################

# Add LSI dimreduc
proj <- addIterativeLSI(
  ArchRProj = proj,
  useMatrix = "TileMatrix",
  name = "IterativeLSI",
  iterations = 4,
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
  force = T,
  seed=6
)

proj <- addClusters(
  input = proj,
  reducedDims = "IterativeLSI",
  method = "Seurat",
  name = "ATAC_clusters",
  resolution = 0.7,
  dimsToUse = 1:50,force = T
)

# Add UMAP based on LSI dims
proj <- addUMAP(proj,nNeighbors = 30,minDist = 0.3,dimsToUse = 1:50,metric = "cosine",force = T,reducedDims="IterativeLSI")
saveRDS(proj,"proj_LSI_AND_UMAP.rds")

###################################################################################################

# Estimate gene activity in ATAC data and perform cell type annotation:

# Add Gene activity matrix using ArchR model
proj <- addGeneScoreMatrix(proj,matrixName = "ArchRGeneScore",force = T)


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


ggplot(atac.emb.all,aes_string(x = "UMAP1",y="UMAP2",color = "Predicted.Group.ArchR"))+
  geom_point(size = .1)+
  theme_classic()+
  ggtitle("scATAC-seq: Predicted Group ArchR")+
  theme(plot.title = element_text(face = "bold"))+
  xlab("UMAP_1")+
  ylab("UMAP_2")+
  theme(legend.key.size = unit(0.2, "cm"))+
  guides(colour = guide_legend(override.aes = list(size=3)))+
  ggsave(paste0("PredictedGroup_ArchR_ATAC.pdf"),width = 8,height = 6)

ggplot(atac.emb.all,aes_string(x = "UMAP1",y="UMAP2",color = "Sample"))+
  geom_point(size = .1)+
  theme_classic()+
  ggtitle(paste0("scATAC-seq: Sample"))+
  theme(plot.title = element_text(face = "bold"))+
  xlab("UMAP_1")+
  ylab("UMAP_2")+
  theme(legend.key.size = unit(0.2, "cm"))+
  scale_color_manual(values = sampleColors)
  #guides(colour = guide_legend(override.aes = list(size=3)))+
  ggsave(paste0("Sample_ATAC.pdf"),width = 8,height = 6)





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
  scale_color_manual(values = sampleColors)+
  guides(colour = guide_legend(override.aes = list(size=3)))+ggsave("RNA_All_sample.pdf")
rna.sample.plot +ggsave("RNA_All_samples.pdf",width = 8,height = 6)

#########################################################################################################


# # PART 2: Marker gene ovarlap RNA/ATAC
#
#
# # Differential pseudo-gene activity analysis:
# ############################
#
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

#
#
#
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
proj.archr <- addBgdPeaks(proj.archr,force = T)
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

saveRDS(proj.archr,"./final_archr_proj_archrGS.rds")
# # RNA heatmap
####################################################################################
topN <-Wilcox.markers %>% group_by(cluster) %>% dplyr::filter(p_val_adj <= 0.01) %>% top_n(30, desc(avg_logFC))


# Downsample cells from each cluster
rna.sub <- subset(rna,downsample =300)
rna.sub <- NormalizeData(rna.sub)
rna.sub <- rna.sub[rownames(rna.sub) %in% topN$gene,]
rna.sub <- ScaleData(rna.sub,features = rownames(rna.sub))

mat <- rna.sub@assays$RNA@scale.data

cluster_anno<- rna.sub@meta.data$cell.type

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



##############
# END OF SCRIPT
##############
