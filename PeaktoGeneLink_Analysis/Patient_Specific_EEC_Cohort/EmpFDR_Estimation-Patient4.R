
library(ArchR)
library(ChIPpeakAnno)
library(stats)
ArchR::addArchRThreads(threads = 64)
source("./Archr_Peak_Null_Permute.R")
source("./Archr_Peak_RawPval.R")
source("./P2G_Heatmap_Distal.R")
dir.create("./Significant_P2G_Outputs_Patient4")
SAMPLE.ID <- "EEC"
key.word.1 <- "pithelia"
key.word.2 <- "-Ciliated"


# Permuated p2g links 
#################################################################################################################
proj <- readRDS("./final_archr_proj_archrGS.rds")
idxSample <- BiocGenerics::which(proj$Sample =="36639L")
cellsSample <- proj$cellNames[idxSample]
proj <- proj[cellsSample, ]

store <- as.numeric(0)
for (i in 1:100){
  proj.null <- addPermPeak2GeneLinks(
    ArchRProj = proj ,
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
  p2geneDF <- metadata(proj.null@peakSet)$Peak2GeneLinks
  p2geneDF$geneName <- mcols(metadata(p2geneDF)$geneSet)$name[p2geneDF$idxRNA]
  p2geneDF$peakName <- (metadata(p2geneDF)$peakSet %>% {paste0(seqnames(.), "_", start(.), "_", end(.))})[p2geneDF$idxATAC]
  p2g.df.null <- as.data.frame(p2geneDF)
  p2g.df.null <- p2g.df.null[complete.cases(p2g.df.null),]
  p2g.null.sub <- dplyr::filter(p2g.df.null,RawPVal <=1e-12)
  store[i] <- nrow(p2g.null.sub)
}
# Compute eFDR for alpha 1e-12 my taking the median of the store vector dividing by # of significant observed tests
saveRDS(store,"./Significant_P2G_Outputs_Patient4/store_null_tests.rds")


# Observed data final run:
proj <- readRDS("./final_archr_proj_archrGS.rds")
idxSample <- BiocGenerics::which(proj$Sample =="36639L")
cellsSample <- proj$cellNames[idxSample]
proj <- proj[cellsSample, ]
# Add p2g links (no restrictions on FDR, Correlation, Variance cutoff) with raw pvalue
##############################################################################################################
proj <- addPeak2GeneLinks(
  ArchRProj = proj ,
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

p2geneDF <- metadata(proj@peakSet)$Peak2GeneLinks
p2geneDF$geneName <- mcols(metadata(p2geneDF)$geneSet)$name[p2geneDF$idxRNA]
p2geneDF$peakName <- (metadata(p2geneDF)$peakSet %>% {paste0(seqnames(.), ":", start(.), "-", end(.))})[p2geneDF$idxATAC]
annot <- readRDS(metadata(p2geneDF)$seATAC)
p2geneDF$peakType <- annot@rowRanges$peakType[p2geneDF$idxATAC]
p2g.df.obs <- as.data.frame(p2geneDF)
p2g.df.obs <- p2g.df.obs[complete.cases(p2g.df.obs),]

pdf("./Significant_P2G_Outputs_Patient4/P2G_Correlation.pdf",width = 5,height = 3.5)
hist(p2g.df.obs$Correlation,col = "lightblue",
     main = paste0("Histogram of\n",nrow(p2g.df.obs)," P2G correlations"),xlab = "Correlation")
dev.off()

pdf("./Significant_P2G_Outputs_Patient4/P2g_Pval.pdf",width = 5,height = 3.5)
hist(p2g.df.obs$RawPVal,col="lightblue",
     main = paste0("Histogram of\n",nrow(p2g.df.obs)," P2G p-values"),xlab = "p-value")
abline(v=0.01,col = "red")
dev.off()

saveRDS(p2g.df.obs,"./Significant_P2G_Outputs_Patient4/All_P2G_Observed.rds")

################################################################################################################

p2g.df.hist <- dplyr::filter(p2g.df.obs,RawPVal <=1e-12 & Correlation >= 0.45)

pdf("./Significant_P2G_Outputs_Patient4/GenesPerPeak_histogram.pdf",width = 5,height = 3.5)
hist(table(p2g.df.hist$idxRNA),main="Distribution of genes per peaks")
dev.off()

pdf("./Significant_P2G_Outputs_Patient4/PeaksPerGene_histogram.pdf",width = 5,height = 3.5)
hist(table(p2g.df.hist$idxATAC),main="Distribution of peaks per gene")
dev.off()


#Subset to postive correlation P2Gs
p2g.df.sub <- dplyr::filter(p2g.df.obs,RawPVal <=1e-12 & Correlation >= 0.45& peakType == "Distal")

# Plot peak2 gene heatmap

p2g.df.sub$idx <- paste0(p2g.df.sub$idxATAC,"-",p2g.df.sub$idxRNA)
p2g.df.sub.plot <- p2g.df.sub

# Make color scheme for heatmap based on original UMAP colors:
proj$cluster.num <- factor(gsub("-.*","",proj$predictedGroup_ArchR))
my_levels <- c("8",
               "21",
               "2",
               "16",
               "5",
               "3","26",
               "0","18",
               "7",
               "25",
               "27")
# Relevel object@ident
proj@cellColData$cluster.new <- factor(x=proj$cluster.num, levels = my_levels)
# Make order of colors:
epithelial.cols <- colorRampPalette(c("#A0E989", "#337B24"))
epithelial.cols <- epithelial.cols(8)
fibro.cols <- c("#f593f1","#f272ed","#ef52e9","#e415dd","#c412bd","#a30f9d","#62095e")
smooth.cols <- c("#c2abd6","#b47fe5","#8948a1")
endo.cols <- c("#93CEFF","#4A99FF","#286ab5")
t.cols <- c("gray40","gray60")
macro.cols <- c("#ff6600","#ff9d5c")
mast.cols <- "gold3"
b.cols <- c("#B22222","#CD5C5C")
cols <- c(epithelial.cols,fibro.cols,smooth.cols,endo.cols,t.cols,macro.cols,mast.cols,b.cols)
cols <- cols[c(2,7,10,13,17,19,21,22,23,25,26,27)]
names(cols) <- levels(factor(proj@cellColData$cluster.new))

source("P2G_Heatmap_Distal.R")
test <- plotPeak2GeneHeatmap.distal(proj,peaks=p2g.df.sub.plot,groupBy = "cluster.new",k=length(levels(factor(proj$cluster.new))),
                                    corCutOff = .45,varCutOffATAC = 0,varCutOffRNA = 0,FDRCutOff = 1,nPlot =100000,palGroup=cols,
                                    palATAC = paletteContinuous("solarExtra"),
                                    palRNA = paletteContinuous("solarExtra"))

pdf("./Significant_P2G_Outputs_Patient4/Peak2Gene_Heatmap_Legend.pdf",width = 8,height = 12)
draw(test, heatmap_legend_side = "bottom")
dev.off()

p2g.heat.df <- plotPeak2GeneHeatmap.distal(proj,peaks=p2g.df.sub.plot,groupBy = "predictedGroup_ArchR",k=length(levels(factor(proj$predictedGroup_ArchR))),
                                           corCutOff = .45,varCutOffATAC = 0,varCutOffRNA = 0,FDRCutOff = 1,returnMatrices = T,nPlot = 100000)


# Save P2G peaknames and Kmeans cluster for the genes of interest (specific to cancer cells)

p2g.df.sub.plot$kmeans <- p2g.heat.df$RNA$kmeansId
saveRDS(p2g.df.sub.plot,"./Significant_P2G_Outputs_Patient4/Distal_P2Gs_Kmeans.rds")

pdf("./Significant_P2G_Outputs_Patient4/GenesPerPeak_histogram-distal.pdf",width = 5,height = 3.5)
hist(table(p2g.df.sub.plot$idxRNA),main="Distribution of genes per distal peaks")
dev.off()

pdf("./Significant_P2G_Outputs_Patient4/PeaksPerGene_histogram-distal.pdf",width = 5,height = 3.5)
hist(table(p2g.df.sub.plot$idxATAC),main="Distribution of distal peaks per gene")
dev.off()


