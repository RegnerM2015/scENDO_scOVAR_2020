
library(ArchR)
library(ChIPpeakAnno)
library(stats)
ArchR::addArchRThreads(threads = 64)
source("./Archr_Peak_Null_Permute.R")
source("./Archr_Peak_RawPval.R")
dir.create("./Significant_P2G_Outputs")
SAMPLE.ID <- "HGSOC"
key.word.1 <- "pithelia"


# Add p2g links (no restrictions on FDR, Correlation, Variance cutoff) with raw pvalue
##############################################################################################################
proj <- readRDS("./final_archr_proj_archrGS.rds")
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
first.quart <- summary(p2g.df.obs$RawPVal)[2]


# Permuated p2g links 
#################################################################################################################
proj <- readRDS("./final_archr_proj_archrGS.rds")
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
#Example histograms
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

store.prop <- numeric(0)
test <- readRDS(paste0("./",SAMPLE.ID,"/Peak2GeneLinks/seATAC-Group-KNN.rds"))
test <- test@metadata$KNNList@listData
for ( i in 1:length(test)){
  
  test[[i]] <- gsub("\\#.*","",test[[i]])
  num <- max(table(test[[i]]))
  store.prop[i] <- num/100
}
saveRDS(store.prop,"./Significant_P2G_Outputs/store_null_KNN_proportions.rds")

pdf("./Significant_P2G_Outputs/PatientPurityPerAgg-null.pdf",width = 5,height = 3.5)
hist(store.prop,main="Distribtion of patient purity per cell aggregate")
dev.off()


p2geneDF <- metadata(proj.null@peakSet)$Peak2GeneLinks
p2geneDF$geneName <- mcols(metadata(p2geneDF)$geneSet)$name[p2geneDF$idxRNA]
p2geneDF$peakName <- (metadata(p2geneDF)$peakSet %>% {paste0(seqnames(.), "_", start(.), "_", end(.))})[p2geneDF$idxATAC]
p2g.df.null <- as.data.frame(p2geneDF)
p2g.df.null <- p2g.df.null[complete.cases(p2g.df.null),]
saveRDS(p2g.df.null,"./Significant_P2G_Outputs/All_P2G_Null_example.rds")
pdf("./Significant_P2G_Outputs/P2G_Correlation-null.pdf",width = 5,height = 3.5)
hist(p2g.df.null$Correlation,col = "lightblue",main = "Histogram of null peak-to-gene correlations",xlab = "Correlation")
dev.off()

pdf("./Significant_P2G_Outputs/P2g_Pval-null.pdf",width = 5,height = 3.5)
hist(p2g.df.null$RawPVal,col="lightblue",main = "Histogram null peak-to-gene p-values",xlab = "p-value")
abline(v=0.01,col = "red")
dev.off()

# Compute eFDR for alpha 1e-12
print(median(store)/nrow(p2g.df.obs[p2g.df.obs$RawPVal <=1e-12,]))
saveRDS(store,"./Significant_P2G_Outputs/store_null_tests.rds")





# Repeat null correlation tests for first quartile alpha threshold
# Permuated p2g links 
#################################################################################################################
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
  p2g.null.sub <- dplyr::filter(p2g.df.null,RawPVal <=first.quart)
  store[i] <- nrow(p2g.null.sub)
}
#Example histograms
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

store.prop <- numeric(0)
test <- readRDS(paste0("./",SAMPLE.ID,"/Peak2GeneLinks/seATAC-Group-KNN.rds"))
test <- test@metadata$KNNList@listData
for ( i in 1:length(test)){
  
  test[[i]] <- gsub("\\#.*","",test[[i]])
  num <- max(table(test[[i]]))
  store.prop[i] <- num/100
}
saveRDS(store.prop,"./Significant_P2G_Outputs/store_null_KNN_proportions-firstquart.rds")

pdf("./Significant_P2G_Outputs/PatientPurityPerAgg-null-firstquart.pdf",width = 5,height = 3.5)
hist(store.prop,main="Distribtion of patient purity per cell aggregate")
dev.off()


p2geneDF <- metadata(proj.null@peakSet)$Peak2GeneLinks
p2geneDF$geneName <- mcols(metadata(p2geneDF)$geneSet)$name[p2geneDF$idxRNA]
p2geneDF$peakName <- (metadata(p2geneDF)$peakSet %>% {paste0(seqnames(.), "_", start(.), "_", end(.))})[p2geneDF$idxATAC]
p2g.df.null <- as.data.frame(p2geneDF)
p2g.df.null <- p2g.df.null[complete.cases(p2g.df.null),]
saveRDS(p2g.df.null,"./Significant_P2G_Outputs/All_P2G_Null_example-firstquart.rds")
pdf("./Significant_P2G_Outputs/P2G_Correlation-null-firstquart.pdf",width = 5,height = 3.5)
hist(p2g.df.null$Correlation,col = "lightblue",main = "Histogram of null peak-to-gene correlations",xlab = "Correlation")
dev.off()

pdf("./Significant_P2G_Outputs/P2g_Pval-null-firstquart.pdf",width = 5,height = 3.5)
hist(p2g.df.null$RawPVal,col="lightblue",main = "Histogram null peak-to-gene p-values",xlab = "p-value")
abline(v=0.01,col = "red")
dev.off()

# Compute eFDR for alpha first.quart
print(median(store)/nrow(p2g.df.obs[p2g.df.obs$RawPVal <=first.quart,]))
saveRDS(store,"./Significant_P2G_Outputs/store_null_tests-firstquart.rds")




# Observed data final run:
proj <- readRDS("./final_archr_proj_archrGS.rds")
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
saveRDS(proj,"./final_archr_proj_archrGS-P2Gs.rds")
store.prop <- numeric(0)
test <- readRDS(paste0("./",SAMPLE.ID,"/Peak2GeneLinks/seATAC-Group-KNN.rds"))
test <- test@metadata$KNNList@listData
for ( i in 1:length(test)){
  
  test[[i]] <- gsub("\\#.*","",test[[i]])
  num <- max(table(test[[i]]))
  store.prop[i] <- num/100
}
saveRDS(store.prop,"./Significant_P2G_Outputs/store_KNN_proportions.rds")

pdf("./Significant_P2G_Outputs/PatientPurityPerAgg.pdf",width = 5,height = 3.5)
hist(store.prop,main="Distribtion of patient purity per cell aggregate")
dev.off()


p2geneDF <- metadata(proj@peakSet)$Peak2GeneLinks
p2geneDF$geneName <- mcols(metadata(p2geneDF)$geneSet)$name[p2geneDF$idxRNA]
p2geneDF$peakName <- (metadata(p2geneDF)$peakSet %>% {paste0(seqnames(.), ":", start(.), "-", end(.))})[p2geneDF$idxATAC]
annot <- readRDS(metadata(p2geneDF)$seATAC)
p2geneDF$peakType <- annot@rowRanges$peakType[p2geneDF$idxATAC]
p2g.df.obs <- as.data.frame(p2geneDF)
p2g.df.obs <- p2g.df.obs[complete.cases(p2g.df.obs),]

pdf("./Significant_P2G_Outputs/P2G_Correlation.pdf",width = 5,height = 3.5)
hist(p2g.df.obs$Correlation,col = "lightblue",
     main = paste0("Histogram of\n",nrow(p2g.df.obs)," P2G correlations"),xlab = "Correlation")
dev.off()

pdf("./Significant_P2G_Outputs/P2g_Pval.pdf",width = 5,height = 3.5)
hist(p2g.df.obs$RawPVal,col="lightblue",
     main = paste0("Histogram of\n",nrow(p2g.df.obs)," P2G p-values"),xlab = "p-value")
abline(v=0.01,col = "red")
dev.off()

saveRDS(p2g.df.obs,"./Significant_P2G_Outputs/All_P2G_Observed.rds")

################################################################################################################

p2g.df.hist <- dplyr::filter(p2g.df.obs,RawPVal <=1e-12 & Correlation >= 0.45)

pdf("./Significant_P2G_Outputs/GenesPerPeak_histogram.pdf",width = 5,height = 3.5)
hist(table(p2g.df.hist$idxRNA),main="Distribution of genes per peaks")
dev.off()

pdf("./Significant_P2G_Outputs/PeaksPerGene_histogram.pdf",width = 5,height = 3.5)
hist(table(p2g.df.hist$idxATAC),main="Distribution of peaks per gene")
dev.off()


#Subset to postive correlation P2Gs
p2g.df.sub <- dplyr::filter(p2g.df.obs,RawPVal <=1e-12 & Correlation >= 0.45& peakType == "Distal")

# Plot peak2 gene heatmap

p2g.df.sub$idx <- paste0(p2g.df.sub$idxATAC,"-",p2g.df.sub$idxRNA)
p2g.df.sub.plot <- p2g.df.sub

# Make color scheme for heatmap based on original UMAP colors:
proj$cluster.num <- factor(gsub("-.*","",proj$predictedGroup_ArchR))
my_levels <- c("0","2","3","7","11","16",
               "1","4","9","10",
               "14",
               "5","12","18","21",
               "6","8","13",
               "17")
# Relevel object@ident
proj@cellColData$cluster.new <- factor(x=proj$cluster.num, levels = my_levels)
# Make order of colors:
epithelial.cols <- colorRampPalette(c("#A0E989", "#337B24"))
epithelial.cols <- epithelial.cols(6)
fibro.cols <- c("#FB8CAB","#E65C9C","#CF268A","#Af1281")
endo.cols <- "#377EB8"
t.cols <- c("#CCCCCC","#333333","#666666","#999999")
macro.cols <- colorRampPalette(c("#FECC51", "#FD6A02"))
macro.cols <- macro.cols(3)
b.cols <- "#E41A1C"
cols <- c(epithelial.cols,fibro.cols,endo.cols,t.cols,macro.cols,b.cols)
names(cols) <- levels(factor(proj@cellColData$cluster.new))

source("P2G_Heatmap_Distal.R")
test <- plotPeak2GeneHeatmap.distal(proj,peaks=p2g.df.sub.plot,groupBy = "cluster.new",k=length(levels(factor(proj$cluster.new))),
                                    corCutOff = .45,varCutOffATAC = 0,varCutOffRNA = 0,FDRCutOff = 1,nPlot =100000,palGroup=cols,
                                    palATAC = paletteContinuous("solarExtra"),
                                    palRNA = paletteContinuous("solarExtra"))

pdf("./Significant_P2G_Outputs/Peak2Gene_Heatmap_Legend.pdf",width = 8,height = 12)
draw(test, heatmap_legend_side = "bottom")
dev.off()

p2g.heat.df <- plotPeak2GeneHeatmap.distal(proj,peaks=p2g.df.sub.plot,groupBy = "predictedGroup_ArchR",k=length(levels(factor(proj$predictedGroup_ArchR))),
                                           corCutOff = .45,varCutOffATAC = 0,varCutOffRNA = 0,FDRCutOff = 1,returnMatrices = T,nPlot = 100000)


# Save P2G peaknames and Kmeans cluster for the genes of interest (specific to cancer cells)

p2g.df.sub.plot$kmeans <- p2g.heat.df$RNA$kmeansId
saveRDS(p2g.df.sub.plot,"./Significant_P2G_Outputs/Distal_P2Gs_Kmeans.rds")

pdf("./Significant_P2G_Outputs/GenesPerPeak_histogram-distal.pdf",width = 5,height = 3.5)
hist(table(p2g.df.sub.plot$idxRNA),main="Distribution of genes per distal peaks")
dev.off()

pdf("./Significant_P2G_Outputs/PeaksPerGene_histogram-distal.pdf",width = 5,height = 3.5)
hist(table(p2g.df.sub.plot$idxATAC),main="Distribution of distal peaks per gene")
dev.off()

store.kmeans <- c("4","5","6","7","8","9","10","11","18","19")
# Extract P2Gs for kmeans clusters of interest
p2g.df.sub.plot.cancer.kmeans <- p2g.df.sub.plot[p2g.df.sub.plot$kmeans %in% store.kmeans,]
saveRDS(p2g.df.sub.plot.cancer.kmeans,"./Significant_P2G_Outputs/Cancer_enriched_P2G_table.rds")
p2g <- GRanges(p2g.df.sub.plot.cancer.kmeans$peakName)

encode.all <- read.delim("./GRCh38-ccREs.bed",header =F)
colnames(encode.all)[1:3] <- c("seqnames","start","end")
widths <- encode.all$end - encode.all$start
encode.all <- GRanges(encode.all)

ft.peaks <- readRDS("Fallopian_Tube_Cell_line_Peaks.rds")
widths.2 <- ft.peaks$end - ft.peaks$start
ft.peaks <- GRanges(ft.peaks)

oe.peaks <- readRDS("Ovarian_Epithelial_Cell_line_Peaks.rds")
widths.3 <- oe.peaks$end - oe.peaks$start
oe.peaks <- GRanges(oe.peaks)

widths <- c(widths,widths.2,widths.3)
tot <- 3.3e+9*0.98/mean(widths)
##
ol <- findOverlapsOfPeaks(unique(p2g), encode.all,unique(oe.peaks),unique(ft.peaks), minoverlap = 1,connectedPeaks = "min")
saveRDS(ol,"./Significant_P2G_Outputs/Ol.rds")
overlappingPeaks <- ol$overlappingPeaks
peaklist <- ol$peaklist
saveRDS(overlappingPeaks,"./Significant_P2G_Outputs/OverlappingPeaks.rds")
names <- names(overlappingPeaks)[4:6]
print(names)
total <- c(overlappingPeaks[[names[1]]]$overlapFeature,
           overlappingPeaks[[names[2]]]$overlapFeature,
           overlappingPeaks[[names[3]]]$overlapFeature)

pdf("./Significant_P2G_Outputs/ChipPeakAnno_Pie_Overlaps.pdf")
pie1(table(total))
dev.off()

ft.overlap <- overlappingPeaks[[names[1]]]
levels(factor(ft.overlap$overlapFeature))
colnames(ft.overlap)[8:12] <- c("seqnames2","start2","end2","width2","strand2")
ft.overlap$overlap <- ifelse(ft.overlap$overlapFeature == "includeFeature",ft.overlap$width2,"fill")
ft.overlap$overlap <- ifelse(ft.overlap$overlapFeature == "inside",ft.overlap$width,ft.overlap$overlap)
ft.overlap$overlap <- ifelse(ft.overlap$overlapFeature == "overlapEnd",ft.overlap$end2 - ft.overlap$start,ft.overlap$overlap)
ft.overlap$overlap <- ifelse(ft.overlap$overlapFeature == "overlapStart",ft.overlap$end - ft.overlap$start2,ft.overlap$overlap)
levels(factor(ft.overlap$overlap))
ft.overlap$overlap <- as.numeric(ft.overlap$overlap)
hist(ft.overlap$overlap)

oe.overlap <- overlappingPeaks[[names[2]]]
levels(factor(oe.overlap$overlapFeature))
colnames(oe.overlap)[8:12] <- c("seqnames2","start2","end2","width2","strand2")
oe.overlap$overlap <- ifelse(oe.overlap$overlapFeature == "includeFeature",oe.overlap$width2,"fill")
oe.overlap$overlap <- ifelse(oe.overlap$overlapFeature == "inside",oe.overlap$width,oe.overlap$overlap)
oe.overlap$overlap <- ifelse(oe.overlap$overlapFeature == "overlapEnd",oe.overlap$end2 - oe.overlap$start,oe.overlap$overlap)
oe.overlap$overlap <- ifelse(oe.overlap$overlapFeature == "overlapStart",oe.overlap$end - oe.overlap$start2,oe.overlap$overlap)
levels(factor(oe.overlap$overlap))
oe.overlap$overlap <- as.numeric(oe.overlap$overlap)
hist(oe.overlap$overlap)

encode.overlap <- overlappingPeaks[[names[3]]]
levels(factor(encode.overlap$overlapFeature))
colnames(encode.overlap)[8:12] <- c("seqnames2","start2","end2","width2","strand2")
encode.overlap$overlap <- ifelse(encode.overlap$overlapFeature == "includeFeature",encode.overlap$width2,"fill")
encode.overlap$overlap <- ifelse(encode.overlap$overlapFeature == "inside",encode.overlap$width,encode.overlap$overlap)
encode.overlap$overlap <- ifelse(encode.overlap$overlapFeature == "overlapEnd",encode.overlap$end2 - encode.overlap$start,encode.overlap$overlap)
encode.overlap$overlap <- ifelse(encode.overlap$overlapFeature == "overlapStart",encode.overlap$end - encode.overlap$start2,encode.overlap$overlap)
levels(factor(encode.overlap$overlap))
encode.overlap$overlap <- as.numeric(encode.overlap$overlap)
hist(encode.overlap$overlap)

head(ft.overlap[1:3,])
head(oe.overlap[1:3,])
head(encode.overlap[1:3,])

ft.overlap <- ft.overlap[,c(1:12,ncol(ft.overlap),grep("overlapFeature",colnames(ft.overlap)))]
oe.overlap <- oe.overlap[,c(1:12,ncol(oe.overlap),grep("overlapFeature",colnames(oe.overlap)))]
encode.overlap <- encode.overlap[,c(1:12,ncol(encode.overlap),grep("overlapFeature",colnames(encode.overlap)))]

all.overlap <- rbind(ft.overlap,oe.overlap,encode.overlap)

pdf("./Significant_P2G_Outputs/Distribtion_of_All_Overlaps.pdf")
hist(all.overlap$overlap)
dev.off()
types <- c("overlapEnd","overlapStart")
all.overlap <- all.overlap[all.overlap$overlapFeature %in% types,]
pdf("./Significant_P2G_Outputs/Distribution_of_Partial_Overlaps.pdf")
hist(all.overlap$overlap)
dev.off()

pdf("./Significant_P2G_Outputs/ChipPeakAnno_Venn_Overlaps.pdf")
makeVennDiagram(ol,totalTest = tot,connectedPeaks = "min")
dev.off()
venn <- makeVennDiagram(ol,totalTest = tot,connectedPeaks = "min")
saveRDS(venn,"./Significant_P2G_Outputs/venn.rds")
peak.names <- paste0(peaklist$unique.p2g.@seqnames,":",peaklist$unique.p2g.@ranges)

# Find cancer specific peak to gene links 
p2g.df.sub.plot.cancer.kmeans <- p2g.df.sub.plot.cancer.kmeans[p2g.df.sub.plot.cancer.kmeans$peakName %in% peak.names,]
saveRDS(p2g.df.sub.plot.cancer.kmeans,"./Significant_P2G_Outputs/Cancer_specific_P2G_table.rds")

# Plot P2G heatmap for cancer specific distal elements
##############################################################################################################

source("P2G_Heatmap_Distal.R")
test <- plotPeak2GeneHeatmap.distal(proj,peaks=p2g.df.sub.plot.cancer.kmeans,groupBy = "cluster.new",k=length(levels(factor(proj$cluster.new))),
                                    corCutOff = .45,varCutOffATAC = 0,varCutOffRNA = 0,FDRCutOff = 1,nPlot =100000,palGroup=cols,
                                    palATAC = paletteContinuous("solarExtra"),
                                    palRNA = paletteContinuous("solarExtra"))

pdf("./Significant_P2G_Outputs/Peak2Gene_Heatmap_Legend-cancerspecific.pdf",width = 8,height = 12)
draw(test, heatmap_legend_side = "bottom")
dev.off()

p2g.heat.df <- plotPeak2GeneHeatmap.distal(proj,peaks=p2g.df.sub.plot.cancer.kmeans,groupBy = "predictedGroup_ArchR",k=length(levels(factor(proj$predictedGroup_ArchR))),
                                           corCutOff = .45,varCutOffATAC = 0,varCutOffRNA = 0,FDRCutOff = 1,returnMatrices = T,nPlot = 100000)

p2g.df.sub.plot.cancer.kmeans$kmeans <- p2g.heat.df$RNA$kmeansId
saveRDS(p2g.df.sub.plot.cancer.kmeans,"./Significant_P2G_Outputs/Distal_P2Gs_Kmeans-cancerspecific.rds")

pdf("./Significant_P2G_Outputs/GenesPerPeak_histogram-distal-cancerspecific.pdf",width = 5,height = 3.5)
hist(table(p2g.df.sub.plot.cancer.kmeans$idxRNA),main="Distribution of genes per distal peaks")
dev.off()

pdf("./Significant_P2G_Outputs/PeaksPerGene_histogram-distal-cancerspecific.pdf",width = 5,height = 3.5)
hist(table(p2g.df.sub.plot.cancer.kmeans$idxATAC),main="Distribution of distal peaks per gene")
dev.off()
##############################################################################################


