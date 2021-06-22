
library(ArchR)
library(ChIPpeakAnno)
library(stats)
ArchR::addArchRThreads(threads = 32)
source("./Archr_Peak_Null_Permute.R")
source("./Archr_Peak_RawPval.R")
source("./getMatrixFromProject.R")
h5disableFileLocking()
SAMPLE.ID <- "All"
key.word.1 <- "pithelia"
key.word.2 <- "-Ciliated"
proj <- readRDS("./final_archr_proj_archrGS.rds")

p2g.df.obs <- readRDS("./All_P2G_Observed.rds")
#Subset to postive correlation P2Gs
p2g.df.sub <- dplyr::filter(p2g.df.obs,RawPVal <=1e-12 & Correlation >= 0.45& peakType == "Distal")


# Marker Peaks for malignant clusters with 100% patient specificity 
markersPeaks <- getMarkerFeatures(
  ArchRProj = proj, 
  useMatrix = "PeakMatrix", 
  groupBy = "predictedGroup_ArchR",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon"
)

markerList <- getMarkers(markersPeaks, cutOff = "FDR <= 0.01 & Log2FC >= 1.25")
saveRDS(markerList,"markerspeaks.rds")


names(markerList)

idx.1 <- grep(key.word.1,names(markerList))
idx.2 <- grep(key.word.2,names(markerList))
idx.3 <- grep("16-Fibroblast",names(markerList))
idx.4 <- grep("0-Fibroblast",names(markerList))
idx.5 <- grep("27-Fibroblast",names(markerList))

idx <- c(idx.1,idx.2,idx.3,idx.4,idx.5)

markerList.sub <- markerList[idx]

names(markerList.sub)

# Only tumor cell clusters that are 100% patient specific 

df <- as.data.frame(proj@cellColData) %>% dplyr::group_by(predictedGroup_ArchR) %>% 
  dplyr::count(Sample) 

idents <- table(df$predictedGroup_ArchR)
idents.sub <- idents[idents ==1 & names(idents) %in% names(markerList.sub)]

markerList.sub <- markerList.sub[names(idents.sub)] 

# Intersect marker peaks for each cluster with P2Gs and write to bed file
for ( i in names(markerList.sub)){
  data <- unlist(markerList.sub[i])
  markers <- paste0(data$seqnames,":",data$start,"-",data$end)
  markers.int <- intersect(p2g.df.sub$peakName,markers)
  markers.int <- stringr::str_replace(string=markers.int,pattern = ":",replacement = "-")
  bed <- data.frame(do.call(rbind, strsplit(markers.int, "-", fixed=TRUE)))
  write.table(bed,paste0("Marker_Enhancers_ArchR_",i,".bed"),row.names = F,col.names = F,quote = F,sep = "\t")
}

# Concatenate marker peaks from each cluster and remove redundant peaks
data <- unlist(markerList.sub)
markers <- paste0(data$seqnames,":",data$start,"-",data$end)
markers.int <- intersect(p2g.df.sub$peakName,markers)
markers.int <- unique(markers.int)
markers.int <- stringr::str_replace(string=markers.int,pattern = ":",replacement = "-")
bed <- data.frame(do.call(rbind, strsplit(markers.int, "-", fixed=TRUE)))
write.table(bed,"Marker_Enhancers_ArchR-nodups.bed",row.names = F,col.names = F,quote = F,sep = "\t")



writeLines(capture.output(sessionInfo()), "sessionInfo_Intersect_DiffPeaks_P2Gs.txt")