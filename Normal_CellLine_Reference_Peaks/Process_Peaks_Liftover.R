library(GenomicRanges)
library(liftOver)
library(rtracklayer)
library(genomation)
library(plyranges)

# Make Ovarian epithelial peakset
files <- list.files(pattern = "iOSE")
grObject.1 <- readNarrowPeak(files[1])
grObject.2 <- readNarrowPeak(files[2])
grObject.3 <- readNarrowPeak(files[3])
grObject.4 <- readNarrowPeak(files[4])

grObject <- unlist(GenomicRanges::reduce(GRangesList(grObject.1,grObject.2,
                                  grObject.3,grObject.4)))

# Liftover:
chainObject = import.chain("hg19ToHg38.over.chain")
results <- as.data.frame(liftOver(grObject, chainObject))
saveRDS(results,"Ovarian_Epithelial_Cell_line_Peaks.rds")


# Make fallopian tube peakset
files <- list.files(pattern = "iFT")
grObject.1 <- readNarrowPeak(files[1])
grObject.2 <- readNarrowPeak(files[2])
grObject.3 <- readNarrowPeak(files[3])
grObject.4 <- readNarrowPeak(files[4])

grObject <- unlist(GenomicRanges::reduce(GRangesList(grObject.1,grObject.2,
                                                     grObject.3,grObject.4)))

# Liftover:
chainObject = import.chain("hg19ToHg38.over.chain")
results <- as.data.frame(liftOver(grObject, chainObject))
saveRDS(results,"Fallopian_Tube_Cell_line_Peaks.rds")



