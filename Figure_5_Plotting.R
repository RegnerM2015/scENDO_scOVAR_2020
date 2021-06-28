library(ArchR)
library(GenomicRanges)
library(liftOver)
library(rtracklayer)
library(genomation)
library(plyranges)

# Read BRD4 chip and liftover
brd4 <- read.table("./GSE77568_CHIP_DMSO_vs_input.txt")
brd4 <- data.frame(chrom=brd4$V2,start=brd4$V3,end=brd4$V4)
brd4 <- GRanges(brd4)
# Liftover:
chainObject = import.chain("hg19ToHg38.over.chain")
results <- liftOver(brd4, chainObject)
brd4 <- as.data.frame(results)
brd4 <- GRanges(brd4)


# Read H3K27ac chip and liftover
h3k27ac <- rtracklayer::import("GSM2702221_H3K27ac_OVCAR-3.narrowPeak.gz")
# Liftover:
chainObject = import.chain("hg19ToHg38.over.chain")
results <- liftOver(h3k27ac, chainObject)
h3k27ac <- as.data.frame(results)
h3k27ac <- GRanges(h3k27ac)

# Read in other annotation features:
encode.all <- read.delim("./GRCh38-ccREs.bed",header =F)
colnames(encode.all)[1:3] <- c("seqnames","start","end")
encode.all <- makeGRangesFromDataFrame(encode.all)

ft.peaks <- readRDS("./Fallopian_Tube_Cell_line_Peaks.rds")
ft.peaks <- ft.peaks[,3:5]
colnames(ft.peaks)[1:3] <- c("seqnames","start","end")
ft.peaks <- makeGRangesFromDataFrame(ft.peaks)

ov.peaks <- readRDS("./Ovarian_Epithelial_Cell_line_Peaks.rds")
ov.peaks <- ov.peaks[,3:5]
colnames(ov.peaks)[1:3] <- c("seqnames","start","end")
ov.peaks <- makeGRangesFromDataFrame(ov.peaks)

# # Guide RNAs +/- 20bp:
# guides <- data.frame(chrom=rep("chr8",4),
#                      start=c("97743051","97743144","97758770","97759373"),
#                      end=c("97743110","97743202","97758829","97759431"))
# guides <- GRanges(guides)


# Plot both patients separately:
atac <- readRDS("./final_archr_proj_archrGS-P2Gs.rds")

levels(factor(atac$predictedGroup_ArchR))

labels <- c("2-Epithelial cell",
            "3-Epithelial cell",
            "0-Epithelial cell",
            "7-Epithelial cell",
            "11-Epithelial cell",
            "16-Epithelial cell")

atac$new.groups <- ifelse(atac$predictedGroup_ArchR %in% labels, "Cancer","Normal")


atac$new.groups <- ifelse(atac$Sample =="3BAE2L" &  atac$new.groups == "Cancer", "Patient8.cancer",atac$new.groups)
atac$new.groups <- ifelse(atac$Sample =="3E5CFL" &  atac$new.groups == "Cancer", "Patient9.cancer",atac$new.groups)


plot <- plotBrowserTrack(atac,geneSymbol = "LAPTM4B", groupBy = "new.groups",
                         pal=c("gray48","springgreen2","springgreen4"),
                         upstream = 41500,
                         downstream = 10000,features=GRangesList(TrackA=brd4,
                                                                 TrackB=h3k27ac)
)


pdf("LAPTM4B_malignant_normal-Patients-OVCAR3.pdf",width =6,height = 3)
grid::grid.draw(plot[[1]])
dev.off()

plot <- plotBrowserTrack(atac,geneSymbol = "LAPTM4B", groupBy = "new.groups",
                         pal=c("gray48","springgreen2","springgreen4"),
                         upstream = 41500,
                         downstream = 10000,features=GRangesList(TrackA=encode.all,
                                                                 TrackB=ft.peaks,
                                                                 TrackC=ov.peaks)
)


pdf("LAPTM4B_malignant_normal-Patients.pdf",width =6,height = 3.5)
grid::grid.draw(plot[[1]])
dev.off()





# Plot both malignant fractions together:

atac <- readRDS("./final_archr_proj_archrGS-P2Gs.rds")

levels(factor(atac$predictedGroup_ArchR))

labels <- c("2-Epithelial cell",
            "3-Epithelial cell",
            "0-Epithelial cell",
            "7-Epithelial cell",
            "11-Epithelial cell",
            "16-Epithelial cell")

atac$new.groups <- ifelse(atac$predictedGroup_ArchR %in% labels, "2-Cancer","1-Normal")


plot <- plotBrowserTrack(atac,geneSymbol = "LAPTM4B", groupBy = "new.groups",
                         pal=c("gray48","springgreen4"),
                         upstream = 41500,
                         downstream = 10000,features=GRangesList(TrackA=brd4,
                                                                 TrackB=h3k27ac)
)


pdf("LAPTM4B_malignant_normal-OVCAR3.pdf",width =6,height = 3)
grid::grid.draw(plot[[1]])
dev.off()

plot <- plotBrowserTrack(atac,geneSymbol = "LAPTM4B", groupBy = "new.groups",
                         pal=c("gray48","springgreen4"),
                         upstream = 41500,
                         downstream = 10000,features=GRangesList(TrackA=encode.all,
                                                                 TrackB=ft.peaks,
                                                                 TrackC=ov.peaks)
)


pdf("LAPTM4B_malignant_normal.pdf",width =6,height = 3.5)
grid::grid.draw(plot[[1]])
dev.off()



# Plot suppl. malignant fractions IGV figure:
atac <- readRDS("./final_archr_proj_archrGS-P2Gs.rds")

levels(factor(atac$predictedGroup_ArchR))

labels <- c("2-Epithelial cell",
            "3-Epithelial cell",
            "0-Epithelial cell",
            "7-Epithelial cell",
            "11-Epithelial cell",
            "16-Epithelial cell")

atac$new.groups <- ifelse(atac$predictedGroup_ArchR %in% labels, "Cancer","Normal")

idxSample <- BiocGenerics::which(atac$new.groups == "Cancer")
cellsSample <- atac$cellNames[idxSample]
atac <- atac[cellsSample, ]

plot <- plotBrowserTrack(atac,geneSymbol = "LAPTM4B", groupBy = "Sample",
                         pal=c("springgreen2","springgreen4"),
                         upstream = 41500,
                         downstream = 10000
)


pdf("LAPTM4B_malignant_suppl_SNP.pdf",width =6,height = 3)
grid::grid.draw(plot[[1]])
dev.off()



writeLines(capture.output(sessionInfo()), "sessionInfo.txt")
