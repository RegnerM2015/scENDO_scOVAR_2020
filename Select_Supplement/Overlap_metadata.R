library(ChIPpeakAnno)

ol <- readRDS("FullCohort_Ol.rds")
overlap <- ol$overlappingPeaks


# Fallopian tube peaks
ft.overlap <- overlap$`unique.p2g.///unique.ft.peaks.`

colnames(ft.overlap)[8:12] <- c("seqnames2","start2","end2","width2","strand2")
ft.overlap$overlap <- ifelse(ft.overlap$overlapFeature == "includeFeature",ft.overlap$width2,"fill")
ft.overlap$overlap <- ifelse(ft.overlap$overlapFeature == "inside",ft.overlap$width,ft.overlap$overlap)
ft.overlap$overlap <- ifelse(ft.overlap$overlapFeature == "overlapEnd",ft.overlap$end2 - ft.overlap$start,ft.overlap$overlap)
ft.overlap$overlap <- ifelse(ft.overlap$overlapFeature == "overlapStart",ft.overlap$end - ft.overlap$start2,ft.overlap$overlap)
type <- c("overlapStart","overlapEnd")
ft.overlap.sub <- ft.overlap[ft.overlap$overlapFeature %in% type,]



# Ovarian peaks
oe.overlap <- overlap$`unique.p2g.///unique.oe.peaks.`

colnames(oe.overlap)[8:12] <- c("seqnames2","start2","end2","width2","strand2")
oe.overlap$overlap <- ifelse(oe.overlap$overlapFeature == "includeFeature",oe.overlap$width2,"fill")
oe.overlap$overlap <- ifelse(oe.overlap$overlapFeature == "inside",oe.overlap$width,oe.overlap$overlap)
oe.overlap$overlap <- ifelse(oe.overlap$overlapFeature == "overlapEnd",oe.overlap$end2 - oe.overlap$start,oe.overlap$overlap)
oe.overlap$overlap <- ifelse(oe.overlap$overlapFeature == "overlapStart",oe.overlap$end - oe.overlap$start2,oe.overlap$overlap)
type <- c("overlapStart","overlapEnd")
oe.overlap.sub <- oe.overlap[oe.overlap$overlapFeature %in% type,]


# ENCODE
encode.overlap <- overlap$`unique.p2g.///encode.all`

colnames(encode.overlap)[8:12] <- c("seqnames2","start2","end2","width2","strand2")
encode.overlap$overlap <- ifelse(encode.overlap$overlapFeature == "includeFeature",encode.overlap$width2,"fill")
encode.overlap$overlap <- ifelse(encode.overlap$overlapFeature == "inside",encode.overlap$width,encode.overlap$overlap)
encode.overlap$overlap <- ifelse(encode.overlap$overlapFeature == "overlapEnd",encode.overlap$end2 - encode.overlap$start,encode.overlap$overlap)
encode.overlap$overlap <- ifelse(encode.overlap$overlapFeature == "overlapStart",encode.overlap$end - encode.overlap$start2,encode.overlap$overlap)
type <- c("overlapStart","overlapEnd")
encode.overlap.sub <- encode.overlap[encode.overlap$overlapFeature %in% type,]



pdf("./Distribution_of_Overlaps-FullCohort.pdf")
par(mfrow=c(3,2)) 
pie1(table(ft.overlap$overlapFeature))
hist(as.numeric(ft.overlap.sub$overlap))
pie1(table(oe.overlap$overlapFeature))
hist(as.numeric(oe.overlap.sub$overlap))
pie1(table(encode.overlap$overlapFeature))
hist(as.numeric(encode.overlap.sub$overlap))
dev.off()






ol <- readRDS("EEC_Ol.rds")
overlap <- ol$overlappingPeaks


# Fallopian tube peaks
ft.overlap <- overlap$`unique.p2g.///unique.ft.peaks.`

colnames(ft.overlap)[8:12] <- c("seqnames2","start2","end2","width2","strand2")
ft.overlap$overlap <- ifelse(ft.overlap$overlapFeature == "includeFeature",ft.overlap$width2,"fill")
ft.overlap$overlap <- ifelse(ft.overlap$overlapFeature == "inside",ft.overlap$width,ft.overlap$overlap)
ft.overlap$overlap <- ifelse(ft.overlap$overlapFeature == "overlapEnd",ft.overlap$end2 - ft.overlap$start,ft.overlap$overlap)
ft.overlap$overlap <- ifelse(ft.overlap$overlapFeature == "overlapStart",ft.overlap$end - ft.overlap$start2,ft.overlap$overlap)
type <- c("overlapStart","overlapEnd")
ft.overlap.sub <- ft.overlap[ft.overlap$overlapFeature %in% type,]



# Ovarian peaks
oe.overlap <- overlap$`unique.p2g.///unique.oe.peaks.`

colnames(oe.overlap)[8:12] <- c("seqnames2","start2","end2","width2","strand2")
oe.overlap$overlap <- ifelse(oe.overlap$overlapFeature == "includeFeature",oe.overlap$width2,"fill")
oe.overlap$overlap <- ifelse(oe.overlap$overlapFeature == "inside",oe.overlap$width,oe.overlap$overlap)
oe.overlap$overlap <- ifelse(oe.overlap$overlapFeature == "overlapEnd",oe.overlap$end2 - oe.overlap$start,oe.overlap$overlap)
oe.overlap$overlap <- ifelse(oe.overlap$overlapFeature == "overlapStart",oe.overlap$end - oe.overlap$start2,oe.overlap$overlap)
type <- c("overlapStart","overlapEnd")
oe.overlap.sub <- oe.overlap[oe.overlap$overlapFeature %in% type,]


# ENCODE
encode.overlap <- overlap$`unique.p2g.///encode.all`

colnames(encode.overlap)[8:12] <- c("seqnames2","start2","end2","width2","strand2")
encode.overlap$overlap <- ifelse(encode.overlap$overlapFeature == "includeFeature",encode.overlap$width2,"fill")
encode.overlap$overlap <- ifelse(encode.overlap$overlapFeature == "inside",encode.overlap$width,encode.overlap$overlap)
encode.overlap$overlap <- ifelse(encode.overlap$overlapFeature == "overlapEnd",encode.overlap$end2 - encode.overlap$start,encode.overlap$overlap)
encode.overlap$overlap <- ifelse(encode.overlap$overlapFeature == "overlapStart",encode.overlap$end - encode.overlap$start2,encode.overlap$overlap)
type <- c("overlapStart","overlapEnd")
encode.overlap.sub <- encode.overlap[encode.overlap$overlapFeature %in% type,]



pdf("./Distribution_of_Overlaps-EEC.pdf")
par(mfrow=c(3,2)) 
pie1(table(ft.overlap$overlapFeature))
hist(as.numeric(ft.overlap.sub$overlap))
pie1(table(oe.overlap$overlapFeature))
hist(as.numeric(oe.overlap.sub$overlap))
pie1(table(encode.overlap$overlapFeature))
hist(as.numeric(encode.overlap.sub$overlap))
dev.off()






ol <- readRDS("HGSOC_Ol.rds")
overlap <- ol$overlappingPeaks


# Fallopian tube peaks
ft.overlap <- overlap$`unique.p2g.///unique.ft.peaks.`

colnames(ft.overlap)[8:12] <- c("seqnames2","start2","end2","width2","strand2")
ft.overlap$overlap <- ifelse(ft.overlap$overlapFeature == "includeFeature",ft.overlap$width2,"fill")
ft.overlap$overlap <- ifelse(ft.overlap$overlapFeature == "inside",ft.overlap$width,ft.overlap$overlap)
ft.overlap$overlap <- ifelse(ft.overlap$overlapFeature == "overlapEnd",ft.overlap$end2 - ft.overlap$start,ft.overlap$overlap)
ft.overlap$overlap <- ifelse(ft.overlap$overlapFeature == "overlapStart",ft.overlap$end - ft.overlap$start2,ft.overlap$overlap)
type <- c("overlapStart","overlapEnd")
ft.overlap.sub <- ft.overlap[ft.overlap$overlapFeature %in% type,]



# Ovarian peaks
oe.overlap <- overlap$`unique.p2g.///unique.oe.peaks.`

colnames(oe.overlap)[8:12] <- c("seqnames2","start2","end2","width2","strand2")
oe.overlap$overlap <- ifelse(oe.overlap$overlapFeature == "includeFeature",oe.overlap$width2,"fill")
oe.overlap$overlap <- ifelse(oe.overlap$overlapFeature == "inside",oe.overlap$width,oe.overlap$overlap)
oe.overlap$overlap <- ifelse(oe.overlap$overlapFeature == "overlapEnd",oe.overlap$end2 - oe.overlap$start,oe.overlap$overlap)
oe.overlap$overlap <- ifelse(oe.overlap$overlapFeature == "overlapStart",oe.overlap$end - oe.overlap$start2,oe.overlap$overlap)
type <- c("overlapStart","overlapEnd")
oe.overlap.sub <- oe.overlap[oe.overlap$overlapFeature %in% type,]


# ENCODE
encode.overlap <- overlap$`unique.p2g.///encode.all`

colnames(encode.overlap)[8:12] <- c("seqnames2","start2","end2","width2","strand2")
encode.overlap$overlap <- ifelse(encode.overlap$overlapFeature == "includeFeature",encode.overlap$width2,"fill")
encode.overlap$overlap <- ifelse(encode.overlap$overlapFeature == "inside",encode.overlap$width,encode.overlap$overlap)
encode.overlap$overlap <- ifelse(encode.overlap$overlapFeature == "overlapEnd",encode.overlap$end2 - encode.overlap$start,encode.overlap$overlap)
encode.overlap$overlap <- ifelse(encode.overlap$overlapFeature == "overlapStart",encode.overlap$end - encode.overlap$start2,encode.overlap$overlap)
type <- c("overlapStart","overlapEnd")
encode.overlap.sub <- encode.overlap[encode.overlap$overlapFeature %in% type,]



pdf("./Distribution_of_Overlaps-HGSOC.pdf")
par(mfrow=c(3,2)) 
pie1(table(ft.overlap$overlapFeature))
hist(as.numeric(ft.overlap.sub$overlap))
pie1(table(oe.overlap$overlapFeature))
hist(as.numeric(oe.overlap.sub$overlap))
pie1(table(encode.overlap$overlapFeature))
hist(as.numeric(encode.overlap.sub$overlap))
dev.off()
