library(ArchR)
library(utils)

cancer.peaks <- readRDS("../Cancer_specific_P2G_table.rds")
distal.peaks <- readRDS("../Distal_P2Gs_Kmeans.rds") 
all.peaks <- readRDS("../All_P2G_Observed.rds")

# Enhancer 1 coords: chr8:97736256-97736756
enh.1 <- data.frame(chrom="chr8",start="97736256",end="97736756")
write.table(enh.1,"enhancer_1.bed",row.names = F,col.names = F,quote = F,sep = "\t")

# Enhancer 2 coords: chr8:97741891-97742391
enh.2 <- data.frame(chrom="chr8",start="97741891",end="97742391")
write.table(enh.2,"enhancer_2.bed",row.names = F,col.names = F,quote = F,sep = "\t")

# Enhancer 3 coords: chr8:97757746-97758246
enh.3 <- data.frame(chrom="chr8",start="97757746",end="97758246")
write.table(enh.3,"enhancer_3.bed",row.names = F,col.names = F,quote = F,sep = "\t")


# Enhancer 4 coords: chr8:97760127-97760627
enh.4 <- data.frame(chrom="chr8",start="97760127",end="97760627")
write.table(enh.4,"enhancer_4.bed",row.names = F,col.names = F,quote = F,sep = "\t")


# Enhancer 5 coords: chr8:97763518-97764018
enh.5 <- data.frame(chrom="chr8",start="97763518",end="97764018")
write.table(enh.5,"enhancer_5.bed",row.names = F,col.names = F,quote = F,sep = "\t")


# Promoter coords: chr8:97773165-97773665 / chr8:97774858-97775358
prom <- data.frame(chrom="chr8",start="97773165",end="97775358")
write.table(prom,"promoter.bed",row.names = F,col.names = F,quote = F,sep = "\t")

