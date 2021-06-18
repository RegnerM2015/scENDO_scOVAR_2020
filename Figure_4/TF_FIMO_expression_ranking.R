library(Seurat)
library(dplyr)
library(ggplot2)
library(forcats)


ovar_HGSOC_scRNA_processed <- readRDS("../ovar_HGSOC_scRNA_processed.rds")

labels <- c("0-Epithelial cell","7-Epithelial cell","11-Epithelial cell","16-Epithelial cell")


hgsoc <- ovar_HGSOC_scRNA_processed[,ovar_HGSOC_scRNA_processed$cell.type %in% labels] 

counts <- hgsoc@assays$RNA@data
counts.bulk <- as.data.frame(rowSums(as.matrix(counts)))

counts.bulk$gene <- rownames(counts.bulk)



fimo_outs <- list.files(pattern = "_fimo")

for (i in fimo_outs){
  fimo <- read.table(paste0(i,"/fimo.txt"))
  counts.bulk.new <- counts.bulk[counts.bulk$gene %in% fimo$V2,]
  
  colnames(fimo)[2] <- "gene"
  
  test <- merge(fimo,counts.bulk.new,by = "gene")
  colnames(test)[11] <- "Expr"
  
  test <- dplyr::filter(test,V9 <= 0.1)
  test <- dplyr::arrange(test,desc(Expr))
  write.table(test,paste0(i,"_w_expr.txt"),sep = "\t",quote = F,row.names = F,col.names = F)
}

# Make boxplot for enhancer 2 TF expression levels:
enh.2.tfs.1 <- data.frame(Data = hgsoc@assays$RNA@data[12582,],
                          TF ="A.KLF6")
print(rownames(hgsoc)[12582])
enh.2.tfs.2 <- data.frame(Data = hgsoc@assays$RNA@data[15548,],
                          TF="B.YY1")
print(rownames(hgsoc)[15548])
enh.2.tfs.3 <- data.frame(Data = hgsoc@assays$RNA@data[6857,],
                          TF="C.TFAP2A")
print(rownames(hgsoc)[6857])


# Make boxplot for enhancer 4 TF expression levels:
enh.4.tfs.1 <- data.frame(Data = hgsoc@assays$RNA@data[9985,],
                          TF ="D.CEBPD")
print(rownames(hgsoc)[9985])
enh.4.tfs.2 <- data.frame(Data = hgsoc@assays$RNA@data[15548,],
                          TF="E.YY1")
print(rownames(hgsoc)[15548])
enh.4.tfs.3 <- data.frame(Data = hgsoc@assays$RNA@data[6857,],
                          TF="F.TFAP2A")
print(rownames(hgsoc)[6857])

# Make boxplot for promoter TF expression levels:
prom.tfs.1 <- data.frame(Data = hgsoc@assays$RNA@data[12582,],
                         TF ="G.KLF6")
print(rownames(hgsoc)[12582])
prom.tfs.2 <- data.frame(Data = hgsoc@assays$RNA@data[15260,],
                         TF="H.HIF1A")
print(rownames(hgsoc)[15260])
prom.tfs.3 <- data.frame(Data = hgsoc@assays$RNA@data[6925,],
                         TF="I.SOX4")
print(rownames(hgsoc)[6925])
prom.tfs.4 <- data.frame(Data = hgsoc@assays$RNA@data[15548,],
                         TF="J.YY1")
print(rownames(hgsoc)[15548])


tfs <- rbind(enh.2.tfs.1,enh.2.tfs.2,enh.2.tfs.3,
                  enh.4.tfs.1,enh.4.tfs.2,enh.4.tfs.3,
                  prom.tfs.1,prom.tfs.2,prom.tfs.3,prom.tfs.4)
ggplot(tfs,aes(x=fct_rev(TF),y=Data))+
  geom_boxplot(lwd=0.45,outlier.size = 0.95,fatten = 0.95)+theme_classic()+ylim(c(0,4))+
  coord_flip()+
  ggsave("Vln.pdf",width =4,height = 8)


writeLines(capture.output(sessionInfo()), "sessionInfo.txt")

