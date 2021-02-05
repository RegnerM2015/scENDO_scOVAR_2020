
ovar_HGSOC_scRNA_processed <- readRDS("../ovar_HGSOC_scRNA_processed.rds")


hgsoc <-subset(ovar_HGSOC_scRNA_processed,cell.type == "11-Epithelial cell")

counts <- hgsoc@assays$RNA@data
counts.bulk <- as.data.frame(rowSums(counts))

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


VlnPlot(ovar_HGSOC_scRNA_processed,features = c("KLF6","KLF6","YY1","IRF1"),
        ncol = 4,pt.size = 0,idents = "11-Epithelial cell",same.y.lims = T)+
  ggsave("Enhancer_1_vln.pdf",width = 8,height = 4)

VlnPlot(ovar_HGSOC_scRNA_processed,features = c("CEBPD","YY1","YY1","ATF4"),ncol = 4,pt.size = 0,
        idents = "11-Epithelial cell",same.y.lims = T)+
  ggsave("Enhancer_2_vln.pdf",width = 8,height = 4)
  
  
  
