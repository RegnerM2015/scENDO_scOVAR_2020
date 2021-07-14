library(Seurat)
library(ComplexHeatmap)
library(viridis)
library(dplyr)
library(ArchR)
rna <- readRDS("./endo_ovar_All_scRNA_processed.rds")

Idents(rna) <- rna$RNA_snn_res.0.7

levels(factor(rna$RNA_snn_res.0.7))
clusters <- c(6,8,12,14,15,18,24,25,26,29)

rna <- rna[,rna$RNA_snn_res.0.7 %in% clusters]
rna <- NormalizeData(rna)
rna <- ScaleData(rna,features = rownames(rna))
markers <- FindAllMarkers(rna,only.pos =T)
saveRDS(markers,"RNA_markers.rds")
markers <- markers[markers$p_val_adj <= 0.01,]

markers.top <- markers %>% group_by(cluster) %>% top_n( n=30,wt=avg_logFC)

# Downsample cells from each cluster
rna.sub <- subset(rna,downsample =457)
rna.sub <- NormalizeData(rna.sub)
rna.sub <- ScaleData(rna.sub,features = rownames(rna.sub))
rna.sub <- FindVariableFeatures(rna.sub,nfeatures = 2000)


rna.pseudobulk <- data.frame(rownames(rna.sub))
for (i in levels(factor(rna.sub$RNA_snn_res.0.7))){
  cells <- rownames(dplyr::filter(as.data.frame(rna.sub@meta.data),RNA_snn_res.0.7 == i))
  
  rna.sub.new <- rna.sub@assays$RNA@data[,colnames(rna.sub@assays$RNA@data) %in% cells]
  
  rna.bulk <- Matrix::rowSums(rna.sub.new)
  
  rna.pseudobulk$i <- rna.bulk
  colnames(rna.pseudobulk)[dim(rna.pseudobulk)[2]] <- i
  
  print("Iteration complete")
}
mat <- rna.pseudobulk
rownames(mat) <- mat$rownames.rna.sub.
mat <- mat[,-1]
head(rownames(mat))

mat <- mat[rownames(mat) %in% markers.top$gene,]
#mat <- mat[VariableFeatures(rna.sub),]

cluster_anno<-levels(factor(rna.sub$RNA_snn_res.0.7))
library(viridis)
library(ComplexHeatmap)
col_fun = circlize::colorRamp2(c(-2, 0, 2),viridis(n = 3))

heatmapRNA <- Heatmap(t(scale(t(mat))), name = "Expression",
                      cluster_rows = T,
                      show_row_dend = F,
                      clustering_method_columns = "ward.D2",
                      clustering_method_rows="ward.D2",
                      col = col_fun)
plotPDF(heatmapRNA, name = "Heatmap_RNA", width = 8, height = 6)




