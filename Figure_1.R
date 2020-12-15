###################################################
# Matt Regner
# Franco Lab 
# Nov 2020
# Description: plot patient UMAPs, cell type UMAPs,
# histology UMAPs, patient proportion per cluster
###################################################

library(ggplot2)
library(Seurat)
library(scales)
library(forcats)

# Make patient sample metadata and color assignments 

sampleColors <- RColorBrewer::brewer.pal(11,"Paired")
sampleColors[11] <- "#8c8b8b"
pie(rep(1,11), col=sampleColors) 


# Color patient tumors to resemble the cancer ribbon color 
sampleColors <- c(sampleColors[5],sampleColors[7],sampleColors[6],sampleColors[8],sampleColors[10],sampleColors[9],sampleColors[4],sampleColors[3],sampleColors[2],sampleColors[11],sampleColors[1])



sampleAnnot <- data.frame(Sample = c("3533EL","3571DL","36186L","36639L",
                                     "366C5L","37EACL","38FE7L","3BAE2L","3CCF1L","3E4D1L","3E5CFL"),
                          Color = sampleColors,
                          Cancer = c("endometrial","endometrial","endometrial","endometrial","endometrial","endometrial",
                                     "ovarian","ovarian","ovarian","ovarian","ovarian"),
                          Histology = c("endometrioid","endometrioid","endometrioid","endometrioid","endometrioid",
                                        "serous","endometrioid","serous","carcinosarcoma","GIST","serous"),
                          BMI = c(39.89,30.5,38.55,55.29,49.44,29.94,34.8,22.13,23.72,33.96,22.37),
                          Age = c(70,70,70,49,62,74,76,61,69,59,59),
                          Race = c("AA","CAU","CAU","CAU","CAU","CAU","CAU","CAU","CAU","CAU","AS"),
                          Stage = c("IA","IA","IA","IA","IA","IIIA","IA","IIB","IVB","IV","IIIC"),
                          Site = c("Endometrium","Endometrium","Endometrium","Endometrium","Endometrium","Ovary","Ovary","Ovary","Ovary","Ovary","Ovary"),
                          Type = c("Endometrial","Endometrial","Endometrial","Endometrial","Endometrial","Endometrial","Ovarian","Ovarian","Ovarian","Gastric","Ovarian"))



# Read in RNA data for full cohort
rna <- readRDS("endo_ovar_All_scRNA_processed.rds")


# Read in matching ATAC data for full cohort (ArchR project after label transfer)
atac <- readRDS("proj_LSI_GeneScores_Annotations_Int.rds")




# Plot Patient UMAPs for RNA/ATAC

# RNA UMAP first
rna.df <- as.data.frame(rna@reductions$umap@cell.embeddings)
length(which(rownames(rna.df)==rownames(rna@meta.data)))
rna.df$Sample <- rna$Sample

rna.sample.plot <-ggplot(rna.df,aes(x = UMAP_1,y=UMAP_2,color = Sample))+
  geom_point(size = .1)+
  theme_classic()+
  theme(plot.title = element_text(face = "bold"))+
  xlab("UMAP_1")+
  ylab("UMAP_2")+ 
  theme(legend.key.size = unit(0.2, "cm"))+
  scale_color_manual(values = sampleColors)+
  guides(colour = guide_legend(override.aes = list(size=6)))
rna.sample.plot +ggsave("Sample_RNA.pdf",width = 8,height = 6)

# ATAC UMAP second:
atac.df <- plotEmbedding(atac,colorBy = "cellColData",name = "Sample",embedding = "UMAP")
atac.df <- as.data.frame(atac.df$data)
atac.df$Sample <- gsub(".*-","",atac.df$color)

ggplot(atac.df,aes_string(x = "x",y="y",color = "Sample"))+
  geom_point(size = .1)+
  theme_classic()+
  theme(plot.title = element_text(face = "bold"))+
  xlab("UMAP_1")+
  ylab("UMAP_2")+
  theme(legend.key.size = unit(0.2, "cm"))+
  scale_color_manual(values = sampleColors)+
  guides(colour = guide_legend(override.aes = list(size=6)))+
ggsave(paste0("Sample_ATAC.pdf"),width = 8,height = 6)




# Plot Histology UMAPs for RNA/ATAC

# Add Histology variable to RNA data
rna.df$Histology <- rna.df$Sample
rna.df$Histology <- str_replace_all(rna.df$Histology, "3533EL", "Endometrioid")
rna.df$Histology <- str_replace_all(rna.df$Histology, "3571DL", "Endometrioid")
rna.df$Histology <- str_replace_all(rna.df$Histology, "36186L", "Endometrioid")
rna.df$Histology <- str_replace_all(rna.df$Histology, "36639L", "Endometrioid")
rna.df$Histology <- str_replace_all(rna.df$Histology, "366C5L", "Endometrioid")
rna.df$Histology <- str_replace_all(rna.df$Histology, "37EACL", "Serous")
rna.df$Histology <- str_replace_all(rna.df$Histology, "38FE7L", "Endometrioid")
rna.df$Histology <- str_replace_all(rna.df$Histology, "3BAE2L", "Serous")
rna.df$Histology <- str_replace_all(rna.df$Histology, "3E5CFL", "Serous")
rna.df$Histology <- str_replace_all(rna.df$Histology, "3E4D1L", "GIST")
rna.df$Histology <- str_replace_all(rna.df$Histology, "3CCF1L", "Carcinosarcoma")


rna.sample.plot <-ggplot(rna.df,aes(x = UMAP_1,y=UMAP_2,color = Histology))+
  geom_point(size = .1)+
  theme_classic()+
  theme(plot.title = element_text(face = "bold"))+
  xlab("UMAP_1")+
  ylab("UMAP_2")+ 
  theme(legend.key.size = unit(0.2, "cm"))+
  scale_color_manual(values = c("gray60","coral2","mediumpurple1","#00CED1"))+
  guides(colour = guide_legend(override.aes = list(size=6)))
rna.sample.plot +ggsave("Histology_RNA.pdf",width = 8,height = 6)


# Add Histology variable to ATAC data
atac.df$Histology <- atac.df$Sample
atac.df$Histology <- str_replace_all(atac.df$Histology, "3533EL", "Endometrioid")
atac.df$Histology <- str_replace_all(atac.df$Histology, "3571DL", "Endometrioid")
atac.df$Histology <- str_replace_all(atac.df$Histology, "36186L", "Endometrioid")
atac.df$Histology <- str_replace_all(atac.df$Histology, "36639L", "Endometrioid")
atac.df$Histology <- str_replace_all(atac.df$Histology, "366C5L", "Endometrioid")
atac.df$Histology <- str_replace_all(atac.df$Histology, "37EACL", "Serous")
atac.df$Histology <- str_replace_all(atac.df$Histology, "38FE7L", "Endometrioid")
atac.df$Histology <- str_replace_all(atac.df$Histology, "3BAE2L", "Serous")
atac.df$Histology <- str_replace_all(atac.df$Histology, "3E5CFL", "Serous")
atac.df$Histology <- str_replace_all(atac.df$Histology, "3E4D1L", "GIST")
atac.df$Histology <- str_replace_all(atac.df$Histology, "3CCF1L", "Carcinosarcoma")


atac.sample.plot <-ggplot(atac.df,aes(x = x,y=y,color = Histology))+
  geom_point(size = .1)+
  theme_classic()+
  theme(plot.title = element_text(face = "bold"))+
  xlab("UMAP_1")+
  ylab("UMAP_2")+ 
  theme(legend.key.size = unit(0.2, "cm"))+
  scale_color_manual(values = c("gray60","coral2","mediumpurple1","#00CED1"))+
  guides(colour = guide_legend(override.aes = list(size=6)))
atac.sample.plot +ggsave("Histology_ATAC.pdf",width = 8,height = 6)





# Plot cell type UMAPs for RNA/ATAC

# RNA
levels(factor(rna$cell.type))
rna.df$cell.type <- rna$cell.type
rna.df$cell.type <- gsub(".*-","",rna.df$cell.type)
levels(factor(rna.df$cell.type))


rna.df$cell.type <- str_replace_all(rna.df$cell.type, "Ciliated", "Epithelial")
rna.df$cell.type <- str_replace_all(rna.df$cell.type, "Endothelia", "Endothelial")
rna.df$cell.type <- str_replace_all(rna.df$cell.type, "Epithelial cell", "Epithelial")
rna.df$cell.type <- str_replace_all(rna.df$cell.type, "Macrophages", "Macrophage")
rna.df$cell.type <- str_replace_all(rna.df$cell.type, "Fibroblast", "Fibroblast/Stromal")
rna.df$cell.type <- str_replace_all(rna.df$cell.type, "Stromal fibroblasts", "Fibroblast/Stromal")
rna.df$cell.type <- str_replace_all(rna.df$cell.type, "Unciliated epithelia 1", "Epithelial")
rna.df$cell.type <- str_replace_all(rna.df$cell.type, "Unciliated epithelia 2", "Epithelial")
rna.df$cell.type <- str_replace_all(rna.df$cell.type, "Lymphocytes", "NK cell")



cols <- brewer.pal(9,"Set1")
cols[6] <- "gold3"
cols[7] <- "gray30"
cols[9] <- "gray60"
ggplot(rna.df,aes(x = UMAP_1,y=UMAP_2,color = cell.type))+
  geom_point(size = .1)+
  theme_classic()+
  theme(plot.title = element_text(face = "bold"))+
  xlab("UMAP_1")+
  ylab("UMAP_2")+ 
  theme(legend.key.size = unit(0.2, "cm"))+
  scale_color_manual(values = cols)+
  guides(colour = guide_legend(override.aes = list(size=6)))+ggsave("Cell_Type_RNA.pdf",width = 8,height = 6)



# ATAC
levels(factor(atac$predictedGroup_ArchR))


atac.df <- plotEmbedding(atac,colorBy = "cellColData",name = "predictedGroup_ArchR",embedding = "UMAP")
atac.df <- as.data.frame(atac.df$data)
atac.df$cell.type<- gsub(".*-","",atac.df$color)


levels(factor(atac.df$cell.type))

atac.df$cell.type <- str_replace_all(atac.df$cell.type, "Ciliated", "Epithelial")
atac.df$cell.type <- str_replace_all(atac.df$cell.type, "Endothelia", "Endothelial")
atac.df$cell.type <- str_replace_all(atac.df$cell.type, "Epithelial cell", "Epithelial")
atac.df$cell.type <- str_replace_all(atac.df$cell.type, "Macrophages", "Macrophage")
atac.df$cell.type <- str_replace_all(atac.df$cell.type, "Fibroblast", "Fibroblast/Stromal")
atac.df$cell.type <- str_replace_all(atac.df$cell.type, "Stromal fibroblasts", "Fibroblast/Stromal")
atac.df$cell.type <- str_replace_all(atac.df$cell.type, "Unciliated epithelia 1", "Epithelial")
atac.df$cell.type <- str_replace_all(atac.df$cell.type, "Unciliated epithelia 2", "Epithelial")
atac.df$cell.type <- str_replace_all(atac.df$cell.type, "Lymphocytes", "NK cell")



cols <- brewer.pal(9,"Set1")
cols[6] <- "gold3"
cols[7] <- "gray30"
cols[9] <- "gray60"
ggplot(atac.df,aes(x =x,y=y,color = cell.type))+
  geom_point(size = .1)+
  theme_classic()+
  theme(plot.title = element_text(face = "bold"))+
  xlab("UMAP_1")+
  ylab("UMAP_2")+ 
  theme(legend.key.size = unit(0.2, "cm"))+
  scale_color_manual(values = cols)+
  guides(colour = guide_legend(override.aes = list(size=6)))+ggsave("Cell_Type_ATAC.pdf",width = 8,height = 6)



#########################################################################################################

# Patient proportion per subcluster in RNA:
meta <- rna@meta.data

df <- meta %>% group_by(RNA_snn_res.0.7) %>% count(Sample)
colnames(df) <- c("Cluster","Sample","Cells")

# Reorder cluster factor levels to group by cell type 
levels(factor(rna$cell.type))
df %>% 
  dplyr::mutate(cell.type = factor(Cluster,levels = c("1","9","10","11","18","19","21","22","23","32","35",
                                                      "0","7","8","12","14","15","16","17","24","25","26","27","29","30",
                                                      "4",
                                                      "6","20","34",
                                                      "2","3","31",
                                                      "5","13",
                                                      "33",
                                                      "28","36"))) %>% 
  ggplot(aes(fill=Sample, y=Cells, x= fct_rev(cell.type))) + 
  geom_bar(position="fill", stat="identity")+
  coord_flip()+theme_classic()+xlab("Clusters")+ylab("# of cells")+
  scale_fill_manual(values = sampleColors)+ggsave("Cell_Type_Prop_RNA.pdf",width = 4,height = 8)


# Patient proportion per subcluster in ATAC:
meta <- as.data.frame(atac@cellColData)
meta$predictedGroup_ArchR <- gsub("-.*", "", meta$predictedGroup_ArchR)

df <- meta %>% group_by(predictedGroup_ArchR) %>% count(Sample)
colnames(df) <- c("Cluster","Sample","Cells")

# Reorder cluster factor levels to group by cell type 
levels(factor(atac$predictedGroup_ArchR))
df %>% 
  dplyr::mutate(cell.type = factor(Cluster,levels = c("1","9","10","11","18","19","21","22","23","32","35",
                                                      "0","7","8","12","14","15","16","17","24","25","26","27","29","30",
                                                      "4",
                                                      "6","20","34",
                                                      "2","3","31",
                                                      "5","13",
                                                      "33",
                                                      "28","36"))) %>% 
  ggplot(aes(fill=Sample, y=Cells, x= fct_rev(cell.type))) + 
  geom_bar(position="fill", stat="identity")+
  coord_flip()+theme_classic()+xlab("Clusters")+ylab("# of cells")+
  scale_fill_manual(values = sampleColors)+ggsave("Cell_Type_Prop_ATAC.pdf",width = 4,height = 8)



# Write cluster total # of cells to output files
total.atac <- as.data.frame(table(atac$predictedGroup_ArchR))
colnames(total.atac) <- c("Cluster","ATAC cells")

total.rna <- as.data.frame(table(rna$cell.type))
colnames(total.rna) <- c("Cluster","RNA cells")


total <- merge(total.rna,total.atac,by = "Cluster")
WriteXLS::WriteXLS(total,"./Total_Cell_Numbers_per_Cluster.xlsx")



VlnPlot(rna,features = "Total_CNVs")

# 
# 
# # Immune cell proportion barchart
# meta.atac <- as.data.frame(proj@cellColData)
# 
# meta.atac$predictedGroup_ArchR <- sub("-", ":", meta.atac$predictedGroup_ArchR)
# meta.atac$predictedGroup_ArchR <- gsub(".*:", "", meta.atac$predictedGroup_ArchR)
# 
# 
# meta.atac$immune <- ifelse(meta.atac$predictedGroup_ArchR == "T cell" | 
#                              meta.atac$predictedGroup_ArchR == "B cell"|
#                              meta.atac$predictedGroup_ArchR == "NK cell" |
#                              meta.atac$predictedGroup_ArchR == "Macrophage","Immune","Non-immune")
# 
# 
# BAE2L.prop <- nrow(meta.atac[meta.atac$immune == "Immune" & meta.atac$Sample == '3BAE2L',])/nrow(meta.atac[meta.atac$Sample == '3BAE2L',])
# E5CFL.prop <- nrow(meta.atac[meta.atac$immune == "Immune" & meta.atac$Sample == '3E5CFL',])/nrow(meta.atac[meta.atac$Sample == '3E5CFL',])
# 
# df <- data.frame(Sample = c("3BAE2L","3E5CFL"),
#                  Proportion = c(BAE2L.prop,E5CFL.prop))
# 
# meta.atac.sub <- dplyr::filter(meta.atac,immune == "Immune")
# library(ggplot2)
# ggplot(data =meta.atac.sub) + 
#   geom_bar(mapping = aes(x =Sample, fill = Proportion))
# 
# 
# 
# meta.atac.1 <- dplyr::filter(meta.atac,Sample == "3BAE2L")
# ggplot(meta.atac.1, aes(immune)) + 
#   geom_bar(aes(y = (..count..)/sum(..count..))) + 
#   scale_y_continuous(labels=scales::percent,limits = c(0,0.6)) +
#   ylab("% immune cells")+ggsave("3BAE2L_atac_immune.pdf",width = 4,height = 3)
# 
# meta.atac.2 <- dplyr::filter(meta.atac,Sample == "3E5CFL")
# ggplot(meta.atac.2, aes(immune)) + 
#   geom_bar(aes(y = (..count..)/sum(..count..))) + 
#   scale_y_continuous(labels=scales::percent,limits = c(0,0.6)) +
#   ylab("% immune cells")+ggsave("3E5CFL_atac_immune.pdf",width = 4,height = 3)
# 
# 
# 
# 
# library(ggplot2)
# ggplot(meta.atac, aes(x=as.factor(Sample), fill=as.factor(immune)))+
#   geom_bar(aes( y=..count../tapply(..count.., ..x.. ,sum)[..x..]), position="dodge" ) +
#   geom_text(aes( y=..count../tapply(..count.., ..x.. ,sum)[..x..], label=scales::percent(..count../tapply(..count.., ..x.. ,sum)[..x..]) ),
#             stat="count", position=position_dodge(0.9), vjust=-0.5)+
#   scale_y_continuous(labels = scales::percent)+ylab("% immune cells")+
#   theme_classic()+ggsave("atac_immune.pdf",width = 4,height = 3)
# 
# 
# 
# library(ggplot2)
# ggplot(meta.rna, aes(x=as.factor(sample), fill=as.factor(immune)))+
#   geom_bar(aes( y=..count../tapply(..count.., ..x.. ,sum)[..x..]), position="dodge" )+
#   geom_text(aes( y=..count../tapply(..count.., ..x.. ,sum)[..x..], label=scales::percent(..count../tapply(..count.., ..x.. ,sum)[..x..]) ),
#           stat="count", position=position_dodge(0.9), vjust=-0.5)+
#   scale_y_continuous(labels = scales::percent)+ylab("% immune cells")+
#   theme_classic()+ggsave("rna_immune.pdf",width = 4,height = 3)
# 
# 
# 
# 
# 
# library(ggplot2)
# ggplot(meta.rna, aes(x=as.factor(sample), fill=as.factor(immune)))+
#   geom_bar(aes( y=..count../tapply(..count.., ..x.. ,sum)[..x..]), position="dodge" ) +
#   geom_text(aes( y=..count../tapply(..count.., ..x.. ,sum)[..x..], label=scales::percent(..count../tapply(..count.., ..x.. ,sum)[..x..]) ),
#             stat="count", position=position_dodge(0.9), vjust=-0.5)+
#   ylab('Percent of Cylinder Group, %') +
#   scale_y_continuous(labels = scales::percent)
# 
# 
