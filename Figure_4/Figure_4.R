###########################################################
# Matt Regner
# Franco Lab
# Date: April-December 2020
#
# Description: Figure 4, HGSOC enhancer discovery 
# RNA and ATAC data
###########################################################
source("./P2G_Heatmap_Distal.R")
source("./Archr_Peak_Null_Permute.R")
source("./Archr_Peak_RawPval.R")
library(ggplot2)
library(Seurat)
library(scales)
library(forcats)
library(tidyverse)
library(tidyr)
library(ArchR)

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

sampleColors[11] <- "#8fb2c4"# Make blue slightly darker for Figure 4

###########################################################
# Read in RNA data for full cohort
rna <- readRDS("ovar_HGSOC_scRNA_processed.rds")


# Read in matching ATAC data for full cohort (ArchR project after label transfer)
atac <- readRDS("proj_LSI_GeneScores_Annotations_Int.rds")
###########################################################

# Plot Patient UMAPs for RNA/ATAC
###########################################################
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
  scale_color_manual(values = sampleColors[c(8,11)])+
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
  scale_color_manual(values = sampleColors[c(8,11)])+
  guides(colour = guide_legend(override.aes = list(size=6)))+
  ggsave(paste0("Sample_ATAC.pdf"),width = 8,height = 6)
###########################################################




# Plot cell type UMAPs for RNA/ATAC
###########################################################
# RNA UMAP first
rna.df <- as.data.frame(rna@reductions$umap@cell.embeddings)
length(which(rownames(rna.df)==rownames(rna@meta.data)))
rna.df$cell.type <- rna$cell.type
rna.df$cluster <- rna$RNA_snn_res.0.7
rna.df$cluster.new <- factor(rna.df$cluster,levels = c("0","2","3","7","11","15","16","19",
                                                       "1","4","9","10","22",
                                                       "14",
                                                       "5","12","18","21",
                                                       "6","8","13","20","23",
                                                       "17"
                                                       ))


epithelial.cols <- colorRampPalette(c("#A0E989", "#337B24"))
epithelial.cols <- epithelial.cols(8)

fibro.cols <- c("#FB8CAB","#E65C9C","#CF268A","#Af1281","#6B0772")

endo.cols <- "#377EB8"

t.cols <- c("#CCCCCC","#333333","#666666","#999999")

macro.cols <- colorRampPalette(c("#FECC51", "#FD6A02"))
macro.cols <- macro.cols(5)


b.cols <- "#E41A1C"


cols <- c(epithelial.cols,fibro.cols,endo.cols,t.cols,macro.cols,b.cols)

rna.sample.plot <-ggplot(rna.df,aes(x = UMAP_1,y=UMAP_2,color =cluster.new))+
  geom_point(size = .1)+
  theme_classic()+
  theme(plot.title = element_text(face = "bold"))+
  xlab("UMAP_1")+
  ylab("UMAP_2")+ 
  theme(legend.key.size = unit(0.2, "cm"))+
  scale_color_manual(values = cols)+
  guides(colour = guide_legend(override.aes = list(size=3)))
rna.sample.plot +ggsave("Cell_type_RNA.pdf",width = 8,height = 6)
p1 <- ggplot(rna.df,aes(x = UMAP_1,y=UMAP_2,color = cluster.new))+
  geom_point(size = .1)+
  theme_classic()+
  theme(plot.title = element_text(face = "bold"))+
  xlab("UMAP_1")+
  ylab("UMAP_2")+ 
  theme(legend.key.size = unit(0.2, "cm"))+
  scale_color_manual(values = cols)+
  guides(colour = guide_legend(override.aes = list(size=6)))+NoLegend()
LabelClusters(p1,id="cluster.new",color="black",repel = T,size=8)+ggsave("Cell_Type_RNA-labels.pdf",width = 8,height = 7)

# ATAC UMAP second:
atac.df <- plotEmbedding(atac,colorBy = "cellColData",name = "predictedGroup_ArchR",embedding = "UMAP")
atac.df <- as.data.frame(atac.df$data)
atac.df$cell.type <- sub(".*?-","",atac.df$color)

atac.df$cluster <- sub("?-.*","",atac.df$cell.type)
atac.df$cluster.new <- factor(atac.df$cluster,levels = c("0","2","3","7","11","15","16","19",
                                                         "1","4","9","10","22",
                                                         "14",
                                                         "5","12","18","21",
                                                         "6","8","13","20","23",
                                                         "17" 
))



epithelial.cols <- colorRampPalette(c("#A0E989", "#337B24"))
epithelial.cols <- epithelial.cols(8)

fibro.cols <- c("#FB8CAB","#E65C9C","#CF268A","#Af1281","#6B0772")

endo.cols <- "#377EB8"

t.cols <- c("#CCCCCC","#333333","#666666","#999999")

macro.cols <- colorRampPalette(c("#FECC51", "#FD6A02"))
macro.cols <- macro.cols(5)


b.cols <- "#E41A1C"


cols <- c(epithelial.cols,fibro.cols,endo.cols,t.cols,macro.cols,b.cols)


ggplot(atac.df,aes_string(x = "x",y="y",color = "cluster.new"))+
  geom_point(size = .1)+
  theme_classic()+
  theme(plot.title = element_text(face = "bold"))+
  xlab("UMAP_1")+
  ylab("UMAP_2")+
  theme(legend.key.size = unit(0.2, "cm"))+
  scale_color_manual(values = cols)+
  guides(colour = guide_legend(override.aes = list(size=3)))+
  ggsave(paste0("Cell_type_ATAC.pdf"),width = 8,height = 6)

p1 <- ggplot(atac.df,aes(x = x,y=y,color = cluster.new))+
  geom_point(size = .1)+
  theme_classic()+
  theme(plot.title = element_text(face = "bold"))+
  xlab("UMAP_1")+
  ylab("UMAP_2")+ 
  theme(legend.key.size = unit(0.2, "cm"))+
  scale_color_manual(values = cols)+
  guides(colour = guide_legend(override.aes = list(size=6)))+NoLegend()
LabelClusters(p1,id="cluster.new",color="black",repel = T,size=8)+ggsave("Cell_Type_ATAC-labels.pdf",width = 8,height = 7)



p1 <- ggplot(atac.df,aes(x = x,y=y,color = cluster.new))+
  geom_point(size = .1)+
  theme_classic()+
  theme(plot.title = element_text(face = "bold"))+
  xlab("UMAP_1")+
  ylab("UMAP_2")+ 
  theme(legend.key.size = unit(0.2, "cm"))+
  scale_color_manual(values = cols)+
  guides(colour = guide_legend(override.aes = list(size=6)))
LabelClusters(p1,id="cluster.new",color="black",repel = T,size=8)+ggsave("Cell_Type_ATAC-labels-legend.pdf",width = 12,height = 7)

###########################################################

library(dplyr)
# Patient proportion per subcluster in RNA:
###########################################################
meta <- rna@meta.data

df <- meta %>% dplyr::group_by(RNA_snn_res.0.7) %>% dplyr::count(Sample)
colnames(df) <- c("Cluster","Sample","Cells")

# Reorder cluster factor levels to group by cell type 
levels(factor(rna$cell.type))
df %>% 
  dplyr::mutate(cell.type = factor(Cluster,levels = c("0","2","3","7","11","15","16","19",
                                                      "1","4","9","10","22",
                                                      "14",
                                                      "5","12","18","21",
                                                      "6","8","13","20","23",
                                                      "17"))) %>% 
  ggplot(aes(fill=Sample, y=Cells, x= fct_rev(cell.type))) + 
  scale_fill_manual(values = sampleColors[c(8,11)])+
  geom_bar(position="fill", stat="identity")+
  coord_flip()+theme_classic()+xlab("Clusters")+ylab("# of cells")+
  ggsave("Cell_Type_Prop_RNA.pdf",width = 4,height = 11)


# Patient proportion per subcluster in ATAC:
meta <- as.data.frame(atac@cellColData)
meta$predictedGroup_ArchR <- gsub("-.*", "", meta$predictedGroup_ArchR)

df <- meta %>% dplyr::group_by(predictedGroup_ArchR) %>% dplyr::count(Sample)
colnames(df) <- c("Cluster","Sample","Cells")

# Reorder cluster factor levels to group by cell type 
levels(factor(atac$predictedGroup_ArchR))
df %>% 
  dplyr::mutate(cell.type = factor(Cluster,levels = c("0","2","3","7","11","15","16","19",
                                                      "1","4","9","10","22",
                                                      "14",
                                                      "5","12","18","21",
                                                      "6","8","13","20","23",
                                                      "17"))) %>% 
  ggplot(aes(fill=Sample, y=Cells, x= fct_rev(cell.type))) + 
  scale_fill_manual(values = sampleColors[c(8,11)])+
  geom_bar(position="fill", stat="identity")+
  coord_flip()+theme_classic()+xlab("Clusters")+ylab("# of cells")+
  ggsave("Cell_Type_Prop_ATAC.pdf",width = 4,height = 11)
###########################################################

# Write cluster total # of cells to output files
###########################################################
total.atac <- as.data.frame(table(atac$predictedGroup_ArchR))
colnames(total.atac) <- c("Cluster","ATAC cells")

total.rna <- as.data.frame(table(rna$cell.type))
colnames(total.rna) <- c("Cluster","RNA cells")
total.rna$Cluster <- gsub("/","_",total.rna$Cluster)

total <- merge(total.rna,total.atac,by = "Cluster")
WriteXLS::WriteXLS(total,"./Total_Cell_Numbers_per_Cluster.xlsx")
###########################################################




# Suppl. Figure CNV plot
# Plot CNV box: 
meta <- rna@meta.data
meta$cluster <- factor(meta$RNA_snn_res.0.7,levels = rev(c("0","2","3","7","11","15","16","19",
                                                           "1","4","9","10","22",
                                                           "14",
                                                           "5","12","18","21",
                                                           "6","8","13","20","23",
                                                           "17")))

ggplot(meta,aes(x=cluster,y=Total_CNVs,fill=cluster))+geom_boxplot()+coord_flip()+
  theme_classic()+scale_fill_manual(values = rev(cols))+NoLegend()+
  ggsave("CNV_BoxPlot.pdf",width = 4,height = 8)


# Plot browser track for LAPTM4B enhancers (cancer-specifc and cancer-enriched):
# 1) Plot browser track
# 2) Plot matching LAPTM4B expression in scRNA-seq
####################################################################

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

atac <- readRDS("final_archr_proj_archrGS-P2Gs.rds")
# ATAC
levels(factor(atac$predictedGroup_ArchR))
my_levels <- as.character(c(0,2,3,7,11,16,
                            1,4,9,10,
                            14,
                            5,12,18,21,
                            6,8,13,
                            17))

for ( i in levels(factor(atac$predictedGroup_ArchR))){
  num <-  gsub("-.*","",i)
  idx <- match(num,my_levels)
  atac$predictedGroup_ArchR <- str_replace(atac$predictedGroup_ArchR,pattern = i,replacement = paste0(idx,"_",atac$predictedGroup_ArchR))
  print("iter complete")
}

#########################################################################
# Make modified getP2G function:
getPeak2GeneLinks.mod <- function(
  ArchRProj = NULL, 
  corCutOff = 0.45, 
  PValCutOff = 0.0001,
  varCutOffATAC = 0.25,
  varCutOffRNA = 0.25,
  resolution = 1, 
  returnLoops = TRUE
){
  
  .validInput(input = ArchRProj, name = "ArchRProj", valid = "ArchRProject")
  .validInput(input = corCutOff, name = "corCutOff", valid = "numeric")
  .validInput(input = PValCutOff, name = "PValCutOff", valid = "numeric")
  .validInput(input = varCutOffATAC, name = "varCutOffATAC", valid = "numeric")
  .validInput(input = varCutOffRNA, name = "varCutOffRNA", valid = "numeric")
  .validInput(input = resolution, name = "resolution", valid = c("integer", "null"))
  .validInput(input = returnLoops, name = "returnLoops", valid = "boolean")
  
  if(is.null(ArchRProj@peakSet)){
    return(NULL)
  }
  
  if(is.null(metadata(ArchRProj@peakSet)$Peak2GeneLinks)){
    
    return(NULL)
    
  }else{
    
    p2g <- metadata(ArchRProj@peakSet)$Peak2GeneLinks
    p2g <- p2g[which(p2g$Correlation >= corCutOff & p2g$RawPVal <= PValCutOff), ,drop=FALSE]
    
    if(!is.null(varCutOffATAC)){
      p2g <- p2g[which(p2g$VarQATAC > varCutOffATAC),]
    }
    
    if(!is.null(varCutOffRNA)){
      p2g <- p2g[which(p2g$VarQRNA > varCutOffRNA),]
    }
    
    if(returnLoops){
      
      peakSummits <- resize(metadata(p2g)$peakSet, 1, "center")
      geneStarts <- resize(metadata(p2g)$geneSet, 1, "start")
      
      if(!is.null(resolution)){
        summitTiles <- floor(start(peakSummits) / resolution) * resolution + floor(resolution / 2)
        geneTiles <- floor(start(geneStarts) / resolution) * resolution + floor(resolution / 2)
      }else{
        summitTiles <- start(peakSummits)
        geneTiles <- start(geneTiles)
      }
      
      loops <- .constructGR(
        seqnames = seqnames(peakSummits[p2g$idxATAC]),
        start = summitTiles[p2g$idxATAC],
        end = geneTiles[p2g$idxRNA]
      )
      mcols(loops)$value <- p2g$Correlation
      mcols(loops)$FDR <- p2g$FDR
      
      loops <- loops[order(mcols(loops)$value, decreasing=TRUE)]
      loops <- unique(loops)
      loops <- loops[width(loops) > 0]
      loops <- sort(sortSeqlevels(loops))
      
      loops <- SimpleList(Peak2GeneLinks = loops)
      
      return(loops)
      
    }else{
      
      return(p2g)
      
    }
    
  }
  
}



#########################################################################
cancer.p2gs <- readRDS("Cancer_specific_P2G_table.rds")

# plot <- plotBrowserTrack(atac,geneSymbol ="LAPTM4B", groupBy = "predictedGroup_ArchR",
#                          features = GRangesList(TrackA = encode.all,TrackB = ft.peaks,TrackC = ov.peaks), 
#                          loops = getPeak2GeneLinks.mod(atac,corCutOff = 0.45,
#                                                        PValCutOff = 1e-12,varCutOffATAC = 0,
#                                                        varCutOffRNA = 0),upstream = 50000,downstream = 50000)
# 
# pdf("LAPTM4B_final.pdf",width = 6,height = 8)
# grid::grid.draw(plot[[1]])
# dev.off()

plot <- plotBrowserTrack(atac,geneSymbol = "LAPTM4B", groupBy = "predictedGroup_ArchR",
                         features = GRangesList(TrackA = encode.all,TrackB = ft.peaks,TrackC = ov.peaks), 
                         loops = getPeak2GeneLinks.mod(atac,corCutOff = 0.45,
                                                       PValCutOff = 1e-12,varCutOffATAC = 0,
                                                       varCutOffRNA = 0),
                         pal=cols[-c(6,8,13,22,23)],
                         upstream = 41500,
                         downstream = 10000
                         )


pdf("LAPTM4B_final.pdf",width = 6,height = 8)
grid::grid.draw(plot[[1]])
dev.off()

# Plot matching scRNA-seq data
rna <- readRDS("./ovar_HGSOC_scRNA_processed.rds")

my_levels <- as.character(c(0,2,3,7,11,16,
                            1,4,9,10,
                            14,
                            5,12,18,21,
                            6,8,13,
                            17))
rna.sub <- rna[,rna$RNA_snn_res.0.7 %in% my_levels]
# Make violin plots for LAPTM4B
# Relevel object@ident
rna.sub@active.ident <- factor(x =rna.sub$RNA_snn_res.0.7, levels = rev(my_levels))

p1 <- VlnPlot(rna.sub,features = "LAPTM4B",pt.size = 0)+coord_flip()+NoLegend()
p1 <- ggplot(p1$data,aes(y=ident,x=LAPTM4B))+geom_boxplot(aes(fill=ident),lwd=0.45,outlier.size = 0.95,fatten = 0.95)+NoLegend()+
  theme_classic()+NoLegend()+ylab("Cluster #")+xlab("LAPTM4B expression")+ggsave("Vln.pdf",width=3,height = 8)


p1 <- p1$data
colnames(p1)
test <- kruskal.test(LAPTM4B~ident,data = p1)
print(test)
print(test$p.value)


# Check LAPTM4B statistical significance:
test <- FindMarkers(rna.sub,ident.1 = as.character(c(0,2,3,7,11,16)),
                    ident.2 = as.character(c(1,4,9,10,
                                14,
                                5,12,18,21,
                                6,8,13,
                                17)),only.pos = T)
test$genes <- rownames(test)
test <- test[test$gene == "LAPTM4B",]
print(test)


writeLines(capture.output(sessionInfo()), "sessionInfo.txt")

