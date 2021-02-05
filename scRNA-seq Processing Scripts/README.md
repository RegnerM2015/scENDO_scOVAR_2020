# scRNA-seq Processing Scripts
### /Individual_Samples/*Patient[1-9]_scRNA-seq.R*
1. Preprocessing & QC
1. Feature selection, dimension reduction and clustering 
1. inferCNV
1. SingleR reference-based annotation of cell types 
### /Full_Cohort/Patients1-11_scRNA-seq.R
1. Preprocessing & QC
1. Feature selection, dimension reduction and clustering 
1. Retain inferred CNVs from individual sample processing 
1. Assign cell type labels to clusters based on the majority label within each cluster
1. Verify SingleR cell type labels with cell type gene signatures from PanglaoDB using Seurat's *AddModuleScore()*
### /EEC_Cohort/Patients1-5_scRNA-seq.R
1. Preprocessing & QC
1. Feature selection, dimension reduction and clustering 
1. Retain inferred CNVs from individual sample processing 
1. Assign cell type labels to clusters based on the majority label within each cluster
1. Verify SingleR cell type labels with cell type gene signatures from PanglaoDB using Seurat's *AddModuleScore()*
### /HGSOC_Cohort/Patients8-9_scRNA-seq.R
1. Preprocessing & QC
1. Feature selection, dimension reduction and clustering 
1. Retain inferred CNVs from individual sample processing 
1. Assign cell type labels to clusters based on the majority label within each cluster
1. Verify SingleR cell type labels with cell type gene signatures from PanglaoDB using Seurat's *AddModuleScore()*


#### The scRNA-seq pipeline is also depicted here: (https://github.com/RegnerM2015/scENDO_scOVAR_2020/blob/main/Flowchart_Pipelines.png)
