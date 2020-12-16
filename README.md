# A multi-omic single-cell atlas of gynecologic malignancies Regner et al., 2020 

# Link: TBD

# Please cite: TBD

![alt text](https://github.com/RegnerM2015/scENDO_scOVAR_2020/blob/main/Fig1_Overview.png)

# scRNA-seq Processing Scripts/
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

# scATAC-seq Processing Scripts/
### /Full_Cohort/Patients1-11_scATAC-seq.R
1. Preprocessing & QC
1. Feature selection, dimension reduction and label transfering using scRNA-seq cell type subcluster labels
1. Peak calling within each inferred cell type subcluster
1. Peak-to-gene linkage analysis 
1. Generate lists of cell type-specific enhancers for input into Total Functional Score of Enhancer Elements (TFSEE)
### /EEC_Cohort/Patients1-5_scATAC-seq.R
1. Preprocessing & QC
1. Feature selection, dimension reduction and label transfering using scRNA-seq cell type subcluster labels
1. Peak calling within each inferred cell type subcluster
1. Peak-to-gene linkage analysis 
### /HGSOC_Cohort/Patients8-9_scATAC-seq.R
1. Preprocessing & QC
1. Feature selection, dimension reduction and label transfering using scRNA-seq cell type subcluster labels
1. Peak calling within each inferred cell type subcluster
1. Peak-to-gene linkage analysis 

# Figure_1/Figure_1.R
1. Plot scRNA-seq and scATAC-seq UMAPs for the Full Cohort (Patients 1-11)
1. Plot stacked bar charts showing the contribution of each patient to each cell type subcluster

# Figure_2/Figure_2.R
1. Plot heatmap of distal-only peak-to-gene links calculated in the Full Cohort (Patients 1-11)
1. Plot scATAC-seq browser track for gene of interest and corresponding scRNA-seq expression in violin plot
1. Plot violin plot for pathway expression signature calculated by Seurat's *AddModuleScore()*

# Figure_3/Figure_3.R
1. Plot scRNA-seq and scATAC-seq UMAPs for the EEC Cohort (Patients 1-5)
1. Plot stacked bar charts showing the contribution of each patient to each cell type subcluster
1. Plot heatmap of distal-only peak-to-gene links calculated in the EEC Cohort (Patients 1-5)
1. Plot scATAC-seq browser track for gene of interest and corresponding scRNA-seq expression in violin plot

# Figure_4/Figure_4.R
1. Plot scRNA-seq and scATAC-seq UMAPs for the HGSOC Cohort (Patients 8-9)
1. Plot stacked bar charts showing the contribution of each patient to each cell type subcluster
1. Plot heatmap of distal-only peak-to-gene links calculated in the HGSOC Cohort (Patients 8-9)
1. Plot scATAC-seq browser track for gene of interest and corresponding scRNA-seq expression in violin plot

# Figure_5/TFSEE_Step0_Reformat_FileNames.sh
1. Remove whitespace from filenames
1. Concatenate all cell type-specific enhancer lists into one BED file 

# Figure_5/TFSEE_Step1_NonZero_Enhancers.R
1. 

# Figure_5/TFSEE_Step1.5_Reformat_CSV.R
1. Edit input csv files for CellRanger bamslice 

# Figure_5/TFSEE_Step2_Motif_Enrich.sh
1. Run CellRanger bamslice to generate BAM file for each malignant cell type cluster 
1. Find variants in each malignant cell type cluster relative to hg38 assembly
1. bedtools getfasta to extract enhancer sequences (

# Figure_5/TFSEE_Step3_Motif_Pred_Matrix.py

# Figure_5/TFSEE_Step4_Matrix_Calculations_Filter.R


# Supplementary Figure Scripts
Please email regnerm@live.unc.edu for inquiries into the supplemental figure scripts

# Data download: TBD
