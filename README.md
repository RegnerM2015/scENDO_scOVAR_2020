# Title: A multi-omic single-cell atlas of human gynecological malignancies 
### Matthew J. Regner1,3,8, Kamila Wisniewska1,8, Susana Garcia-Recio1, Aatish Thennavan1,4, Raul Mendez-Giraldez1, Venkat S. Malladi5, Gabrielle Hawkins6, Joel S. Parker1,2, Charles M. Perou1,2,7, Victoria L. Bae-Jump1,6, and Hector L. Franco1,2,3,*
### 1.	The Lineberger Comprehensive Cancer Center, University of North Carolina at Chapel Hill, Chapel Hill, NC, 27599, USA
### 2.	The Department of Genetics, School of Medicine, University of North Carolina at Chapel Hill, Chapel Hill, NC, 27599, USA
### 3.	Bioinformatics and Computational Biology Graduate Program, University of North Carolina at Chapel Hill, Chapel Hill, NC, 27599, USA
### 4.	Oral and Craniofacial Biomedicine Program, University of North Carolina at Chapel Hill, Chapel Hill, NC, 27599, USA
### 5.	Department of Bioinformatics, University of Texas Southwestern Medical Center, Dallas, Texas 75390, USA
### 6.	Division of Gynecology Oncology, University of North Carolina at Chapel Hill, Chapel Hill, NC, 27599, USA
### 7.	Department of Pathology and Laboratory Medicine, University of North Carolina at Chapel Hill, Chapel Hill, NC, 27599, USA
### 8.	These authors contributed equally.


### * Corresponding Author 


### Please cite: TBD, submitted to Nature Cancer

![alt text](https://github.com/RegnerM2015/scENDO_scOVAR_2020/blob/main/Cartoon_Fig_1.png)

# Summary of scRNA-seq and scATAC-seq pipelines 
![alt text](https://github.com/RegnerM2015/scENDO_scOVAR_2020/blob/main/Flowchart_Pipelines.png)


# Figure_3
### /Figure_3.R
1. Plot scRNA-seq and scATAC-seq UMAPs for the EEC Cohort (Patients 1-5)
1. Plot stacked bar charts showing the contribution of each patient to each cell type subcluster
1. Plot heatmap of distal-only peak-to-gene links calculated in the EEC Cohort (Patients 1-5)
1. Plot scATAC-seq browser track for gene of interest and corresponding scRNA-seq expression in violin plot

# Figure_4
### /Figure_4.R
1. Plot scRNA-seq and scATAC-seq UMAPs for the HGSOC Cohort (Patients 8-9)
1. Plot stacked bar charts showing the contribution of each patient to each cell type subcluster
1. Plot heatmap of distal-only peak-to-gene links calculated in the HGSOC Cohort (Patients 8-9)
1. Plot scATAC-seq browser track for gene of interest and corresponding scRNA-seq expression in violin plot

# Figure_5
### /TFSEE_Step0_Reformat_FileNames.sh
1. Remove whitespace from filenames
1. Gather marker enhancer bed files for select cell type clusters 
1. Concatenate all cell type-specific enhancer lists into one BED file 
### /TFSEE_Step1_NonZero_Enhancers.R
1. Screen for enhancers that have signal across all malignant cell types 
1. Subset cell type-specific enhancers BED file after screen
### /TFSEE_Step2_Motif_Enrich.sh
1. bedtools getfasta to extract enhancer sequences (cell type-specific enhancers BED file)
1. MEME motif searching and matching (meme/tomtom)
### /TFSEE_Step3_Motif_Pred_Matrix.py
1. Parse MEME/TOMTOM outputs to generate a matrix of TF motif prediction p-values 
### /TFSEE_Step4_Matrix_Calculations_Filter.R
1. Generate enhancer activity matrix, TF expression matrix, and read in TF motif prediction matrix
2. Rlog transform enhancer activity matrix and TF expression matrix
1. Enrich for robustly expressed TFs (mean signal >1st quartile) and for robustly accessible enhancers (mean signal >1st quartile)
1. Scale enhancer activity matrix, TF expression matrix, and -log10(TF motif prediction matrix + 0.5) 
1. Perform TFSEE matrix operations (matrix multiplication followed by element-wise multiplication)
1. Plot TFSEE unsupervised clustering heatmap with TFs marked for druggability status as determined from the canSAR database 
### /Scatter_Sarcoma_Epi_Ovarian.R
1. Compute log2FC in TFSEE score between sarcoma cluster and epithelial cluster 
1. Compute log2FC in TF expression (rlog transformed) between sarcoma and epithelial cluster 
1. Make scatter plot of log2FC in TFSEE score versus log2FC in TF expression 



# Supplementary Figure Scripts
Please email regnerm@live.unc.edu for inquiries into the supplemental figure scripts

# Data download: dbGaP Study Accession: phs002340.v1.p1
