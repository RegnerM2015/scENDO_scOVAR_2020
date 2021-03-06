# Figure_5
### /TFSEE_Step0_Reformat_FileNames.sh
1. Remove whitespace from filenames
1. Gather marker enhancer bed files for cell type subclusters that are CNV postive and/or epithelial
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


#### *An alternate version of the analysis was performed using only marker enhancers enriched in the 11 cell type subclusters with 100% patient specificity: (https://github.com/RegnerM2015/scENDO_scOVAR_2020/tree/main/Figure_5/Figure_5_Alternate)
