# Figure_4
### /Figure_4.R
1. Plot scRNA-seq and scATAC-seq UMAPs for the HGSOC Cohort (Patients 8-9)
1. Plot stacked bar charts showing the contribution of each patient to each cell type subcluster
1. Plot heatmap of distal-only peak-to-gene links calculated in the HGSOC Cohort (Patients 8-9)
1. Plot scATAC-seq browser track for gene of interest and corresponding scRNA-seq expression in violin plot

### /HGSOC_LAPTM4B_TFBS.sh
1. Call variants in the Patient 9 11-Epithelial cluster-specific bam file relative to hg38 reference genome
1. Generate Patient 9 11-Epithelial cluster-specific genome sequence
1. Extract enhancer sequences from Patient 9 11-Epithelial cluster-specific genome sequence
1. Perform FIMO motif scanning using Patient-specfic enhancer sequences 

### /TF_expression_FIMO.R
1. Read in FIMO output
1. Sort FIMO TF motif hits by TF expression and q-values
1. Plot violin plots for matching TF expression in scRNA-seq for top ranking TF motifs

