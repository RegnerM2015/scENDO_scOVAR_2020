
## Figure 2 plotting
The starting inputs for this script are 1) the multi-sample (full cohort) Seurat object and 2) the multi-Sample ArchR project with all peak-to-gene associations, and 3) the peak-to-gene link metadata table. The outputs are the ATAC browser track per cluster for the gene locus of interest, the expression boxplot per cluster, the gene signature boxplot per cluster, the histograms showing the distributions of correlation and p-values, and the histograms showing the distribution of genes per peak and vice versa. 

### [/Figure_2/Figure_2.R](https://github.com/RegnerM2015/scENDO_scOVAR_2020/tree/main/Figure_2)

1. Plot histograms of correlation values and p-values for both observed and null condition
2. Plot bar charts for peaks per gene and genes per peak information
3. Plot scATAC-seq browser track for gene of interest and corresponding scRNA-seq expression in box plot
4. Plot box plot for pathway expression signature calculated by Seurat's *AddModuleScore()*

Interested in more exciting research in cancer genomics? Visit https://www.thefrancolab.org/ to learn more!
