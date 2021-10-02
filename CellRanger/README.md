# CellRanger (version 3.1.0)

### /Example_Patients1-2_scRNA-seq_mkfastq.sh
Demultiplex Illumina base calls into FASTQ files with CellRanger mkfastq:
```
#!/usr/bin/env bash

#SBATCH --job-name HMGVCBGX9
#SBATCH -c 16
#SBATCH --mem 80g
#SBATCH --partition allnodes
#SBATCH --output HMGVCBGX9_demultiplex.job.out
#SBATCH --error HMGVCBGX9_demultiplex.job.err

DATA=/datastore/nextgenout5/share/labs/bioinformatics/seqware/francolab_10x_copy
OUT=/datastore/nextgenout5/share/labs/francolab/scRNA-seq_Endometrial.05.15.2019

cellranger mkfastq --id=HMGVCBGX9 \
                   --run=${DATA}/190514_NS500270_0297_AHMGVCBGX9 \
                   --csv=${OUT}/Samplesheet.csv \
                   --qc \
                   --localcores=16
```

### /Example_Patient1_scRNA-seq_CellRanger-count.sh
Run CellRanger count pipeline for Patient 1 with refdata-cellranger-GRCh38-3.0.0:
```
#!/usr/bin/env bash

#SBATCH --job-name 3533EL-RNA_F6
#SBATCH -c 16
#SBATCH --mem 80g
#SBATCH --partition allnodes
#SBATCH --output 3533EL-RNA_F6.cellranger-count.job.update.out 
#SBATCH --error 3533EL-RNA_F6.cellranger-count.job.update.err 

DATA=/datastore/nextgenout5/share/labs/francolab/Data

cellranger count  --id=3533EL-RNA_F6_Update \
                      --fastqs=./fastq_path/HMGVCBGX9/3533EL-RNA_F6 \
                      --transcriptome=${DATA}/refdata-cellranger-GRCh38-3.0.0 \
                      --sample=3533EL-RNA_F6 \
                      --localcores=16 \
                      --localmem=80
```


### /Example_Patient2_scRNA-seq_CellRanger-count.sh
Run CellRanger count pipeline for Patient 2 with refdata-cellranger-GRCh38-3.0.0:
```
#!/usr/bin/env bash

#SBATCH --job-name 3571DL-RNA_G1
#SBATCH -c 16
#SBATCH --mem 80g
#SBATCH --partition allnodes
#SBATCH --output 3571DL-RNA_G1.cellranger-count.job.update.out 
#SBATCH --error 3571DL-RNA_G1.cellranger-count.job.update.err 

DATA=/datastore/nextgenout5/share/labs/francolab/Data

cellranger count  --id=3571DL-RNA_G1_Update \
                      --fastqs=./fastq_path/HMGVCBGX9/3571DL-RNA_G1 \
                      --transcriptome=${DATA}/refdata-cellranger-GRCh38-3.0.0 \
                      --sample=3571DL-RNA_G1 \
                      --localcores=16 \
                      --localmem=80

```

### Visit https://support.10xgenomics.com/single-cell-gene-expression/software/overview/welcome for more details on the CellRanger workflow

# CellRanger ATAC (version 1.2.0)
### /Example_Patients1-2_scATAC-seq_mkfastq.sh
Demultiplex Illumina base calls into FASTQ files with CellRanger-atac mkfastq:
```
#!/usr/bin/env bash

#SBATCH --job-name H333JBGXB_2
#SBATCH -c 16
#SBATCH --mem 80g
#SBATCH --partition allnodes
#SBATCH --output H333JBGXB_demultiplex.2.job.out
#SBATCH --error H333JBGXB_demultiplex.2.job.err

DATA=/datastore/nextgenout5/share/labs/bioinformatics/seqware/francolab_10x_copy
OUT=/datastore/nextgenout5/share/labs/francolab/scATAC-seq_Endometrial.05.16.2019

cellranger-atac mkfastq --id=H333JBGXB_2 \
                        --run=${DATA}/190515_NS500270_0298_AH333JBGXB \
                        --csv=${OUT}/Samplesheet.2.csv \
                        --qc \
                        --localcores=16
```
### /Example_Patient1_scATAC-seq_CellRanger-count.sh
Run CellRanger-atac count pipeline for Patient 1 with refdata-cellranger-atac-GRCh38-1.2.0:
```
#!/usr/bin/env bash

#SBATCH --job-name 3533EL-ATAC_A3
#SBATCH -c 16
#SBATCH --mem 80g
#SBATCH --partition allnodes
#SBATCH --output 3533EL-ATAC_A3.cellranger-count.job.update.out 
#SBATCH --error 3533EL-ATAC_A3.cellranger-count.job.update.err 

DATA=/datastore/nextgenout5/share/labs/francolab/Data

cellranger-atac count --id=3533EL-ATAC_A3_Update \
                      --fastqs=./fastq_path2/H333JBGXB/3533EL-ATAC_A3 \
                      --reference=${DATA}/refdata-cellranger-atac-GRCh38-1.2.0 \
                      --sample=3533EL-ATAC_A3 \
                      --localcores=16
```
### /Example_Patient2_scATAC-seq_CellRanger-count.sh
Run CellRanger-atac count pipeline for Patient 2 with refdata-cellranger-atac-GRCh38-1.2.0:
```
#!/usr/bin/env bash

#SBATCH --job-name 3571DL-ATAC_A4
#SBATCH -c 16
#SBATCH --mem 80g
#SBATCH --partition allnodes
#SBATCH --output 3571DL-ATAC_A4.cellranger-count.job.update.out 
#SBATCH --error 3571DL-ATAC_A4.cellranger-count.job.update.err 

DATA=/datastore/nextgenout5/share/labs/francolab/Data

cellranger-atac count --id=3571DL-ATAC_A4_Update \
                      --fastqs=./fastq_path2/H333JBGXB/3571DL-ATAC_A4 \
                      --reference=${DATA}/refdata-cellranger-atac-GRCh38-1.2.0 \
                      --sample=3571DL-ATAC_A4 \
                      --localcores=16
```

### Visit https://support.10xgenomics.com/single-cell-atac/software/overview/welcome for more details on the CellRanger ATAC workflow



## Important note about the provided data formats used in this study: 
The processed scRNA-seq data for each patient sample is provided as a filtered feature barcode matrix generated by the Cell Ranger software pipeline. We note that this format is widely used in the single-cell genomics community. For more information about this file format, please visit the Cell Ranger website: https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/output/matrices.

The processed scATAC-seq data for each patient sample is provided as a fragments file, or “BED-like tabular file where each line represents a unique ATAC-seq fragment captured by the assay.” For more information about this file format, please visit the Cell Ranger website: https://support.10xgenomics.com/single-cell-atac/software/pipelines/latest/output/fragments. We provide the fragments files instead of the filtered peak-barcode matrices generated by cellranger-atac, because we used the fragments files as the starting input into our scATAC-seq analysis performed in the ArchR R package(Granja et al., 2021). We did not use the filtered peak-barcode matrices generated by cellranger-atac because this algorithm calls peaks in a pseudo-bulk fashion (i.e. using all signal from all cells present in the sample). This pseudo-bulk approach effectively drowns out the cell-type specific patterns in chromatin accessibility and would harm the contribution of ATAC signal from rare cell types(Granja et al., 2021).

We highlight that our provided data formats are consistent with the data formats provided by Granja et al. (2019) who presented a leukemia dataset of matched scATAC-seq and CITE-seq (combined single-cell antibody-derived tag [scADT] and RNA sequencing [scRNA]): https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE139369. For their processed scATAC-seq data, they provide the fragments for each patient sample. For their processed CITE-seq data for each patient sample, they provide sparse count matrices for the matched scRNA and scADT data. 



Interested in more exciting research in cancer genomics? Visit https://www.thefrancolab.org/ to learn more!
