import Bio.motifs
import pandas as pd
import numpy as np
import csv
import math
import re
import string
from sklearn import preprocessing
from sklearn.decomposition import PCA
import matplotlib.pyplot as plt
import seaborn as sns
import scipy
from scipy.stats import pearsonr, spearmanr, cumfreq




### Parse MEME and TOMTOM Motif data

# Loop through meme output
meme_cell_dict = {"./enhancers_Marker_Enhancers_ArchR_0-Fibroblast-updated.bed_meme":"./enhancers_Marker_Enhancers_ArchR_0-Fibroblast-updated.bed_tomtom",
"./enhancers_Marker_Enhancers_ArchR_10-Epithelial_cell-updated.bed_meme":"./enhancers_Marker_Enhancers_ArchR_10-Epithelial_cell-updated.bed_tomtom",
"./enhancers_Marker_Enhancers_ArchR_16-Fibroblast-updated.bed_meme":"./enhancers_Marker_Enhancers_ArchR_16-Fibroblast-updated.bed_tomtom",
"./enhancers_Marker_Enhancers_ArchR_17-Epithelial_cell-updated.bed_meme":"./enhancers_Marker_Enhancers_ArchR_17-Epithelial_cell-updated.bed_tomtom",
"./enhancers_Marker_Enhancers_ArchR_19-Epithelial_cell-updated.bed_meme":"./enhancers_Marker_Enhancers_ArchR_19-Epithelial_cell-updated.bed_tomtom",
"./enhancers_Marker_Enhancers_ArchR_21-Unciliated_epithelia_1-updated.bed_meme":"./enhancers_Marker_Enhancers_ArchR_21-Unciliated_epithelia_1-updated.bed_tomtom",
"./enhancers_Marker_Enhancers_ArchR_27-Fibroblast-updated.bed_meme":"./enhancers_Marker_Enhancers_ArchR_27-Fibroblast-updated.bed_tomtom",
"./enhancers_Marker_Enhancers_ArchR_3-Epithelial_cell-updated.bed_meme":"enhancers_Marker_Enhancers_ArchR_3-Epithelial_cell-updated.bed_tomtom",
"./enhancers_Marker_Enhancers_ArchR_31-Unciliated_epithelia_1-updated.bed_meme":"./enhancers_Marker_Enhancers_ArchR_31-Unciliated_epithelia_1-updated.bed_tomtom",
"./enhancers_Marker_Enhancers_ArchR_34-Epithelial_cell-updated.bed_meme":"./enhancers_Marker_Enhancers_ArchR_34-Epithelial_cell-updated.bed_tomtom",
"./enhancers_Marker_Enhancers_ArchR_9-Epithelial_cell-updated.bed_meme":"./enhancers_Marker_Enhancers_ArchR_9-Epithelial_cell-updated.bed_tomtom"}


# Read Target ID to Motif into dictionary
motif_id_dict = {}
with open("./JASPAR_2020_Homo_Sapiens.txt", "r") as data:
    motif_ids = csv.DictReader(data, delimiter="\t")
    for line in motif_ids:
        motif_id_dict[line['ID']] = line['NAME']

meme_tomtom = pd.DataFrame()
for meme,tom in meme_cell_dict.items():
    # load meme output
    meme_file = '%s/meme.xml' % (meme)
    record = Bio.motifs.parse(open(meme_file),'meme')
    # Loop through all motifs and make dataframe
    meme_positions = pd.DataFrame()
    for motif in record:
        name = motif.name
        ones = [1] * len(motif.instances)
        names = []
        for instance in motif.instances:
            names.append(instance.sequence_name)
        new = pd.DataFrame({name: ones},index = names)
        temp = pd.concat([meme_positions, new], axis=1).fillna(0)
        meme_positions = temp
    # Read tomtom file
    tomtom_file = "%s/tomtom.txt" % (tom)
    tomtom_dict = {}
    with open(tomtom_file, "r") as data:
        tomtom = csv.DictReader(data, delimiter="\t")
        for line in tomtom:
            target = line['Target ID']
            motif = line['#Query ID']
            pval = float(line['p-value'])
            tfs = motif_id_dict[target].upper()
            motif_pvalue = { motif:  [pval]}
            # JASPAR :: means that any TF can either protein, split the protein
            tf_list = tfs.split("::")
            for tf in tf_list:
                # Reduce split form splice to single value [ID]_#
                single_isoform = tf.split("_")[0]
                if single_isoform in tomtom_dict.keys():
                    if motif in tomtom_dict[single_isoform].keys():
                        tomtom_dict[single_isoform][motif].append(pval)
                    else:
                        tomtom_dict[single_isoform].update(motif_pvalue)
                else:
                    tomtom_dict[single_isoform] = motif_pvalue
    # Make dataframe
    tomtom_motif = pd.DataFrame()
    for key,motif in tomtom_dict.items():
        pvalue_dict = {}
        # Loop through motifs to see if length greater than 1, if so do pvalue scaling
        for m,p in motif.items():
            if len(p) > 1:
                stouffer_statistic, stouffer_pval = scipy.stats.combine_pvalues(p,method = 'stouffer', weights = None)
                pvalue_dict[m] = stouffer_pval
            else:
                pvalue_dict[m] = p[0]
        pvalues = np.array(list(pvalue_dict.values()))
        new = pd.DataFrame({key: pvalues},index = pvalue_dict.keys())
        temp = pd.concat([tomtom_motif, new], axis=1).fillna(0).sort_index(level=int)
        tomtom_motif = temp
    # Reorder
    tomtom_motif_reorder = tomtom_motif.reindex( list(meme_positions.columns.values)).fillna(0)
    # dot product
    meme_tomtom_cell = meme_positions.dot(tomtom_motif_reorder)
    # Scale and add
    vals = meme_tomtom_cell.max(axis=0)
    print(vals)
    meme_tomtom_cell = np.log10(meme_tomtom_cell+0.01)
    meme_tomtom_cell = np.negative(meme_tomtom_cell)
    scaler = preprocessing.MinMaxScaler()
    meme_tomtom_cell_transform = meme_tomtom_cell.T
    print(meme_tomtom_cell_transform.shape)
    norm = scaler.fit_transform(meme_tomtom_cell_transform.values) # norm across enhancers for each enhancer
    meme_tomtom_cell_std = pd.DataFrame(data=norm.T, columns=list(meme_tomtom_cell.columns.values), index = meme_tomtom_cell.index )
    # Add to previous data
    temp = meme_tomtom.add(meme_tomtom_cell_std, fill_value=0).fillna(0).sort_index(level=int)
    meme_tomtom = temp



# Transform meme tom_tom
motif_enhancers = meme_tomtom.T

# Rename column headers
#motif_enhancers.rename(columns=lambda x: x.split('-')[0], inplace=True)
#motif_enhancers.rename(columns=lambda x: x.replace(':', "_"), inplace=True)

# Standardize to range 0-1
scaler = preprocessing.MinMaxScaler()
motif_enhancers_transform = motif_enhancers.T
print(motif_enhancers_transform.shape)
norm = scaler.fit_transform(motif_enhancers_transform.values) # norm across enhancers for each enhancer
motif_enhancers_scaled = pd.DataFrame(data=norm.T, columns=list(motif_enhancers.columns.values), index = motif_enhancers.index)



print(motif_enhancers.index)
motif_enhancers.to_csv("motif_enhancers.txt",sep = "\t")
motif_enhancers_scaled.to_csv("motif_enhancers_scaled.txt",sep = "\t")
