
###########################################################
# Matt Regner
# Franco Lab
# Date: May-December 2020
# 
# Sample: All 
# Description: This script performs the following tasks  
#         1) Edit CSV files to patient specific clusters
#         2) Write updated csv files for input into Cell Ranger
#            Bamslice 
###########################################################

path = '/datastore/nextgenout5/share/labs/francolab/scENDO_scOVAR_Proj/TFSEE_Methods/'

csv.files <- list.files(pattern = "_config.csv")

csv.3533EL <- read.csv(csv.files[1])
csv.3533EL$library_id <- gsub(" ","_",csv.3533EL$library_id)
csv.3533EL$barcodes_csv <- paste0(path,csv.3533EL$library_id,"_barcodes.csv")
csv.3533EL <- csv.3533EL[3,]
write.csv(csv.3533EL,"All_3533EL_config-updated.csv",quote=FALSE,row.names = F)


# Sample has no specific clusters 
csv.3571DL <- read.csv(csv.files[2])


# Sample has no specific clusters 
csv.36186L <- read.csv(csv.files[3])



csv.36639L <- read.csv(csv.files[4])
csv.36639L$library_id <- gsub(" ","_",csv.36639L$library_id)
csv.36639L$barcodes_csv <- paste0(path,csv.36639L$library_id,"_barcodes.csv")
write.csv(csv.36639L,"All_36639L_config-updated.csv",quote=FALSE,row.names = F)


# Sample has no specific clusters 
csv.366C5L <- read.csv(csv.files[5])




csv.37EACL <- read.csv(csv.files[6])
csv.37EACL$library_id <- gsub(" ","_",csv.37EACL$library_id)
csv.37EACL$barcodes_csv <- paste0(path,csv.37EACL$library_id,"_barcodes.csv")
write.csv(csv.37EACL,"All_37EACL_config-updated.csv",quote=FALSE,row.names = F)



csv.38FE7L <- read.csv(csv.files[7])
csv.38FE7L$library_id <- gsub(" ","_",csv.38FE7L$library_id)
csv.38FE7L$barcodes_csv <- paste0(path,csv.38FE7L$library_id,"_barcodes.csv")
write.csv(csv.38FE7L,"All_38FE7L_config-updated.csv",quote=FALSE,row.names = F)




csv.3BAE2L <- read.csv(csv.files[8])
csv.3BAE2L$library_id <- gsub(" ","_",csv.3BAE2L$library_id)
csv.3BAE2L$barcodes_csv <- paste0(path,csv.3BAE2L$library_id,"_barcodes.csv")
csv.3BAE2L <- csv.3BAE2L[1,]
write.csv(csv.3BAE2L,"All_3BAE2L_config-updated.csv",quote=FALSE,row.names = F)




csv.3CCF1L <- read.csv(csv.files[9])
csv.3CCF1L$library_id <- gsub(" ","_",csv.3CCF1L$library_id)
csv.3CCF1L$barcodes_csv <- paste0(path,csv.3CCF1L$library_id,"_barcodes.csv")
write.csv(csv.3CCF1L,"All_3CCF1L_config-updated.csv",quote=FALSE,row.names = F)



csv.3E4D1L <- read.csv(csv.files[10])
csv.3E4D1L$library_id <- gsub(" ","_",csv.3E4D1L$library_id)
csv.3E4D1L$barcodes_csv <- paste0(path,csv.3E4D1L$library_id,"_barcodes.csv")
write.csv(csv.3E4D1L,"All_3E4D1L_config-updated.csv",quote=FALSE,row.names = F)




csv.3E5CFL <- read.csv(csv.files[11])
csv.3E5CFL$library_id <- gsub(" ","_",csv.3E5CFL$library_id)
csv.3E5CFL$barcodes_csv <- paste0(path,csv.3E5CFL$library_id,"_barcodes.csv")
write.csv(csv.3E5CFL,"All_3E5CFL_config-updated.csv",quote=FALSE,row.names = F)

