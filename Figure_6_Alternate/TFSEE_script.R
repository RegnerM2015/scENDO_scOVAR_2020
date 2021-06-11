TFSEE <- function(peak.mat,activity.mat,
                       expr.mat,motif.p.val,genome,motif.list,return.names){
  if ( is.null(expr.mat) == TRUE){
    
    # Get motif matches for example motifs in peaks
    
    
    if (return.names == TRUE){
      
      motif_ix <- matchMotifs(motif.list,GRanges(rownames(peak.mat)), genome = genome,out = "matches",p.cutoff = motif.p.val)
      print("Finding motifs in peaks...")
      motif.mat <- motifMatches(motif_ix)
      head(motif.mat[,1:5])
      
      colnames(motif.mat) <- motif_ix$name
      rownames(motif.mat) <- rownames(peak.mat)
      
      length(colnames(motif.mat))
      length(unique(colnames(motif.mat)))
      
      tf.names <- gsub("\\s*\\([^\\)]+\\)","",as.character(colnames(motif.mat)))
      tf.names <- str_split(tf.names,pattern = "::")
      tf.names <- unlist(tf.names)
      
      remove <- base::setdiff(tf.names,rownames(activity.mat))
      tf.names <- base::setdiff(tf.names,remove)
      
      length(tf.names)
      length(intersect(tf.names,rownames(activity.mat)))
      
      dim(motif.mat)
      
      expression.matrix <- matrix(, nrow=0, ncol = ncol(activity.mat) )
      for ( i in 1:length(colnames(motif.mat))){
        
        if (as.character(colnames(motif.mat)[i]) == "EWSR1-FLI1"){
          
          tf.feature <- str_split(tf.name,pattern = "-")
        }else{
          tf.name <- gsub("\\s*\\([^\\)]+\\)","",as.character(colnames(motif.mat)[i]))
          tf.feature <- str_split(tf.name,pattern = "::")
        }
        
        
        if (length(tf.feature[[1]]) > 1){
          print("Creating TF dimer metagene")
          
          expr.meta <-  matrix(, nrow=0, ncol = ncol(activity.mat) )
          for ( k in 1:length(tf.feature[[1]])){
            
            if ( tf.feature[[1]][k] %in% remove){
              
              if (tf.feature[[1]][k] %in% rownames(activity.mat)){
                print("TF not found in activity matrix!")
              }else{
                print("TF not found in activity matrix!")
              }
              
              
            }else{
              expression.vector <- activity.mat[tf.feature[[1]][k],]
              expr.meta <- rbind(expr.meta,expression.vector)
              
            }
            
          }
          if ( nrow(expr.meta) > 0){
            expression.meta.min <- colMins(expr.meta)
            expression.matrix <- rbind(expression.matrix,expression.meta.min)
            rownames(expression.matrix)[length(rownames(expression.matrix))] <- colnames(motif.mat)[i]
            print("TF dimer metagene made")
          }else{
            print("Skipping TF dimer")
          }
          
          
          
        }else{
          if ( tf.feature[[1]][1] %in% remove){
            
            if (tf.feature[[1]][1] %in% rownames(activity.mat)){
              print("TF not found in activity matrix!")
            }else{
              print("TF not found in activity matrix!")
            }
            
          }else{
            expression.vector <- activity.mat[tf.feature[[1]][1],]
            expression.matrix <- rbind(expression.matrix,expression.vector)
            rownames(expression.matrix)[length(rownames(expression.matrix))] <- colnames(motif.mat)[i]
          }
          
        }
        
      }
      expr.mat <- expression.matrix
      
      motif.mat <- motif.mat[,colnames(motif.mat) %in% rownames(expr.mat)]
      
      
      length(which(colnames(motif.mat)==rownames(expr.mat)))
      
      
      dir.create("./TFSEE_Outputs")
      print("Writing intermediate outputs to TFSEE_Outputs...")
      saveRDS(peak.mat,"./TFSEE_Outputs/peak_matrix.rds")
      
      saveRDS(expr.mat,"./TFSEE_Outputs/expression_matrix.rds")
      
      saveRDS(motif.mat,"./TFSEE_Outputs/motif_matrix.rds")
      ##########################################################
      peak.mat <- readRDS("./TFSEE_Outputs/peak_matrix.rds")
      expression.mat <- readRDS("./TFSEE_Outputs/expression_matrix.rds")
      motif.mat <- readRDS("./TFSEE_Outputs/motif_matrix.rds")
      
      
      dim(motif.mat)
      dim(peak.mat)
      dim(expression.mat)
      
      
      peak.mat <- t(peak.mat)
      expression.mat <- t(expression.mat)
      
      
      head(motif.mat[,1:3])
      head(peak.mat[,1:3])
      head(expression.mat[,1:3])
      
      print("Rescaling peak accessibility matrix...")
      peak.mat <- apply(peak.mat, 2, rescale,to=c(0,1))
      print("Rescaling TF expression matrix...")
      expression.mat <- apply(expression.mat, 2, rescale,to=c(0,1))
      
      
      dim(motif.mat)
      dim(peak.mat)
      dim(expression.mat)
      
      motif.mat <- t(motif.mat)
      peak.mat <- t(peak.mat)
      expression.mat <- t(expression.mat)
      
      dim(motif.mat)
      dim(peak.mat)
      dim(expression.mat)
      
      
      length(which(rownames(peak.mat)==colnames(motif.mat)))
      
      print("Generating raw motif activity matrix (this may take a few minutes)")
      int <- motif.mat%*%peak.mat
      
      dim(int)
      dim(expression.mat)
      
      print("Saving raw motif activity matrix")
      saveRDS(int,"./TFSEE_Outputs/raw_motif_activity_intermediate_matrix.rds")
      
      length(which(rownames(expression.mat)==rownames(int)))
      length(which(colnames(expression.mat)==colnames(int)))
      
      
      expression.mat <- expression.mat[order(match(rownames(expression.mat), rownames(int))), ]
      
      dim(int)
      dim(expression.mat)
      int <- int[order(match(rownames(int), rownames(expression.mat))), ]
      
      length(which(rownames(expression.mat)==rownames(int)))
      length(which(colnames(expression.mat)==colnames(int)))
      
      print("Generating final TFSEE activity matrix")
      res <- int*expression.mat
      
      print("Saving final TFSEE activity matrix")
      saveRDS(res,"./TFSEE_Outputs/TFSEE_final.rds")
      print("TFSEE finished successfully!")
      return(res)
      
      
    }else{
      
      motif_ix <- matchMotifs(motif.list,GRanges(rownames(peak.mat)), genome = genome,out = "matches",p.cutoff = motif.p.val)
      print("Finding motifs in peaks...")
      motif.mat <- motifMatches(motif_ix)
      head(motif.mat[,1:5])
      
      motif.id <- colnames(motif.mat)
      names(motif.id) <- motif_ix$name
      colnames(motif.mat) <- motif_ix$name
      rownames(motif.mat) <- rownames(peak.mat)
      
      length(colnames(motif.mat))
      length(unique(colnames(motif.mat)))
      
      tf.names <- gsub("\\s*\\([^\\)]+\\)","",as.character(colnames(motif.mat)))
      tf.names <- str_split(tf.names,pattern = "::")
      tf.names <- unlist(tf.names)
      
      remove <- base::setdiff(tf.names,rownames(activity.mat))
      tf.names <- base::setdiff(tf.names,remove)
      
      # All
      length(tf.names)
      length(intersect(tf.names,rownames(activity.mat)))
      
      dim(motif.mat)
      
      expression.matrix <- matrix(, nrow=0, ncol = ncol(activity.mat) )
      for ( i in 1:length(colnames(motif.mat))){
        
        if (as.character(colnames(motif.mat)[i]) == "EWSR1-FLI1"){
          
          tf.feature <- str_split(tf.name,pattern = "-")
        }else{
          tf.name <- gsub("\\s*\\([^\\)]+\\)","",as.character(colnames(motif.mat)[i]))
          tf.feature <- str_split(tf.name,pattern = "::")
        }
        
        
        if (length(tf.feature[[1]]) > 1){
          print("Creating TF dimer metagene")
          
          expr.meta <-  matrix(, nrow=0, ncol = ncol(activity.mat) )
          for ( k in 1:length(tf.feature[[1]])){
            
            if ( tf.feature[[1]][k] %in% remove){
              
              if (tf.feature[[1]][k] %in% rownames(activity.mat)){
                print("TF not found in activity matrix!")
              }else{
                print("TF not found in activity matrix!")
              }
              
              
            }else{
              expression.vector <- activity.mat[tf.feature[[1]][k],]
              expr.meta <- rbind(expr.meta,expression.vector)
              
            }
            
          }
          if ( nrow(expr.meta) > 0){
            expression.meta.min <- colMins(expr.meta)
            expression.matrix <- rbind(expression.matrix,expression.meta.min)
            rownames(expression.matrix)[length(rownames(expression.matrix))] <- colnames(motif.mat)[i]
            print("TF dimer metagene made")
          }else{
            print("Skipping TF dimer")
          }
          
          
          
        }else{
          if ( tf.feature[[1]][1] %in% remove){
            
            if (tf.feature[[1]][1] %in% rownames(activity.mat)){
              print("TF not found in activity matrix!")
            }else{
              print("TF not found in activity matrix!")
            }
            
          }else{
            expression.vector <- activity.mat[tf.feature[[1]][1],]
            expression.matrix <- rbind(expression.matrix,expression.vector)
            rownames(expression.matrix)[length(rownames(expression.matrix))] <- colnames(motif.mat)[i]
          }
          
        }
        
      }
      expr.mat <- expression.matrix
      
      motif.mat <- motif.mat[,colnames(motif.mat) %in% rownames(expr.mat)]
      
      
      length(which(colnames(motif.mat)==rownames(expr.mat)))
      
      
      
      dir.create("./TFSEE_Outputs")
      print("Writing intermediate outputs to TFSEE_Outputs...")
      saveRDS(peak.mat,"./TFSEE_Outputs/peak_matrix.rds")
      
      saveRDS(expr.mat,"./TFSEE_Outputs/expression_matrix.rds")
      
      saveRDS(motif.mat,"./TFSEE_Outputs/motif_matrix.rds")
      ###########################################################
      
      peak.mat <- readRDS("./TFSEE_Outputs/peak_matrix.rds")
      expression.mat <- readRDS("./TFSEE_Outputs/expression_matrix.rds")
      motif.mat <- readRDS("./TFSEE_Outputs/motif_matrix.rds")
      
      
      dim(motif.mat)
      dim(peak.mat)
      dim(expression.mat)
      
      
      peak.mat <- t(peak.mat)
      expression.mat <- t(expression.mat)
      
      
      head(motif.mat[,1:3])
      head(peak.mat[,1:3])
      head(expression.mat[,1:3])
      
      print("Rescaling peak accessibility matrix...")
      peak.mat <- apply(peak.mat, 2, rescale,to=c(0,1))
      print("Rescaling TF expression matrix...")
      expression.mat <- apply(expression.mat, 2, rescale,to=c(0,1))
      
      
      dim(motif.mat)
      dim(peak.mat)
      dim(expression.mat)
      
      motif.mat <- t(motif.mat)
      peak.mat <- t(peak.mat)
      expression.mat <- t(expression.mat)
      
      dim(motif.mat)
      dim(peak.mat)
      dim(expression.mat)
      
      
      length(which(rownames(peak.mat)==colnames(motif.mat)))
      
      print("Generating raw motif activity matrix (this may take a few minutes)")
      int <- motif.mat%*%peak.mat
      
      dim(int)
      dim(expression.mat)
      
      print("Saving raw motif activity matrix")
      saveRDS(int,"./TFSEE_Outputs/raw_motif_activity_intermediate_matrix.rds")
      
      length(which(rownames(expression.mat)==rownames(int)))
      length(which(colnames(expression.mat)==colnames(int)))
      
      
      expression.mat <- expression.mat[order(match(rownames(expression.mat), rownames(int))), ]
      
      dim(int)
      dim(expression.mat)
      int <- int[order(match(rownames(int), rownames(expression.mat))), ]
      
      length(which(rownames(expression.mat)==rownames(int)))
      length(which(colnames(expression.mat)==colnames(int)))
      
      print("Generating final TFSEE activity matrix")
      res <- int*expression.mat
      
      length(rownames(res))
      length(motif.id)
      motif.idx <- names(motif.id)
      motif.idx <- motif.idx[motif.idx %in% rownames(res)]
      length(rownames(res))
      length(motif.idx)
      length(which(motif.idx==rownames(res)))
      motif.idx <- motif.idx[order(match(motif.idx, rownames(res)))]
      length(which(motif.idx==rownames(res)))
      
      
      motif.id <- motif.id[motif.idx]
      length(which(names(motif.id)==rownames(res)))
      
      rownames(res) <- motif.id
      
      print("Saving final TFSEE activity matrix")
      saveRDS(res,"./TFSEE_Outputs/TFSEE_final.rds")
      print("TFSEE finished successfully!")
      return(res)
    }
    
  }else{
    
    if(return.names == TRUE){
      # Get motif matches for example motifs in peaks
      
      motif_ix <- matchMotifs(motif.list,GRanges(rownames(peak.mat)), genome = genome,out = "matches",p.cutoff = motif.p.val)
      print("Finding motifs in peaks...")
      motif.mat <- motifMatches(motif_ix)
      head(motif.mat[,1:5])
      
      colnames(motif.mat) <- motif_ix$name
      rownames(motif.mat) <- rownames(peak.mat)
      
      length(colnames(motif.mat))
      length(unique(colnames(motif.mat)))
      
      tf.names <- gsub("\\s*\\([^\\)]+\\)","",as.character(colnames(motif.mat)))
      tf.names <- str_split(tf.names,pattern = "::")
      tf.names <- unlist(tf.names)
      
      remove <- base::setdiff(tf.names,rownames(expr.mat))
      tf.names <- base::setdiff(tf.names,remove)
      
      # All
      length(tf.names)
      length(intersect(tf.names,rownames(expr.mat)))
      
      dim(motif.mat)
      
      expression.matrix <- matrix(, nrow=0, ncol = ncol(expr.mat) )
      for ( i in 1:length(colnames(motif.mat))){
        
        if (as.character(colnames(motif.mat)[i]) == "EWSR1-FLI1"){
          
          tf.feature <- str_split(tf.name,pattern = "-")
        }else{
          tf.name <- gsub("\\s*\\([^\\)]+\\)","",as.character(colnames(motif.mat)[i]))
          tf.feature <- str_split(tf.name,pattern = "::")
        }
        
        
        if (length(tf.feature[[1]]) > 1){
          print("Creating TF dimer metagene")
          
          expr.meta <-  matrix(, nrow=0, ncol = ncol(expr.mat) )
          for ( k in 1:length(tf.feature[[1]])){
            
            if ( tf.feature[[1]][k] %in% remove){
              
              if (tf.feature[[1]][k] %in% rownames(activity.mat)){
                expression.vector <- activity.mat[tf.feature[[1]][k],]
                expr.meta <- rbind(expr.meta,expression.vector)
              }else{
                print("TF not found in expression matrix or activity matrix!")
              }
              
              
            }else{
              expression.vector <- expr.mat[tf.feature[[1]][k],]
              expr.meta <- rbind(expr.meta,expression.vector)
            }
            
          }
          if ( nrow(expr.meta) > 0){
            expression.meta.min <- colMins(as.matrix(expr.meta))
            expression.matrix <- rbind(expression.matrix,expression.meta.min)
            rownames(expression.matrix)[length(rownames(expression.matrix))] <- colnames(motif.mat)[i]
            print("TF dimer metagene made")
          }else{
            print("Both TFs not found in expression matrix or activity matrix!")
          }
          
          
        }else{
          if ( tf.feature[[1]][1] %in% remove){
            
            if (tf.feature[[1]][1] %in% rownames(activity.mat)){
              expression.vector <- activity.mat[tf.feature[[1]][1],]
              expression.matrix <- rbind(expression.matrix,expression.vector)
              rownames(expression.matrix)[length(rownames(expression.matrix))] <- colnames(motif.mat)[i]
              
            }else{
              print("TF not found in expression matrix or activity matrix!")
            }
            
          }else{
            expression.vector <- expr.mat[tf.feature[[1]][1],]
            expression.matrix <- rbind(expression.matrix,expression.vector)
            rownames(expression.matrix)[length(rownames(expression.matrix))] <- colnames(motif.mat)[i]
          }
          
        }
        
      }
      expr.mat <- expression.matrix
      
      motif.mat <- motif.mat[,colnames(motif.mat) %in% rownames(expr.mat)]
      
      
      length(which(colnames(motif.mat)==rownames(expr.mat)))
      
      
      dir.create("./TFSEE_Outputs")
      print("Writing intermediate outputs to TFSEE_Outputs...")
      saveRDS(peak.mat,"./TFSEE_Outputs/peak_matrix-prescaling.rds")
      
      saveRDS(expr.mat,"./TFSEE_Outputs/expression_matrix-prescaling.rds")
      
      saveRDS(motif.mat,"./TFSEE_Outputs/motif_matrix-prescaling.rds")
      ###########################################################
      
      
      
      peak.mat <- readRDS("./TFSEE_Outputs/peak_matrix-prescaling.rds")
      expression.mat <- readRDS("./TFSEE_Outputs/expression_matrix-prescaling.rds")
      motif.mat <- readRDS("./TFSEE_Outputs/motif_matrix-prescaling.rds")
      
      
      dim(motif.mat)
      dim(peak.mat)
      dim(expression.mat)
      
      
      peak.mat <- t(as.matrix(peak.mat))
      expression.mat <- t(as.matrix(expression.mat))
      motif.mat <- as.matrix(motif.mat)
      
      head(motif.mat[,1:3])
      head(peak.mat[,1:3])
      head(expression.mat[,1:3])
      
      print("Rescaling peak accessibility matrix...")
      peak.mat <- apply(peak.mat, 2, rescale,to=c(0,1))
      print("Rescaling TF expression matrix...")
      expression.mat <- apply(expression.mat, 2, rescale,to=c(0,1))
      
      
      dim(motif.mat)
      dim(peak.mat)
      dim(expression.mat)
      
      motif.mat <- t(motif.mat)
      peak.mat <- t(peak.mat)
      expression.mat <- t(expression.mat)
      
      dim(motif.mat)
      dim(peak.mat)
      dim(expression.mat)
      
      
      length(which(rownames(peak.mat)==colnames(motif.mat)))
      
      print("Generating raw motif activity matrix (this may take a few minutes)")
      int <- motif.mat%*%peak.mat
      
      dim(int)
      dim(expression.mat)
      
      print("Saving raw motif activity matrix")
      saveRDS(int,"./TFSEE_Outputs/raw_motif_activity_intermediate_matrix.rds")
      
      length(which(rownames(expression.mat)==rownames(int)))
      length(which(colnames(expression.mat)==colnames(int)))
      
      
      expression.mat <- expression.mat[order(match(rownames(expression.mat), rownames(int))), ]
      
      dim(int)
      dim(expression.mat)
      int <- int[order(match(rownames(int), rownames(expression.mat))), ]
      
      length(which(rownames(expression.mat)==rownames(int)))
      length(which(colnames(expression.mat)==colnames(int)))
      
      print("Generating final TFSEE activity matrix")
      res <- int*expression.mat
      
      print("Saving final TFSEE activity matrix")
      saveRDS(res,"./TFSEE_Outputs/TFSEE_final.rds")
      
      print("Saving post-scaling intermediates")
      saveRDS(peak.mat,"./TFSEE_Outputs/peak_matrix-postscaling.rds")
      
      saveRDS(expression.mat,"./TFSEE_Outputs/expression_matrix-postscaling.rds")
      
      saveRDS(motif.mat,"./TFSEE_Outputs/motif_matrix-postscaling.rds")
      
      print("TFSEE finished successfully!")
      return(res)
    }else{
      # Get motif matches for example motifs in peaks
      
      motif_ix <- matchMotifs(motif.list,GRanges(rownames(peak.mat)), genome = genome,out = "matches",p.cutoff = motif.p.val)
      print("Finding motifs in peaks...")
      motif.mat <- motifMatches(motif_ix)
      head(motif.mat[,1:5])
      
      motif.id <- colnames(motif.mat)
      names(motif.id) <- motif_ix$name
      colnames(motif.mat) <- motif_ix$name
      rownames(motif.mat) <- rownames(peak.mat)
      
      length(colnames(motif.mat))
      length(unique(colnames(motif.mat)))
      
      tf.names <- gsub("\\s*\\([^\\)]+\\)","",as.character(colnames(motif.mat)))
      tf.names <- str_split(tf.names,pattern = "::")
      tf.names <- unlist(tf.names)
      
      remove <- base::setdiff(tf.names,rownames(expr.mat))
      tf.names <- base::setdiff(tf.names,remove)
      
      # All
      length(tf.names)
      length(intersect(tf.names,rownames(expr.mat)))
      
      dim(motif.mat)
      
      expression.matrix <- matrix(, nrow=0, ncol = ncol(expr.mat) )
      for ( i in 1:length(colnames(motif.mat))){
        
        if (as.character(colnames(motif.mat)[i]) == "EWSR1-FLI1"){
          
          tf.feature <- str_split(tf.name,pattern = "-")
        }else{
          tf.name <- gsub("\\s*\\([^\\)]+\\)","",as.character(colnames(motif.mat)[i]))
          tf.feature <- str_split(tf.name,pattern = "::")
        }
        
        
        if (length(tf.feature[[1]]) > 1){
          print("Creating TF dimer metagene")
          
          expr.meta <-  matrix(, nrow=0, ncol = ncol(expr.mat) )
          for ( k in 1:length(tf.feature[[1]])){
            
            if ( tf.feature[[1]][k] %in% remove){
              
              if (tf.feature[[1]][k] %in% rownames(activity.mat)){
                expression.vector <- activity.mat[tf.feature[[1]][k],]
                expr.meta <- rbind(expr.meta,expression.vector)
              }else{
                print("TF not found in expression matrix or activity matrix!")
              }
              
              
            }else{
              expression.vector <- expr.mat[tf.feature[[1]][k],]
              expr.meta <- rbind(expr.meta,expression.vector)
            }
            
          }
          if ( nrow(expr.meta) > 0){
            expression.meta.min <- colMins(as.matrix(expr.meta))
            expression.matrix <- rbind(expression.matrix,expression.meta.min)
            rownames(expression.matrix)[length(rownames(expression.matrix))] <- colnames(motif.mat)[i]
            print("TF dimer metagene made")
          }else{
            print("Both TFs not found in expression matrix or activity matrix!")
          }
          
          
        }else{
          if ( tf.feature[[1]][1] %in% remove){
            
            if (tf.feature[[1]][1] %in% rownames(activity.mat)){
              expression.vector <- activity.mat[tf.feature[[1]][1],]
              expression.matrix <- rbind(expression.matrix,expression.vector)
              rownames(expression.matrix)[length(rownames(expression.matrix))] <- colnames(motif.mat)[i]
              
            }else{
              print("TF not found in expression matrix or activity matrix!")
            }
            
          }else{
            expression.vector <- expr.mat[tf.feature[[1]][1],]
            expression.matrix <- rbind(expression.matrix,expression.vector)
            rownames(expression.matrix)[length(rownames(expression.matrix))] <- colnames(motif.mat)[i]
          }
          
        }
        
      }
      expr.mat <- expression.matrix
      
      motif.mat <- motif.mat[,colnames(motif.mat) %in% rownames(expr.mat)]
      
      
      length(which(colnames(motif.mat)==rownames(expr.mat)))
      
      
      dir.create("./TFSEE_Outputs")
      print("Writing intermediate outputs to TFSEE_Outputs...")
      saveRDS(peak.mat,"./TFSEE_Outputs/peak_matrix.rds")
      
      saveRDS(expr.mat,"./TFSEE_Outputs/expression_matrix.rds")
      
      saveRDS(motif.mat,"./TFSEE_Outputs/motif_matrix.rds")
      ###########################################################
      
      
      
      peak.mat <- readRDS("./TFSEE_Outputs/peak_matrix.rds")
      expression.mat <- readRDS("./TFSEE_Outputs/expression_matrix.rds")
      motif.mat <- readRDS("./TFSEE_Outputs/motif_matrix.rds")
      
      
      dim(motif.mat)
      dim(peak.mat)
      dim(expression.mat)
      
      
      peak.mat <- t(as.matrix(peak.mat))
      expression.mat <- t(as.matrix(expression.mat))
      motif.mat <- as.matrix(motif.mat)
      
      head(motif.mat[,1:3])
      head(peak.mat[,1:3])
      head(expression.mat[,1:3])
      
      print("Rescaling peak accessibility matrix...")
      peak.mat <- apply(peak.mat, 2, rescale,to=c(0,1))
      print("Rescaling TF expression matrix...")
      expression.mat <- apply(expression.mat, 2, rescale,to=c(0,1))
      
      
      dim(motif.mat)
      dim(peak.mat)
      dim(expression.mat)
      
      motif.mat <- t(motif.mat)
      peak.mat <- t(peak.mat)
      expression.mat <- t(expression.mat)
      
      dim(motif.mat)
      dim(peak.mat)
      dim(expression.mat)
      
      
      length(which(rownames(peak.mat)==colnames(motif.mat)))
      
      print("Generating raw motif activity matrix (this may take a few minutes)")
      int <- motif.mat%*%peak.mat
      
      dim(int)
      dim(expression.mat)
      
      print("Saving raw motif activity matrix")
      saveRDS(int,"./TFSEE_Outputs/raw_motif_activity_intermediate_matrix.rds")
      
      length(which(rownames(expression.mat)==rownames(int)))
      length(which(colnames(expression.mat)==colnames(int)))
      
      
      expression.mat <- expression.mat[order(match(rownames(expression.mat), rownames(int))), ]
      
      dim(int)
      dim(expression.mat)
      int <- int[order(match(rownames(int), rownames(expression.mat))), ]
      
      length(which(rownames(expression.mat)==rownames(int)))
      length(which(colnames(expression.mat)==colnames(int)))
      
      print("Generating final TFSEE activity matrix")
      res <- int*expression.mat
      
      length(rownames(res))
      length(motif.id)
      motif.idx <- names(motif.id)
      motif.idx <- motif.idx[motif.idx %in% rownames(res)]
      length(rownames(res))
      length(motif.idx)
      length(which(motif.idx==rownames(res)))
      motif.idx <- motif.idx[order(match(motif.idx, rownames(res)))]
      length(which(motif.idx==rownames(res)))
      
      
      motif.id <- motif.id[motif.idx]
      length(which(names(motif.id)==rownames(res)))
      
      rownames(res) <- motif.id
      
      print("Saving final TFSEE activity matrix")
      saveRDS(res,"./TFSEE_Outputs/TFSEE_final.rds")
      print("TFSEE finished successfully!")
      return(res)
    }
    
  }
  
}
