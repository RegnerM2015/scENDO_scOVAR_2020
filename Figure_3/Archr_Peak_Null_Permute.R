##########################################################################################
# Peak2Gene Links Methods
##########################################################################################

#' Add Peak2GeneLinks to an ArchRProject
#' 
#' This function will add peak-to-gene links to a given ArchRProject
#' 
#' @param ArchRProj An `ArchRProject` object.
#' @param reducedDims The name of the `reducedDims` object (i.e. "IterativeLSI") to retrieve from the designated `ArchRProject`.
#' @param dimsToUse A vector containing the dimensions from the `reducedDims` object to use in clustering.
#' @param scaleDims A boolean value that indicates whether to z-score the reduced dimensions for each cell. This is useful for minimizing
#' the contribution of strong biases (dominating early PCs) and lowly abundant populations. However, this may lead to stronger sample-specific
#' biases since it is over-weighting latent PCs. If set to `NULL` this will scale the dimensions based on the value of `scaleDims` when the
#' `reducedDims` were originally created during dimensionality reduction. This idea was introduced by Timothy Stuart.
#' @param corCutOff A numeric cutoff for the correlation of each dimension to the sequencing depth. If the dimension has a
#' correlation to sequencing depth that is greater than the `corCutOff`, it will be excluded from analysis.
#' @param k The number of k-nearest neighbors to use for creating single-cell groups for correlation analyses.
#' @param knnIteration The number of k-nearest neighbor groupings to test for passing the supplied `overlapCutoff`.
#' @param overlapCutoff The maximum allowable overlap between the current group and all previous groups to permit the current
#' group be added to the group list during k-nearest neighbor calculations.
#' @param maxDist The maximum allowable distance in basepairs between two peaks to consider for co-accessibility.
#' @param scaleTo The total insertion counts from the designated group of single cells is summed across all relevant peak regions
#' from the `peakSet` of the `ArchRProject` and normalized to the total depth provided by `scaleTo`.
#' @param log2Norm A boolean value indicating whether to log2 transform the single-cell groups prior to computing co-accessibility correlations.
#' @param predictionCutoff A numeric describing the cutoff for RNA integration to use when picking cells for groupings.
#' @param seed A number to be used as the seed for random number generation required in knn determination. It is recommended
#' to keep track of the seed used so that you can reproduce results downstream.
#' @param threads The number of threads to be used for parallel computing.
#' @param verbose A boolean value that determines whether standard output should be printed.
#' @param logFile The path to a file to be used for logging ArchR output.
#' @export
addPermPeak2GeneLinks <- function(
  ArchRProj = NULL,
  reducedDims = "IterativeLSI",
  useMatrix = "GeneIntegrationMatrix",
  dimsToUse = 1:30,
  scaleDims = NULL,
  corCutOff = 0.75,
  cellsToUse = NULL,
  k = 100, 
  knnIteration = 500, 
  overlapCutoff = 0.8, 
  maxDist = 250000,
  scaleTo = 10^4,
  log2Norm = TRUE,
  predictionCutoff = 0.4,
  seed = 1, 
  threads = max(floor(getArchRThreads() / 2), 1),
  verbose = TRUE,
  logFile = createLogFile("addPeak2GeneLinks")
){
  
  .validInput(input = ArchRProj, name = "ArchRProj", valid = c("ArchRProj"))
  .validInput(input = reducedDims, name = "reducedDims", valid = c("character"))
  .validInput(input = dimsToUse, name = "dimsToUse", valid = c("numeric", "null"))
  .validInput(input = scaleDims, name = "scaleDims", valid = c("boolean", "null"))
  .validInput(input = corCutOff, name = "corCutOff", valid = c("numeric", "null"))
  .validInput(input = cellsToUse, name = "cellsToUse", valid = c("character", "null"))
  .validInput(input = k, name = "k", valid = c("integer"))
  .validInput(input = knnIteration, name = "knnIteration", valid = c("integer"))
  .validInput(input = overlapCutoff, name = "overlapCutoff", valid = c("numeric"))
  .validInput(input = maxDist, name = "maxDist", valid = c("integer"))
  .validInput(input = scaleTo, name = "scaleTo", valid = c("numeric"))
  .validInput(input = log2Norm, name = "log2Norm", valid = c("boolean"))
  .validInput(input = threads, name = "threads", valid = c("integer"))
  .validInput(input = verbose, name = "verbose", valid = c("boolean"))
  .validInput(input = logFile, name = "logFile", valid = c("character"))
  
  tstart <- Sys.time()
  .startLogging(logFile = logFile)
  .logThis(mget(names(formals()),sys.frame(sys.nframe())), "addPeak2GeneLinks Input-Parameters", logFile = logFile)
  
  .logDiffTime(main="Getting Available Matrices", t1=tstart, verbose=verbose, logFile=logFile)
  AvailableMatrices <- getAvailableMatrices(ArchRProj)
  
  if("PeakMatrix" %ni% AvailableMatrices){
    stop("PeakMatrix not in AvailableMatrices")
  }
  
  if(useMatrix %ni% AvailableMatrices){
    stop(paste0(useMatrix, " not in AvailableMatrices"))
  }
  
  ArrowFiles <- getArrowFiles(ArchRProj)
  
  tstart <- Sys.time()
  
  dfAll <- .safelapply(seq_along(ArrowFiles), function(x){
    DataFrame(
      cellNames = paste0(names(ArrowFiles)[x], "#", h5read(ArrowFiles[x], paste0(useMatrix, "/Info/CellNames"))),
      predictionScore = h5read(ArrowFiles[x], paste0(useMatrix, "/Info/predictionScore"))
    )
  }, threads = threads) %>% Reduce("rbind", .)
  
  .logDiffTime(
    sprintf("Filtered Low Prediction Score Cells (%s of %s, %s)", 
            sum(dfAll[,2] < predictionCutoff), 
            nrow(dfAll), 
            round(sum(dfAll[,2] < predictionCutoff) / nrow(dfAll), 3)
    ), t1=tstart, verbose=verbose, logFile=logFile)
  
  keep <- sum(dfAll[,2] >= predictionCutoff) / nrow(dfAll)
  dfAll <- dfAll[which(dfAll[,2] > predictionCutoff),]
  
  set.seed(seed)
  
  #Get Peak Set
  peakSet <- getPeakSet(ArchRProj)
  .logThis(peakSet, "peakSet", logFile = logFile)
  
  #Gene Info
  geneSet <- .getFeatureDF(ArrowFiles, useMatrix, threads = threads)
  geneStart <- GRanges(geneSet$seqnames, IRanges(geneSet$start, width = 1), name = geneSet$name, idx = geneSet$idx)
  .logThis(geneStart, "geneStart", logFile = logFile)
  
  #Get Reduced Dims
  rD <- getReducedDims(ArchRProj, reducedDims = reducedDims, corCutOff = corCutOff, dimsToUse = dimsToUse)
  if(!is.null(cellsToUse)){
    rD <- rD[cellsToUse, ,drop=FALSE]
  }
  #Subsample
  idx <- sample(seq_len(nrow(rD)), knnIteration, replace = !nrow(rD) >= knnIteration)
  
  #KNN Matrix
  .logDiffTime(main="Computing KNN", t1=tstart, verbose=verbose, logFile=logFile)
  knnObj <- .computeKNN(data = rD, query = rD[idx,], k = k)
  
  #Determin Overlap
  .logDiffTime(main="Identifying Non-Overlapping KNN pairs", t1=tstart, verbose=verbose, logFile=logFile)
  keepKnn <- determineOverlapCpp(knnObj, floor(overlapCutoff * k))
  
  #Keep Above Cutoff
  knnObj <- knnObj[keepKnn==0,]
  .logDiffTime(paste0("Identified ", nrow(knnObj), " Groupings!"), t1=tstart, verbose=verbose, logFile=logFile)
  
  #Convert To Names List
  knnObj <- lapply(seq_len(nrow(knnObj)), function(x){
    rownames(rD)[knnObj[x, ]]
  }) %>% SimpleList
  
  #Check Chromosomes
  chri <- gtools::mixedsort(unique(paste0(seqnames(peakSet))))
  chrj <- gtools::mixedsort(unique(paste0(seqnames(geneStart))))
  chrij <- intersect(chri, chrj)
  
  #Features
  geneDF <- mcols(geneStart)
  peakDF <- mcols(peakSet)
  geneDF$seqnames <- seqnames(geneStart)
  peakDF$seqnames <- seqnames(peakSet)
  
  #Group Matrix RNA
  .logDiffTime(main="Getting Group RNA Matrix", t1=tstart, verbose=verbose, logFile=logFile)
  groupMatRNA <- .getGroupMatrix(
    ArrowFiles = getArrowFiles(ArchRProj), 
    featureDF = geneDF, 
    groupList = knnObj, 
    useMatrix = useMatrix,
    threads = threads,
    verbose = FALSE
  )
  .logThis(groupMatRNA, "groupMatRNA", logFile = logFile)
  
  #Group Matrix ATAC
  .logDiffTime(main="Getting Group ATAC Matrix", t1=tstart, verbose=verbose, logFile=logFile)
  groupMatATAC <- .getGroupMatrix(
    ArrowFiles = getArrowFiles(ArchRProj), 
    featureDF = peakDF, 
    groupList = knnObj, 
    useMatrix = "PeakMatrix",
    threads = threads,
    verbose = FALSE
  )
  .logThis(groupMatATAC, "groupMatATAC", logFile = logFile)
  
  .logDiffTime(main="Normalizing Group Matrices", t1=tstart, verbose=verbose, logFile=logFile)
  
  groupMatRNA <- t(t(groupMatRNA) / colSums(groupMatRNA)) * scaleTo

  groupMatATAC <- t(t(groupMatATAC) / colSums(groupMatATAC)) * scaleTo
  
  if(log2Norm){
    groupMatRNA  <- log2(groupMatRNA + 1)
    groupMatATAC <- log2(groupMatATAC + 1)    
  }
  
  names(geneStart) <- NULL
  
  seRNA <- SummarizedExperiment(
    assays = SimpleList(RNA = groupMatRNA), 
    rowRanges = geneStart
  )
  metadata(seRNA)$KNNList <- knnObj
  .logThis(seRNA, "seRNA", logFile = logFile)
  names(peakSet) <- NULL
  
  seATAC <- SummarizedExperiment(
    assays = SimpleList(ATAC = groupMatATAC), 
    rowRanges = peakSet
  )
  metadata(seATAC)$KNNList <- knnObj
  .logThis(seATAC, "seATAC", logFile = logFile)
  
  rm(groupMatRNA, groupMatATAC)
  gc()
  
  #Overlaps
  .logDiffTime(main="Finding Peak Gene Pairings", t1=tstart, verbose=verbose, logFile=logFile)
  o <- DataFrame(
    findOverlaps(
      .suppressAll(resize(seRNA, 2 * maxDist + 1, "center")), 
      resize(rowRanges(seATAC), 1, "center"), 
      ignore.strand = TRUE
    )
  )
  colnames(seRNA@assays@data[[1]])
  #Get Distance from Fixed point A B 
  o$distance <- distance(rowRanges(seRNA)[o[,1]] , rowRanges(seATAC)[o[,2]] )
  colnames(o) <- c("B", "A", "distance")

  # #Null Correlations
  # .logDiffTime(main="Computing Background Correlations", t1=tstart, verbose=verbose, logFile=logFile)
  # nullCor <- .getNullCorrelations(seATAC, seRNA, o, 1000)
  # saveRDS(nullCor,"nullCor.rds")
  # .logDiffTime(main="Computing Correlations", t1=tstart, verbose=verbose, logFile=logFile)
  # 
  seATAC.mat <- assay(seATAC)
  seATAC.mat <- seATAC.mat[,sample(ncol(seATAC.mat))]#Permute ATAC aggregates
  
  o$Correlation <- rowCorCpp(as.integer(o$A), as.integer(o$B), seATAC.mat,assay(seRNA))
  o$VarAssayA <- .getQuantiles(matrixStats::rowVars(seATAC.mat))[o$A]
  o$VarAssayB <- .getQuantiles(matrixStats::rowVars(assay(seRNA)))[o$B]
  o$TStat <- (o$Correlation / sqrt((1-o$Correlation^2)/(ncol(seATAC.mat)-2))) #T-statistic P-value
  o$Pval <- 2*pt(-abs(o$TStat), ncol(seATAC.mat) - 2)
  o$FDR <- p.adjust(o$Pval, method = "fdr")
  # o$EmpPval <- 2*pnorm(-abs(((o$Correlation - mean(nullCor[[2]])) / sd(nullCor[[2]]))))
  # o$EmpFDR <- p.adjust(o$EmpPval, method = "fdr")
  
  out <- o[, c("A", "B", "Correlation","Pval", "FDR", "VarAssayA", "VarAssayB")]
  colnames(out) <- c("idxATAC", "idxRNA", "Correlation", "RawPVal","FDR", "VarQATAC", "VarQRNA")
  mcols(peakSet) <- NULL
  names(peakSet) <- NULL
  metadata(out)$peakSet <- peakSet
  metadata(out)$geneSet <- geneStart
  
  #Save Group Matrices
  dir.create(file.path(getOutputDirectory(ArchRProj), "Peak2GeneLinks"), showWarnings = FALSE)
  outATAC <- file.path(getOutputDirectory(ArchRProj), "Peak2GeneLinks", "seATAC-Group-KNN.rds")
  saveRDS(seATAC, outATAC, compress = FALSE)
  outRNA <- file.path(getOutputDirectory(ArchRProj), "Peak2GeneLinks", "seRNA-Group-KNN.rds")
  saveRDS(seRNA, outRNA, compress = FALSE)
  metadata(out)$seATAC <- outATAC
  metadata(out)$seRNA <- outRNA
  
  metadata(ArchRProj@peakSet)$Peak2GeneLinks <- out
  
  .logDiffTime(main="Completed Peak2Gene Correlations!", t1=tstart, verbose=verbose, logFile=logFile)
  .endLogging(logFile = logFile)
  
  ArchRProj
  
}

 .getNullCorrelations <- function(seA, seB, o, n){

   o$seq <- seqnames(seA)[o$A]

   nullCor <- lapply(seq_along(unique(o$seq)), function(i){

     #Get chr from olist
     chri <- unique(o$seq)[i]
     message(chri, " ", appendLF = FALSE)

     #Randomly get n seA
     id <- which(as.character(seqnames(seA)) != chri)
     if(length(id) > n){
       transAidx <- sample(id, n)
     }else{
       transAidx <- id
     }

     #Calculate Correlations
     grid <- expand.grid(transAidx, unique(o[o$seq==chri,]$B))

     idxA <- unique(grid[,1])
     idxB <- unique(grid[,2])

     seSubA <- seA[idxA]
     seSubB <- seB[idxB]

     grid[,3] <- match(grid[,1], idxA)
     grid[,4] <- match(grid[,2], idxB)

     colnames(grid) <- c("A", "B")
     out <- rowCorCpp(grid[,3], grid[,4], assay(seSubA), assay(seSubB))
     out <- na.omit(out)

     return(out)

   }) %>% SimpleList
   message("")

   summaryDF <- lapply(nullCor, function(x){
     data.frame(mean = mean(x), sd = sd(x), median = median(x), n = length(x))
   }) %>% Reduce("rbind",.)

   return(list(summaryDF, unlist(nullCor)))

 }

#' Get the peak-to-gene links from an ArchRProject
#' 
#' This function obtains peak-to-gene links from an ArchRProject.
#' 
#' @param ArchRProj An `ArchRProject` object.
#' @param corCutOff A numeric describing the minimum numeric peak-to-gene correlation to return.
#' @param FDRCutOff A numeric describing the maximum numeric peak-to-gene false discovery rate to return.
#' @param varCutOffATAC A numeric describing the minimum variance quantile of the ATAC peak accessibility when selecting links.
#' @param varCutOffRNA A numeric describing the minimum variance quantile of the RNA gene expression when selecting links.
#' @param resolution A numeric describing the bp resolution to return loops as. This helps with overplotting of correlated regions.
#' @param returnLoops A boolean indicating to return the peak-to-gene links as a `GRanges` "loops" object designed for use with
#' the `ArchRBrowser()` or as an `ArchRBrowserTrack()`.
#' @export
getPeak2GeneLinks <- function(
  ArchRProj = NULL, 
  corCutOff = 0.45, 
  FDRCutOff = 0.0001,
  varCutOffATAC = 0.25,
  varCutOffRNA = 0.25,
  resolution = 1, 
  returnLoops = TRUE
){
  
  .validInput(input = ArchRProj, name = "ArchRProj", valid = "ArchRProject")
  .validInput(input = corCutOff, name = "corCutOff", valid = "numeric")
  .validInput(input = FDRCutOff, name = "FDRCutOff", valid = "numeric")
  .validInput(input = varCutOffATAC, name = "varCutOffATAC", valid = "numeric")
  .validInput(input = varCutOffRNA, name = "varCutOffRNA", valid = "numeric")
  .validInput(input = resolution, name = "resolution", valid = c("integer", "null"))
  .validInput(input = returnLoops, name = "returnLoops", valid = "boolean")
  
  if(is.null(ArchRProj@peakSet)){
    return(NULL)
  }
  
  if(is.null(metadata(ArchRProj@peakSet)$Peak2GeneLinks)){
    
    return(NULL)
    
  }else{
    
    p2g <- metadata(ArchRProj@peakSet)$Peak2GeneLinks
    p2g <- p2g[which(p2g$Correlation >= corCutOff & p2g$FDR <= FDRCutOff), ,drop=FALSE]
    
    if(!is.null(varCutOffATAC)){
      p2g <- p2g[which(p2g$VarQATAC > varCutOffATAC),]
    }
    
    if(!is.null(varCutOffRNA)){
      p2g <- p2g[which(p2g$VarQRNA > varCutOffRNA),]
    }
    
    if(returnLoops){
      
      peakSummits <- resize(metadata(p2g)$peakSet, 1, "center")
      geneStarts <- resize(metadata(p2g)$geneSet, 1, "start")
      
      if(!is.null(resolution)){
        summitTiles <- floor(start(peakSummits) / resolution) * resolution + floor(resolution / 2)
        geneTiles <- floor(start(geneStarts) / resolution) * resolution + floor(resolution / 2)
      }else{
        summitTiles <- start(peakSummits)
        geneTiles <- start(geneTiles)
      }
      
      loops <- .constructGR(
        seqnames = seqnames(peakSummits[p2g$idxATAC]),
        start = summitTiles[p2g$idxATAC],
        end = geneTiles[p2g$idxRNA]
      )
      mcols(loops)$value <- p2g$Correlation
      mcols(loops)$FDR <- p2g$FDR
      
      loops <- loops[order(mcols(loops)$value, decreasing=TRUE)]
      loops <- unique(loops)
      loops <- loops[width(loops) > 0]
      loops <- sort(sortSeqlevels(loops))
      
      loops <- SimpleList(Peak2GeneLinks = loops)
      
      return(loops)
      
    }else{
      
      return(p2g)
      
    }
    
  }
  
}


# Generated by using Rcpp::compileAttributes() -> do not edit by hand
# Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

rowCorCpp <- function(idxX, idxY, X, Y) {
  .Call('_ArchR_rowCorCpp', PACKAGE = 'ArchR', idxX, idxY, X, Y)
}

rleSumsStrandedChr <- function(rle, x, strand, width) {
  .Call('_ArchR_rleSumsStrandedChr', PACKAGE = 'ArchR', rle, x, strand, width)
}

rleSumsStranded <- function(rleList, grList, width, as_integer) {
  .Call('_ArchR_rleSumsStranded', PACKAGE = 'ArchR', rleList, grList, width, as_integer)
}

tabulate2dCpp <- function(x, xmin, xmax, y, ymin, ymax) {
  .Call('_ArchR_tabulate2dCpp', PACKAGE = 'ArchR', x, xmin, xmax, y, ymin, ymax)
}

computeSparseRowVariances <- function(j, val, rm, n) {
  .Call('_ArchR_computeSparseRowVariances', PACKAGE = 'ArchR', j, val, rm, n)
}

determineOverlapCpp <- function(m, overlapCut) {
  .Call('_ArchR_determineOverlapCpp', PACKAGE = 'ArchR', m, overlapCut)
}

kmerIdxCpp <- function(str, window, n, kmer) {
  .Call('_ArchR_kmerIdxCpp', PACKAGE = 'ArchR', str, window, n, kmer)
}

kmerPositionFrequencyCpp <- function(string_vector, strand_vector, window, w, kmer) {
  .Call('_ArchR_kmerPositionFrequencyCpp', PACKAGE = 'ArchR', string_vector, strand_vector, window, w, kmer)
}

kmerIDFrequencyCpp <- function(string_vector, id_vector, n_id, window, w, kmer) {
  .Call('_ArchR_kmerIDFrequencyCpp', PACKAGE = 'ArchR', string_vector, id_vector, n_id, window, w, kmer)
}




##########################################################################################
# Validation Methods
##########################################################################################

.validInput <- function(input = NULL, name = NULL, valid = NULL){
  
  valid <- unique(valid)
  
  if(is.character(valid)){
    valid <- tolower(valid)
  }else{
    stop("Validator must be a character!")
  }
  
  if(!is.character(name)){
    stop("name must be a character!")
  }
  
  if("null" %in% tolower(valid)){
    valid <- c("null", valid[which(tolower(valid) != "null")])
  }
  
  av <- FALSE
  
  for(i in seq_along(valid)){
    
    vi <- valid[i]
    
    if(vi == "integer" | vi == "wholenumber"){
      
      if(all(is.numeric(input))){
        #https://stackoverflow.com/questions/3476782/check-if-the-number-is-integer
        cv <- min(abs(c(input%%1, input%%1-1))) < .Machine$double.eps^0.5
      }else{
        cv <- FALSE
      }
      
    }else if(vi == "null"){
      
      cv <- is.null(input)
      
    }else if(vi == "bool" | vi == "boolean" | vi == "logical"){
      
      cv <- is.logical(input)
      
    }else if(vi == "numeric"){
      
      cv <- is.numeric(input)
      
    }else if(vi == "vector"){
      
      cv <- is.vector(input)
      
    }else if(vi == "matrix"){
      
      cv <- is.matrix(input)
      
    }else if(vi == "sparsematrix"){
      
      cv <- is(input, "dgCMatrix")
      
    }else if(vi == "character"){
      
      cv <- is.character(input)
      
    }else if(vi == "factor"){
      
      cv <- is.factor(input)
      
    }else if(vi == "rlecharacter"){
      
      cv1 <- is(input, "Rle")
      if(cv1){
        cv <- is(input@values, "factor") || is(input@values, "character")
      }else{
        cv <- FALSE
      }
      
    }else if(vi == "palette"){
      
      cv <- all(.isColor(input))
      
    }else if(vi == "timestamp"){
      
      cv <- is(input, "POSIXct")
      
    }else if(vi == "dataframe" | vi == "data.frame" | vi == "df"){
      
      cv1 <- is.data.frame(input)
      cv2 <- is(input, "DataFrame")
      cv <- any(cv1, cv2)
      
    }else if(vi == "fileexists"){
      
      cv <- all(file.exists(input))
      
    }else if(vi == "direxists"){
      
      cv <- all(dir.exists(input))
      
    }else if(vi == "granges" | vi == "gr"){
      
      cv <- is(input, "GRanges")
      
    }else if(vi == "grangeslist" | vi == "grlist"){
      
      cv <- .isGRList(input)
      
    }else if(vi == "list" | vi == "simplelist"){
      
      cv1 <- is.list(input)
      cv2 <- is(input, "SimpleList")
      cv <- any(cv1, cv2)
      
    }else if(vi == "bsgenome"){
      
      cv1 <- is(input, "BSgenome")
      cv2 <- tryCatch({
        library(input)
        eval(parse(text=input))
      }, error = function(e){
        FALSE
      })
      cv <- any(cv1, cv2)
      
    }else if(vi == "se" | vi == "summarizedexperiment"){
      
      cv <- is(input, "SummarizedExperiment")
      
    }else if(vi == "seurat" | vi == "seuratobject"){
      
      cv <- is(input, "Seurat")
      
    }else if(vi == "txdb"){
      
      cv <- is(input, "TxDb")
      
    }else if(vi == "orgdb"){
      
      cv <- is(input, "OrgDb")
      
    }else if(vi == "bsgenome"){
      
      cv <- is(input, "BSgenome")
      
    }else if(vi == "parallelparam"){
      
      cv <- is(input, "BatchtoolsParam")
      
    }else if(vi == "archrproj" | vi == "archrproject"){
      
      cv <- is(input, "ArchRProject")
      ###validObject(input) check this doesnt break anything if we
      ###add it. Useful to make sure all ArrowFiles exist! QQQ
      
    }else{
      
      stop("Validator is not currently supported by ArchR!")
      
    }
    
    if(cv){
      av <- TRUE
      break
    }   
    
  }
  
  if(av){
    
    return(invisible(TRUE))
    
  }else{
    
    stop("Input value for '", name,"' is not a ", paste(valid, collapse="," ), ", (",name," = ",class(input),") please supply valid input!")
    
  }
  
}

#https://stackoverflow.com/questions/3476782/check-if-the-number-is-integer
.isWholenumber <- function(x, tol = .Machine$double.eps^0.5){
  abs(x - round(x)) < tol
}

#https://stackoverflow.com/questions/13289009/check-if-character-string-is-a-valid-color-representation
.isColor <- function(x = NULL){
  unlist(lapply(x, function(y) tryCatch(is.matrix(col2rgb(y)), error = function(e) FALSE)))
}

.isDiscrete <- function(x = NULL){
  is.factor(x) || is.character(x) || is.logical(x)
}

.isGRList <- function(x){
  isList <- grepl("list", class(x), ignore.case=TRUE)
  if(!isList){
    FALSE
  }else{
    allGR <- all(unlist(lapply(x, function(x) is(x, "GRanges") )))
    if(allGR){
      TRUE
    }else{
      FALSE
    }
  }
}

#' Get/Validate BSgenome
#' 
#' This function will attempt to get or validate an input as a BSgenome.
#' 
#' @param genome This option must be one of the following: (i) the name of a valid ArchR-supported genome ("hg38", "hg19", or "mm10"),
#' (ii) the name of a `BSgenome` package (for ex. "BSgenome.Hsapiens.UCSC.hg19"), or (iii) a `BSgenome` object.
#' @param masked A boolean describing whether or not to access the masked version of the selected genome. See `BSgenome::getBSgenome()`.
#' @export
validBSgenome <- function(genome = NULL, masked = FALSE){
  
  .validInput(input = genome, name = "genome", valid = c("character", "bsgenome"))
  .validInput(input = masked, name = "masked", valid = c("boolean"))
  
  stopifnot(!is.null(genome))
  if(inherits(genome, "BSgenome")){
    return(genome)
  }else if(is.character(genome)){
    genome <- tryCatch({
      .requirePackage(genome)
      bsg <- eval(parse(text = genome))
      if(inherits(bsg, "BSgenome")){
        return(bsg)
      }else{
        stop("genome is not a BSgenome valid class!")
      }
    }, error = function(x){
      BSgenome::getBSgenome(genome, masked = masked)
    })  
    return(genome)
  }else{
    stop("Cannot validate BSgenome options are a valid BSgenome or character for getBSgenome")
  }  
}

.validTxDb <- function(TxDb = NULL){
  stopifnot(!is.null(TxDb))
  if(inherits(TxDb, "TxDb")){
    return(TxDb)
  }else if(is.character(TxDb)){
    return(getTxDb(TxDb)) #change
  }else{
    stop("Cannot validate TxDb options are a valid TxDb or character for getTxDb")
  }
}

.validOrgDb <- function(OrgDb = NULL){
  stopifnot(!is.null(OrgDb))
  if(inherits(OrgDb, "OrgDb")){
    return(OrgDb)
  }else if(is.character(OrgDb)){
    return(getOrgDb(OrgDb)) #change
  }else{
    stop("Cannot validate OrgDb options are a valid OrgDb or character for getOrgDb")
  }
}

.validGRanges <- function(gr = NULL){
  stopifnot(!is.null(gr))
  if(inherits(gr, "GRanges")){
    return(gr)
  }else{
    stop("Error cannot validate genomic range!")
  }
}

.validGeneAnnotation <- function(geneAnnotation = NULL){
  
  if(!inherits(geneAnnotation, "SimpleList")){
    if(inherits(geneAnnotation, "list")){
      geneAnnotation <- as(geneAnnotation, "SimpleList")
    }else{
      stop("geneAnnotation must be a list/SimpleList of 3 GRanges for : Genes GRanges, Exons GRanges and TSS GRanges!")
    }
  }
  if(identical(sort(tolower(names(geneAnnotation))), c("exons", "genes", "tss"))){
    
    gA <- SimpleList()
    gA$genes <- .validGRanges(geneAnnotation[[grep("genes", names(geneAnnotation), ignore.case = TRUE)]])
    gA$exons <- .validGRanges(geneAnnotation[[grep("exons", names(geneAnnotation), ignore.case = TRUE)]])
    gA$TSS <- .validGRanges(geneAnnotation[[grep("TSS", names(geneAnnotation), ignore.case = TRUE)]])
    
  }else{
    stop("geneAnnotation must be a list/SimpleList of 3 GRanges for : Genes GRanges, Exons GRanges and TSS GRanges!")
  }
  
  gA
  
}

.validGenomeAnnotation <- function(genomeAnnotation = NULL){
  
  if(!inherits(genomeAnnotation, "SimpleList")){
    if(inherits(genomeAnnotation, "list")){
      genomeAnnotation <- as(genomeAnnotation, "SimpleList")
    }else{
      stop("genomeAnnotation must be a list/SimpleList of 3 GRanges for : blacklist GRanges, chromSizes GRanges and genome BSgenome package string (ie hg38 or BSgenome.Hsapiens.UCSC.hg38)!")
    }
  }
  
  if(identical(sort(tolower(names(genomeAnnotation))), c("blacklist", "chromsizes", "genome"))){
    
    gA <- SimpleList()
    gA$blacklist <- .validGRanges(genomeAnnotation[[grep("blacklist", names(genomeAnnotation), ignore.case = TRUE)]])
    if(genomeAnnotation[[grep("genome", names(genomeAnnotation), ignore.case = TRUE)]]=="nullGenome"){
      gA$genome <- "nullGenome"
    }else{
      bsg <- validBSgenome(genomeAnnotation[[grep("genome", names(genomeAnnotation), ignore.case = TRUE)]])
      gA$genome <- bsg@pkgname
    }
    gA$chromSizes <- .validGRanges(genomeAnnotation[[grep("chromsizes", names(genomeAnnotation), ignore.case = TRUE)]])
    
  }else{
    
    stop("genomeAnnotation must be a list/SimpleList of 3 GRanges for : blacklist GRanges, chromSizes GRanges and genome BSgenome package string (ie hg38 or BSgenome.Hsapiens.UCSC.hg38)!")
    
  }
  
  gA
  
}

.validGeneAnnoByGenomeAnno <- function(geneAnnotation, genomeAnnotation){
  
  allSeqs <- unique(paste0(seqnames(genomeAnnotation$chromSizes)))
  
  geneSeqs <- unique(paste0(seqnames(geneAnnotation$genes)))
  if(!all(geneSeqs %in% allSeqs)){
    geneNotIn <- geneSeqs[which(geneSeqs %ni% allSeqs)]
    message("Found Gene Seqnames not in GenomeAnnotation chromSizes, Removing : ", paste0(geneNotIn, collapse=","))
    geneAnnotation$genes <- .subsetSeqnamesGR(geneAnnotation$genes, names = allSeqs)
  }
  
  exonSeqs <- unique(paste0(seqnames(geneAnnotation$exons)))
  if(!all(exonSeqs %in% allSeqs)){
    exonNotIn <- exonSeqs[which(exonSeqs %ni% allSeqs)]
    message("Found Exon Seqnames not in GenomeAnnotation chromSizes, Removing : ", paste0(exonNotIn, collapse=","))
    geneAnnotation$exons <- .subsetSeqnamesGR(geneAnnotation$exons, names = allSeqs)
  }
  
  TSSSeqs <- unique(paste0(seqnames(geneAnnotation$TSS)))
  if(!all(TSSSeqs %in% allSeqs)){
    TSSNotIn <- TSSSeqs[which(TSSSeqs %ni% allSeqs)]
    message("Found TSS Seqnames not in GenomeAnnotation chromSizes, Removing : ", paste0(TSSNotIn, collapse=","))
    geneAnnotation$TSS <- .subsetSeqnamesGR(geneAnnotation$TSS, names = allSeqs)
  }
  
  geneAnnotation
  
}


.validArchRProject <- function(ArchRProj = NULL){
  if(!inherits(ArchRProj, "ArchRProject")){
    stop("Not a valid ArchRProject as input!")
  }else{
    ArchRProj
  }
}


####################################
# Log Tools
####################################

#' Set ArchR Logging
#' 
#' This function will set ArchR logging
#'
#' @param useLogs A boolean describing whether to use logging with ArchR.
#' @export
addArchRLogging <- function(useLogs = TRUE){
  .validInput(input = useLogs, name = "useLogs", valid = "boolean")
  message("Setting ArchRLogging = ", useLogs)
  options(ArchR.logging = useLogs)
  return(invisible(0))
}

#' Get ArchR Logging
#' 
#' This function will get ArchR logging
#'
#' @export
getArchRLogging <- function(){
  ArchRLogging <- options()[["ArchR.logging"]]
  if(!is.logical(ArchRLogging)){
    options(ArchR.logging = TRUE)
    return(TRUE)
  }
  ArchRLogging
}

#' Set ArchR Debugging
#' 
#' This function will set ArchR Debugging which will save an RDS if an error is encountered.
#'
#' @param debug A boolean describing whether to use logging with ArchR.
#' @export
addArchRDebugging <- function(debug = FALSE){
  .validInput(input = debug, name = "debug", valid = "boolean")
  message("Setting ArchRDebugging = ", debug)
  options(ArchR.logging = debug)
  return(invisible(0))
}

#' Get ArchR Debugging
#' 
#' This function will get ArchR Debugging which will save an RDS if an error is encountered.
#'
#' @export
getArchRDebugging <- function(){
  ArchRDebugging <- options()[["ArchR.debugging"]]
  if(!is.logical(ArchRDebugging)){
    options(ArchR.debugging = FALSE)
    return(FALSE)
  }
  ArchRDebugging
}

#' Create a Log File for ArchR
#' 
#' This function will create a log file for ArchR functions. If ArchRLogging is not TRUE
#' this function will return NULL.
#'
#' @param name A character string to add a more descriptive name in log file.
#' @param logDir The path to a directory where log files should be written.
#' @export
createLogFile <- function(
  name = NULL,
  logDir = "ArchRLogs",
  useLogs = getArchRLogging()
){
  
  .validInput(input = name, name = "name", valid = "character")
  .validInput(input = logDir, name = "logDir", valid = "character")
  .validInput(input = useLogs, name = "useLogs", valid = "boolean")
  
  if(!useLogs){
    return(NULL)
  }
  dir.create(logDir, showWarnings = FALSE)
  if(is.null(name)){
    logFile <- .tempfile(pattern = "ArchR", fileext = ".log", tmpdir = logDir)
  }else{
    logFile <- .tempfile(pattern = paste0("ArchR-", name), fileext = ".log", tmpdir = logDir)
  }
  logFile
}

.messageDiffTime <- function(...){ #Deprecated
  .logDiffTime(...)
}

.logDiffTime <- function(
  main = "",
  t1 = NULL,
  verbose = TRUE,
  addHeader = FALSE,
  t2 = Sys.time(),
  units = "mins",
  header = "###########",
  tail = "elapsed.",
  precision = 3,
  logFile = NULL,
  useLogs = getArchRLogging()
){
  
  if(verbose){
    
    timeStamp <- tryCatch({
      dt <- abs(round(difftime(t2, t1, units = units),precision))
      if(addHeader){
        msg <- sprintf("%s\n%s : %s, %s %s %s\n%s", header, Sys.time(), main, dt, units, tail, header)
      }else{
        msg <- sprintf("%s : %s, %s %s %s", Sys.time(), main, dt, units, tail)
      }
      message(msg)
    }, error = function(x){
      message("Time Error : ", x)
    })
    
  }
  
  if(!useLogs){
    return(invisible(0))
  }
  
  if(!is.null(logFile)){
    if(file.exists(logFile)){
      logStamp <- tryCatch({
        dt <- abs(round(difftime(t2, t1, units = units),precision))
        if(addHeader){
          msg <- sprintf("%s\n%s : %s, %s %s %s\n%s", header, Sys.time(), main, dt, units, tail, header)
        }else{
          msg <- sprintf("%s : %s, %s %s %s", Sys.time(), main, dt, units, tail)
        }
        cat(paste0(msg,"\n"), file = logFile, append = TRUE)
      }, error = function(x){
        0
      })
    }
  }
  
  return(invisible(0))
  
}

.startLogging <- function(
  logFile = NULL, 
  useLogs = getArchRLogging()
){
  
  if(!useLogs){
    return(invisible(0))
  }
  
  if(is.null(logFile)){
    return(invisible(0))
  }
  
  if(file.exists(logFile)){
    return(invisible(0))
  }
  
  .getRam <- function(OS = .Platform$OS.type){
    if(grepl("linux", OS, ignore.case = TRUE)){
      ram <- paste0("Linux : ", as.numeric(system("awk '/MemTotal/ {print $2}' /proc/meminfo", intern = TRUE)))
    }else if(grepl("unix", OS, ignore.case = TRUE)){
      ram <- system("/usr/sbin/system_profiler SPHardwareDataType", intern = TRUE)
      ram <- paste0("MAC : ", gsub("Memory:","",gsub(" ","", grep("Memory", ram, value = TRUE))))
    }else{
      ram <- NA
    }
  }
  
  message("ArchR logging to : ", logFile, 
          "\nIf there is an issue, please report to github with logFile!")
  
  #Begin With
  cat(.ArchRLogo(ascii = "Package", messageLogo = FALSE), file = logFile, append = FALSE) 
  cat("\nLogging With ArchR!\n\n", file = logFile, append = TRUE) 
  cat(paste0("Start Time : ",Sys.time(),"\n\n"), file = logFile, append = TRUE)
  
  #ArchR Info
  cat("------- ArchR Info\n\n", file = logFile, append = TRUE)
  cat(paste0("ArchRThreads = ", getArchRThreads()), file = logFile, append = TRUE)
  tryCatch({
    if(!is.null(getArchRGenome())){
      cat(paste0("\nArchRGenome = ", getArchRGenome()), file = logFile, append = TRUE)
    }
  }, error = function(x){
  })
  cat("\n\n", file = logFile, append = TRUE)
  
  #Add Info
  cat("------- System Info\n\n", file = logFile, append = TRUE)
  cat(paste0("Computer OS = ", .Platform$OS.type), file = logFile, append = TRUE)
  tryCatch({
    cat(paste0("\nTotal Cores = ", detectCores()), file = logFile, append = TRUE)
  }, error = function(x){
  })
  # tryCatch({
  #     cat(paste0("\nTotal RAM = ", .getRam()), file = logFile, append = TRUE)
  # }, error = function(x){
  # })
  cat("\n\n", file = logFile, append = TRUE)
  
  #Session Info
  cat("------- Session Info\n\n", file = logFile, append = TRUE)
  utils::capture.output(sessionInfo(), file = logFile, append = TRUE)
  cat("\n\n------- Log Info\n\n", file = logFile, append = TRUE)
  
  return(invisible(0))
  
}

.logMessage <- function(
  ..., 
  logFile = NULL,
  verbose = FALSE,   
  useLogs = getArchRLogging()
){
  
  msg <- utils::capture.output(message(...), type = "message")
  msg <- paste0(msg, collapse = "\n")
  
  if(is.null(msg)){
    stop("Message must be provided when logging!")
  }
  
  if(verbose){
    message(sprintf("%s : %s", Sys.time(), msg))
  }
  
  if(!useLogs){
    return(invisible(0))
  }
  
  if(is.null(logFile)){
    return(invisible(0))
  }
  
  cat(sprintf("\n%s : %s\n", Sys.time(), msg), file = logFile, append = TRUE)
  
  return(invisible(0))
  
}

.logHeader <- function(
  name = NULL, 
  logFile = NULL,   
  useLogs = getArchRLogging()
){
  
  if(!useLogs){
    return(invisible(0))
  }
  
  if(is.null(logFile)){
    return(invisible(0))
  }
  
  if(is.null(name)){
    stop("Name must be provided when logging!")
  }
  
  header <- "###########"
  cat(sprintf("\n%s\n%s : %s\n%s\n\n", header, Sys.time(), name, header), file = logFile, append = TRUE)
  
  return(invisible(0))
}

.logStop <- function(
  ..., 
  logFile = NULL,
  useLogs = getArchRLogging()
){
  
  msg <- utils::capture.output(message(...), type = "message")
  msg <- paste0(msg, collapse = "\n")
  
  if(is.null(msg)){
    stop("Message must be provided when logging!")
  }
  
  if(useLogs){
    if(!is.null(logFile)){
      cat(sprintf("\n%s : %s\n", Sys.time(), msg), file = logFile, append = TRUE)
    }
  }
  
  stop(sprintf("%s\n", msg), call. = FALSE)
  
  return(invisible(0))
  
}

.logError <- function(
  e = NULL,
  fn = NULL,
  info = NULL, 
  errorList = NULL,
  logFile = NULL,   
  useLogs = getArchRLogging(),
  throwError = TRUE,
  debug = getArchRDebugging()
){
  
  header <- "************************************************************"
  
  if(is.null(logFile)){
    useLogs <- FALSE
  }
  
  if(useLogs){
    #To Log File
    cat(sprintf("\n%s\n%s : ERROR Found in %s for %s \nLogFile = %s\n\n", header, Sys.time(), fn, info, logFile), file = logFile, append = TRUE)
    
    utils::capture.output(print(e), file = logFile, append = TRUE)
    
    if(!is.null(errorList)){
      tryCatch({
        #saveRDS(errorList, "Save-Error.rds")
        .logThis(errorList, name = "errorList", logFile)
      }, error = function(e){
        cat("Error recording errorList", file = logFile, append = TRUE)
      })
    }
    
    cat(sprintf("\n%s\n\n", header), file = logFile, append = TRUE)
  }
  
  #To Console
  cat(sprintf("\n%s\n%s : ERROR Found in %s for %s \nLogFile = %s\n\n", header, Sys.time(), fn, info, logFile))
  
  if(debug){
    if(!is.null(errorList)){
      debugFile <- paste0(gsub("\\.log", "", logFile), "-debug.rds")
      cat(sprintf("\n%s : ArchRDebugging is set to TRUE, DebugFile = %s\n", Sys.time(), debugFile))
      saveRDS(errorList, debugFile)
    }
  }
  
  print(e)
  
  cat(sprintf("\n%s\n\n", header))
  
  if(throwError) stop("Exiting See Error Above")
  
  return(invisible(0))
  
}


.logThis <- function(
  x = NULL, 
  name = NULL, 
  logFile = NULL, 
  useLogs = getArchRLogging()
){
  
  if(!useLogs){
    return(invisible(0))
  }
  
  if(is.null(logFile)){
    return(invisible(0))
  }
  
  if(!file.exists(logFile)){
    stop("logFile does not exist! Something may have deleted this file! Exiting...")
  }
  if(is.null(name)){
    stop("Name must be provided when logging!")
  }
  cat(paste0("\n", Sys.time(), " : ", name, ", Class = ", class(x), "\n"), file = logFile, append = TRUE)
  
  if(missing(x)){
    cat("Data is Missing\n\n", file = logFile, append = TRUE)
    return(invisible(0))
  }
  
  if(is.matrix(x)){
    
    px <- x[head(seq_len(nrow(x)), 5), head(seq_len(ncol(x)), 5), drop = FALSE]
    suppressMessages(utils::capture.output(print(px), file = logFile, append = TRUE))
    cat("\n", file = logFile, append = TRUE)
    cat(paste0(name, ": nRows = ", nrow(x), ", nCols = ", ncol(x), "\n"), file = logFile, append = TRUE)
    
  }else if(is.data.frame(x)){
    
    cat(paste0(name, ": nRows = ", nrow(x), ", nCols = ", ncol(x), "\n"), file = logFile, append = TRUE)
    suppressMessages(utils::capture.output(print(head(x)), file = logFile, append = TRUE))
    cat("\n", file = logFile, append = TRUE)
    
  }else if(is(x, "dgCMatrix") | is(x, "dgeMatrix")){
    
    cat(paste0(name, ": nRows = ", nrow(x), ", nCols = ", ncol(x), "\n"), file = logFile, append = TRUE)
    cat(paste0(name, ": NonZeroEntries = ", length(x@x), ", EntryRange = [ ", paste0(range(x@x), collapse=" , "), " ]\n"), file = logFile, append = TRUE)    
    px <- x[head(seq_len(nrow(x)), 5), head(seq_len(ncol(x)), 5), drop = FALSE]
    suppressMessages(utils::capture.output(print(px), file = logFile, append = TRUE))
    cat("\n", file = logFile, append = TRUE)
    
  }else if(is(x, "GRanges")){
    
    suppressMessages(utils::capture.output(print(x), file = logFile, append = TRUE))
    
  }else if(is(x, "SummarizedExperiment")){
    
    suppressMessages(utils::capture.output(print(x), file = logFile, append = TRUE))
    
  }else if(is(x, "DataFrame")){
    
    suppressMessages(utils::capture.output(print(x), file = logFile, append = TRUE))
    
  }else if(is(x, "ArchRProj")){
    
    suppressMessages(utils::capture.output(print(proj), file = logFile, append = TRUE))
    
  }else if(is(x, "SimpleList") | is(x, "list")){
    
    for(i in seq_along(x)){
      
      y <- x[[i]]
      
      if(missing(y)){
        next
      }
      
      if(is.matrix(y)){
        
        cat("\n", file = logFile, append = TRUE)
        cat(paste0(paste0(name,"$", names(x[i])), ": nRows = ", nrow(y), ", nCols = ", ncol(y), "\n"), file = logFile, append = TRUE)        
        px <- y[head(seq_len(nrow(y)), 5), head(seq_len(ncol(y)), 5), drop = FALSE]
        suppressMessages(utils::capture.output(print(px), file = logFile, append = TRUE))
        cat("\n", file = logFile, append = TRUE)
        
      }else if(is.data.frame(y)){
        
        cat("\n", file = logFile, append = TRUE)
        cat(paste0(paste0(name,"$", names(x[i])), ": nRows = ", nrow(y), ", nCols = ", ncol(y), "\n"), file = logFile, append = TRUE)
        suppressMessages(utils::capture.output(print(head(y)), file = logFile, append = TRUE))
        cat("\n", file = logFile, append = TRUE)
        
      }else if(is(y, "dgCMatrix")){
        
        cat("\n", file = logFile, append = TRUE)
        cat(paste0(paste0(name,"$", names(x[i])), ": nRows = ", nrow(y), ", nCols = ", ncol(y), "\n"), file = logFile, append = TRUE)
        cat(paste0(paste0(name,"$", names(x[i])), ": NonZeroEntries = ", length(y@x), ", EntryRange = [ ", paste0(range(y@x), collapse=" , "), " ]\n"), file = logFile, append = TRUE)        
        px <- y[head(seq_len(nrow(y)), 5), head(seq_len(ncol(y)), 5), drop = FALSE]
        suppressMessages(utils::capture.output(print(px), file = logFile, append = TRUE))
        cat("\n", file = logFile, append = TRUE)
        
      }else if(is(y, "SimpleList") | is(y, "list")){
        
        for(j in seq_along(y)){
          
          z <- y[[j]]
          
          if(missing(z)){
            next
          }
          
          if(is.matrix(z)){
            
            cat("\n", file = logFile, append = TRUE)
            cat(paste0(paste0(name,"$", names(y[j])), ": nRows = ", nrow(z), ", nCols = ", ncol(z), "\n"), file = logFile, append = TRUE)            
            px <- z[head(seq_len(nrow(z)), 5), head(seq_len(ncol(z)), 5), drop = FALSE]
            suppressMessages(utils::capture.output(print(px), file = logFile, append = TRUE))
            cat("\n", file = logFile, append = TRUE)
            
          }else if(is.data.frame(z)){
            
            cat("\n", file = logFile, append = TRUE)
            cat(paste0(paste0(name,"$", names(y[j])), ": nRows = ", nrow(z), ", nCols = ", ncol(z), "\n"), file = logFile, append = TRUE)
            suppressMessages(utils::capture.output(print(head(z)), file = logFile, append = TRUE))
            cat("\n", file = logFile, append = TRUE)
            
          }else if(is(z, "dgCMatrix")){
            
            cat("\n", file = logFile, append = TRUE)
            cat(paste0(paste0(name,"$", names(y[j])), ": nRows = ", nrow(z), ", nCols = ", ncol(z), "\n"), file = logFile, append = TRUE)
            cat(paste0(paste0(name,"$", names(y[j])), ": NonZeroEntries = ", length(z@x), ", EntrzRange = [ ", paste0(range(z@x), collapse=" , "), " ]\n"), file = logFile, append = TRUE)            
            px <- z[head(seq_len(nrow(z)), 5), head(seq_len(ncol(z)), 5), drop = FALSE]
            suppressMessages(utils::capture.output(print(px), file = logFile, append = TRUE))
            cat("\n", file = logFile, append = TRUE)
            
          }else if(is(y, "SimpleList") | is(y, "list")){
            
            #Only print 2x nested lists
            
          }else{
            
            tryCatch({
              cat("\n", file = logFile, append = TRUE)
              cat(paste0(paste0(name,"$", names(y[j])), ": length = ", length(z), "\n"), file = logFile, append = TRUE)
              suppressMessages(utils::capture.output(print(head(z)), file = logFile, append = TRUE))
              cat("\n", file = logFile, append = TRUE)
            }, error = function(q){
            })
            
          }
          
        }
        
      }else{
        
        tryCatch({
          cat("\n", file = logFile, append = TRUE)
          cat(paste0(paste0(name,"$", names(x[i])), ": length = ", length(y), "\n"), file = logFile, append = TRUE)
          suppressMessages(utils::capture.output(print(head(y)), file = logFile, append = TRUE))
          cat("\n", file = logFile, append = TRUE)
        }, error = function(q){
        })
        
      }
      
    }
    
  }else{
    
    tryCatch({
      cat("\n", file = logFile, append = TRUE)
      cat(paste0(name, ": length = ", length(x), "\n"), file = logFile, append = TRUE)
      suppressMessages(utils::capture.output(print(head(x)), file = logFile, append = TRUE))
      cat("\n", file = logFile, append = TRUE)
    }, error = function(q){
    })
    
  }
  
  cat("\n", file = logFile, append = TRUE)
  return(invisible(0))
  
}

.endLogging <- function(
  logFile = NULL, 
  useLogs = getArchRLogging()
){
  
  if(!useLogs){
    return(invisible(0))
  }
  
  if(is.null(logFile)){
    return(invisible(0))
  }
  
  rL <- readLines(logFile)
  o <- tryCatch({
    t1 <- gsub("Start Time : ","", grep("Start Time", rL, ignore.case = TRUE, value = TRUE))
    mn <- as.numeric(difftime(Sys.time(), as.POSIXct(t1), units = "mins"))
    hr <- as.numeric(difftime(Sys.time(), as.POSIXct(t1), units = "hours"))
    cat("\n------- Completed\n\n", file = logFile, append = TRUE)
    cat(paste0("End Time : ",Sys.time(),"\n"), file = logFile, append = TRUE)
    cat(paste0("Elapsed Time Minutes = ", mn), file = logFile, append = TRUE)
    cat(paste0("\nElapsed Time Hours = ", hr), file = logFile, append = TRUE)
    cat("\n\n-------\n\n\n\n", file = logFile, append = TRUE)
    message("ArchR logging successful to : ", logFile)
  }, error = function(x){
  })
  
  # tryCatch({
  #   R.utils::gzip(logFile, paste0(logFile, ".gz"))
  #   message("ArchR logging successful to : ", paste0(logFile, ".gz"))
  # }, error = function(x){
  # })
  
  return(invisible(0))
  
}


##########################################################################################
# Helper Intermediate Methods
##########################################################################################

.mergeParams <- function(paramInput = NULL, paramDefault = NULL){
  for(i in seq_along(paramDefault)){
    if(!(names(paramDefault)[i] %in% names(paramInput))){
      paramInput[[names(paramDefault)[i]]] <- paramDefault[[i]]
    }
  }
  return(paramInput)
}

.requirePackage <- function(x = NULL, load = TRUE, installInfo = NULL, source = NULL){
  if(x %in% rownames(installed.packages())){
    if(load){
      suppressPackageStartupMessages(require(x, character.only = TRUE))
    }else{
      return(0)
    }
  }else{
    if(!is.null(source) & is.null(installInfo)){
      if(tolower(source) == "cran"){
        installInfo <- paste0('install.packages("',x,'")')
      }else if(tolower(source) == "bioc"){
        installInfo <- paste0('BiocManager::install("',x,'")')
      }else{
        stop("Unrecognized package source, available are cran/bioc!")
      }
    }
    if(!is.null(installInfo)){
      stop(paste0("Required package : ", x, " is not installed/found!\n  Package Can Be Installed : ", installInfo))
    }else{
      stop(paste0("Required package : ", x, " is not installed/found!"))
    }
  }
}

##########################################################################################
# Stat/Summary Methods
##########################################################################################

.computeKNN <- function(
  data = NULL,
  query = NULL,
  k = 50,
  includeSelf = FALSE,
  ...
){
  .validInput(input = data, name = "data", valid = c("dataframe", "matrix"))
  .validInput(input = query, name = "query", valid = c("dataframe", "matrix"))
  .validInput(input = k, name = "k", valid = c("integer"))
  .validInput(input = includeSelf, name = "includeSelf", valid = c("boolean"))
  if(is.null(query)){
    query <- data
    searchSelf <- TRUE
  }else{
    searchSelf <- FALSE
  }
  .requirePackage("nabor", source = "cran")
  if(searchSelf & !includeSelf){
    knnIdx <- nabor::knn(data = data, query = query, k = k + 1, ...)$nn.idx
    knnIdx <- knnIdx[,-1,drop=FALSE]
  }else{
    knnIdx <- nabor::knn(data = data, query = query, k = k, ...)$nn.idx
  }
  knnIdx
}

.rowZscores <- function(m = NULL, min = -2, max = 2, limit = FALSE){
  z <- sweep(m - rowMeans(m), 1, matrixStats::rowSds(m),`/`)
  if(limit){
    z[z > max] <- max
    z[z < min] <- min
  }
  return(z)
}

.computeROC <- function(labels = NULL, scores = NULL, name="ROC"){
  .calcAUC <- function(TPR = NULL, FPR = NULL){
    # http://blog.revolutionanalytics.com/2016/11/calculating-auc.html
    dFPR <- c(diff(FPR), 0)
    dTPR <- c(diff(TPR), 0)
    out <- sum(TPR * dFPR) + sum(dTPR * dFPR)/2
    return(out)
  }
  labels <- labels[order(scores, decreasing=TRUE)]
  df <- data.frame(
    False_Positive_Rate = cumsum(!labels)/sum(!labels),
    True_Positive_Rate =  cumsum(labels)/sum(labels)
  )
  df$AUC <- round(.calcAUC(df$True_Positive_Rate,df$False_Positive_Rate),3)
  df$name <- name
  return(df)
}

.getQuantiles <- function(v = NULL, len = length(v)){
  if(length(v) < len){
    v2 <- rep(0, len)
    v2[seq_along(v)] <- v
  }else{
    v2 <- v
  }
  p <- trunc(rank(v2))/length(v2)
  if(length(v) < len){
    p <- p[seq_along(v)]
  }
  return(p)
}

.rowScale <- function(mat = NULL, min = NULL, max = NULL){
  if(!is.null(min)){
    rMin <- min
  }else{
    rMin <- matrixStats::rowMins(mat)
  }
  if(!is.null(max)){
    rMax <- max
  }else{
    rMax <- matrixStats::rowMaxs(mat)
  }
  rScale <- rMax - rMin
  matDiff <- mat - rMin
  matScale <- matDiff/rScale
  out <- list(mat=matScale, min=rMin, max=rMax)
  return(out)
}

.quantileCut <- function(x = NULL, lo = 0.025, hi = 0.975, maxIf0 = TRUE){
  q <- quantile(x, probs = c(lo,hi))
  if(q[2] == 0){
    if(maxIf0){
      q[2] <- max(x)
    }
  }
  x[x < q[1]] <- q[1]
  x[x > q[2]] <- q[2]
  return(x)
}

.normalizeCols <- function(mat = NULL, colSm = NULL, scaleTo = NULL){
  if(is.null(colSm)){
    colSm <- Matrix::colSums(mat)
  }
  if(!is.null(scaleTo)){
    mat@x <- scaleTo * mat@x / rep.int(colSm, Matrix::diff(mat@p))
  }else{
    mat@x <- mat@x / rep.int(colSm, Matrix::diff(mat@p))
  }
  return(mat)
}

.safeSubset <- function(mat = NULL, subsetRows = NULL, subsetCols = NULL){
  
  if(!is.null(subsetRows)){
    idxNotIn <- which(subsetRows %ni% rownames(mat))
    if(length(idxNotIn) > 0){
      subsetNamesNotIn <- subsetRows[idxNotIn]
      matNotIn <- Matrix::sparseMatrix(i=1,j=1,x=0,dims=c(length(idxNotIn), ncol = ncol(mat)))
      rownames(matNotIn) <- subsetNamesNotIn
      mat <- rbind(mat, matNotIn)
    }
    mat <- mat[subsetRows,]
  }
  
  if(!is.null(subsetCols)){
    idxNotIn <- which(subsetCols %ni% colnames(mat))
    if(length(idxNotIn) > 0){
      subsetNamesNotIn <- subsetCols[idxNotIn]
      matNotIn <- Matrix::sparseMatrix(i=1,j=1,x=0,dims=c(nrow(mat), ncol = length(idxNotIn)))
      colnames(matNotIn) <- subsetNamesNotIn
      mat <- cbind(mat, matNotIn)
    }
    mat <- mat[,subsetCols]
  }
  
  mat
  
}

.groupMeans <- function(mat = NULL, groups=NULL, na.rm = TRUE, sparse = FALSE){
  stopifnot(!is.null(groups))
  stopifnot(length(groups)==ncol(mat))
  gm <- lapply(unique(groups), function(x){
    if(sparse){
      Matrix::rowMeans(mat[,which(groups==x),drop=F], na.rm=na.rm)
    }else{
      rowMeans(mat[,which(groups==x),drop=F], na.rm=na.rm)
    }
  }) %>% Reduce("cbind",.)
  colnames(gm) <- unique(groups)
  return(gm)
}

.groupSums <- function(mat = NULL, groups=NULL, na.rm = TRUE, sparse = FALSE){
  stopifnot(!is.null(groups))
  stopifnot(length(groups)==ncol(mat))
  gm <- lapply(unique(groups), function(x){
    if(sparse){
      Matrix::rowSums(mat[,which(groups==x),drop=F], na.rm=na.rm)
    }else{
      rowSums(mat[,which(groups==x),drop=F], na.rm=na.rm)
    }
  }) %>% Reduce("cbind",.)
  colnames(gm) <- unique(groups)
  return(gm)
}

.groupSds <- function(mat = NULL, groups = NULL, na.rm = TRUE, sparse = FALSE){
  stopifnot(!is.null(groups))
  stopifnot(length(groups)==ncol(mat))
  gs <- lapply(unique(groups), function(x){
    if (sparse){
      matrixStats::rowSds(as.matrix(mat[, which(groups == x), drop = F]), na.rm = na.rm)
    }else{
      matrixStats::rowSds(mat[, which(groups == x), drop = F], na.rm = na.rm)
    }
  }) %>% Reduce("cbind",.)
  colnames(gs) <- unique(groups)
  return(gs)
}

.centerRollMean <- function(v = NULL, k = NULL){
  o1 <- data.table::frollmean(v, k, align = "right", na.rm = FALSE)
  if(k%%2==0){
    o2 <- c(rep(o1[k], floor(k/2)-1), o1[-seq_len(k-1)], rep(o1[length(o1)], floor(k/2)))
  }else if(k%%2==1){
    o2 <- c(rep(o1[k], floor(k/2)), o1[-seq_len(k-1)], rep(o1[length(o1)], floor(k/2)))
  }else{
    stop("Error!")
  }
  o2
}

##########################################################################################
# Miscellaneous Methods
##########################################################################################

.splitEvery <- function(x = NULL, n = NULL){
  #https://stackoverflow.com/questions/3318333/split-a-vector-into-chunks-in-r
  if(is.atomic(x)){
    split(x, ceiling(seq_along(x) / n))
  }else{
    split(x, ceiling(seq_len(nrow(x)) / n))
  }
}

.suppressAll <- function(expr = NULL){
  suppressPackageStartupMessages(suppressMessages(suppressWarnings(expr)))
}

.getAssay <- function(se = NULL, assayName = NULL){
  .assayNames <- function(se){
    names(SummarizedExperiment::assays(se))
  }
  if(is.null(assayName)){
    o <- SummarizedExperiment::assay(se)
  }else if(assayName %in% .assayNames(se)){
    o <- SummarizedExperiment::assays(se)[[assayName]]
  }else{
    stop(sprintf("assayName '%s' is not in assayNames of se : %s", assayName, paste(.assayNames(se),collapse=", ")))
  }
  return(o)
}

.fileExtension <- function (x = NULL){
  pos <- regexpr("\\.([[:alnum:]]+)$", x)
  ifelse(pos > -1L, substring(x, pos + 1L), "")
}

.checkPath <- function(u = NULL, path = NULL, throwError = TRUE){
  if(is.null(u)){
    out <- TRUE
  }
  out <- lapply(u, function(x, error = TRUE){
    if (Sys.which(x) == "") {
      if(!is.null(path) && file.exists(file.path(path, x))){
        o <- TRUE
      }else{
        if(throwError){
          stop(x, " not found in path, please add ", x, " to path!")
        }else{
          o <- FALSE
        }
      }
    }else{
      o <- TRUE
    }
    return(o)
  }) %>% unlist %>% all
  return(out)
}

.tempfile <- function(pattern = "tmp", tmpdir = "tmp", fileext = "", addDOC = TRUE){
  
  dir.create(tmpdir, showWarnings = FALSE)
  
  if(addDOC){
    doc <- paste0("-Date-", Sys.Date(), "_Time-", gsub(":","-", stringr::str_split(Sys.time(), pattern=" ",simplify=TRUE)[1,2]))
  }else{
    doc <- ""
  }
  
  tempfile(pattern = paste0(pattern, "-"), tmpdir = tmpdir, fileext = paste0(doc, fileext))
  
}

.ArchRLogo <- function(ascii = "Logo", messageLogo = TRUE){
  Ascii <- list(
    Package = c("
                ___      .______        ______  __    __  .______      
                /   \\\     |   _  \\\      /      ||  |  |  | |   _  \\\     
                /  ^  \\\    |  |_)  |    |  ,----'|  |__|  | |  |_)  |    
                /  /_\\\  \\\   |      /     |  |     |   __   | |      /     
                /  _____  \\\  |  |\\\  \\\\___ |  `----.|  |  |  | |  |\\\  \\\\___.
                /__/     \\__\\ | _| `._____| \\______||__|  |__| | _| `._____|
                "),
    
    #modified from cyu@athena.mit.edu
    Logo = c("
             / |
             /    \\\
             .                                  /      |.
             \\\\\\                              /        |.
             \\\\\\                          /           `|.
             \\\\\\                      /              |.
             \\\                    /                |\\\
             \\\\#####\\\           /                  ||
             ==###########>      /                   ||
             \\\\##==......\\\    /                     ||
             ______ =       =|__ /__                     ||      \\\\\\\
             ,--' ,----`-,__ ___/'  --,-`-===================##========>
             \\\               '        ##_______ _____ ,--,__,=##,__   ///
             ,    __==    ___,-,__,--'#'  ==='      `-'    | ##,-/
             -,____,---'       \\\\####\\\\________________,--\\\\_##,/
             ___      .______        ______  __    __  .______      
             /   \\\     |   _  \\\      /      ||  |  |  | |   _  \\\     
             /  ^  \\\    |  |_)  |    |  ,----'|  |__|  | |  |_)  |    
             /  /_\\\  \\\   |      /     |  |     |   __   | |      /     
             /  _____  \\\  |  |\\\  \\\\___ |  `----.|  |  |  | |  |\\\  \\\\___.
             /__/     \\__\\ | _| `._____| \\______||__|  |__| | _| `._____|
             ")
  )
  
  if(messageLogo){
    message(Ascii[[ascii]])
  }else{
    Ascii[[ascii]]
  }
  
}

##########################################################################################
# Batch Methods
##########################################################################################

.safelapply <- function(..., threads = 1, preschedule = FALSE){
  
  if(tolower(.Platform$OS.type) == "windows"){
    threads <- 1
  }
  
  if(threads > 1){
    
    o <- mclapply(..., mc.cores = threads, mc.preschedule = preschedule)
    
    errorMsg <- list()
    
    for(i in seq_along(o)){ #Make Sure this doesnt explode!
      if(inherits(o[[i]], "try-error")){
        capOut <- utils::capture.output(o[[i]])
        capOut <- capOut[!grepl("attr\\(\\,|try-error", capOut)]
        capOut <- head(capOut, 10)
        capOut <- unlist(lapply(capOut, function(x) substr(x, 1, 250)))
        capOut <- paste0("\t", capOut)
        errorMsg[[length(errorMsg) + 1]] <- paste0(c(paste0("Error Found Iteration ", i, " : "), capOut), "\n")
      }
    }
    
    if(length(errorMsg) != 0){
      
      errorMsg <- unlist(errorMsg)
      errorMsg <- head(errorMsg, 50)
      errorMsg[1] <- paste0("\n", errorMsg[1])
      stop(errorMsg)
      
    }
    
  }else{
    
    o <- lapply(...)
    
  }
  
  o
  
}

.batchlapply <- function(args = NULL, sequential = FALSE){
  
  if(is.null(args$tstart)){
    args$tstart <- Sys.time()
  }
  
  #Determine Parallel Backend
  if(inherits(args$parallelParam, "BatchtoolsParam")){
    
    .logStop("Batchtools not yet fully supported please use local parallel threading!", logFile = args$logFile)
    
    .logDiffTime("Batch Execution w/ BatchTools through BiocParallel!", t1 = args$tstart, verbose = TRUE, logFile = args$logFile)
    
    require(BiocParallel)
    
    args$parallelParam <- btParam
    #Unlink registry Directory
    if(dir.exists(args$registryDir)){
      #Clean Up Registry
      unlink(args$registryDir, recursive = TRUE)# Delete registry directory
    }
    
    #Set Up Registry For Runnning
    args$parallelParam$registryargs <- batchtoolsRegistryargs(
      file.dir = args$registryDir,
      work.dir = getwd(),
      packages = character(0L),
      namespaces = character(0L),
      source = character(0L),
      load = character(0L)
    )
    
    #Register
    BPPARAM <- args$parallelParam
    register(BPPARAM)
    
    #Add To Args
    args$BPPARAM <- BPPARAM
    
    if("..." %in% names(args)){
      args["..."] <- NULL
    }
    
    #Run
    args <- args[names(args) %ni% c("threads", "parallelParam", "subThreading")]
    outlist <- do.call(bplapply, args)
    
  }else{
    
    .logDiffTime("Batch Execution w/ safelapply!", t1 = args$tstart, verbose = TRUE, logFile = args$logFile)
    if(sequential){
      args$subThreads <- args$threads
      args$threads <- 1
    }else{
      if(args$threads > length(args$X)){
        args$subThreads <- floor( args$threads / length(args$X) )
        args$threads <- length(args$X)
      }else{
        args$subThreads <- 1
      }
    }
    
    args <- args[names(args) %ni% c("registryDir", "parallelParam", "subThreading")]
    outlist <- do.call(.safelapply, args)
    
  }
  
  return(outlist)
  
}

.retryCatch <- function(expr, ..., maxAttempts = 3, warnAttempts = FALSE, nameFN = "FN", printInfo = NULL, logFile = NULL){
  currentAttempt <- 0
  completed <- FALSE
  while(!completed & currentAttempt <= maxAttempts){
    currentAttempt <- currentAttempt + 1
    if(currentAttempt > 1){
      .logMessage(nameFN, " : Error occured, attempting again (", currentAttempt - 1, " of ", maxAttempts, ")", logFile = logFile)
    }
    ###########################################################
    tryResult <- tryCatch({
      #########################################################
      #Try Catch Statement Here
      if(warnAttempts){
        out <- return(expr)
      }else{
        out <- suppressWarnings(return(expr))
      }
      #########################################################
      list(out = out, completed = TRUE)
    }, error = function(e){
      list(out = e, completed = FALSE)
    }, ...)
    ###########################################################
    completed <- tryResult$completed
  }
  if(!completed){
    .logMessage(nameFN, " : Error occured and could not be resolved after ", maxAttempts, " additional attempts!", logFile = logFile)
    if(!is.null(printInfo)){
      .logMessage("Error occured at ", printInfo, logFile = logFile)
    }
    print(tryResult[[1]])
    stop()
  }
  
  tryResult[[1]]
  
}


##########################################################################################
# Developer Utils
##########################################################################################

.devMode <- function(package = "ArchR"){
  # fn <- unclass(lsf.str(envir = asNamespace(package), all = TRUE))
  # for(i in seq_along(fn)){
  #   tryCatch({
  #     assign(fn[i], paste0(package,':::', fn[i]), envir=globalenv())
  #     #eval(parse(text=paste0(fn[i], paste0('<<-',package,':::'), fn[i])))
  #   }, error = function(x){
  #   })
  # }
  fn <- unclass(lsf.str(envir = asNamespace("ArchR"), all = TRUE))
  for(i in seq_along(fn)){
    tryCatch({
      eval(parse(text=paste0(fn[i], '<-ArchR:::', fn[i])))
    }, error = function(x){
    })
  }
}

.convertToPNG <- function(
  ArchRProj = NULL,
  paths = c("QualityControl"),
  recursive = TRUE,
  outDir = "Figures",
  command = "mv"
){
  
  #If error try
  #brew install fontconfig
  
  .requirePackage("pdftools", source = "cran")
  
  if(!is.null(ArchRProj)){
    paths <- c(paths, file.path(getOutputDirectory(ArchRProj), "Plots"))
  }
  
  pdfFiles <- lapply(seq_along(paths), function(i){
    if(recursive){
      dirs <- list.dirs(paths[i], recursive = FALSE, full.names = FALSE)
      if(length(dirs) > 0){
        pdfs <- lapply(seq_along(dirs), function(j){
          list.files(file.path(paths[i], dirs[j]), full.names = TRUE, pattern = "\\.pdf")
        }) %>% unlist
      }else{
        pdfs <- c()
      }
      pdfs <- c(list.files(paths[i], full.names = TRUE, pattern = "\\.pdf"), pdfs)
    }else{
      pdfs <- list.files(paths[i], full.names = TRUE, pattern = "\\.pdf")
    }
    pdfs
  }) %>% unlist
  
  dir.create(outDir, showWarnings = FALSE)
  
  for(i in seq_along(pdfFiles)){
    print(i)
    tryCatch({
      pdf_convert(
        pdfFiles[i], 
        format = "png", 
        pages = NULL, 
        filenames = file.path(outDir, gsub("\\.pdf", "_%d.png",basename(pdfFiles[i]))),
        dpi = 300, 
        opw = "", 
        upw = "", 
        verbose = TRUE
      )
      system(paste0(command, " ", pdfFiles[i], " ", file.path(outDir, basename(pdfFiles[i]))))
    },error=function(x){
      0
    })
  }
  
}

####################################################################
# Hidden Helper Utils for Arrow Files
####################################################################

.validArrow <- function(ArrowFile = NULL){
  o <- h5closeAll()
  if(h5read(ArrowFile,"Class")!="Arrow"){
    warning(
      "This file is not a valid ArrowFile, this most likely is a bug with previous function where the class was not added.\n",
      "To fix your ArrowFiles :\n",
      "\tlapply(getArrowFiles(ArchRProj), function(x) h5write(obj = 'Arrow', file = x, name = 'Class'))",
      "\nThis will be an error in future versions."
    )
  }
  o <- h5closeAll()
  return(ArrowFile)
}

.isProtectedArray <- function(matrixName = NULL, exclude = NULL){
  protectedArrays <- tolower(c("peakmatrix", "tilematrix", "genescorematrix"))
  if(!is.null(exclude)){
    protectedArrays <- protectedArrays[protectedArrays %ni% tolower(exclude)]
  }
  if(tolower(matrixName) %in% protectedArrays){
    stop(sprintf("Error %s cannot be used as this conflicts with another predefined matrix function!", matrixName))
  }
  matrixName
}

.availableArrays <- function(ArrowFiles = NULL, threads = getArchRThreads()){
  threads <- min(threads, length(ArrowFiles))
  o <- h5closeAll()
  availableArrays <- .safelapply(seq_along(ArrowFiles), function(x){
    groups <- h5ls(ArrowFiles[x]) %>% {.[.$group=="/" & .$otype=="H5I_GROUP","name"]}
    groups <- groups[!grepl("Fragments|Metadata", groups)]
    groups
  }, threads = threads) %>% Reduce("intersect", .)
  o <- h5closeAll()
  return(availableArrays)
}

.availableSeqnames <- function(ArrowFiles = NULL, subGroup = "Fragments", threads = getArchRThreads()){
  threads <- min(threads, length(ArrowFiles))
  o <- h5closeAll()
  seqList <- .safelapply(seq_along(ArrowFiles), function(x){
    seqnames <- h5ls(ArrowFiles[x]) %>% {.[.$group==paste0("/",subGroup),]$name}
    seqnames <- seqnames[!grepl("Info", seqnames)]
    seqnames
  }, threads = threads)
  if(!all(unlist(lapply(seq_along(seqList), function(x) identical(seqList[[x]],seqList[[1]]))))){
    stop("Not All Seqnames Identical!")
  }
  o <- h5closeAll()
  return(paste0(seqList[[1]]))
}

.availableChr <- function(ArrowFiles = NULL, subGroup = "Fragments"){
  seqnames <- .availableSeqnames(ArrowFiles, subGroup)
  # if(getArchRChrPrefix()){
  #   seqnames <- seqnames[grep("chr", seqnames, ignore.case = TRUE)]
  # }
  if(length(seqnames) == 0){
    stop("No Chr Found in ArrowFiles!")
  }
  return(seqnames)
}

.availableCells <- function(ArrowFile = NULL, subGroup = NULL, passQC = TRUE){
  if(is.null(subGroup)){
    o <- h5closeAll()
    cellNames <- h5read(ArrowFile, "Metadata/CellNames")
    if(passQC){
      passQC <- tryCatch({
        h5read(ArrowFile, "Metadata/PassQC")
      }, error = function(x){
        rep(1, length(cellNames))
      })
      cellNames <- cellNames[which(passQC==1)]
    }
    sampleName <- h5read(ArrowFile, paste0("Metadata/Sample"))
    o <- h5closeAll()
  }else{
    o <- h5closeAll()
    cellNames <- h5read(ArrowFile, paste0(subGroup, "/Info/CellNames"))
    sampleName <- h5read(ArrowFile, paste0("Metadata/Sample"))
    o <- h5closeAll()
  }
  return(paste0(sampleName,"#",cellNames))
}

.sampleName <- function(ArrowFile = NULL){
  o <- h5closeAll()
  sampleName <- h5read(ArrowFile, paste0("Metadata/Sample"))
  o <- h5closeAll()
  return(sampleName)
}

.summarizeArrowContent <- function(ArrowFile = NULL){
  
  o <- h5closeAll()
  
  #Get Contents of ArrowFile
  h5DF <- h5ls(ArrowFile)
  
  #Re-Organize Content Info
  h5DF <- h5DF[-which(h5DF$group == "/"),]
  groups <- stringr::str_split(h5DF$group, pattern = "/", simplify=TRUE)[,2]
  groupList <- split(h5DF, groups)
  
  #Split Nested Lists
  groupList2 <- lapply(seq_along(groupList), function(x){
    groupDFx <- groupList[[x]]
    groupx <- gsub(paste0("/", names(groupList)[x]),"",groupDFx$group)
    if(all(groupx=="")){
      groupDFx
    }else{
      subDF <- groupDFx[-which(groupx == ""),]
      split(subDF, stringr::str_split(subDF$group, pattern = "/", simplify=TRUE)[,3])
    }
  })
  names(groupList2) <- names(groupList)
  
  
  o <- h5closeAll()
  
  return(groupList2)
  
}

.getMetadata <- function(ArrowFile = NULL){
  
  o <- h5closeAll()
  
  #Get Contents of ArrowFile
  sampleName <- h5read(ArrowFile, paste0("Metadata/Sample"))
  arrowMD <- .summarizeArrowContent(ArrowFile)$Metadata
  
  #Which are same dimensions as cell names
  arrowMD <- arrowMD[which(arrowMD$dim == arrowMD$dim[arrowMD$name=="CellNames"]),]
  
  #Load these into a S4 DataFrame
  md <- lapply(seq_len(nrow(arrowMD)), function(x){
    dfx <- DataFrame(h5read(ArrowFile, paste0(arrowMD$group[x],"/",arrowMD$name[x])))
    colnames(dfx) <- arrowMD$name[x]
    dfx
  }) %>% Reduce("cbind", .)
  
  #Correct CellNames
  md$CellNames <- paste0(sampleName,"#",md$CellNames)
  md$Sample <- Rle(sampleName, nrow(md))
  rownames(md) <- md$CellNames
  md <- md[, -which(colnames(md)=="CellNames")]
  md <- md[,order(colnames(md))]
  
  o <- h5closeAll()
  
  return(md)
}

.getFeatureDF <- function(ArrowFiles = NULL, subGroup = "TileMatrix", threads = getArchRThreads()){
  
  threads <- min(threads, length(ArrowFiles))
  
  .helpFeatureDF <- function(ArrowFile = NULL, subGroup = NULL){
    o <- h5closeAll()
    featureDF <- DataFrame(h5read(ArrowFile, paste0(subGroup,"/Info/FeatureDF")))
    featureDF$seqnames <- Rle(as.character(featureDF$seqnames))
    o <- h5closeAll()
    return(featureDF)
  }
  
  fdf <- .helpFeatureDF(ArrowFiles[1], subGroup = subGroup)
  
  if(length(ArrowFiles) > 1){
    ArrowFiles <- ArrowFiles[-1]
    checkIdentical <- .safelapply(seq_along(ArrowFiles), function(x){
      fdfx <- .helpFeatureDF(ArrowFiles[x], subGroup = subGroup)
      identical(fdfx, fdf)
    }, threads = threads) %>% unlist %>% all
    if(!checkIdentical){
      stop("Error not all FeatureDF for asssay is the same!")
    }
  }
  
  #Re-Order for Split Check!
  newOrder <- split(seq_len(nrow(fdf)), fdf$seqnames) %>% {lapply(seq_along(.), function(x) .[[x]])} %>% Reduce("c", .)
  fdf[newOrder,]
  
}

#####################################################################
# Dropping Group From Hdf5 File
#####################################################################
.createArrowGroup <- function(
  ArrowFile = NULL, 
  group = "GeneScoreMatrix", 
  force = FALSE, 
  verbose = FALSE,
  logFile = NULL
){
  
  ArrowInfo <- .summarizeArrowContent(ArrowFile)
  if(group == "Fragments"){ #This shouldnt happen but just in case
    .logMessage(".createArrowGroup : Cannot create Group over Fragments in Arrow!", logFile = logFile)
    stop("Cannot create Group over Fragments in Arrow!")
  }
  
  if(group %in% names(ArrowInfo)){
    #We Should Check How Big it is if it exists
    ArrowGroup <- ArrowInfo[[group]]
    ArrowGroup <- ArrowGroup[names(ArrowGroup) %ni% c("Info")]
    if(length(ArrowGroup) > 0){
      if(!force){
        .logMessage(".createArrowGroup : Arrow Group already exists! Set force = TRUE to continue!", logFile = logFile)
        stop("Arrow Group already exists! Set force = TRUE to continue!")
      }else{
        .logMessage(".createArrowGroup : Arrow Group already exists! Dropping Group from ArrowFile! This will take ~10-30 seconds!", logFile = logFile)
        if(verbose) message("Arrow Group already exists! Dropping Group from ArrowFile! This will take ~10-30 seconds!")
        o <- .dropGroupsFromArrow(ArrowFile = ArrowFile, dropGroups = group, verbose = verbose, logFile = logFile)
        tryCatch({h5createGroup(ArrowFile , group)}, error=function(e){})
        invisible(return(0))
      }
    }
  }else{
    tryCatch({h5createGroup(ArrowFile , group)}, error=function(e){})
    invisible(return(0))
  }
  
}

.dropGroupsFromArrow <- function(
  ArrowFile = NULL, 
  dropGroups = NULL,
  level = 0,
  verbose = FALSE,
  logFile = NULL
){
  
  tstart <- Sys.time()
  
  #Summarize Arrow Content
  ArrowInfo <- .summarizeArrowContent(ArrowFile)
  
  .logMessage(".dropGroupsFromArrow : Initializing Temp ArrowFile", logFile = logFile)
  
  #We need to transfer first
  outArrow <- .tempfile(fileext = ".arrow")
  o <- h5closeAll()
  o <- h5createFile(outArrow)
  o <- h5write(obj = "Arrow", file = outArrow, name = "Class")
  o <- h5write(obj = paste0(packageVersion("ArchR")), file = outArrow, name = "ArchRVersion")
  
  #1. Metadata First
  .logMessage(".dropGroupsFromArrow : Adding Metadata to Temp ArrowFile", logFile = logFile)
  groupName <- "Metadata"
  o <- h5createGroup(outArrow, groupName)
  
  mData <- ArrowInfo[[groupName]]
  
  for(i in seq_len(nrow(mData))){
    h5name <- paste0(groupName, "/", mData$name[i])
    h5write(.h5read(ArrowFile, h5name), file = outArrow, name = h5name)
  }
  
  #2. Other Groups
  .logMessage(".dropGroupsFromArrow : Adding SubGroups to Temp ArrowFile", logFile = logFile)
  
  groupsToTransfer <- names(ArrowInfo)
  groupsToTransfer <- groupsToTransfer[groupsToTransfer %ni% "Metadata"]
  if(!is.null(dropGroups)){
    groupsToTransfer <- groupsToTransfer[tolower(groupsToTransfer) %ni% tolower(dropGroups)]
  }
  
  for(k in seq_along(groupsToTransfer)){
    
    .logDiffTime(paste0("Transferring ", groupsToTransfer[k]), tstart, verbose = verbose, logFile = logFile)
    
    #Create Group
    groupName <- groupsToTransfer[k]
    o <- h5createGroup(outArrow, groupName)
    
    #Sub Data
    mData <- ArrowInfo[[groupName]]
    
    #Get Order Of Sub Groups (Mostly Related to Seqnames)
    seqOrder <- sort(names(mData))
    if(any(grepl("chr", seqOrder))){
      seqOrder <- c(seqOrder[!grepl("chr", seqOrder)], seqOrder[grepl("chr", seqOrder)])
    }
    
    for(j in seq_along(seqOrder)){
      
      if(verbose) message(j, " ", appendLF = FALSE)
      
      #Create Group
      groupJ <- paste0(groupName, "/", seqOrder[j])
      o <- h5createGroup(outArrow, groupJ)
      
      #Sub mData
      mDataj <- mData[[seqOrder[j]]]
      
      #Transfer Components
      for(i in seq_len(nrow(mDataj))){
        h5name <- paste0(groupJ, "/", mDataj$name[i])
        .suppressAll(h5write(.h5read(ArrowFile, h5name), file = outArrow, name = h5name, level = level))
      }
      
    }
    
    gc()
    
    if(verbose) message("")
    
  }
  
  .logMessage(".dropGroupsFromArrow : Move Temp ArrowFile to ArrowFile", logFile = logFile)
  
  rmf <- file.remove(ArrowFile)
  out <- .fileRename(from = outArrow, to = ArrowFile)
  
  .logDiffTime("Completed Dropping of Group(s)", tstart, logFile = logFile, verbose = verbose)
  
  ArrowFile
  
}

.copyArrows <- function(
  inArrows = NULL,
  outArrows = NULL,
  cellsKeep = NULL,
  level = 0,
  verbose = FALSE,
  logFile = NULL,
  threads = 1
){
  
  stopifnot(length(inArrows) == length(outArrows))
  
  unlist(.safelapply(seq_along(inArrows), function(x){
    .copyArrowSingle(
      inArrow = inArrows[x],
      outArrow = outArrows[x],
      cellsKeep = cellsKeep,
      level = level,
      verbose = verbose,
      logFile = logFile
    )
  }, threads = threads))
  
}

.copyArrowSingle <- function(
  inArrow = NULL,
  outArrow = NULL,
  cellsKeep = NULL,
  level = 0,
  verbose = FALSE,
  logFile = NULL
){
  
  tstart <- Sys.time()
  
  #Summarize Arrow Content
  ArrowInfo <- .summarizeArrowContent(inArrow)
  sampleName <- .sampleName(inArrow)
  
  .logMessage(".copyArrow : Initializing Out ArrowFile", logFile = logFile)
  
  #We need to transfer first
  o <- .suppressAll(file.remove(outArrow))
  o <- h5closeAll()
  o <- h5createFile(outArrow)
  o <- h5write(obj = "Arrow", file = outArrow, name = "Class")
  o <- h5write(obj = paste0(packageVersion("ArchR")), file = outArrow, name = "ArchRVersion")
  
  #1. Metadata First
  .logMessage(".copyArrow : Adding Metadata to Out ArrowFile", logFile = logFile)
  groupName <- "Metadata"
  o <- h5createGroup(outArrow, groupName)
  
  mData <- ArrowInfo[[groupName]]
  
  for(i in seq_len(nrow(mData))){
    h5name <- paste0(groupName, "/", mData$name[i])
    h5write(.h5read(inArrow, h5name), file = outArrow, name = h5name)
  }
  
  #2. scATAC-Fragments
  .logDiffTime(paste0("Transferring Fragments"), tstart, verbose = verbose, logFile = logFile)
  
  #Create Group
  groupName <- "Fragments"
  o <- h5createGroup(outArrow, groupName)
  
  #Sub Data
  mData <- ArrowInfo[[groupName]]
  
  #Get Order Of Sub Groups (Mostly Related to Seqnames)
  seqOrder <- sort(names(mData))
  if(any(grepl("chr", seqOrder))){
    seqOrder <- c(seqOrder[!grepl("chr", seqOrder)], seqOrder[grepl("chr", seqOrder)])
  }
  
  for(j in seq_along(seqOrder)){
    
    if(verbose) message(j, " ", appendLF = FALSE)
    
    #Create Group
    groupJ <- paste0(groupName, "/", seqOrder[j])
    o <- h5createGroup(outArrow, groupJ)
    
    #Sub mData
    mDataj <- mData[[seqOrder[j]]]
    
    #Read In Fragments
    RGLengths <- .h5read(inArrow, paste0(groupJ, "/RGLengths"))
    RGValues <- .h5read(inArrow, paste0(groupJ, "/RGValues"))
    RGRle <- Rle(paste0(sampleName, "#", RGValues), RGLengths)
    
    #Determine Which to Keep
    idx <- BiocGenerics::which(RGRle %bcin% cellsKeep)
    RGRle <- RGRle[idx]
    RGLengths <- RGRle@lengths
    
    #print(head(RGRle@values))
    RGValues <- stringr::str_split(RGRle@values, pattern = "#", simplify = TRUE)[,2]
    
    #Create Data Sets
    # o <- .suppressAll(h5createDataset(outArrow, paste0(groupJ, "/Ranges"), storage.mode = "integer", dims = c(length(RGRle), 2), level = level))
    # o <- .suppressAll(h5createDataset(outArrow, paste0(groupJ, "/RGLengths"), storage.mode = "integer", dims = c(length(RGRle), 1), level = level))
    # o <- .suppressAll(h5createDataset(outArrow, paste0(groupJ, "/RGValues"), storage.mode = "character", dims = c(length(RGRle), 1), level = level, 
    #         size = max(nchar(RGValues) + 1)))
    
    #Write Barcodes
    o <- .suppressAll(h5write(RGLengths, file = outArrow, name = paste0(groupJ, "/RGLengths"), level = level))
    o <- .suppressAll(h5write(RGValues, file = outArrow, name = paste0(groupJ, "/RGValues"), level = level))
    
    #Write Ranges
    o <- .suppressAll(
      h5write(
        obj = .h5read(inArrow, paste0(groupJ, "/Ranges"))[idx, ], 
        file = outArrow, 
        name = paste0(groupJ, "/Ranges"), 
        level = level
      )
    )
    
  }
  
  if(verbose) message("")
  
  #3. Other Matrices
  .logMessage(".copyArrow : Adding SubMatrices to Out ArrowFile", logFile = logFile)
  groupsToTransfer <- names(ArrowInfo)
  groupsToTransfer <- groupsToTransfer[groupsToTransfer %ni% c("Metadata", "Fragments")]
  
  for(k in seq_along(groupsToTransfer)){
    
    .logDiffTime(paste0("Transferring ", groupsToTransfer[k]), tstart, verbose = verbose, logFile = logFile)
    
    #Create Group
    groupName <- groupsToTransfer[k]
    o <- h5createGroup(outArrow, groupName)
    
    #Sub Data
    mData <- ArrowInfo[[groupName]]
    
    #Get Order Of Sub Groups (Mostly Related to Seqnames)
    seqOrder <- sort(names(mData))
    if(any(grepl("chr", seqOrder))){
      seqOrder <- c(seqOrder[!grepl("chr", seqOrder)], seqOrder[grepl("chr", seqOrder)])
    }
    
    cellNames <- paste0(sampleName, "#", .h5read(inArrow, paste0(groupName, "/Info/CellNames")))
    featureDF <- .getFeatureDF(ArrowFile = inArrow, subGroup = groupName)
    
    for(j in seq_along(seqOrder)){
      
      if(verbose) message(j, " ", appendLF = FALSE)
      
      #Create Group
      groupJ <- paste0(groupName, "/", seqOrder[j])
      
      if(seqOrder[j] == "Info"){
        
        o <- h5createGroup(outArrow, groupJ)
        
        #Sub mData
        mDataj <- mData[[seqOrder[j]]]
        idxCL <- which(mDataj$dim == mDataj$dim[mDataj$name=="CellNames"])
        idxCL <- idxCL[mDataj$name[idxCL] %ni% "FeatureDF"]
        idxKeep <- which(cellNames %in% cellsKeep)
        
        #Transfer Components
        for(i in seq_len(nrow(mDataj))){
          
          h5name <- paste0(groupJ, "/", mDataj$name[i])
          
          if(i %in% idxCL){
            .suppressAll(h5write(.h5read(inArrow, h5name)[idxKeep], file = outArrow, name = h5name))
          }else{
            .suppressAll(h5write(.h5read(inArrow, h5name), file = outArrow, name = h5name))
          }
          
        }
        
      }else{
        
        #Sub mData
        mDataj <- mData[[seqOrder[j]]]
        addAnalysis <- mDataj[mDataj$name %ni% c("i", "jLengths", "jValues", "x"), "name"]
        
        mat <- .getMatFromArrow(
          ArrowFile = inArrow,
          featureDF = featureDF[BiocGenerics::which(featureDF$seqnames %bcin% seqOrder[j]),],
          useMatrix = groupName,
          cellNames = cellNames[cellNames %in% cellsKeep]
        )
        
        o <- .addMatToArrow(
          mat = mat, 
          ArrowFile = outArrow, 
          Group = paste0(groupName, "/", seqOrder[j]), 
          binarize = all(mat@x == 1),
          addColSums = "colSums" %in% addAnalysis,
          addRowSums = "rowSums" %in% addAnalysis,
          addRowMeans = "rowMeans" %in% addAnalysis,
          addRowVars = "rowVars" %in% addAnalysis,
          addRowVarsLog2 = "rowVarsLog2" %in% addAnalysis
        )
        
        rm(mat)
        
      }
      
    }
    
    gc()
    
    if(verbose) message("")
    
  }
  
  .logDiffTime("Completed Copying ArrowFile", tstart, logFile = logFile, verbose = verbose)
  
  outArrow
  
}


####################################################################
# Reading fragments from Arrow Files
####################################################################

#' Get the fragments from an ArrowFile 
#' 
#' This function retrieves the fragments from a given ArrowFile as a GRanges object.
#'
#' @param ArrowFile The path to the ArrowFile from which fragments should be obtained.
#' @param chr A name of a chromosome to be used to subset the fragments `GRanges` object to a specific chromsome if desired.
#' @param cellNames A character vector indicating the cell names of a subset of cells from which fragments whould be extracted.
#' This allows for extraction of fragments from only a subset of selected cells. By default, this function will extract all cells
#' from the provided ArrowFile using `getCellNames()`.
#' @param verbose A boolean value indicating whether to use verbose output during execution of this function. Can be set to `FALSE` for a cleaner output.
#' @param logFile The path to a file to be used for logging ArchR output.
#' @export
getFragmentsFromArrow <- function(
  ArrowFile = NULL, 
  chr = NULL, 
  cellNames = NULL, 
  verbose = TRUE,
  logFile = createLogFile("getFragmentsFromArrow")
){
  
  .validInput(input = ArrowFile, name = "ArrowFile", valid = "character")
  .validInput(input = chr, name = "chr", valid = c("character","null"))
  .validInput(input = cellNames, name = "cellNames", valid = c("character","null"))
  .validInput(input = verbose, name = "verbose", valid = c("boolean"))
  
  tstart <- Sys.time()
  .startLogging(logFile = logFile)
  .logThis(mget(names(formals()),sys.frame(sys.nframe())), "getFragmentsFromArrow Input-Parameters", logFile = logFile)
  
  ArrowFile <- .validArrow(ArrowFile)
  
  if(is.null(chr)){
    chr <- .availableSeqnames(ArrowFile, subGroup = "Fragments")
  }
  
  if(any(chr %ni% .availableSeqnames(ArrowFile, subGroup = "Fragments"))){
    stop("Error Chromosome not in ArrowFile!")
  }
  
  out <- lapply(seq_along(chr), function(x){
    .logDiffTime(sprintf("Reading Chr %s of %s", x, length(chr)), t1 = tstart, verbose = verbose, logFile = logFile)
    .getFragsFromArrow(
      ArrowFile = ArrowFile, 
      chr = chr[x], 
      out = "GRanges", 
      cellNames = cellNames, 
      method = "fast"
    )
  })
  
  
  .logDiffTime("Merging", tstart, t1 = tstart, verbose = verbose, logFile = logFile)
  
  out <- tryCatch({
    
    o <- .suppressAll(unlist(GRangesList(out, compress = FALSE)))
    
    if(.isGRList(o)){
      stop("Still a GRangesList")
    }
    
    o
    
  }, error = function(x){
    
    o <- c()
    
    for(i in seq_along(out)){
      if(!is.null(out[[i]])){
        if(i == 1){
          o <- out[[i]]
        }else{
          o <- c(o, out[[i]])
        }
      }
    }
    
    o
    
  })
  
  out
  
}

.getFragsFromArrow <- function(
  ArrowFile = NULL, 
  chr = NULL, 
  out = "GRanges", 
  cellNames = NULL, 
  method = "fast"
){
  
  if(is.null(chr)){
    stop("Need to provide chromosome to read!")
  }
  
  o <- h5closeAll()
  ArrowFile <- .validArrow(ArrowFile)
  
  avSeq <- .availableSeqnames(ArrowFile)
  if(chr %ni% avSeq){
    stop(paste0("Chromosome ", chr ," not in ArrowFile! Available Chromosomes are : ", paste0(avSeq, collapse=",")))
  }
  
  #Get Sample Name
  sampleName <- .h5read(ArrowFile, paste0("Metadata/Sample"), method = method)
  
  o <- h5closeAll()
  nFrags <- sum(.h5read(ArrowFile, paste0("Fragments/",chr,"/RGLengths"), method = method))
  
  if(nFrags==0){
    if(tolower(out)=="granges"){
      output <- GRanges(seqnames = chr, IRanges(start = 1, end = 1), RG = "tmp")
      output <- output[-1,]
    }else{
      output <- IRanges(start = 1, end = 1)
      mcols(output)$RG <- c("tmp")
      output <- output[-1,]
    }
    return(output)
  }
  
  if(is.null(cellNames) | tolower(method) == "fast"){
    
    output <- .h5read(ArrowFile, paste0("Fragments/",chr,"/Ranges"), method = method) %>% 
    {IRanges(start = .[,1], width = .[,2])}
    mcols(output)$RG <- Rle(
      values = paste0(sampleName, "#", .h5read(ArrowFile, paste0("Fragments/",chr,"/RGValues"), method = method)), 
      lengths = .h5read(ArrowFile, paste0("Fragments/",chr,"/RGLengths"), method = method)
    )
    if(!is.null(cellNames)){
      output <- output[BiocGenerics::which(mcols(output)$RG %bcin% cellNames)]
    }
    
  }else{
    
    if(!any(cellNames %in% .availableCells(ArrowFile))){
      
      stop("None of input cellNames are in ArrowFile availableCells!")
      
    }else{
      
      barRle <- Rle(h5read(ArrowFile, paste0("Fragments/",chr,"/RGValues")), h5read(ArrowFile, paste0("Fragments/",chr,"/RGLengths")))
      barRle@values <- paste0(sampleName, "#", barRle@values)
      idx <- BiocGenerics::which(barRle %bcin% cellNames)
      if(length(idx) > 0){
        output <- h5read(ArrowFile, paste0("Fragments/",chr,"/Ranges"), index = list(idx, 1:2)) %>% 
        {IRanges(start = .[,1], width = .[,2])}
        mcols(output)$RG <- barRle[idx]
      }else{
        output <- IRanges(start = 1, end = 1)
        mcols(output)$RG <- c("tmp")
        output <- output[-1,]
      }
    }
    
  }
  
  o <- h5closeAll()
  
  if(tolower(out)=="granges"){
    if(length(output) > 0){
      output <- GRanges(seqnames = chr, ranges(output), RG = mcols(output)$RG)    
    }else{
      output <- IRanges(start = 1, end = 1)
      mcols(output)$RG <- c("tmp")
      output <- GRanges(seqnames = chr, ranges(output), RG = mcols(output)$RG)
      output <- output[-1,]
    }
  }
  
  return(output)
}

####################################################################
# Reading Matrices/Arrays from Arrow Files
####################################################################

#' Get a data matrix stored in an ArchRProject
#' 
#' This function gets a given data matrix from an `ArchRProject`.
#'
#' @param ArchRProj An `ArchRProject` object to get data matrix from.
#' @param useMatrix The name of the data matrix to retrieve from the given ArrowFile. Options include "TileMatrix", "GeneScoreMatrix", etc.
#' @param useSeqnames A character vector of chromosome names to be used to subset the data matrix being obtained.
#' @param verbose A boolean value indicating whether to use verbose output during execution of  this function. Can be set to FALSE for a cleaner output.
#' @param binarize A boolean value indicating whether the matrix should be binarized before return. This is often desired when working with insertion counts.
#' @param logFile The path to a file to be used for logging ArchR output.
#' @export
getMatrixFromProject <- function(
  ArchRProj = NULL,
  useMatrix = "GeneScoreMatrix",
  useSeqnames = NULL,
  verbose = TRUE,
  binarize = FALSE,
  threads = getArchRThreads(),
  logFile = createLogFile("getMatrixFromProject")
){
  
  .validInput(input = ArchRProj, name = "ArchRProj", valid = c("ArchRProj"))
  .validInput(input = useMatrix, name = "useMatrix", valid = c("character"))
  .validInput(input = useSeqnames, name = "useSeqnames", valid = c("character","null"))
  .validInput(input = verbose, name = "verbose", valid = c("boolean"))
  .validInput(input = binarize, name = "binarize", valid = c("boolean"))
  .validInput(input = threads, name = "threads", valid = c("integer"))
  
  tstart <- Sys.time()
  .startLogging(logFile = logFile)
  .logThis(mget(names(formals()),sys.frame(sys.nframe())), "getMatrixFromProject Input-Parameters", logFile = logFile)
  
  ArrowFiles <- getArrowFiles(ArchRProj)
  
  cellNames <- ArchRProj$cellNames
  
  
  seL <- .safelapply(seq_along(ArrowFiles), function(x){
    
    .logDiffTime(paste0("Reading ", useMatrix," : ", names(ArrowFiles)[x], "(",x," of ",length(ArrowFiles),")"), 
                 t1 = tstart, verbose = FALSE, logFile = logFile)
    
    allCells <- .availableCells(ArrowFile = ArrowFiles[x], subGroup = useMatrix)
    allCells <- allCells[allCells %in% cellNames]
    
    if(length(allCells) != 0){
      
      o <- getMatrixFromArrow(
        ArrowFile = ArrowFiles[x],
        useMatrix = useMatrix,
        useSeqnames = useSeqnames,
        cellNames = allCells, 
        ArchRProj = ArchRProj,
        verbose = FALSE,
        binarize = binarize,
        logFile = logFile
      )
      
      .logDiffTime(paste0("Completed ", useMatrix," : ", names(ArrowFiles)[x], "(",x," of ",length(ArrowFiles),")"), 
                   t1 = tstart, verbose = FALSE, logFile = logFile)
      
      o
      
    }else{
      
      NULL
      
    }
    
  }, threads = threads) 
  
  #ColData
  .logDiffTime("Organizing colData", t1 = tstart, verbose = verbose, logFile = logFile)
  cD <- lapply(seq_along(seL), function(x){
    colData(seL[[x]])
  }) %>% Reduce("rbind", .)
  
  #RowData
  .logDiffTime("Organizing rowData", t1 = tstart, verbose = verbose, logFile = logFile)
  rD1 <- rowData(seL[[1]])
  rD <- lapply(seq_along(seL), function(x){
    identical(rowData(seL[[x]]), rD1)
  }) %>% unlist %>% all
  if(!rD){
    stop("Error with rowData being equal for every sample!")
  }
  
  #Assays
  nAssays <- names(assays(seL[[1]]))
  asy <- lapply(seq_along(nAssays), function(i){
    .logDiffTime(sprintf("Organizing Assays (%s of %s)", i, length(nAssays)), t1 = tstart, verbose = verbose, logFile = logFile)
    m <- lapply(seq_along(seL), function(j){
      assays(seL[[j]])[[nAssays[i]]]
    }) %>% Reduce("cbind", .)
    m
  }) %>% SimpleList()
  names(asy) <- nAssays
  
  .logDiffTime("Constructing SummarizedExperiment", t1 = tstart, verbose = verbose, logFile = logFile)
  se <- SummarizedExperiment(assays = asy, colData = cD, rowData = rD1)  
  rm(seL)
  gc()
  
  .logDiffTime("Finished Matrix Creation", t1 = tstart, verbose = verbose, logFile = logFile)
  
  se
  
}

#' Get a data matrix stored in an ArrowFile
#' 
#' This function gets a given data matrix from an individual ArrowFile.
#'
#' @param ArrowFile The path to an ArrowFile from which the selected data matrix should be obtained.
#' @param useMatrix The name of the data matrix to retrieve from the given ArrowFile. Options include "TileMatrix", "GeneScoreMatrix", etc.
#' @param useSeqnames A character vector of chromosome names to be used to subset the data matrix being obtained.
#' @param cellNames A character vector indicating the cell names of a subset of cells from which fragments whould be extracted.
#' This allows for extraction of fragments from only a subset of selected cells. By default, this function will extract all cells from
#' the provided ArrowFile using `getCellNames()`.
#' @param ArchRProj An `ArchRProject` object to be used for getting additional information for cells in `cellColData`.
#' In some cases, data exists within the `ArchRProject` object that does not exist within the ArrowFiles. To access this data, you can
#' provide the `ArchRProject` object here.
#' @param verbose A boolean value indicating whether to use verbose output during execution of  this function. Can be set to FALSE for a cleaner output.
#' @param binarize A boolean value indicating whether the matrix should be binarized before return. This is often desired when working with insertion counts.
#' @param logFile The path to a file to be used for logging ArchR output.
#' @export
getMatrixFromArrow <- function(
  ArrowFile = NULL, 
  useMatrix = "GeneScoreMatrix",
  useSeqnames = NULL,
  cellNames = NULL, 
  ArchRProj = NULL,
  verbose = TRUE,
  binarize = FALSE,
  logFile = createLogFile("getMatrixFromArrow")
){
  
  .validInput(input = ArrowFile, name = "ArrowFile", valid = "character")
  .validInput(input = useMatrix, name = "useMatrix", valid = "character")
  .validInput(input = useSeqnames, name = "useSeqnames", valid = c("character","null"))
  .validInput(input = cellNames, name = "cellNames", valid = c("character","null"))
  .validInput(input = ArchRProj, name = "ArchRProj", valid = c("ArchRProj","null"))
  .validInput(input = verbose, name = "verbose", valid = c("boolean"))
  .validInput(input = binarize, name = "binarize", valid = c("boolean"))
  
  tstart <- Sys.time()
  .startLogging(logFile = logFile)
  .logThis(mget(names(formals()),sys.frame(sys.nframe())), "getMatrixFromArrow Input-Parameters", logFile = logFile)
  
  ArrowFile <- .validArrow(ArrowFile)
  sampleName <- .sampleName(ArrowFile)
  
  seqnames <- .availableSeqnames(ArrowFile, subGroup = useMatrix)
  featureDF <- .getFeatureDF(ArrowFile, subGroup = useMatrix)
  .logThis(featureDF, paste0("featureDF ", sampleName), logFile = logFile)
  
  if(!is.null(useSeqnames)){
    seqnames <- seqnames[seqnames %in% useSeqnames]
  }
  
  if(length(seqnames) == 0){
    stop("No seqnames available!")
  }
  
  featureDF <- featureDF[BiocGenerics::which(featureDF$seqnames %bcin% seqnames), ]
  
  .logDiffTime(paste0("Getting ",useMatrix," from ArrowFile : ", basename(ArrowFile)), 
               t1 = tstart, verbose = verbose, logFile = logFile)
  
  if(!is.null(cellNames)){
    allCells <- .availableCells(ArrowFile = ArrowFile, subGroup = useMatrix)
    if(!all(cellNames %in% allCells)){
      stop("cellNames must all be within the ArrowFile!!!!")
    }
  }
  
  mat <- .getMatFromArrow(
    ArrowFile = ArrowFile, 
    featureDF = featureDF, 
    cellNames = cellNames, 
    useMatrix = useMatrix,
    binarize = binarize,
    useIndex = FALSE
  )
  .logThis(mat, paste0("mat ", sampleName), logFile = logFile)
  
  .logDiffTime(paste0("Organizing SE ",useMatrix," from ArrowFile : ", basename(ArrowFile)), 
               t1 = tstart, verbose = verbose, logFile = logFile)
  matrixClass <- h5read(ArrowFile, paste0(useMatrix, "/Info/Class"))
  
  if(matrixClass == "Sparse.Assays.Matrix"){
    rownames(mat) <- paste0(featureDF$name)
    splitIdx <- split(seq_len(nrow(mat)), featureDF$seqnames)
    mat <- lapply(seq_along(splitIdx), function(x){
      mat[splitIdx[[x]], , drop = FALSE]
    }) %>% SimpleList
    names(mat) <- names(splitIdx)
    featureDF <- featureDF[!duplicated(paste0(featureDF$name)), ,drop = FALSE]
    featureDF <- featureDF[,which(colnames(featureDF) %ni% "seqnames"), drop=FALSE]
    rownames(featureDF) <- paste0(featureDF$name)
  }else{
    mat <- SimpleList(mat)
    names(mat) <- useMatrix    
  }
  
  colData <- .getMetadata(ArrowFile)
  colData <- colData[colnames(mat[[1]]),,drop=FALSE]
  
  if(!is.null(ArchRProj)){
    projColData <- getCellColData(ArchRProj)[rownames(colData), ]
    colData <- cbind(colData, projColData[ ,colnames(projColData) %ni% colnames(colData)])
  }
  
  rowData <- tryCatch({
    makeGRangesFromDataFrame(featureDF, keep.extra.columns = TRUE)
  }, error = function(x){
    featureDF
  })
  
  se <- SummarizedExperiment(
    assays = mat,
    rowData = rowData,
    colData = colData
  )
  .logThis(se, paste0("se ", sampleName), logFile = logFile)
  
  se
  
}

.getMatFromArrow <- function(
  ArrowFile = NULL, 
  featureDF = NULL, 
  binarize = NULL, 
  cellNames = NULL,
  useMatrix = "TileMatrix", 
  useIndex = FALSE,
  threads = 1
){
  
  if(is.null(featureDF)){
    featureDF <- .getFeatureDF(ArrowFile, useMatrix)
  }
  
  if(any(c("seqnames","idx") %ni% colnames(featureDF))){
    stop("Need to provide featureDF with columns seqnames and idx!")
  }
  
  #Add RowNames for Check at the end
  rownames(featureDF) <- paste0("f", seq_len(nrow(featureDF)))
  
  o <- h5closeAll()
  
  matClass <- h5read(ArrowFile, paste0(useMatrix,"/Info/Class"))
  if(matClass %ni% c("Sparse.Binary.Matrix", "Sparse.Integer.Matrix", "Sparse.Double.Matrix", "Sparse.Assays.Matrix")){
    stop("Arrow Mat is not a valid Sparse Matrix!")
  }
  if(is.null(binarize)){
    if(matClass == "Sparse.Binary.Matrix"){
      binarize <- TRUE
    }else{
      binarize <- FALSE
    }
  }
  if(matClass == "Sparse.Binary.Matrix"){
    if(!binarize){
      stop("Sparse Matrix in Arrow is Binarized! Set binarize = TRUE to use matrix!")
    }
  }
  
  matColNames <- paste0(.sampleName(ArrowFile), "#", h5read(ArrowFile, paste0(useMatrix,"/Info/CellNames")))
  if(!is.null(cellNames)){
    idxCols <- which(matColNames %in% cellNames)
  }else{
    idxCols <- seq_along(matColNames)
  }
  
  seqnames <- unique(featureDF$seqnames)
  
  mat <- .safelapply(seq_along(seqnames), function(x){
    
    seqnamex <- seqnames[x]
    featureDFx <- featureDF[BiocGenerics::which(featureDF$seqnames %bcin% seqnamex),]
    idxRows <- featureDFx$idx
    
    j <- Rle(
      values = h5read(ArrowFile, paste0(useMatrix,"/",seqnamex,"/jValues")), 
      lengths = h5read(ArrowFile, paste0(useMatrix,"/",seqnamex,"/jLengths"))
    )
    
    #Match J
    matchJ <- S4Vectors::match(j, idxCols, nomatch = 0)
    idxJ <- BiocGenerics::which(matchJ > 0)
    if(useIndex){
      i <- h5read(ArrowFile, paste0(useMatrix,"/",seqnamex,"/i"), index = list(idxJ, 1))
    }else{
      i <- h5read(ArrowFile, paste0(useMatrix,"/",seqnamex,"/i"))[idxJ]
    }
    j <- matchJ[idxJ]
    
    #Match I
    matchI <- match(i, idxRows, nomatch = 0)
    idxI <- which(matchI > 0)
    i <- i[idxI]
    j <- j[idxI]
    i <- matchI[idxI]
    
    if(!binarize){
      x <- h5read(ArrowFile, paste0(useMatrix,"/",seqnamex,"/x"))[idxJ][idxI]
    }else{
      x <- rep(1, length(j))
    }
    
    mat <- Matrix::sparseMatrix(
      i=as.vector(i),
      j=j,
      x=x,
      dims = c(length(idxRows), length(idxCols))
    )
    rownames(mat) <- rownames(featureDFx)
    
    rm(matchI, idxI, matchJ, idxJ, featureDFx, idxRows)
    
    return(mat)
    
  }, threads = threads) %>% Reduce("rbind", .)
  
  o <- h5closeAll()
  
  colnames(mat) <- matColNames[idxCols]
  
  #Double Check Order!
  mat <- mat[rownames(featureDF), , drop = FALSE]
  rownames(mat) <- NULL
  
  if(!is.null(cellNames)){
    mat <- mat[,cellNames,drop=FALSE]
  }
  
  return(mat)
  
}

####################################################################
# Helper read functioning
####################################################################
.getGroupMatrix <- function(
  ArrowFiles = NULL, 
  featureDF = NULL, 
  groupList = NULL,
  threads = 1, 
  useIndex = FALSE, 
  verbose = TRUE, 
  useMatrix = "TileMatrix",
  asSparse = FALSE,
  tstart = NULL
){
  
  #########################################
  # Time Info
  #########################################
  if(is.null(tstart)){
    tstart <- Sys.time()
  }
  
  #########################################
  # Construct Matrix
  #########################################
  seqnames <- unique(featureDF$seqnames)
  rownames(featureDF) <- paste0("f", seq_len(nrow(featureDF)))
  cellNames <- unlist(groupList, use.names = FALSE) ### UNIQUE here? doublet check QQQ
  
  allCellsList <- lapply(seq_along(ArrowFiles), function(x){
    allCells <- .availableCells(ArrowFile = ArrowFiles[x], subGroup = useMatrix)
    allCells <- allCells[allCells %in% cellNames]
    if(length(allCells) != 0){
      allCells
    }else{
      NULL
    }
  })
  
  mat <- .safelapply(seq_along(seqnames), function(x){
    
    .logDiffTime(sprintf("Constructing Group Matrix %s of %s", x, length(seqnames)), tstart, verbose = verbose)
    
    #Construct Matrix
    seqnamex <- seqnames[x]
    featureDFx <- featureDF[BiocGenerics::which(featureDF$seqnames %bcin% seqnamex), ]
    
    matChr <- matrix(0, nrow = nrow(featureDFx), ncol = length(groupList))
    colnames(matChr) <- names(groupList)
    rownames(matChr) <- rownames(featureDFx)
    
    for(y in seq_along(ArrowFiles)){
      
      allCells <- allCellsList[[y]]
      
      if(!is.null(allCells)){
        
        maty <- .getMatFromArrow(
          ArrowFile = ArrowFiles[y], 
          useMatrix = useMatrix,
          featureDF = featureDFx, 
          cellNames = allCells, 
          useIndex = useIndex
        )
        
        for(z in seq_along(groupList)){
          
          #Check Cells In Group
          cellsGroupz <- groupList[[z]]
          idx <- BiocGenerics::which(colnames(maty) %in% cellsGroupz)
          
          #If In Group RowSums
          if(length(idx) > 0){
            matChr[,z] <- matChr[,z] + Matrix::rowSums(maty[,idx,drop=FALSE])
          }
          
        }
        
        rm(maty)
        
      }
      
      
      if(y %% 20 == 0 | y %% length(ArrowFiles) == 0){
        gc()
      } 
      
    }
    
    if(asSparse){
      matChr <- as(matChr, "dgCMatrix")
    }
    
    .logDiffTime(sprintf("Finished Group Matrix %s of %s", x, length(seqnames)), tstart, verbose = verbose)
    
    matChr
    
  }, threads = threads) %>% Reduce("rbind", .)
  
  mat <- mat[rownames(featureDF), , drop = FALSE]
  
  .logDiffTime("Successfully Created Group Matrix", tstart, verbose = verbose)
  
  gc()
  
  return(mat)
  
}

.getPartialMatrix <- function(
  ArrowFiles = NULL, 
  featureDF = NULL, 
  cellNames = NULL, 
  progress = TRUE, 
  threads = 1, 
  useMatrix = "TileMatrix",
  doSampleCells = FALSE, 
  sampledCellNames = NULL, 
  tmpPath = .tempfile(pattern = paste0("tmp-partial-mat")), 
  useIndex = FALSE,
  tstart = NULL,
  verbose = TRUE
){
  
  #########################################
  # Time Info
  #########################################
  if(is.null(tstart)){
    tstart <- Sys.time()
  }
  
  #########################################
  # Construct Matrix
  #########################################
  
  mat <- .safelapply(seq_along(ArrowFiles), function(x){
    
    .logDiffTime(sprintf("Getting Partial Matrix %s of %s", x, length(ArrowFiles)), tstart, verbose = verbose)
    
    allCells <- .availableCells(ArrowFile = ArrowFiles[x], subGroup = useMatrix)
    allCells <- allCells[allCells %in% cellNames]
    
    if(length(allCells) == 0){
      if(doSampleCells){
        return(list(mat = NULL, out = NULL))
      }else{
        return(NULL)
      }
    }
    
    o <- h5closeAll()
    matx <- .getMatFromArrow(
      ArrowFile = ArrowFiles[x], 
      featureDF = featureDF, 
      cellNames = allCells,
      useMatrix = useMatrix, 
      useIndex = useIndex
    )
    
    if(doSampleCells){
      
      #Save Temporary Matrix
      outx <- paste0(tmpPath, "-", .sampleName(ArrowFiles[x]), ".rds")
      saveRDS(matx, outx, compress = FALSE)     
      
      #Sample Matrix 
      matx <- matx[, which(colnames(matx) %in% sampledCellNames),drop = FALSE]
      
      return(list(mat = matx, out = outx))
      
    }else{
      
      return(matx)
      
    }
    
  }, threads = threads)
  
  gc()
  
  
  if(doSampleCells){
    
    matFiles <- lapply(mat, function(x) x[[2]]) %>% Reduce("c", .)
    mat <- lapply(mat, function(x) x[[1]]) %>% Reduce("cbind", .)
    mat <- mat[,sampledCellNames]
    
    .logDiffTime("Successfully Created Partial Matrix", tstart, verbose = verbose)
    
    return(list(mat = mat, matFiles = matFiles))
    
  }else{
    
    mat <- Reduce("cbind", mat)
    mat <- mat[,cellNames]
    
    .logDiffTime("Successfully Created Partial Matrix", tstart, verbose = verbose)
    
    return(mat)
    
  }
  
  
}

########################################################################
# Compute Summary Statistics!
########################################################################

.getRowSums <- function(
  ArrowFiles = NULL,
  useMatrix = NULL,
  seqnames = NULL,
  verbose = TRUE,
  tstart = NULL,
  filter0 = FALSE,
  threads = 1,
  addInfo = FALSE
){
  
  if(is.null(tstart)){
    tstart <- Sys.time()
  }
  
  if(is.null(seqnames)){
    seqnames <- .availableSeqnames(ArrowFiles, useMatrix)
  }
  
  #Compute RowSums
  summaryDF <- .safelapply(seq_along(seqnames), function(x){
    o <- h5closeAll()
    for(y in seq_along(ArrowFiles)){
      if(y == 1){
        sumy <- h5read(ArrowFiles[y], paste0(useMatrix, "/", seqnames[x], "/rowSums"))
      }else{
        sumy1 <- h5read(ArrowFiles[y], paste0(useMatrix, "/", seqnames[x], "/rowSums"))
        if(length(sumy1) != length(sumy)){
          stop("rowSums lengths do not match in ArrowFiles for a seqname!")
        }else{
          sumy <- sumy + sumy1
        }
      }
    }
    #Return Setup In Feature DF Format (seqnames, idx columns)
    DataFrame(seqnames = Rle(seqnames[x], lengths = length(sumy)), idx = seq_along(sumy), rowSums = as.vector(sumy))
  }, threads = threads) %>% Reduce("rbind", .)
  
  if(addInfo){
    featureDF <- .getFeatureDF(ArrowFiles, useMatrix)
    rownames(featureDF) <- paste0(featureDF$seqnames, "_", featureDF$idx)
    rownames(summaryDF) <- paste0(summaryDF$seqnames, "_", summaryDF$idx)
    featureDF <- featureDF[rownames(summaryDF), , drop = FALSE]
    featureDF$rowSums <- summaryDF[rownames(featureDF), "rowSums"]
    summaryDF <- featureDF
    rownames(summaryDF) <- NULL
    remove(featureDF)
  }
  
  if(filter0){
    summaryDF <- summaryDF[which(summaryDF$rowSums > 0), ,drop = FALSE]
  }
  
  return(summaryDF)
  
}

.getRowVars <- function(
  ArrowFiles = NULL,
  seqnames = NULL,
  useMatrix = NULL,
  threads = 1
){
  
  .combineVariances <- function(dfMeans = NULL, dfVars = NULL, ns = NULL){
    
    #https://rdrr.io/cran/fishmethods/src/R/combinevar.R
    
    if(ncol(dfMeans) != ncol(dfVars) | ncol(dfMeans) != length(ns)){
      stop("Means Variances and Ns lengths not identical")
    }
    
    combinedMeans <- rowSums(t(t(dfMeans) * ns)) / sum(ns)
    summedVars <- rowSums(t(t(dfVars) * (ns - 1)) + t(t(dfMeans^2) * ns))
    combinedVars <- (summedVars - sum(ns)*combinedMeans^2)/(sum(ns)-1)
    
    data.frame(combinedVars = combinedVars, combinedMeans = combinedMeans)
    
  }
  
  featureDF <- .getFeatureDF(ArrowFiles, useMatrix)
  
  if(!is.null(seqnames)){
    featureDF <- featureDF[BiocGenerics::which(featureDF$seqnames %bcin% seqnames),]
  }
  
  rownames(featureDF) <- paste0("f", seq_len(nrow(featureDF)))
  fnames <- rownames(featureDF)
  
  featureDF <- S4Vectors::split(featureDF, as.character(featureDF$seqnames))
  
  ns <- lapply(seq_along(ArrowFiles), function(y){
    length(.availableCells(ArrowFiles[y], useMatrix))
  }) %>% unlist
  
  #Compute RowVars
  summaryDF <- .safelapply(seq_along(featureDF), function(x){
    
    o <- h5closeAll()
    seqx <- names(featureDF)[x]
    meanx <- matrix(NA, ncol = length(ArrowFiles), nrow = nrow(featureDF[[x]]))
    varx <- matrix(NA, ncol = length(ArrowFiles), nrow = nrow(featureDF[[x]]))
    
    for(y in seq_along(ArrowFiles)){
      meanx[, y] <- h5read(ArrowFiles[y], paste0(useMatrix, "/", seqx, "/rowMeans"))
      varx[, y] <- h5read(ArrowFiles[y], paste0(useMatrix, "/", seqx, "/rowVars"))
    }
    
    cbind(featureDF[[x]], DataFrame(.combineVariances(meanx, varx, ns)))
    
  }, threads = threads) %>% Reduce("rbind", .)
  
  summaryDF <- summaryDF[fnames, , drop = FALSE]
  
  return(summaryDF)
  
}

.getColSums <- function(
  ArrowFiles = NULL,
  seqnames = NULL,
  useMatrix = NULL,
  verbose = TRUE,
  tstart = NULL,
  threads = 1
){
  
  if(is.null(tstart)){
    tstart <- Sys.time()
  }
  
  #Compute ColSums
  cS <- .safelapply(seq_along(seqnames), function(x){
    
    lapply(seq_along(ArrowFiles), function(y){
      
      o <- h5closeAll()
      cSy <- h5read(ArrowFiles[y], paste0(useMatrix, "/", seqnames[x], "/colSums"))
      rownames(cSy) <- .availableCells(ArrowFiles[y], useMatrix)
      cSy
      
    }) %>% Reduce("rbind", .)
    
  }, threads = threads) %>% Reduce("cbind", .) %>% rowSums
  
  .logDiffTime("Successfully Computed colSums", tstart, verbose = verbose)
  
  return(cS)
  
}

# h5read implementation for optimal reading
.h5read <- function(
  file = NULL,
  name = NULL,
  method = "fast",
  index = NULL,
  start = NULL,
  block = NULL,
  count = NULL
){
  
  if(tolower(method) == "fast" & is.null(index) & is.null(start) & is.null(block) & is.null(count)){
    fid <- H5Fopen(file)
    dapl <- H5Pcreate("H5P_DATASET_ACCESS")
    did <- .Call("_H5Dopen", fid@ID, name, dapl@ID, PACKAGE='rhdf5')
    res <- .Call("_H5Dread", did, NULL, NULL, NULL, TRUE, 0L, FALSE, fid@native, PACKAGE='rhdf5')
    invisible(.Call("_H5Dclose", did, PACKAGE='rhdf5'))   
  }else{
    res <- h5read(file = file, name = name, index = index, start = start, block = block, count = count)
  }
  o <- h5closeAll()
  return(res)
}

.getMatrixClass <- function(
  ArrowFiles = NULL, 
  useMatrix = NULL,
  threads = getArchRThreads()
){
  
  threads <- min(length(ArrowFiles), threads)
  
  matrixClass <- .safelapply(seq_along(ArrowFiles), function(i){
    h5read(ArrowFiles[i], paste0(useMatrix, "/Info/Class"))
  }, threads = threads) %>% unlist %>% unique
  
  if(length(matrixClass) != 1){
    stop("Not all matrix classes are the same!")
  }
  
  matrixClass
  
}

.getMatrixUnits <- function(
  ArrowFiles = NULL, 
  useMatrix = NULL,
  threads = getArchRThreads()
){
  
  threads <- min(length(ArrowFiles), threads)
  
  matrixUnits <- .safelapply(seq_along(ArrowFiles), function(i){
    tryCatch({ #This handles backwards compatibility!
      h5read(ArrowFiles[i], paste0(useMatrix, "/Info/Units"))
    }, error = function(x){
      "None"
    })
  }, threads = threads) %>% unlist %>% unique
  
  if(length(matrixUnits) != 1){
    stop("Not all matrix units are the same!")
  }
  
  matrixUnits
  
}

##########################################################################################
# Assay Correlation Methods
##########################################################################################

#' Correlate Matrices within an ArchRProject
#' 
#' This function will correlate 2 matrices within an ArchRProject by name matching.
#' 
#' @param ArchRProj An `ArchRProject` object.
#' @param useMatrix1 A character describing the first matrix to use. See `getAvailableMatrices` for valid options.
#' @param useMatrix2 A character describing the second matrix to use. See `getAvailableMatrices` for valid options.
#' @param useSeqnames1 A character vector describing which seqnames to use in matrix 1.
#' @param useSeqnames2 A character vector describing which seqnames to use in matrix 2.
#' @param removeFromName1 A character vector describing how to filter names in matrix 1. 
#' Options include "underscore", "dash", "numeric" and "dot". The string portion prior to these will be kept.
#' @param removeFromName2 A character vector describing how to filter names in matrix 2. 
#' Options include "underscore", "dash", "numeric" and "dot". The string portion prior to these will be kept.
#' @param log2Norm1 A boolean describing whether to log2 normalize matrix 1.
#' @param log2Norm2 A boolean describing whether to log2 normalize matrix 2.
#' @param reducedDims The name of the `reducedDims` object (i.e. "IterativeLSI") to use from the designated `ArchRProject`.
#' @param dimsToUse A vector containing the dimensions from the `reducedDims` object to use in computing the embedding.
#' @param scaleDims A boolean value that indicates whether to z-score the reduced dimensions for each cell. This is useful for minimizing
#' the contribution of strong biases (dominating early PCs) and lowly abundant populations. However, this may lead to stronger sample-specific
#' biases since it is over-weighting latent PCs. If set to `NULL` this will scale the dimensions based on the value of `scaleDims` when the
#' `reducedDims` were originally created during dimensionality reduction. This idea was introduced by Timothy Stuart.
#' @param corCutOff A numeric cutoff for the correlation of each dimension to the sequencing depth. If the dimension has a correlation to
#' sequencing depth that is greater than the `corCutOff`, it will be excluded from analysis.
#' @param k The number of k-nearest neighbors to use for creating single-cell groups for correlation analyses.
#' @param knnIteration The number of k-nearest neighbor groupings to test for passing the supplied `overlapCutoff`.
#' @param overlapCutoff The maximum allowable overlap between the current group and all previous groups to permit the current group be
#' added to the group list during k-nearest neighbor calculations.
#' @param seed A number to be used as the seed for random number generation required in knn determination. It is recommended to keep track
#' of the seed used so that you can reproduce results downstream.
#' @param threads The number of threads to be used for parallel computing.
#' @param verbose A boolean value that determines whether standard output should be printed.
#' @param logFile The path to a file to be used for logging ArchR output.
#' @export
correlateMatrices <- function(
  ArchRProj = NULL,
  useMatrix1 = NULL,
  useMatrix2 = NULL,
  useSeqnames1 = NULL,
  useSeqnames2 = NULL,
  removeFromName1 = c("underscore", "dash"),
  removeFromName2 = c("underscore", "dash"),
  log2Norm1 = TRUE,
  log2Norm2 = TRUE,
  reducedDims = "IterativeLSI",
  dimsToUse = 1:30,
  scaleDims = NULL,
  corCutOff = 0.75,
  k = 100, 
  knnIteration = 500, 
  overlapCutoff = 0.8, 
  seed = 1, 
  threads = getArchRThreads(),
  verbose = TRUE,
  logFile = createLogFile("correlateMatrices")
){
  
  .validInput(input = ArchRProj, name = "ArchRProj", valid = c("ArchRProj"))
  .validInput(input = useMatrix1, name = "useMatrix1", valid = c("character"))
  .validInput(input = useMatrix2, name = "useMatrix2", valid = c("character"))
  .validInput(input = useSeqnames1, name = "useSeqnames1", valid = c("character", "null"))
  .validInput(input = useSeqnames2, name = "useSeqnames2", valid = c("character", "null"))
  .validInput(input = removeFromName1, name = "removeFromName1", valid = c("character", "null"))
  .validInput(input = removeFromName2, name = "removeFromName2", valid = c("character", "null"))
  .validInput(input = log2Norm1, name = "log2Norm1", valid = c("boolean"))
  .validInput(input = log2Norm2, name = "log2Norm2", valid = c("boolean"))
  .validInput(input = reducedDims, name = "reducedDims", valid = c("character"))
  .validInput(input = dimsToUse, name = "dimsToUse", valid = c("integer", "null"))
  .validInput(input = scaleDims, name = "scaleDims", valid = c("boolean", "null"))
  .validInput(input = corCutOff, name = "corCutOff", valid = c("numeric", "null"))
  .validInput(input = k, name = "k", valid = c("integer"))
  .validInput(input = knnIteration, name = "knnIteration", valid = c("integer"))
  .validInput(input = overlapCutoff, name = "overlapCutoff", valid = c("numeric"))
  .validInput(input = seed, name = "seed", valid = c("integer"))
  .validInput(input = threads, name = "threads", valid = c("integer"))
  .validInput(input = verbose, name = "verbose", valid = c("boolean"))
  .validInput(input = logFile, name = "logFile", valid = c("character"))
  
  tstart <- Sys.time()
  .startLogging(logFile = logFile)
  .logThis(mget(names(formals()),sys.frame(sys.nframe())), "correlateMatrices Input-Parameters", logFile = logFile)
  
  set.seed(seed)
  
  #Get Available Matrices
  matrixNames <- getAvailableMatrices(ArchRProj)
  
  if(useMatrix1 %ni% matrixNames){
    .logStop(paste0("useMatrix1 (",useMatrix1,") not in availableMatrices :\n", paste(matrixNames, collapse = ", ")), logFile = logFile)
  }
  
  if(useMatrix2 %ni% matrixNames){
    .logStop(paste0("useMatrix2 (",useMatrix2,") not in availableMatrices :\n", paste(matrixNames, collapse = ", ")), logFile = logFile)
  }
  
  #Get Matrix Classes
  matrixClass1 <- as.character(h5read(getArrowFiles(ArchRProj)[1], paste0(useMatrix1, "/Info/Class")))
  matrixClass2 <- as.character(h5read(getArrowFiles(ArchRProj)[1], paste0(useMatrix2, "/Info/Class")))
  .logThis(matrixClass1, name = "matrixClass1", logFile = logFile)
  .logThis(matrixClass2, name = "matrixClass2", logFile = logFile)
  
  #Get Feature DFs
  featureDF1 <- .getFeatureDF(getArrowFiles(ArchRProj), useMatrix1)
  featureDF2 <- .getFeatureDF(getArrowFiles(ArchRProj), useMatrix2)
  .logThis(featureDF1, name = "featureDF1", logFile = logFile)
  .logThis(featureDF2, name = "featureDF2", logFile = logFile)
  
  #Check Seqnames
  featureDF1 <- .checkSeqnames(featureDF1, useMatrix1, useSeqnames1, matrixClass1, logFile)
  featureDF2 <- .checkSeqnames(featureDF2, useMatrix2, useSeqnames2, matrixClass2, logFile)
  
  #Create Match Names
  featureDF1$matchName <- featureDF1$name
  if("underscore" %in% tolower(removeFromName1)){
    featureDF1$matchName <- gsub("\\_.*","",featureDF1$matchName)
  }
  if("dash" %in% tolower(removeFromName1)){
    featureDF1$matchName <- gsub("\\-.*","",featureDF1$matchName)
  }
  if("numeric" %in% tolower(removeFromName1)){
    featureDF1$matchName <- gsub("[0-9]+","",featureDF1$matchName)
  }
  if("dot" %in% tolower(removeFromName1)){
    featureDF1$matchName <- gsub("\\..*","",featureDF1$matchName)
  }
  
  featureDF2$matchName <- featureDF2$name
  if("underscore" %in% tolower(removeFromName2)){
    featureDF2$matchName <- gsub("\\_.*","",featureDF2$matchName)
  }
  if("dash" %in% tolower(removeFromName2)){
    featureDF2$matchName <- gsub("\\-.*","",featureDF2$matchName)
  }
  if("numeric" %in% tolower(removeFromName2)){
    featureDF2$matchName <- gsub("[0-9]+","",featureDF2$matchName)
  }
  if("dot" %in% tolower(removeFromName2)){
    featureDF2$matchName <- gsub("\\..*","",featureDF2$matchName)
  }
  
  .logThis(featureDF1, name = "featureDF1", logFile = logFile)
  .logThis(featureDF2, name = "featureDF2", logFile = logFile)
  
  #Now Lets see how many matched pairings
  matchP1 <- sum(featureDF1$matchName %in% featureDF2$matchName) / nrow(featureDF1)
  matchP2 <- sum(featureDF2$matchName %in% featureDF1$matchName) / nrow(featureDF2)
  matchP <- max(matchP1, matchP2)
  
  .logThis(featureDF1$matchName, "featureDF1$matchName", logFile)
  .logThis(featureDF2$matchName, "featureDF2$matchName", logFile)
  
  if(sum(featureDF1$matchName %in% featureDF2$matchName) == 0){
    .logStop("Matching of useMatrix1 and useMatrix2 resulted in no mappings!", logFile = logFile)
  }
  if(matchP < 0.05){
    if(force){
      .logStop("Matching of useMatrix1 and useMatrix2 resulted in less than 5% mappings! Set force = TRUE to continue!", logFile = logFile)
    }else{
      .logMessage("Matching of useMatrix1 and useMatrix2 resulted in less than 5% mappings! Continuing since force = TRUE.", verbose = TRUE, logFile = logFile)
    }
  }
  matchedNames <- intersect(featureDF1$matchName, featureDF2$matchName)
  featureDF1m <- featureDF1[featureDF1$matchName %in% matchedNames, ]
  featureDF2m <- featureDF2[featureDF2$matchName %in% matchedNames, ]
  
  #Create Mappings
  mappingDF <- lapply(seq_len(nrow(featureDF1m)), function(x){
    expand.grid(x, which(paste0(featureDF2m$matchName) %in% paste0(featureDF1m$matchName[x])))
  }) %>% Reduce("rbind", .)
  
  #Test Mappings
  .logDiffTime(main=paste0("Testing ", nrow(mappingDF), " Mappings!"), t1=tstart, verbose=verbose, logFile=logFile)
  
  #Get Reduced Dims
  rD <- getReducedDims(ArchRProj, reducedDims = reducedDims, corCutOff = corCutOff, dimsToUse = dimsToUse)
  
  #Subsample
  idx <- sample(seq_len(nrow(rD)), knnIteration, replace = !nrow(rD) >= knnIteration)
  
  #KNN Matrix
  .logDiffTime(main="Computing KNN", t1=tstart, verbose=verbose, logFile=logFile)
  knnObj <- .computeKNN(data = rD, query = rD[idx,], k = k)
  
  #Determin Overlap
  .logDiffTime(main="Identifying Non-Overlapping KNN pairs", t1=tstart, verbose=verbose, logFile=logFile)
  keepKnn <- determineOverlapCpp(knnObj, floor(overlapCutoff * k))
  
  #Keep Above Cutoff
  knnObj <- knnObj[keepKnn==0,]
  .logDiffTime(main=paste0("Identified ", nrow(knnObj), " Groupings!"), t1=tstart, verbose=verbose, logFile=logFile)
  .logThis(knnObj, name = "knnObj", logFile = logFile)
  
  #Convert To Names List
  knnObj <- lapply(seq_len(nrow(knnObj)), function(x){
    rownames(rD)[knnObj[x, ]]
  }) %>% SimpleList
  .logThis(knnObj, name = "knnObjList", logFile = logFile)
  
  #Get Group Matrices
  .logDiffTime(main="Getting Group Matrix 1", t1=tstart, verbose=verbose, logFile=logFile)
  groupMat1 <- .getGroupMatrix(
    ArrowFiles = getArrowFiles(ArchRProj), 
    featureDF = featureDF1m, 
    groupList = knnObj, 
    useMatrix = useMatrix1,
    threads = threads,
    verbose = FALSE
  )
  
  .logDiffTime(main="Getting Group Matrix 2", t1=tstart, verbose=verbose, logFile=logFile)
  groupMat2 <- .getGroupMatrix(
    ArrowFiles = getArrowFiles(ArchRProj), 
    featureDF = featureDF2m, 
    groupList = knnObj, 
    useMatrix = useMatrix2,
    threads = threads,
    verbose = FALSE
  )
  
  .logThis(groupMat1, name = "groupMat1", logFile = logFile)
  .logThis(groupMat2, name = "groupMat2", logFile = logFile)
  
  #We need to divide by number of cells for the mean
  groupMat1 <- t(t(groupMat1) / k)
  groupMat2 <- t(t(groupMat2) / k)
  
  .logThis(groupMat1, name = "groupMat1", logFile = logFile)
  .logThis(groupMat2, name = "groupMat2", logFile = logFile)
  
  #Now we can normalize
  if(log2Norm1){
    if(any(groupMat1 < 0)){
      .logMessage("Some entries in groupMat1 are less than 0, continuing without Log2 Normalization.\nMost likely this assay is a deviations matrix.", logFile=logFile)
    }else{
      groupMat1 <- log2(groupMat1 + 1)
    }
  }
  if(log2Norm2){
    if(any(groupMat2 < 0)){
      .logMessage("Some entries in groupMat2 are less than 0, continuing without Log2 Normalization.\nMost likely this assay is a deviations matrix.", logFile=logFile)
    }else{
      groupMat2 <- log2(groupMat2 + 1)
    }
  }
  
  .logThis(groupMat1, name = "groupMat1", logFile = logFile)
  .logThis(groupMat2, name = "groupMat2", logFile = logFile)
  
  #Row Correlate
  rowTest <- .rowCorTest(
    X = groupMat1,
    Y = groupMat2,
    idxX = mappingDF[,1],
    idxY = mappingDF[,2],
    verbose = verbose,
    padjMethod = "bonferroni",
    logFile = logFile
  )
  .logThis(rowTest, name = "rowTest", logFile = logFile)
  
  #Output DF
  colnames(featureDF1m) <- paste0(useMatrix1, "_", colnames(featureDF1m))
  colnames(featureDF2m) <- paste0(useMatrix2, "_", colnames(featureDF2m))
  
  df <- DataFrame(
    cbind(rowTest, featureDF1m[mappingDF[,1],,drop=FALSE], featureDF2m[mappingDF[,2],,drop=FALSE])
  )
  
  frontOrder <- c(paste0(useMatrix1, "_name"), paste0(useMatrix2, "_name"), "cor", "padj", "pval")
  df <- df[, c(frontOrder, colnames(df)[colnames(df) %ni% frontOrder])]
  
  .endLogging(logFile = logFile)
  
  return(df)
  
}

.rowCorTest <- function(
  X = NULL, 
  Y = NULL, 
  idxX = seq_row(X), 
  idxY = seq_row(Y), 
  padjMethod = "BH", 
  min = 10, 
  use="complete", 
  verbose = TRUE,
  threads = 1,
  logFile = NULL
){
  .logMessage("Getting Correlations...", verbose = verbose, logFile=logFile)
  corTestList <- .safelapply(seq_along(idxX), function(i){
    if(i %% 250 == 0){
      .logMessage("Computing Correlation (",i," of ", length(idxX), ")", logFile=logFile)
    }
    if(length(which(!is.na(X[idxX[i],]))) > min && length(which(!is.na(Y[idxY[i],]))) > min){
      corx <- .suppressAll(cor.test(X[idxX[i],],Y[idxY[i],],use=use))
      return(list(cor=corx$estimate[1], pval=corx$p.value[1]))
    }else{
      return(list(cor=NA, pval=NA))
    }
  }, threads = threads)
  corTest <- data.frame(
    cor = lapply(corTestList, function(x) x[[1]]) %>% unlist,
    pval = lapply(corTestList, function(x) x[[2]]) %>% unlist
  )
  rm(corTestList)
  corTest$padj <- p.adjust(corTest$pval, method=padjMethod)
  return(corTest)
}

.checkSeqnames <- function(
  featureDF = NULL, 
  useMatrix = NULL, 
  useSeqnames = NULL, 
  matrixClass = NULL,
  logFile = NULL
){
  
  seqnames <- unique(as.vector(featureDF$seqnames))
  useSeqnames <- useSeqnames[useSeqnames %in% seqnames]
  if(length(useSeqnames)==0){
    useSeqnames <- NULL
  }
  
  if(!is.null(useSeqnames)){
    if(length(useSeqnames) == 1){
      featureDF <- featureDF[BiocGenerics::which(featureDF$seqnames %bcin% useSeqnames),]
    }else{
      .logMessage("When accessing features from a matrix of class Sparse.Assays.Matrix it requires 1 seqname!\n",
                  "Continuing with first seqname '", seqnames[1], "'!\n",
                  "If confused, try getFeatures(ArchRProj, '", useMatrix,"') to list out available seqnames for input!", logFile=logFile)
      useSeqnames <- seqnames[1]
      featureDF <- featureDF[BiocGenerics::which(featureDF$seqnames %bcin% useSeqnames),]
    }
  }else{
    if(matrixClass == "Sparse.Assays.Matrix"){
      if(all(seqnames %in% c("deviations", "z"))){
        seqnames <- c("z", "deviations")
      }
      .logMessage("When accessing features from a matrix of class Sparse.Assays.Matrix it requires 1 seqname!\n",
                  "Continuing with first seqname '", seqnames[1], "'!\n",
                  "If confused, try getFeatures(ArchRProj, '", useMatrix,"') to list out available seqnames for input!", logFile=logFile)
      useSeqnames <- seqnames[1]
      featureDF <- featureDF[BiocGenerics::which(featureDF$seqnames %bcin% useSeqnames),]
    }
  }
  if(!(nrow(featureDF) > 1)){
    .logStop("Less than 1 feature is remaining in featureDF please check input!", logFile=logFile)
  }
  
  featureDF
  
}


#' Correlate Trajectories
#' 
#' This function will correlate 2 trajectory matrices from getTrajectory.
#' 
#' @param seTrajectory1 A `SummarizedExperiment` object that results from calling `getTrajectory()`.
#' @param seTrajectory2 A `SummarizedExperiment` object that results from calling `getTrajectory()`.
#' @param corCutOff A numeric describing the cutoff for determining correlated features.
#' @param varCutOff1 The "Variance Quantile Cutoff" to be used for identifying the top variable features across `seTrajectory1`.
#' Only features with a variance above the provided quantile will be retained.
#' @param varCutOff2 The "Variance Quantile Cutoff" to be used for identifying the top variable features across `seTrajectory2`.
#' Only features with a variance above the provided quantile will be retained.
#' @param removeFromName1 A character vector describing how to filter names in matrix 1. 
#' Options include "underscore", "dash", "numeric" and "dot". The string portion prior to these will be kept.
#' @param removeFromName2 A character vector describing how to filter names in matrix 2. 
#' Options include "underscore", "dash", "numeric" and "dot". The string portion prior to these will be kept.
#' @param useRanges A boolean describing whether to use range overlap matching for correlation analysis.
#' @param fix1 A character describing where to resize the coordinates of `seTrajectory1`. Options include "start", "center", "end".
#' @param fix2 A character describing where to resize the coordinates of `seTrajectory2`. Options include "start", "center", "end".
#' @param maxDist A integer specifying the maximum distance between the coordinates of `seTrajectory1` and `seTrajectory2` for 
#' computing correlations.
#' @param log2Norm1 A boolean describing whether to log2 normalize `seTrajectory1`.
#' @param log2Norm2 A boolean describing whether to log2 normalize `seTrajectory2`.
#' @param force A boolean value that determines whether analysis should continue if resizing coordinates in `seTrajectory1` or 
#' `seTrajectory2` does not align with the strandedness. Only when `useRanges = TRUE`.
#' @param logFile The path to a file to be used for logging ArchR output.
#' @export
correlateTrajectories <- function(
  seTrajectory1 = NULL,
  seTrajectory2 = NULL,
  corCutOff = 0.5,
  varCutOff1 = 0.8,
  varCutOff2 = 0.8,
  removeFromName1 = c("underscore", "dash"),
  removeFromName2 = c("underscore", "dash"),
  useRanges = FALSE,
  fix1 = "center",
  fix2 = "start",
  maxDist = 250000,
  log2Norm1 = TRUE,
  log2Norm2 = TRUE,
  force = FALSE,
  logFile = createLogFile("correlateTrajectories")
){
  
  .validInput(input = seTrajectory1, name = "seTrajectory1", valid = c("SummarizedExperiment"))
  .validInput(input = seTrajectory2, name = "seTrajectory2", valid = c("SummarizedExperiment"))
  .validInput(input = corCutOff, name = "corCutOff", valid = c("numeric"))
  .validInput(input = varCutOff1, name = "varCutOff1", valid = c("numeric"))
  .validInput(input = varCutOff2, name = "varCutOff2", valid = c("numeric")) 
  .validInput(input = removeFromName1, name = "removeFromName1", valid = c("character", "null"))
  .validInput(input = removeFromName2, name = "removeFromName2", valid = c("character", "null"))
  .validInput(input = useRanges, name = "useRanges", valid = c("boolean"))
  .validInput(input = fix1, name = "fix1", valid = c("character"))
  .validInput(input = fix2, name = "fix2", valid = c("character"))
  .validInput(input = maxDist, name = "maxDist", valid = c("integer"))
  .validInput(input = log2Norm1, name = "log2Norm1", valid = c("boolean"))
  .validInput(input = log2Norm2, name = "log2Norm2", valid = c("boolean"))
  .validInput(input = force, name = "force", valid = c("boolean"))
  .validInput(input = logFile, name = "logFile", valid = c("character"))
  
  .startLogging(logFile = logFile)
  .logThis(mget(names(formals()),sys.frame(sys.nframe())), "correlateTrajectories Input-Parameters", logFile=logFile)
  
  featureDF1 <- rowData(seTrajectory1)
  featureDF2 <- rowData(seTrajectory2)
  
  .logThis(featureDF1, "featureDF1", logFile = logFile)
  .logThis(featureDF2, "featureDF2", logFile = logFile)
  
  if("name" %in% colnames(featureDF1)){
    rownames(featureDF1) <- paste0(featureDF1$seqnames, ":", featureDF1$name)
    rownames(seTrajectory1) <- paste0(featureDF1$seqnames, ":", featureDF1$name)
  }else{
    if(!useRanges){
      .logStop("seTrajectory1 does not have a name column in rowData. This means most likely the matching format needs useRanges = TRUE!", logFile = logFile)
    }
    rownames(featureDF1) <- paste0(featureDF1$seqnames, ":", featureDF1$start, "_", featureDF1$end)
    rownames(seTrajectory1) <- paste0(featureDF1$seqnames, ":", featureDF1$start, "_", featureDF1$end)
  }
  
  if("name" %in% colnames(featureDF2)){
    rownames(featureDF2) <- paste0(featureDF2$seqnames, ":", featureDF2$name)
    rownames(seTrajectory2) <- paste0(featureDF2$seqnames, ":", featureDF2$name)
  }else{
    if(!useRanges){
      .logStop("seTrajectory2 does not have a name column in rowData. This means most likely the matching format needs useRanges = TRUE!", logFile = logFile)
    }
    rownames(featureDF2) <- paste0(featureDF2$seqnames, ":", featureDF2$start, "_", featureDF2$end)
    rownames(seTrajectory2) <- paste0(featureDF2$seqnames, ":", featureDF2$start, "_", featureDF2$end)
  }
  
  .logThis(rownames(featureDF1), "rownames(featureDF1)", logFile = logFile)
  .logThis(rownames(featureDF2), "rownames(featureDF2)", logFile = logFile)
  
  if(useRanges){
    
    if("start" %ni% colnames(featureDF1)){
      .logStop("start is not in seTrajectory1, this is not a ranges object. Please set useRanges = FALSE", logFile = logFile)
    }
    
    if("start" %ni% colnames(featureDF2)){
      .logStop("start is not in seTrajectory2, this is not a ranges object. Please set useRanges = FALSE", logFile = logFile)
    }
    
    if("strand" %in% colnames(featureDF1)){
      ranges1 <- GRanges(
        seqnames = featureDF1$seqnames, 
        IRanges(
          ifelse(featureDF1$strand == 2, featureDF1$end, featureDF1$start),
          ifelse(featureDF1$strand == 2, featureDF1$start, featureDF1$end)
        ),
        strand = ifelse(featureDF1$strand == 2, "-", "+")
      )
    }else{
      ranges1 <- GRanges(featureDF1$seqnames, IRanges(featureDF1$start, featureDF1$end))
    }
    #mcols(ranges1) <- featureDF1
    names(ranges1) <- rownames(featureDF1)
    rowRanges(seTrajectory1) <- ranges1
    rm(ranges1)
    
    if("strand" %in% colnames(featureDF2)){
      ranges2 <- GRanges(
        seqnames = featureDF2$seqnames, 
        IRanges(
          ifelse(featureDF2$strand == 2, featureDF2$end, featureDF2$start),
          ifelse(featureDF2$strand == 2, featureDF2$start, featureDF2$end)
        ),
        strand = ifelse(featureDF2$strand == 2, "-", "+")
      )
    }else{
      ranges2 <- GRanges(featureDF2$seqnames, IRanges(featureDF2$start, featureDF2$end))
    }
    #mcols(ranges2) <- featureDF2
    names(ranges2) <- rownames(featureDF2)
    rowRanges(seTrajectory2) <- ranges2
    rm(ranges2)
    
    .logThis(ranges1, "ranges1", logFile = logFile)
    .logThis(ranges2, "ranges2", logFile = logFile)
    
    #Find Associations to test
    isStranded1 <- any(as.integer(strand(seTrajectory1)) == 2)
    isStranded2 <- any(as.integer(strand(seTrajectory2)) == 2)
    
    if(fix1 == "center" & isStranded1){
      if(!force){
        .logStop("fix1 equals center when there is strandedness. Most likely you want this as fix1='start' or fix1='end'. Set force = TRUE to bypass this.", logFile = logFile)
      }else{
        .logMessage("fix1 equals center when there is strandedness. Most likely you want this as fix1='start' or fix1='end'. Continuing since force = TRUE", verbose = TRUE, logFile=logFile)
      }         
    }
    
    if(fix1 != "center" & !isStranded1){
      if(!force){
        .logStop("fix1 does not equal center when there is no strandedness. Most likely you want this as fix1='center'. Set force = TRUE to bypass this.", logFile = logFile)
      }else{
        .logMessage("fix1 does not equal center when there is no strandedness. Most likely you want this as fix1='center'. Continuing since force = TRUE", verbose = TRUE, logFile=logFile)
      }         
    }
    
    if(fix2 == "center" & isStranded2){
      if(!force){
        .logStop("fix2 equals center when there is strandedness. Most likely you want this as fix1='start' or fix1='end'. Set force = TRUE to bypass this.", logFile = logFile)
      }else{
        .logMessage("fix2 equals center when there is strandedness. Most likely you want this as fix1='start' or fix1='end'. Continuing since force = TRUE", verbose = TRUE, logFile=logFile)
      }         
    }
    
    if(fix2 != "center" & !isStranded2){
      if(!force){
        .logStop("fix2 does not equal center when there is no strandedness. Most likely you want this as fix1='center'. Set force = TRUE to bypass this.", logFile = logFile)
      }else{
        .logMessage("fix2 does not equal center when there is no strandedness. Most likely you want this as fix1='center'. Continuing since force = TRUE", verbose = TRUE, logFile=logFile)
      }         
    }
    
    #Overlaps
    mappingDF <- DataFrame(
      findOverlaps( 
        resize(rowRanges(seTrajectory1), 1, fix1), 
        .suppressAll(resize(resize(seTrajectory2, 1, fix2), 2 * maxDist + 1, "center")),
        ignore.strand = TRUE
      )
    )
    
    #Get Distance 
    mappingDF$distance <- distance(
      x = ranges(rowRanges(seTrajectory1)[mappingDF[,1]]), 
      y = ranges(rowRanges(resize(seTrajectory2, 1, fix2))[mappingDF[,2]])
    )
    
  }else{
    
    #Create Match Names
    featureDF1$matchName <- featureDF1$name
    if("underscore" %in% tolower(removeFromName1)){
      featureDF1$matchName <- gsub("\\_.*","",featureDF1$matchName)
    }
    if("dash" %in% tolower(removeFromName1)){
      featureDF1$matchName <- gsub("\\-.*","",featureDF1$matchName)
    }
    if("numeric" %in% tolower(removeFromName1)){
      featureDF1$matchName <- gsub("[0-9]+","",featureDF1$matchName)
    }
    if("dot" %in% tolower(removeFromName1)){
      featureDF1$matchName <- gsub("\\..*","",featureDF1$matchName)
    }
    
    featureDF2$matchName <- featureDF2$name
    if("underscore" %in% tolower(removeFromName2)){
      featureDF2$matchName <- gsub("\\_.*","",featureDF2$matchName)
    }
    if("dash" %in% tolower(removeFromName2)){
      featureDF2$matchName <- gsub("\\-.*","",featureDF2$matchName)
    }
    if("numeric" %in% tolower(removeFromName2)){
      featureDF2$matchName <- gsub("[0-9]+","",featureDF2$matchName)
    }
    if("dot" %in% tolower(removeFromName2)){
      featureDF2$matchName <- gsub("\\..*","",featureDF2$matchName)
    }  
    
    #Now Lets see how many matched pairings
    matchP1 <- sum(featureDF1$matchName %in% featureDF2$matchName) / nrow(featureDF1)
    matchP2 <- sum(featureDF2$matchName %in% featureDF1$matchName) / nrow(featureDF2)
    matchP <- max(matchP1, matchP2)
    
    .logThis(featureDF1$matchName, "featureDF1$matchName", logFile)
    .logThis(featureDF2$matchName, "featureDF2$matchName", logFile)
    
    if(sum(featureDF1$matchName %in% featureDF2$matchName) == 0){
      .logMessage("Matching of seTrajectory1 and seTrajectory2 resulted in no mappings!", logFile = logFile)
      stop("Matching of seTrajectory1 and seTrajectory2 resulted in no mappings!")
    }
    if(matchP < 0.05){
      if(!force){
        .logStop("Matching of seTrajectory1 and seTrajectory2 resulted in less than 5% mappings! Set force = TRUE to continue!", logFile=logFile)
      }else{
        .logMessage("Matching of seTrajectory1 and seTrajectory2 resulted in less than 5% mappings! Continuing since force = TRUE.", verbose=TRUE, logFile=logFile)
      }
    }
    
    #Create Mappings
    mappingDF <- lapply(seq_len(nrow(featureDF1)), function(x){
      idx <- which(paste0(featureDF2$matchName) %in% paste0(featureDF1$matchName[x]))
      if(length(idx) > 0){
        expand.grid(x, idx)
      }else{
        NULL
      }
    }) %>% Reduce("rbind", .)
    
  }
  
  colnames(mappingDF)[1:2] <- c("idx1", "idx2")
  mappingDF <- DataFrame(mappingDF)
  
  .logThis(mappingDF, "mappingDF", logFile = logFile)
  
  if(!useRanges){
    mappingDF$matchname1 <- featureDF1$matchName[mappingDF$idx1]
    mappingDF$matchname2 <- featureDF2$matchName[mappingDF$idx2]
    mappingDF$name1 <- rownames(featureDF1)[mappingDF$idx1]
    mappingDF$name2 <- rownames(featureDF2)[mappingDF$idx2]
  }
  
  mappingDF$Correlation <- rowCorCpp(
    idxX = as.integer(mappingDF[,1]), 
    idxY = as.integer(mappingDF[,2]), 
    X = assays(seTrajectory1)[["mat"]], 
    Y = assays(seTrajectory2)[["mat"]]
  )
  mappingDF$VarAssay1 <- .getQuantiles(matrixStats::rowVars(assays(seTrajectory1)[["mat"]]))[as.integer(mappingDF[,1])]
  mappingDF$VarAssay2 <- .getQuantiles(matrixStats::rowVars(assays(seTrajectory2)[["mat"]]))[as.integer(mappingDF[,2])]
  mappingDF$TStat <- (mappingDF$Correlation / sqrt((1-mappingDF$Correlation^2)/(ncol(seTrajectory1)-2))) #T-statistic P-value
  mappingDF$Pval <- 2 * pt(-abs(mappingDF$TStat), ncol(seTrajectory1) - 2)
  mappingDF$FDR <- p.adjust(mappingDF$Pval, method = "fdr")
  
  idxPF <- which(mappingDF$Correlation > corCutOff & mappingDF$VarAssay1 > varCutOff1 & mappingDF$VarAssay2 > varCutOff2)
  .logMessage("Found ", length(idxPF), " Correlated Pairings!", logFile=logFile, verbose=TRUE)
  
  .logThis(mappingDF[idxPF,], "mappingDF-PF", logFile = logFile)
  
  out <- SimpleList(
    correlatedMappings = mappingDF[idxPF,],
    allMappings = mappingDF,
    seTrajectory1 = seTrajectory1,
    seTrajectory2 = seTrajectory2
  )
  
  out
  
}

##########################################################################################
# Co-accessibility Methods
##########################################################################################

#' Add Peak Co-Accessibility to an ArchRProject
#' 
#' This function will add co-accessibility scores to peaks in a given ArchRProject
#'
#' @param ArchRProj An `ArchRProject` object.
#' @param reducedDims The name of the `reducedDims` object (i.e. "IterativeLSI") to retrieve from the designated `ArchRProject`.
#' @param dimsToUse A vector containing the dimensions from the `reducedDims` object to use in clustering.
#' @param scaleDims A boolean value that indicates whether to z-score the reduced dimensions for each cell. This is useful for minimizing
#' the contribution of strong biases (dominating early PCs) and lowly abundant populations. However, this may lead to stronger sample-specific
#' biases since it is over-weighting latent PCs. If set to `NULL` this will scale the dimensions based on the value of `scaleDims` when the
#' `reducedDims` were originally created during dimensionality reduction. This idea was introduced by Timothy Stuart.
#' @param corCutOff A numeric cutoff for the correlation of each dimension to the sequencing depth. If the dimension has a correlation to
#' sequencing depth that is greater than the `corCutOff`, it will be excluded from analysis.
#' @param cellsToUse A character vector of cellNames to compute coAccessibility on if desired to run on a subset of the total cells.
#' @param k The number of k-nearest neighbors to use for creating single-cell groups for correlation analyses.
#' @param knnIteration The number of k-nearest neighbor groupings to test for passing the supplied `overlapCutoff`.
#' @param overlapCutoff The maximum allowable overlap between the current group and all previous groups to permit the current group be
#' added to the group list during k-nearest neighbor calculations.
#' @param maxDist The maximum allowable distance in basepairs between two peaks to consider for co-accessibility.
#' @param scaleTo The total insertion counts from the designated group of single cells is summed across all relevant peak regions from
#' the `peakSet` of the `ArchRProject` and normalized to the total depth provided by `scaleTo`.
#' @param log2Norm A boolean value indicating whether to log2 transform the single-cell groups prior to computing co-accessibility correlations.
#' @param seed A number to be used as the seed for random number generation required in knn determination. It is recommended to keep track
#' of the seed used so that you can reproduce results downstream.
#' @param threads The number of threads to be used for parallel computing.
#' @param verbose A boolean value that determines whether standard output should be printed.
#' @param logFile The path to a file to be used for logging ArchR output.
#' @export
addCoAccessibility <- function(
  ArchRProj = NULL,
  reducedDims = "IterativeLSI",
  dimsToUse = 1:30,
  scaleDims = NULL,
  corCutOff = 0.75,
  cellsToUse = NULL,
  k = 100, 
  knnIteration = 500, 
  overlapCutoff = 0.8, 
  maxDist = 100000,
  scaleTo = 10^4,
  log2Norm = TRUE,
  seed = 1, 
  threads = getArchRThreads(),
  verbose = TRUE,
  logFile = createLogFile("addCoAccessibility")
){
  
  .validInput(input = ArchRProj, name = "ArchRProj", valid = c("ArchRProj"))
  .validInput(input = reducedDims, name = "reducedDims", valid = c("character"))
  .validInput(input = dimsToUse, name = "dimsToUse", valid = c("numeric", "null"))
  .validInput(input = scaleDims, name = "scaleDims", valid = c("boolean", "null"))
  .validInput(input = corCutOff, name = "corCutOff", valid = c("numeric", "null"))
  .validInput(input = cellsToUse, name = "cellsToUse", valid = c("character", "null"))
  .validInput(input = k, name = "k", valid = c("integer"))
  .validInput(input = knnIteration, name = "knnIteration", valid = c("integer"))
  .validInput(input = overlapCutoff, name = "overlapCutoff", valid = c("numeric"))
  .validInput(input = maxDist, name = "maxDist", valid = c("integer"))
  .validInput(input = scaleTo, name = "scaleTo", valid = c("numeric"))
  .validInput(input = log2Norm, name = "log2Norm", valid = c("boolean"))
  .validInput(input = threads, name = "threads", valid = c("integer"))
  .validInput(input = verbose, name = "verbose", valid = c("boolean"))
  .validInput(input = logFile, name = "logFile", valid = c("character"))
  
  tstart <- Sys.time()
  .startLogging(logFile = logFile)
  .logThis(mget(names(formals()),sys.frame(sys.nframe())), "addCoAccessibility Input-Parameters", logFile = logFile)
  
  set.seed(seed)
  
  #Get Peak Set
  peakSet <- getPeakSet(ArchRProj)
  
  #Get Reduced Dims
  rD <- getReducedDims(ArchRProj, reducedDims = reducedDims, corCutOff = corCutOff, dimsToUse = dimsToUse)
  if(!is.null(cellsToUse)){
    rD <- rD[cellsToUse, ,drop=FALSE]
  }
  
  #Subsample
  idx <- sample(seq_len(nrow(rD)), knnIteration, replace = !nrow(rD) >= knnIteration)
  
  #KNN Matrix
  .logDiffTime(main="Computing KNN", t1=tstart, verbose=verbose, logFile=logFile)
  knnObj <- .computeKNN(data = rD, query = rD[idx,], k = k)
  
  #Determin Overlap
  .logDiffTime(main="Identifying Non-Overlapping KNN pairs", t1=tstart, verbose=verbose, logFile=logFile)
  keepKnn <- determineOverlapCpp(knnObj, floor(overlapCutoff * k))
  
  #Keep Above Cutoff
  knnObj <- knnObj[keepKnn==0,]
  .logDiffTime(paste0("Identified ", nrow(knnObj), " Groupings!"), t1=tstart, verbose=verbose, logFile=logFile)
  
  #Convert To Names List
  knnObj <- lapply(seq_len(nrow(knnObj)), function(x){
    rownames(rD)[knnObj[x, ]]
  }) %>% SimpleList
  
  #Check Chromosomes
  chri <- gtools::mixedsort(.availableChr(getArrowFiles(ArchRProj), subGroup = "PeakMatrix"))
  chrj <- gtools::mixedsort(unique(paste0(seqnames(getPeakSet(ArchRProj)))))
  stopifnot(identical(chri,chrj))
  
  #Create Ranges
  peakSummits <- resize(peakSet, 1, "center")
  peakWindows <- resize(peakSummits, maxDist, "center")
  
  #Create Pairwise Things to Test
  o <- DataFrame(findOverlaps(peakSummits, peakWindows, ignore.strand = TRUE))
  o <- o[o[,1] != o[,2],]
  o$seqnames <- seqnames(peakSet)[o[,1]]
  o$idx1 <- peakSet$idx[o[,1]]
  o$idx2 <- peakSet$idx[o[,2]]
  o$correlation <- -999
  
  #Peak Matrix ColSums
  cS <- .getColSums(getArrowFiles(ArchRProj), chri, verbose = FALSE, useMatrix = "PeakMatrix")
  gS <- unlist(lapply(seq_along(knnObj), function(x) sum(cS[knnObj[[x]]], na.rm=TRUE)))
  
  for(x in seq_along(chri)){
    
    .logDiffTime(sprintf("Computing Co-Accessibility %s (%s of %s)", chri[x], x, length(chri)), t1=tstart, verbose=verbose, logFile=logFile)
    
    #Features
    featureDF <- mcols(peakSet)[BiocGenerics::which(seqnames(peakSet) == chri[x]),]
    featureDF$seqnames <- chri[x]
    
    #Group Matrix
    groupMat <- .getGroupMatrix(
      ArrowFiles = getArrowFiles(ArchRProj), 
      featureDF = featureDF, 
      groupList = knnObj, 
      useMatrix = "PeakMatrix",
      threads = threads,
      verbose = FALSE
    )
    
    #Scale
    groupMat <- t(t(groupMat) / gS) * scaleTo
    
    if(log2Norm){
      groupMat <- log2(groupMat + 1)
    }
    
    #Correlations
    idx <- BiocGenerics::which(o$seqnames==chri[x])
    corVals <- rowCorCpp(idxX = o[idx,]$idx1, idxY = o[idx,]$idx2, X = as.matrix(groupMat), Y = as.matrix(groupMat))
    .logThis(head(corVals), paste0("SubsetCorVals-", x), logFile = logFile)
    
    o[idx,]$correlation <- as.numeric(corVals)
    
    .logThis(groupMat, paste0("SubsetGroupMat-", x), logFile = logFile)
    .logThis(o[idx,], paste0("SubsetCoA-", x), logFile = logFile)
    
  }
  
  o$idx1 <- NULL
  o$idx2 <- NULL
  o <- o[!is.na(o$correlation),]
  mcols(peakSet) <- NULL
  o@metadata$peakSet <- peakSet
  
  metadata(ArchRProj@peakSet)$CoAccessibility <- o
  
  .endLogging(logFile = logFile)
  
  ArchRProj
  
}

#' Get the peak co-accessibility from an ArchRProject
#' 
#' This function obtains co-accessibility data from an ArchRProject.
#' 
#' @param ArchRProj An `ArchRProject` object.
#' @param corCutOff A numeric describing the minimum numeric peak-to-peak correlation to return.
#' @param resolution A numeric describing the bp resolution to return loops as. This helps with overplotting of correlated regions.
#' @param returnLoops A boolean indicating to return the co-accessibility signal as a `GRanges` "loops" object designed for use with
#' the `ArchRBrowser()` or as an `ArchRBrowserTrack()`.
#' @export
getCoAccessibility <- function(
  ArchRProj = NULL, 
  corCutOff = 0.5, 
  resolution = 1, 
  returnLoops = TRUE
){
  
  .validInput(input = ArchRProj, name = "ArchRProj", valid = "ArchRProject")
  .validInput(input = corCutOff, name = "corCutOff", valid = "numeric")
  .validInput(input = resolution, name = "resolution", valid = c("integer", "null"))
  .validInput(input = returnLoops, name = "returnLoops", valid = "boolean")
  
  if(is.null(ArchRProj@peakSet)){
    
    return(NULL)
  }
  
  
  if(is.null(metadata(ArchRProj@peakSet)$CoAccessibility)){
    
    return(NULL)
    
  }else{
    
    coA <- metadata(ArchRProj@peakSet)$CoAccessibility
    coA <- coA[coA$correlation >= corCutOff,,drop=FALSE]
    
    if(returnLoops){
      
      peakSummits <- resize(metadata(coA)$peakSet, 1, "center")
      
      if(!is.null(resolution)){
        summitTiles <- floor(start(peakSummits) / resolution) * resolution + floor(resolution / 2)
      }else{
        summitTiles <- start(peakSummits)
      }
      
      loops <- .constructGR(
        seqnames = seqnames(peakSummits[coA[,1]]),
        start = summitTiles[coA[,1]],
        end = summitTiles[coA[,2]]
      )
      mcols(loops)$value <- coA$correlation
      
      loops <- loops[order(mcols(loops)$value, decreasing=TRUE)]
      loops <- unique(loops)
      loops <- loops[width(loops) > 0]
      loops <- sort(sortSeqlevels(loops))
      
      loops <- SimpleList(CoAccessibility = loops)
      
      return(loops)
      
    }else{
      
      return(coA)
      
    }
    
  }
  
}

.constructGR <- function(
  seqnames = NULL, 
  start = NULL, 
  end = NULL, 
  ignoreStrand = TRUE
){
  .validInput(input = seqnames, name = "seqnames", valid = c("character", "rleCharacter"))
  .validInput(input = start, name = "start", valid = c("integer"))
  .validInput(input = end, name = "end", valid = c("integer"))
  .validInput(input = ignoreStrand, name = "ignoreStrand", valid = c("boolean"))
  df <- data.frame(seqnames, start, end)
  idx <- which(df[,2] > df[,3])
  df[idx,2:3] <-  df[idx,3:2]
  if(!ignoreStrand){
    strand <- rep("+",nrow(df))
    strand[idx] <- "-" 
  }else{
    strand <- rep("*",nrow(df))
  }
  gr <- GRanges(df[,1], IRanges(df[,2],df[,3]), strand = strand)
  return(gr)
}


#' @export
peak2GeneHeatmap <- function(...){
  .Deprecated("plotPeak2GeneHeatmap")
  plotPeak2GeneHeatmap(...)
}

#' Plot Peak2Gene Heatmap from an ArchRProject
#' 
#' This function plots side by side heatmaps of linked ATAC and Gene regions from `addPeak2GeneLinks`.
#' 
#' @param ArchRProj An `ArchRProject` object.
#' @param corCutOff A numeric describing the minimum numeric peak-to-gene correlation to return.
#' @param FDRCutOff A numeric describing the maximum numeric peak-to-gene false discovery rate to return.
#' @param varCutOffATAC A numeric describing the minimum variance quantile of the ATAC peak accessibility when selecting links.
#' @param varCutOffRNA A numeric describing the minimum variance quantile of the RNA gene expression when selecting links.
#' @param k An integer describing the number of k-means clusters to group peak-to-gene links prior to plotting heatmaps.
#' @param nPlot An integer describing the maximum number of peak-to-gene links to plot in heatmap.
#' @param limitsATAC An integer describing the maximum number of peak-to-gene links to plot in heatmap.
#' @param limitsRNA An integer describing the maximum number of peak-to-gene links to plot in heatmap.
#' @param groupBy The name of the column in `cellColData` to use for labeling KNN groupings. The maximum group appeared in the KNN groupings is used.
#' @param palGroup A color palette describing the colors in `groupBy`. For example, if groupBy = "Clusters" try paletteDiscrete(ArchRProj$Clusters) for a color palette.
#' @param palATAC A color palette describing the colors to be used for the ATAC heatmap. For example, paletteContinuous("solarExtra").
#' @param palRNA A color palette describing the colors to be used for the RNA heatmap. For example, paletteContinuous("blueYellow").
#' @param verbose A boolean value that determines whether standard output should be printed.
#' @param returnMatrices A boolean value that determines whether the matrices should be returned with kmeans id versus plotting.
#' @param seed A number to be used as the seed for random number generation. It is recommended to keep track of the seed used so that you can
#' reproduce results downstream.
#' @param logFile The path to a file to be used for logging ArchR output.
#' @export
plotPeak2GeneHeatmap <- function(
  ArchRProj = NULL, 
  corCutOff = 0.45, 
  FDRCutOff = 0.0001,
  varCutOffATAC = 0.25,
  varCutOffRNA = 0.25,
  k = 25,
  nPlot = 25000,
  limitsATAC = c(-2, 2),
  limitsRNA = c(-2, 2),
  groupBy = "Clusters",
  palGroup = NULL,
  palATAC = paletteContinuous("solarExtra"),
  palRNA = paletteContinuous("blueYellow"),
  verbose = TRUE,
  returnMatrices = FALSE,
  seed = 1,
  logFile = createLogFile("plotPeak2GeneHeatmap")
){
  
  .validInput(input = ArchRProj, name = "ArchRProj", valid = c("ArchRProj"))
  .validInput(input = corCutOff, name = "corCutOff", valid = c("numeric"))
  .validInput(input = FDRCutOff, name = "FDRCutOff", valid = c("numeric"))
  .validInput(input = varCutOffATAC, name = "varCutOffATAC", valid = "numeric")
  .validInput(input = varCutOffRNA, name = "varCutOffRNA", valid = "numeric")
  .validInput(input = k, name = "k", valid = c("integer"))
  .validInput(input = nPlot, name = "nPlot", valid = c("integer"))
  .validInput(input = limitsATAC, name = "limitsATAC", valid = c("numeric"))
  .validInput(input = limitsRNA, name = "limitsRNA", valid = c("numeric"))
  .validInput(input = groupBy, name = "groupBy", valid = c("character"))
  .validInput(input = palGroup, name = "palGroup", valid = c("palette", "null"))
  .validInput(input = palATAC, name = "palATAC", valid = c("palette", "null"))
  .validInput(input = palRNA, name = "palRNA", valid = c("palette", "null"))
  .validInput(input = verbose, name = "verbose", valid = c("boolean"))
  .validInput(input = returnMatrices, name = "returnMatrices", valid = c("boolean"))
  .validInput(input = seed, name = "seed", valid = c("integer"))
  .validInput(input = logFile, name = "logFile", valid = c("character"))
  
  tstart <- Sys.time()
  .startLogging(logFile = logFile)
  .logThis(mget(names(formals()),sys.frame(sys.nframe())), "peak2GeneHeatmap Input-Parameters", logFile = logFile)
  
  if(is.null(metadata(ArchRProj@peakSet)$Peak2GeneLinks)){
    stop("No Peak2GeneLinks Found! Try addPeak2GeneLinks!")
  }
  
  #########################################
  # Get Inputs
  #########################################
  ccd <- getCellColData(ArchRProj, select = groupBy)
  p2g <- metadata(ArchRProj@peakSet)$Peak2GeneLinks
  p2g <- p2g[which(p2g$Correlation >= corCutOff & p2g$FDR <= FDRCutOff), ,drop=FALSE]
  
  if(!is.null(varCutOffATAC)){
    p2g <- p2g[which(p2g$VarQATAC > varCutOffATAC),]
  }
  
  if(!is.null(varCutOffRNA)){
    p2g <- p2g[which(p2g$VarQRNA > varCutOffRNA),]
  }
  
  if(nrow(p2g) == 0){
    stop("No peak2genelinks found with your cutoffs!")
  }
  
  if(!file.exists(metadata(p2g)$seATAC)){
    stop("seATAC does not exist! Did you change paths? If this does not work, please try re-running addPeak2GeneLinks!")
  }
  if(!file.exists(metadata(p2g)$seRNA)){
    stop("seRNA does not exist! Did you change paths? If this does not work, please try re-running addPeak2GeneLinks!")
  }
  mATAC <- readRDS(metadata(p2g)$seATAC)[p2g$idxATAC, ]
  mRNA <- readRDS(metadata(p2g)$seRNA)[p2g$idxRNA, ]
  p2g$peak <- paste0(rowRanges(mATAC))
  p2g$gene <- rowData(mRNA)$name
  gc()
  
  mATAC <- assay(mATAC)
  mRNA <- assay(mRNA)
  
  #########################################
  # Determine Groups from KNN
  #########################################
  .logDiffTime(main="Determining KNN Groups!", t1=tstart, verbose=verbose, logFile=logFile)
  KNNList <- as(metadata(readRDS(metadata(p2g)$seRNA))$KNNList, "list")
  KNNGroups <- lapply(seq_along(KNNList), function(x){
    KNNx <- KNNList[[x]]
    names(sort(table(ccd[KNNx, 1, drop = TRUE]), decreasing = TRUE))[1]
  }) %>% unlist
  cD <- DataFrame(row.names=paste0("K_", seq_len(ncol(mATAC))), groupBy = KNNGroups)
  pal <- paletteDiscrete(values=gtools::mixedsort(unique(ccd[,1])))
  if(!is.null(palGroup)){
    pal[names(palGroup)[names(palGroup) %in% names(pal)]] <- palGroup[names(palGroup) %in% names(pal)]
  }
  colorMap <- list(groupBy = pal)
  attr(colorMap[[1]], "discrete") <- TRUE
  
  #########################################
  # Organize Matrices
  #########################################
  mATAC <- .rowZscores(mATAC)
  mRNA <- .rowZscores(mRNA)
  rownames(mATAC) <- NULL
  rownames(mRNA) <- NULL
  colnames(mATAC) <- paste0("K_", seq_len(ncol(mATAC)))
  colnames(mRNA) <- paste0("K_", seq_len(ncol(mRNA)))
  rownames(mATAC) <- paste0("P2G_", seq_len(nrow(mATAC)))
  rownames(mRNA) <- paste0("P2G_", seq_len(nrow(mRNA)))
  rownames(p2g) <- paste0("P2G_", seq_len(nrow(p2g)))
  
  .logDiffTime(main="Ordering Peak2Gene Links!", t1=tstart, verbose=verbose, logFile=logFile)
  if(!is.null(seed)){
    set.seed(seed)
  }
  k1 <- kmeans(mATAC, k)
  if(nrow(mATAC) > nPlot){
    nPK <- nPlot * table(k1$cluster) / length(k1$cluster) 
    splitK <- split(seq_len(nrow(mATAC)), k1$cluster)
    kDF <- lapply(seq_along(splitK), function(x){
      idx <- sample(splitK[[x]], floor(nPK[x]))
      k <- rep(x, length(idx))
      DataFrame(k = k, idx = idx)
    }) %>% Reduce("rbind", .)
  }else{
    kDF <- DataFrame(k = k1$cluster, idx = seq_len(nrow(mATAC)))
  }
  bS <- .binarySort(t(.groupMeans(t(mATAC[kDF[,2],]), kDF[,1])),  clusterCols = TRUE, cutOff = 1)
  rowOrder <- rownames(bS[[1]])
  colOrder <- colnames(bS[[1]])
  kDF[,3] <- as.integer(mapLabels(paste0(kDF[,1]), newLabels = paste0(seq_along(rowOrder)), oldLabels = rowOrder))
  
  if(returnMatrices){
    
    out <- SimpleList(
      ATAC = SimpleList(
        matrix = mATAC[kDF[,2],colOrder],
        kmeansId = kDF[,3],
        colData = cD[colOrder,,drop=FALSE]
      ),
      RNA = SimpleList(
        matrix = mRNA[kDF[,2],colOrder],
        kmeansId = kDF[,3],
        colData = cD[colOrder,,drop=FALSE]
      ),
      Peak2GeneLinks = p2g[kDF[,2],]
    )
    
    return(out)
    
  }
  
  
  #########################################
  # Plot Heatmaps
  #########################################
  .logDiffTime(main="Constructing ATAC Heatmap!", t1=tstart, verbose=verbose, logFile=logFile)
  htATAC <- .ArchRHeatmap(
    mat = mATAC[kDF[,2],colOrder],
    scale = FALSE,
    limits = limitsATAC,
    color = palATAC, 
    colData = cD[colOrder,,drop=FALSE],
    colorMap = colorMap,
    clusterCols = FALSE,
    clusterRows = FALSE,
    split = kDF[,3],
    labelRows = FALSE,
    labelCols = FALSE,
    draw = FALSE,
    name = paste0("ATAC Z-Scores\n", nrow(mATAC), " P2GLinks")
  )
  
  .logDiffTime(main = "Constructing RNA Heatmap!", t1 = tstart, verbose = verbose, logFile = logFile)
  htRNA <- .ArchRHeatmap(
    mat = mRNA[kDF[,2],colOrder], 
    scale = FALSE,
    limits = limitsRNA,
    color = palRNA, 
    colData = cD[colOrder,,drop=FALSE],
    colorMap = colorMap,
    clusterCols = FALSE,
    clusterRows = FALSE,
    split = kDF[,3],
    labelRows = FALSE,
    labelCols = FALSE,
    draw = FALSE,
    name = paste0("RNA Z-Scores\n", nrow(mRNA), " P2GLinks")
  )
  
  .endLogging(logFile = logFile)
  
  htATAC + htRNA
  
}

########################################################################################################
# Helpers for Nice Heatmap with Bioconductors ComplexHeamtap
########################################################################################################
.ArchRHeatmap <- function(
  mat = NULL, 
  scale = FALSE,
  limits = c(min(mat), max(mat)),
  colData = NULL, 
  color = paletteContinuous(set = "solarExtra", n = 100),
  clusterCols = TRUE,
  clusterRows = FALSE,
  labelCols = FALSE,
  labelRows = FALSE,
  colorMap = NULL,
  useRaster = TRUE,
  rasterQuality = 5,
  split = NULL,
  fontSizeRows = 10,
  fontSizeCols = 10,
  fontSizeLabels = 8,
  colAnnoPerRow = 4,
  showRowDendrogram = FALSE,
  showColDendrogram = FALSE,
  customRowLabel = NULL,
  customRowLabelIDs = NULL,
  customColLabel = NULL,
  customColLabelIDs = NULL,
  customLabelWidth = 0.75,
  rasterDevice = "png",
  padding = 45,
  borderColor = NA,
  draw = TRUE,
  name = "Heatmap"
){
  
  #Packages
  .requirePackage("ComplexHeatmap", source = "bioc")
  .requirePackage("circlize", source = "cran")
  
  #Z-score
  if (scale) {
    message("Scaling Matrix..")
    mat <- .rowZscores(mat, limit = FALSE)
    name <- paste0(name," Z-Scores")
  }
  
  #Get A Color map if null
  if (is.null(colorMap)) {
    colorMap <- .colorMapAnno(colData)
  }
  
  #Prepare ColorMap format for Complex Heatmap
  if (!is.null(colData)){
    colData = data.frame(colData)
    colorMap <- .colorMapForCH(colorMap, colData) #change
    showLegend <- .checkShowLegend(colorMap[match(names(colorMap), colnames(colData))]) #change
  }else {
    colorMap <- NULL
    showLegend <- NULL
  }
  
  #Prepare Limits if needed
  breaks <- NULL
  if (!is.null(limits)) {
    mat[mat > max(limits)] <- max(limits)
    mat[mat < min(limits)] <- min(limits)
  }else{
    limits <- c(round(min(mat),2), round(max(mat),2))
  }
  
  #Scale Values 0 - 1
  mat <- (mat - min(limits)) / (max(limits) - min(limits))
  breaks <- seq(0, 1, length.out = length(color))
  color <- circlize::colorRamp2(breaks, color)
  
  if(exists('anno_mark', where='package:ComplexHeatmap', mode='function')){
    anno_check_version_rows <- ComplexHeatmap::anno_mark
    anno_check_version_cols <- ComplexHeatmap::anno_mark
  }else{
    anno_check_version_rows <- ComplexHeatmap::row_anno_link
    anno_check_version_cols <- ComplexHeatmap::column_anno_link
  }
  
  #Annotation Heatmap
  if(!is.null(colData) & !is.null(customColLabel)){
    message("Adding Annotations..")
    if(is.null(customColLabelIDs)){
      customColLabelIDs <- colnames(mat)[customColLabel]
    }
    ht1Anno <- HeatmapAnnotation(
      df = colData,
      col = colorMap, 
      show_legend = showLegend,
      show_annotation_name = TRUE,
      gp = gpar(col = "NA"),
      annotation_legend_param =
        list(
          nrow = min(colAnnoPerRow, max(round(nrow(colData)/colAnnoPerRow), 1))
        ),
      foo = anno_check_version_cols(at = customColLabel, labels = customColLabelIDs, labels_gp = gpar(fontsize = fontSizeLabels))
    )
    
  }else if(!is.null(colData)){
    message("Adding Annotations..")
    ht1Anno <- HeatmapAnnotation(
      df = colData,
      col = colorMap, 
      show_legend = showLegend,
      show_annotation_name = TRUE,
      gp = gpar(col = "NA"),
      annotation_legend_param =
        list(
          nrow = min(colAnnoPerRow, max(round(nrow(colData)/colAnnoPerRow), 1))
        )
    )
  }else if(is.null(colData) & !is.null(customColLabel)){
    if(is.null(customColLabelIDs)){
      customColLabelIDs <- colnames(mat)[customColLabel]
    }
    message("Adding Annotations..")
    #print(customColLabel)
    #print(customColLabelIDs)
    #ht1Anno <- columnAnnotation(foo = anno_check_version_cols(
    #   at = customColLabel, labels = customColLabelIDs),
    #   width = unit(customLabelWidth, "cm") + max_text_width(customColLabelIDs))
    #ht1Anno <- HeatmapAnnotation(foo = anno_mark(at = c(1:4, 20, 60, 1097:1100), labels = month.name[1:10]))
    ht1Anno <- HeatmapAnnotation(foo = anno_check_version_cols(at = customColLabel, labels = customColLabelIDs, labels_gp = gpar(fontsize = fontSizeLabels)))
  }else{
    ht1Anno <- NULL
  }
  
  message("Preparing Main Heatmap..")
  ht1 <- Heatmap(
    
    #Main Stuff
    matrix = as.matrix(mat),
    name = name,
    col = color, 
    
    #Heatmap Legend
    heatmap_legend_param = 
      list(
        at = c(0, 1),
        labels = c(round(min(limits),2), round(max(limits),2)),
        color_bar = "continuous", 
        legend_direction = "horizontal",
        legend_width = unit(3, "cm")
      ), 
    rect_gp = gpar(col = borderColor), 
    
    #Column Options
    show_column_names = labelCols,
    cluster_columns = clusterCols, 
    show_column_dend = showColDendrogram,
    clustering_method_columns = "ward.D2",
    column_names_gp = gpar(fontsize = fontSizeCols), 
    column_names_max_height = unit(100, "mm"),
    
    #Row Options
    show_row_names = labelRows,
    row_names_gp = gpar(fontsize = fontSizeRows), 
    cluster_rows = clusterRows, 
    show_row_dend = showRowDendrogram, 
    clustering_method_rows = "ward.D2",
    split = split, 
    
    #Annotation
    top_annotation = ht1Anno, 
    
    #Raster Info
    use_raster = useRaster, 
    raster_device = rasterDevice, 
    raster_quality = rasterQuality
  )
  
  if(!is.null(customRowLabel)){
    if(is.null(customRowLabelIDs)){
      customRowLabelIDs <- rownames(mat)[customRowLabel]
    }
    ht1 <- ht1 + rowAnnotation(link = 
                                 anno_check_version_rows(at = customRowLabel, labels = customRowLabelIDs, labels_gp = gpar(fontsize = fontSizeLabels)),
                               width = unit(customLabelWidth, "cm") + max_text_width(customRowLabelIDs))
  }
  
  if(draw){
    draw(ht1, 
         padding = unit(c(padding, padding, padding, padding), "mm"), 
         heatmap_legend_side = "bot", 
         annotation_legend_side = "bot")
  }else{
    ht1
  }
  
}


.colorMapForCH <- function(colorMap = NULL, colData = NULL){
  colorMap <- colorMap[which(names(colorMap) %in% colnames(colData))]
  colorMapCH <- lapply(seq_along(colorMap), function(x){
    if(attr(colorMap[[x]],"discrete")){
      colorx <- colorMap[[x]]
    }else{
      vals <- colData[[names(colorMap)[x]]][!is.na(colData[[names(colorMap)[x]]])]
      s <-  seq(min(vals), max(vals), length.out = length(colorMap[[x]]))
      colorx <- circlize::colorRamp2(s, colorMap[[x]])
    }
    if(any(is.na(names(colorx)))){
      names(colorx)[is.na(names(colorx))] <- paste0("NA",seq_along(names(colorx)[is.na(names(colorx))]))
    }
    return(colorx)
  })
  names(colorMapCH) <- names(colorMap)
  return(colorMapCH)
}

.checkShowLegend <- function(colorMap = NULL, max_discrete = 30){
  show <- lapply(seq_along(colorMap), function(x){
    if(attr(colorMap[[x]],"discrete") && length(unique(colorMap[[x]])) > max_discrete){
      sl <- FALSE
    }else{
      sl <- TRUE
    }
    return(sl)
  }) %>% unlist
  names(show) <- names(colorMap)
  return(show)
}

.colorMapAnno <- function(colData = NULL, customAnno = NULL, discreteSet = "stallion", continuousSet = "solarExtra"){
  discreteCols <- sapply(colData,function(x) !is.numeric(x))
  if(!is.null(customAnno)){
    colorMap <- lapply(seq_along(discreteCols),function(x){
      if(discreteCols[x]){
        colors <- paletteDiscrete(values = colData[[names(discreteCols[x])]], set = discreteSet)
        names(colors) <- unique(colData[[names(discreteCols[x])]])
        attr(colors, "discrete") <- TRUE
      }else{
        colors <- paletteContinuous(set = continuousSet)
        attr(colors, "discrete") <- FALSE
      }
      if(length(which(customAnno[,1] %in% names(discreteCols[x]))) > 0){
        if(length(which(customAnno[,2] %in% names(colors))) > 0){
          customAnnox <- customAnno[which(customAnno[,2] %in% names(colors)),]
          colors[which(names(colors) %in% customAnnox[,2])] <- paste0(customAnnox[match(names(colors),customAnnox[,2]),3])
        }
      }
      return(colors)
    })
    names(colorMap) <- colnames(colData)
    return(colorMap)
  }else{
    colorMap <- lapply(seq_along(discreteCols), function(x){
      if(discreteCols[x]){
        colors <- paletteDiscrete(values = colData[[names(discreteCols[x])]], set = discreteSet)
        names(colors) <- unique(colData[[names(discreteCols[x])]])
        attr(colors, "discrete") <- TRUE
      }else{
        colors <- paletteContinuous(set = continuousSet)
        attr(colors, "discrete") <- FALSE
      }
      return(colors)
    })
    names(colorMap) <- colnames(colData)
    return(colorMap)
  }
  
}

.binarySort <- function(m = NULL, scale = FALSE, cutOff = 1, lmat = NULL, clusterCols = TRUE){
  
  if(is.null(lmat)){
    #Compute Row-Zscores
    if(scale){
      lmat <- sweep(m - rowMeans(m), 1, matrixStats::rowSds(m), `/`)
    }else{
      lmat <- m
    }
    lmat <- lmat >= cutOff
  }
  
  #Transpose
  m <- t(m)
  lmat <- t(lmat)
  
  #Identify Column Ordering
  if(clusterCols){
    hc <- hclust(dist(m))
    colIdx <- hc$order
    m <- t(m[colIdx,])
    lmat <- t(lmat[colIdx,])
  }else{
    m <- t(m)
    lmat <- t(lmat)
    hc <- NULL
  }
  
  #Identify Row Ordering
  rowIdx <- do.call("order", c(as.data.frame(lmat)[seq_len(ncol(lmat))], list(decreasing = TRUE)))
  m <- t(m[rowIdx,])
  lmat <- t(lmat[rowIdx,])
  
  #Transpose
  m <- t(m)
  lmat <- t(lmat)
  
  return(list(mat = m, hclust = hc))
  
}




