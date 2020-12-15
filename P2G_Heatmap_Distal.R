library(ArchR)
plotPeak2GeneHeatmap.distal <- function(
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
  
  # if(!is.null(varCutOffATAC)){
  #   p2g <- p2g[which(p2g$VarQATAC > varCutOffATAC),]
  # }
  # 
  # if(!is.null(varCutOffRNA)){
  #   p2g <- p2g[which(p2g$VarQRNA > varCutOffRNA),]
  # }
  # 
  # if(nrow(p2g) == 0){
  #   stop("No peak2genelinks found with your cutoffs!")
  # }
  # 
  # if(!file.exists(metadata(p2g)$seATAC)){
  #   stop("seATAC does not exist! Did you change paths? If this does not work, please try re-running addPeak2GeneLinks!")
  # }
  # if(!file.exists(metadata(p2g)$seRNA)){
  #   stop("seRNA does not exist! Did you change paths? If this does not work, please try re-running addPeak2GeneLinks!")
  # }
  mATAC <- readRDS("/datastore/nextgenout5/share/labs/francolab/scENDO_scOVAR_Proj/scATAC-seq_Processing/All_ATAC_v2/All/Peak2GeneLinks/seATAC-Group-KNN.rds")[p2g$idxATAC, ]
  p2g$peak <- paste0(rowRanges(mATAC))
  meta.atac <- rowRanges(mATAC)
  p2g$type <- paste0(meta.atac$peakType)
  
  p2g <- p2g[which(p2g$type == "Distal"), ,drop=FALSE]
  
  mATAC <- SummarizedExperiment::subset(mATAC,peakType == "Distal")
  
  mRNA <- readRDS("/datastore/nextgenout5/share/labs/francolab/scENDO_scOVAR_Proj/scATAC-seq_Processing/All_ATAC_v2/All/Peak2GeneLinks/seRNA-Group-KNN.rds")[p2g$idxRNA, ]
  p2g$gene <- rowData(mRNA)$name
  gc()
  
  mATAC <- assay(mATAC)
  mRNA <- assay(mRNA)
  
  #########################################
  # Determine Groups from KNN
  #########################################
  .logDiffTime(main="Determining KNN Groups!", t1=tstart, verbose=verbose, logFile=logFile)
  KNNList <- as(metadata(readRDS("/datastore/nextgenout5/share/labs/francolab/scENDO_scOVAR_Proj/scATAC-seq_Processing/All_ATAC_v2/All/Peak2GeneLinks/seRNA-Group-KNN.rds"))$KNNList, "list")
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
  
  #Log Info
  .logThis(colorMap, "colorMap", logFile = logFile)
  .logThis(colOrder, "colOrder", logFile = logFile)
  .logThis(kDF, "kDF", logFile = logFile)
  .logThis(mATAC, "mATAC", logFile = logFile)
  .logThis(mRNA, "mRNA", logFile = logFile)
  .logThis(cD[colOrder,,drop=FALSE], "cD", logFile = logFile)
  .logThis(mATAC[kDF[,2],colOrder], "mATAC2", logFile = logFile)
  .logThis(mRNA[kDF[,2],colOrder], "mRNA2", logFile = logFile)
  
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
        cv <- min(abs(c(input%%1, input%%1-1)), na.rm = TRUE) < .Machine$double.eps^0.5
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
        #.safeSaveRDS(errorList, "Save-Error.rds")
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
      .safeSaveRDS(errorList, debugFile)
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
# Safe saveRDS check
##########################################################################################

.safeSaveRDS <- function(
  object = NULL, 
  file = "", 
  ascii = FALSE, 
  version = NULL, 
  compress = TRUE, 
  refhook = NULL
){
  #Try to save a test data.frame in location
  testDF <- data.frame(a=1,b=2)
  canSave <- suppressWarnings(tryCatch({
    saveRDS(object = testDF, file = file, ascii = ascii, version = version, compress = compress, refhook = refhook)
    TRUE
  }, error = function(x){
    FALSE
  }))
  if(!canSave){
    dirExists <- dir.exists(dirname(file))
    if(dirExists){
      stop("Cannot saveRDS. File Path : ", file)
    }else{
      stop("Cannot saveRDS because directory does not exist (",dirname(file),"). File Path : ", file)
    }
  }else{
    saveRDS(object = object, file = file, ascii = ascii, version = version, compress = compress, refhook = refhook)
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

