#' @name twoclasslimma
#' @description perform limma differential analysis for ONLY two conditions
#' @param subtype a data frame with a condition' columns and rownames of samples
#' @param featmat a matrix of features vs samples; FPKM or TPM or other profiling matrix without log2 transformation is recommended.
#' @param treatVar a string value to indicate the name of treatment group
#' @param ctrlVar a string value to indicate the name of control group
#' @param prefix a string value to indicate the prefix of output file
#' @param overwt a logic value to indicate if to overwrite existing results; FALSE by default
#' @param sort.p a logic value to indicate if to sort adjusted p value for output table; TRUE by default
#' @param verbose a logic value to indicate if to only output id, log2fc, pvalue, and padj; TRUE by default.
#' @param res.path a string value to indicate the path to store limma result
#' @return
#' @export
#' @import limma
#' @examples
twoclasslimma <- function(subtype   = NULL,
                          featmat   = NULL,
                          treatVar  = NULL,
                          ctrlVar   = NULL,
                          prefix    = NULL,
                          overwt    = FALSE,
                          sort.p    = TRUE,
                          verbose   = TRUE,
                          res.path  = getwd()) {
  
  library(limma)
  if(!is.element("condition", colnames(subtype))) {
    stop("argument of subtype must contain a column named with 'condition'!")
  }
  # Create comparison list for differential analysis between two classes.
  createList  <- function(subtype = NULL) {
    
    tumorsam <- rownames(subtype)
    sampleList <- list()
    treatsamList <- list()
    treatnameList <- c()
    ctrlnameList <- c()
    
    sampleList[[1]] <- tumorsam
    treatsamList[[1]] <- intersect(tumorsam, rownames(subtype[which(subtype$condition == treatVar),,drop = F]))
    treatnameList[1] <- treatVar
    ctrlnameList[1] <- ctrlVar
    
    return(list(sampleList, treatsamList, treatnameList, ctrlnameList))
  }
  
  complist <- createList(subtype = subtype)
  
  sampleList <- complist[[1]]
  treatsamList <- complist[[2]]
  treatnameList <- complist[[3]]
  ctrlnameList <- complist[[4]]
  allsamples <- colnames(featmat)
  
  # log transformation
  if(max(featmat) < 25 | (max(featmat) >= 25 & min(featmat) < 0)) {
    message("--feature matrix seems to have been standardised (z-score or log transformation), no more action will be performed.")
    gset <- featmat
  }
  if(max(featmat) >= 25 & min(featmat) >= 0){
    message("--log2 transformation done for feature matrix.")
    gset <- log2(featmat + 1)
  }
  
  options(warn = 1)
  for (k in 1:length(sampleList)) {
    samples <- sampleList[[k]]
    treatsam <- treatsamList[[k]]
    treatname <- treatnameList[k]
    ctrlname <- ctrlnameList[k]
    
    compname <- paste(treatname, "_vs_", ctrlname, sep="")
    tmp <- rep("others", times = length(allsamples))
    names(tmp) <- allsamples
    tmp[samples] <- "control"
    tmp[treatsam] <- "treatment"
    
    if(!is.null(prefix)) {
      outfile <- file.path(res.path, paste(prefix, "_limma_test_result.", compname, ".txt", sep = ""))
    } else {
      outfile <- file.path(res.path, paste("limma_test_result.", compname, ".txt", sep = ""))
    }
    
    if (file.exists(outfile) & (overwt == FALSE)) {
      cat(paste0("limma of ",compname, " exists and skipped...\n"))
      next
    }
    
    pd <- data.frame(Samples = names(tmp),
                     Group = as.character(tmp),
                     stringsAsFactors = FALSE)
    
    design <-model.matrix(~ -1 + factor(pd$Group, levels = c("treatment","control")))
    colnames(design) <- c("treatment","control")
    
    fit <- limma::lmFit(gset, design = design);
    contrastsMatrix <- limma::makeContrasts(treatment - control, levels = c("treatment", "control"))
    fit2 <- limma::contrasts.fit(fit, contrasts = contrastsMatrix)
    fit2 <- limma::eBayes(fit2, 0.01)
    resData <- limma::topTable(fit2, adjust = "fdr", sort.by = "B", number = 100000)
    resData <- as.data.frame(subset(resData, select=c("logFC","t","B","P.Value","adj.P.Val")))
    resData$id <- rownames(resData)
    colnames(resData) <- c("log2fc","t","B","pvalue","padj","id")
    resData$fc <- 2^resData$log2fc
    
    if(sort.p) {
      resData <- resData[order(resData$padj),]
    } else {
      resData <- as.data.frame(resData)
    }
    if(verbose) {
      resData <- resData[,c("id","fc","log2fc","pvalue","padj")]
    } else {
      resData <- resData[,c("id","fc","log2fc","t","B","pvalue","padj")]
    }
    write.table(resData, file = outfile, row.names = FALSE, sep = "\t", quote = FALSE)
    cat(paste0("limma of ",compname, " done...\n"))
  }
  options(warn = 0)
}