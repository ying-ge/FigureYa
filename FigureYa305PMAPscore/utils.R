neo_get_km_data <- function(mut_sam, gene_Ucox, symbol_Entrez, 
                               sur, path_Ucox_mul, sig, data.dir = NULL, 
                               organism = "hsa", TRAIN = FALSE){
  a <- colnames(mut_sam)
  rownames(mut_sam) <- gsub(pattern = "-", replacement = ".", 
                            rownames(mut_sam))
  mut_sam <- as.matrix(mut_sam[match(rownames(gene_Ucox)[which(gene_Ucox$HR < 1)], rownames(mut_sam)), ])
  mut_sam <- get_Entrez_ID(mut_sam, symbol_Entrez, Entrez_ID = TRUE)
  mut_sample <- as.matrix(mut_sam)
  colnames(mut_sample) <- a
  res = list()
  for (i in colnames(mut_sample)){
    path_matrix <- as.numeric(mut_sample[, i])
    names(path_matrix) <- rownames(mut_sample)
    path_matrix <- path_matrix[which(path_matrix != 0)]
    if (length(path_matrix[which(path_matrix != 0)]) > 0) {
      res[[i]] <- newspia(de = path_matrix[which(path_matrix != 0)], 
                          all = rownames(mut_sample), organism = "hsa", 
                          beta = NULL, verbose = FALSE, data.dir = data.dir)
    }
  }
  .myDataEnv <- new.env(parent = emptyenv())
  datload <- paste(organism, "SPIA", sep = "")
  if (is.null(data.dir)) {
    if (!paste(datload, ".RData", sep = "") %in% dir(system.file("extdata", 
                                                                 package = "PMAPscore"))) {
      cat("The KEGG pathway data for your organism is not present in the extdata folder of the SPIA package!!!")
      cat("\n")
      cat("Please generate one first using makeSPIAdata and specify its location using data.dir argument or copy it in the extdata folder of the SPIA package!")
    }
    else {
      load(file = paste(system.file("extdata", package = "PMAPscore"), 
                        paste("/", organism, "SPIA", sep = ""), ".RData", 
                        sep = ""), envir = .myDataEnv)
    }
  }
  if (!is.null(data.dir)) {
    if (!paste(datload, ".RData", sep = "") %in% dir(data.dir)) {
      cat(paste(data.dir, " does not contin a file called ", 
                paste(datload, ".RData", sep = "")))
    }
    else {
      load(file = paste(data.dir, paste(datload, ".RData", sep = ""), sep = ""), envir = .myDataEnv)
    }
  }
  path.info = .myDataEnv[["path.info"]]
  pathname <- c()
  for (j in 1:length(path.info)) {
    pathname <- c(pathname, path.info[[j]]$title)
  }
  pfs_score <- matrix(data = 0, nrow = length(pathname), ncol = dim(mut_sample)[2])
  rownames(pfs_score) <- pathname
  colnames(pfs_score) <- a
  for (i in colnames(mut_sample)){
    loc2 <- match(res[[i]][, 1], rownames(pfs_score))
    pfs_score[loc2, i] <- res[[i]][, 3]
  }
  km_data <- get_risk_score(sig, pfs_score, path_Ucox_mul, 
                            sur, TRAIN = TRAIN)
  km_data <- as.data.frame(km_data)
  return(km_data)
}


neo_get_sam_cla <- function(km_data, cut_off = -0.986){
  result <- c()
  for (i in 1:dim(km_data)[1]) {
    if (km_data[i, dim(km_data)[2]] > cut_off) {
      newdata <- cbind(rownames(km_data)[i], km_data[i, dim(km_data)[2]], "high")
    }
    else {
      newdata <- cbind(rownames(km_data)[i], km_data[i, dim(km_data)[2]], "low")
    }
    result <- rbind(result, newdata)
  }
  colnames(result) <- c("sample", "risk_score", "class")
  rownames(result) <- result[, 1]
  return(result)
}

getPathGene <- function(path.info, gene_symbol_Entrez){
  pathname <- unlist(lapply(path.info, function(x) x$title))
  path_gene <- lapply(path.info, function(x){
    gene_symbol_Entrez$Hugo_Symbol[match(x$nodes, gene_symbol_Entrez$Entrez_Gene_Id)]
  })
  names(path_gene) = pathname
  return(path_gene)
}