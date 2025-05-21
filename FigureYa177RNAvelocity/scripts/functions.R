Read10X_V3 <- function (data.dir = NULL) 
{
  full.data <- list()
  for (i in seq_along(data.dir)) {
    run <- data.dir[i]
    if (!dir.exists(run)) {
      stop("Directory provided does not exist")
    }
    if (!grepl("\\/$", run)) {
      run <- paste(run, "/", sep = "")
    }
    barcode.loc <- paste0(run, "barcodes.tsv")
    gene.loc <- paste0(run, "features.tsv")
    matrix.loc <- paste0(run, "matrix.mtx")
    if (!file.exists(barcode.loc)) {
      barcode.loc <- paste0(run, "barcodes.tsv.gz")
      if (!file.exists(barcode.loc)) {
        stop("Barcode file missing")
      }
    }
    if (!file.exists(gene.loc)) {
      gene.loc <- paste0(run, "features.tsv.gz")
      if (!file.exists(gene.loc)) {
        stop("Gene name file missing")  
      }
    }
    if (!file.exists(matrix.loc)) {
      matrix.loc <- paste0(run, "matrix.mtx.gz")
      if (!file.exists(matrix.loc)) {
        stop("Expression matrix file missing")  
      }
    }
    data <- readMM(file = matrix.loc)
    cell.names <- readLines(barcode.loc)
    gene.names <- readLines(gene.loc)
    if (all(grepl(pattern = "\\-1$", x = cell.names))) {
      cell.names <- as.vector(x = as.character(x = sapply(X = cell.names, 
                                                          FUN = ExtractField, field = 1, delim = "-")))
    }
    rownames(x = data) <- make.unique(names = as.character(x = sapply(X = gene.names, 
                                                                      FUN = ExtractField, field = 2, delim = "\\t")))
    if (is.null(x = names(x = data.dir))) {
      if (i < 2) {
        colnames(x = data) <- cell.names
      }
      else {
        colnames(x = data) <- paste0(i, "_", cell.names)
      }
    }
    else {
      colnames(x = data) <- paste0(names(x = data.dir)[i], 
                                   "_", cell.names)
    }
    full.data <- append(x = full.data, values = data)
  }
  full.data <- do.call(cbind, full.data)
  return(full.data)
}
