MagicImpute_LowMemory <- function(seu.list, tmp.dir="tmp.data2/magic", ...) {
  if (!dir.exists(tmp.dir)) {
    dir.create(tmp.dir, recursive = T)
  }
  message("MAGIC imputation: low memory mode.")
  for (i in seq_along(seu.list)) {
    seu <- seu.list[[i]]
    seu <- Rmagic::magic(seu, ...)
    imputed.matrix <- t(as.data.frame(seu[["MAGIC_RNA"]]@data))
    imputed.matrix <- as_tibble(imputed.matrix)
    target.file <- file.path(tmp.dir, sprintf("%s.MAGIC_RNA.feather", names(seu.list)[i]))
    arrow::write_feather(imputed.matrix, target.file)
    cells.file <- file.path(tmp.dir, sprintf("%s.MAGIC_RNA.feather.cellnames.txt", names(seu.list)[i]))
    writeLines(colnames(seu), cells.file)
    message(sprintf("[%s|%s] %s done!", i, length(seu.list), names(seu.list)[i]))
  }
}

SplitBucketsByFeatures <- function(data.dir, target.dir, n.buckets=20) {
  if (!dir.exists(target.dir)) {
    dir.create(target.dir, recursive = T)
  }
  message("Generating cell names ...")
  feather.files <- list.files(data.dir, pattern = ".feather$", full.names = T)
  cellname.files <- sapply(feather.files, function(xx) paste0(xx, ".cellnames.txt"))
  cellnames <- pbapply::pblapply(cellname.files, readLines) %>% base::Reduce(c, .)
  writeLines(cellnames, file.path(target.dir, "cellnames.txt"))

  gene.names <- colnames(arrow::read_feather(feather.files[1]))
  genes.buckets <- data.frame(
    row.names = gene.names,
    bucket = sort((1:length(gene.names)-1) %% n.buckets)
  )
  saveRDS(genes.buckets, file.path(target.dir, "genes.buckets.rds"))

  for (i in unique(genes.buckets$bucket)) {
    message(sprintf("Generating bucket %s ...", i))
    genes.used <- rownames(subset(genes.buckets, bucket == i))
    imputed.matrix <- pbapply::pblapply(feather.files, function(xx) {
      tmp.matrix <- arrow::read_feather(xx)
      tmp.matrix <- tmp.matrix[, genes.used]
    }) %>% base::Reduce(rbind, .)
    out.file <- file.path(target.dir, sprintf("bucket_%s.feather", i))
    arrow::write_feather(imputed.matrix, out.file)
  }
}


CalculateKnnDREMI_LowMemory <- function(data.dir, source.gene, target.genes, ...) {
  genes.buckets <- readRDS(file.path(data.dir, "genes.buckets.rds"))
  if (!exists("source.gene")) {
    stop("The parameter `source.gene` must be provided.")
  }
  if (length(source.gene) > 1) {
    source.gene <- source.gene[0]
    warning(sprintf("Too many source genes were provided. Only use `%s` for calculation.", source.gene))
  }
  bg.genes <- rownames(genes.buckets)
  if (!source.gene %in% bg.genes) {
    stop(sprintf("`%s` was not recorded in the imputed.matrix", source.gene))
  }
  if (is.null(target.genes)) {
    target.genes <- setdiff(bg.genes, source.gene)
  } else if (length(setdiff(target.genes, bg.genes)) > 0) {
    warning(sprintf("%s genes not in `imputed.matrix` were removed.", length(setdiff(target.genes, bg.genes))))
    target.genes <- intersect(target.genes, bg.genes)
  }

  source.bucket <- genes.buckets[source.gene, , drop = F]$bucket
  in.file <- file.path(data.dir, sprintf("bucket_%s.feather", source.bucket))
  imputed.matrix <- arrow::read_feather(in.file)
  source.matrix <- imputed.matrix[, source.gene, drop = F]
  genes.buckets <- genes.buckets[target.genes, , drop = F]

  res.list <- lapply(
    X = sort(unique(genes.buckets$bucket)),
    FUN = function(i) {
      in.file <- file.path(data.dir, sprintf("bucket_%s.feather", i))
      genes.used <- rownames(subset(genes.buckets, bucket == i))
      imputed.matrix <- arrow::read_feather(in.file)
      imputed.matrix <- imputed.matrix[, genes.used, drop = F]
      imputed.matrix <- as_tibble(cbind(imputed.matrix, source.matrix))
      message(sprintf("bucket: %s", i))
      message(sprintf("shape: (%s, %s)", dim(imputed.matrix)[1], dim(imputed.matrix)[2]))
      CalculateKnnDREMI(imputed.matrix = imputed.matrix,
                        source.gene = source.gene,
                        target.genes = NULL,
                        ...)
    })
  if (length(res.list[[1]]) == 1) {
    return(list("dremi" = lapply(res.list, function(xx) xx[[1]]) %>% base::Reduce(rbind, .)))
  } else {
    return(list(
      "dremi" = lapply(res.list, function(xx) xx[[1]]) %>% base::Reduce(rbind, .),
      "drevi" = lapply(res.list, function(xx) xx[[2]]) %>% base::Reduce(c, .)
    ))
  }
}
