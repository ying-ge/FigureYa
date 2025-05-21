#' @title Compute kNN conditional Density Resampled Estimate of Mutual Information.
#' @details  Calculates k-Nearest Neighbor conditional Density Resampled Estimate
#' of Mutual Information as defined in Van Dijk et al, 2018. Given a source gene
#' and many target genes, this function will calculate DREMI scores for each
#' source-target pairs.
#' @param imputed.matrix input imputed matrix, rows are cells and columns are genes.
#' @param source.gene gene name.
#' @param target.genes a vector of gene names. If is NULL, it will regard all
#' genes as targets. Default: NULL
#' @param k int, number of neighbors. Default: 10L
#' @param n.bins int, number of bins for density resampling. Default: 20L
#' @param n.mesh int, In each bin, density will be calculated around (mesh ** 2)
#' points. Default: 3L
#' @param n.cores int, number of cores used for calculation. Default: 1
#' @param return.drevi bool, If True, return the DREVI normalized density matrix
#' in addition to the DREMI score. Default: FALSE
CalculateKnnDREMI <- function(imputed.matrix, source.gene, target.genes=NULL, k=10L, n.bins=20L, n.mesh=3L, return.drevi=FALSE, n.cores=1) {
  # check parameters
  if (!exists("imputed.matrix")) {
    stop("The parameter `impute.matrix` must be provided.")
  }
  if (!exists("source.gene")) {
    stop("The parameter `source.gene` must be provided.")
  }
  if (k < 0 || n.bins < 0 || n.mesh < 0) {
    stop("Inappropriate parameter(s) in `k`, `n.bins` or `n.mesh`.")
  }
  if (length(source.gene) > 1) {
    source.gene <- source.gene[0]
    warning(sprintf("Too many source genes were provided. Only use `%s` for calculation.", source.gene))
  }
  bg.genes <- colnames(imputed.matrix)
  if (!source.gene %in% bg.genes) {
    stop(sprintf("`%s` was not recorded in the imputed.matrix", source.gene))
  }
  if (is.null(target.genes)) {
    target.genes <- setdiff(bg.genes, source.gene)
  } else if (length(setdiff(target.genes, bg.genes)) > 0) {
    warning(sprintf("%s genes not in `imputed.matrix` were removed.", length(setdiff(target.genes, bg.genes))))
    target.genes <- intersect(target.genes, bg.genes)
  }

  scprep <- reticulate::import("scprep")
  k <- as.integer(k)
  n.bins <- as.integer(n.bins)
  n.mesh <- as.integer(n.mesh)
  expr.1 <- imputed.matrix[, source.gene]

  if (n.cores == 1) {
    dremi.list <- pbapply::pblapply(
      X = target.genes,
      FUN = function(gn) {
        scprep$stats$knnDREMI(x=expr.1, y=imputed.matrix[, gn], plot=FALSE, k=k, n_bins=n.bins, n_mesh=n.mesh, return_drevi=return.drevi)
      })
  } else {
    dremi.list <- parallel::mclapply(
      X = target.genes,
      FUN = function(gn) {
        scprep$stats$knnDREMI(x=expr.1, y=imputed.matrix[, gn], plot=FALSE, k=k, n_bins=n.bins, n_mesh=n.mesh, return_drevi=return.drevi)
      },
      mc.cores = n.cores
    )
  }

  if (!return.drevi) {
    dremi.score <- data.frame(
      gene_name = target.genes,
      dremi_score = unlist(dremi.list)
    )
    return(list(
      "dremi" = dremi.score
    ))
  } else {
    dremi.score <- data.frame(
      gene_name = target.genes,
      dremi_score = sapply(dremi.list, function(xx) xx[[1]])
    )
    drevi.matirx <- lapply(dremi.list, function(xx) xx[[2]])
    names(drevi.matirx) <- target.genes
    return(list(
      "dremi" = dremi.score,
      "drevi" = drevi.matirx
    ))
  }
}

DREVIPlot <- function(drevi.vector, title = "", bins=20) {
  mm <- matrix(drevi.vector, nrow = bins, byrow = F) %>% as.data.frame()
  colnames(mm) <- 1:bins
  mm <- mm %>% mutate(Y = 1:nrow(mm), .before = 1)
  mm <- mm %>% pivot_longer(cols = 2:ncol(.), names_to = "X", values_to = "density")
  mm$X <- as.integer(mm$X)

  nf <- lm(Y ~ splines::ns(X, df=6), data=mm, weights = mm$density)

  dat <- mutate(mm, smooth=fitted(nf))

  ggplot(dat, aes(x=X, y=Y)) +
    geom_tile(aes(fill=ifelse(density > quantile(density, 0.95), quantile(density, 0.95), density))) +
    geom_line(aes(X, smooth), color="red", lwd=1.5) +
    scale_fill_viridis_c(option="inferno") +
    ggtitle(title) +
    theme_bw(base_size = 15) +
    theme(legend.position = "none",
          plot.title = element_text(hjust = .5, face = "bold"))
}
