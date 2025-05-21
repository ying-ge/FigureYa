normalizeGeneCounts <- function(counts, TxDb, method) {
  require("GenomicFeatures")
  length <- width(genes(TxDb))/10 ^ 3 #gene length per kb
  if (method == "CPM") {
    normalized <-
      as.data.frame(apply(counts,2, function(x) (x/sum(x)) * 10 ^ 6))
  } else if (method == "RPKM") {
    normalized <-
      as.data.frame(t(apply(counts, 1, "/", colSums(counts) / 10 ^ 6)) / length)
  } else if (method == "TPM") {
    normalized <-
      as.data.frame(t(apply(
        counts / length, 1, "/", colSums(counts / length) / 10 ^ 6
        )))
  } else {
    stop("method must be CPM, RPKM or TPM")
  }
    return(normalized)
}
