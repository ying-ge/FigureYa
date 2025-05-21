## libraries
library(optparse)
## arguments
option_list <- list(
  make_option(c('-q', '--quiet'), action='store_false', default=TRUE, help='Suppresses status updates.', dest='print'),
  make_option(c('-r', '--rds'), type="character", help='Input data in .rds format (features X samples).'),
  make_option(c('-o', '--out'), type="character", help='Name of output file; an ARACNe compatible .tsv.')
)
opt <- parse_args(OptionParser(option_list = option_list))

#' Saves a matrix in a format for input to ARACNe
#'
#' @param dat.mat Matrix of data (genes X samples).
#' @param out.file Output file where matrix will be saved.
#' @param subset Switch for subsetting the matrix to 500 samples. Default TRUE.
ARACNeTable <- function(dat.mat, out.file, subset = TRUE) {
  dat.mat <- dat.mat[!duplicated(rownames(dat.mat)), ]
  if (subset) {
    dat.mat <- dat.mat[, sample(colnames(dat.mat), min(ncol(dat.mat), 500)) ]
  }
  sample.names <- colnames(dat.mat)
  gene.ids <- rownames(dat.mat)
  m <- dat.mat
  mm <- rbind( c("gene", sample.names), cbind(gene.ids, m))
  write.table( x = mm , file = out.file , 
               sep="\t", quote = F , row.names = F , col.names = F )
}

## read in data
dat.mat <- readRDS(opt$rds)
ARACNeTable(dat.mat, opt$out)
