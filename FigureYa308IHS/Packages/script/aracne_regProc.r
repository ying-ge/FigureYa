library(optparse)
library(viper)
## load arguments
option_list <- list(
  make_option(c('-a', '--a_file'), type="character", help='ARACNe .tsv.'),
  make_option(c('-e', '--exp_file'), type="character", help='EXP .tsv.'),
  make_option(c('-o', '--out_dir'), type="character", help='Output .tsv.'),
  make_option(c('-n', '--out_name'), type="character", help='Output .tsv.')
)
opt <- parse_args(OptionParser(option_list = option_list))
## load expresison
eset <- read.table(opt$exp_file, header = TRUE, sep = '\t', row.names = 1, stringsAsFactors = FALSE)
eset <- eset[match(unique(rownames(eset)), rownames(eset)),]
eset <- as.matrix(eset)
print(dim(eset))
## process interactome
processed.reg <- aracne2regulon(afile=opt$a_file, eset=eset, format='3col')
saveRDS(processed.reg, file = paste(opt$out_dir, opt$out_name, '_unPruned.rds', sep=''))
pruned.reg <- pruneRegulon(processed.reg, 50, adaptive=FALSE, eliminate=TRUE)
saveRDS(pruned.reg, file = paste(opt$out_dir, opt$out_name, '_pruned.rds', sep=''))
