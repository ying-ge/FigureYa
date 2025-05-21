# Consolidate bootstrap runs

collapse.key = '--aracne-edge--'
tmp <- commandArgs(TRUE)
wd <- tmp[1] # Working directory
expfn <- tmp[2] # expression data
regfn <- tmp[3] # regulators
mht <- tmp[4] # Correction (none,bonferroni, fdr)
alpha <- as.numeric(tmp[5]) # p-value

print (paste("Looking for bootstrapped networks in: ", wd))
fn <- list.files(path=wd, pattern="bootstrapNetwork") # Get bootstrap filenames
res <- resmi <- edges <- NULL
for (fn1 in fn) {
	print (fn1)
    tmp <- strsplit(readLines(file.path(wd, fn1)), "\t")
    mi <- as.numeric(sapply(tmp, function(x) x[3]))
    names(mi) <- sapply(tmp, function(x) paste(x[1:2], collapse=collapse.key))
    edges <- c(edges, length(mi))
    ed <- unique(c(names(res), names(mi)))
    pos1 <- match(ed, names(res))
    pos1na <- which(is.na(pos1))
    pos2 <- match(ed, names(mi))
    pos2na <- which(is.na(pos2))
    res <- res[pos1]
    resmi <- resmi[pos1]
    res[pos1na] <- resmi[pos1na] <- 0
    names(res)[pos1na] <- names(resmi)[pos1na] <- ed[pos1na]
    mi <- mi[pos2]
    mi[pos2na] <- 0
    names(mi)[pos2na] <- ed[pos2na]
    res <- res+sign(mi)
    resmi <- resmi + mi
}
resmi <- resmi/res
regs <- readLines(regfn)
genes <- sapply(strsplit(readLines(expfn)[-1], "\t"), function(x) x[1])
tedges <- length(genes)*length(which(genes %in% regs))
lambda <- sum(edges)/tedges
res <- p.adjust(ppois(res, lambda, lower.tail=F), mht)
res <- sort(res[res < alpha])
## FIXME: possible bug, if more than one '-' in the name (i.e. a single gene has a dash in it)
## this strsplit will break:
## workaround is remove these, but fix by 
tmp <- cbind(t(sapply(strsplit(names(res), collapse.key), function(x) x)), signif(resmi[match(names(res), names(resmi))], 5), signif(res, 5))
output.file = file.path(wd, "finalNetwork_4col.tsv")
print (paste("Writing to output file:", output.file))
cat(unlist(apply(tmp, 1, paste, collapse="\t"), use.names=FALSE), sep="\n", file=output.file)