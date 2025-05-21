eXtremeLogFC <- function(logFC_Matrix, N = 500) {
	logFC_Matrix <- as.data.frame(logFC_Matrix)
	numOfGenes <- dim(logFC_Matrix)[1]
	XLogFC <- sapply(logFC_Matrix, function(logFC) {
		logFCOrder <- order(logFC)
		FC0 <- logFCOrder[(N+1):(numOfGenes-N)]
		logFC[FC0] <- 0L
		return(logFC)
	})
	rownames(XLogFC) <- rownames(logFC_Matrix)
	colnames(XLogFC) <- colnames(logFC_Matrix)
	XLogFC
}



XSum <- function(logFC_Matrix, upGenes, dnGenes) {
	upGeneLogFC <- logFC_Matrix[match(upGenes, rownames(logFC_Matrix), nomatch = 0L), ]
	dnGeneLogFC <- logFC_Matrix[match(dnGenes, rownames(logFC_Matrix), nomatch = 0L), ]
	colSums(upGeneLogFC) - colSums(dnGeneLogFC)
}

