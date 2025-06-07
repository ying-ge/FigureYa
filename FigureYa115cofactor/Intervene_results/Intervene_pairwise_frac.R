#!/usr/bin/env Rscript
if (suppressMessages(!require("corrplot"))) suppressMessages(install.packages("corrplot", repos="http://cran.us.r-project.org"))
library("corrplot")
pdf("/Users/air/Documents/FigureYa/FigureYa115cofactor/Intervene_results/Intervene_pairwise_frac.pdf", width=8, height=8)

intersection_matrix <- as.matrix(read.table("/Users/air/Documents/FigureYa/FigureYa115cofactor/Intervene_results/Intervene_pairwise_frac_matrix.txt"))
intersection_matrix <- cor(intersection_matrix, method="pearson")
corrplot(intersection_matrix, method ="color", title="Pairwise intersection-pearson correlation(Fraction of overlap)", tl.col="black", tl.cex=0.8, is.corr = TRUE, diag=FALSE, addrect=1, mar=c(0,0,2,1), rect.col = "black")
invisible(dev.off())
