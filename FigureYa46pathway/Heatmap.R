library(pheatmap)

Relative_RPKM <- read.table("Relative_RPKM.txt", header = TRUE, row.names = 1, sep = "\t")

pheatmap(Relative_RPKM, 
         scale = "none",
         color = colorRampPalette(c("#eee814", "#3d4c9e"))(100),
         cluster_row = FALSE, cluster_col = FALSE,
         cellwidth = 9, cellheight = 7,
         border_color = "black")

Absolute_RPKM <- read.table("Absolute_RPKM.txt", header = TRUE, row.names = 1, sep = "\t")

pheatmap(Absolute_RPKM, 
         scale = "none",
         color = colorRampPalette(c("white", "red"))(100),
         cluster_row = FALSE, cluster_col = FALSE,
         cellwidth = 9, cellheight = 7,
         border_color = "black")
