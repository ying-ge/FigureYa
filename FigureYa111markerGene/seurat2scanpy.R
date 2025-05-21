#
#
seurat2scanpy <- function(x, 
                          type="data",
                          dir_name = "seurat"){
  if ( ! dir.exists(dir_name)){
    dir.create(dir_name)
  }
  
  # write matrix
  mt <- slot(x@assays$RNA, type)
  Matrix::writeMM(mt, file = file.path(dir_name, "matrix.mtx"))
  # write barcode and cell info
  genes <- data.frame(x=row.names(mt),y=row.names(mt))
  write.table(genes, file = file.path(dir_name, "genes.tsv"),
              row.names = F, col.names = F, sep="\t",
              quote = F)
  cells <- colnames(mt)
  write.table(cells, file = file.path(dir_name, "barcodes.tsv"),
              row.names = F, col.names = F,
              quote = F)
  # write meta data
  medata <- as.matrix(slot(x, "meta.data"))
  medata <- as.data.frame(medata, stringsAsFactors = F)
  medata$bulk_labels <- as.character(x@active.ident)
  write.csv(medata, file = file.path(dir_name, "metadata.csv"))
  
}