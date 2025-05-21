Seurat.pipeline <- function(RdsName, outputname, mingene = 200, maxgene = 7500, MTfraction =0.2, pc, res) {
  seu_object <- RdsName
  # filter
  mito.features <- grep(pattern = "^MT-", x = rownames(x = seu_object), value = TRUE)
  percent.mito <- Matrix::colSums(x = GetAssayData(object = seu_object, slot = 'counts')[mito.features, ]) / Matrix::colSums(x = GetAssayData(object = seu_object, slot = 'counts'))
  seu_object[['percent.mito']] <- percent.mito
  
  # Save plots as pdf file
  pdf(file = paste(outputname, "_pc", pc, "_res", res, ".pdf", sep = ""), width = 7, height = 5.5)
  
  print(VlnPlot(object = seu_object, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3))
  print(FeatureScatter(object = seu_object, feature1 = "nCount_RNA", feature2 = "percent.mito"))
  print(FeatureScatter(object = seu_object, feature1 = "nCount_RNA", feature2 = "nFeature_RNA"))
  
  seu_object <- subset(x = seu_object, subset = nFeature_RNA > 200 & nFeature_RNA < 7500 & percent.mito < 0.25)
  ## 2 Normalizing the data
  seu_object <- NormalizeData(object = seu_object, normalization.method = "LogNormalize", scale.factor = 1e4)
  ## 3 Detection of variable features across the single cell
  seu_object <- FindVariableFeatures(object = seu_object, selection.method = 'mean.var.plot', mean.cutoff = c(0.0125, 3), 
                                     dispersion.cutoff = c(0.5, Inf))
  length(x = VariableFeatures(object = seu_object))
  
  ### 4 Scaling the data and removing unwanted sources of variation
  seu_object <- ScaleData(object = seu_object, features = rownames(x = seu_object), 
                          vars.to.regress = c("nCount_RNA", "percent.mito"))
  
  ### 5 Perform linear dimensional reduction
  seu_object <- RunPCA(object = seu_object, features = VariableFeatures(object = seu_object), verbose = FALSE)
  
  ##### Determine statistically significant principal components
  print(ElbowPlot(object = seu_object,ndims = 50))
  
  #### cluster cells
  seu_object <- FindNeighbors(object = seu_object, dims = 1:pc)
  seu_object <- FindClusters(object = seu_object, resolution = res)
  
  ### tSNE
  seu_object <- RunTSNE(object = seu_object, dims = 1:pc)
  # note that you can set `label = TRUE` or use the LabelClusters function to help label individual clusters
  print(DimPlot(object = seu_object, reduction = 'tsne', label = TRUE))
  dev.off()
  
  # Save RDS files
  saveRDS(object = seu_object, file = paste(outputname, "_pc", pc, "_res", res, ".rds", sep = ""))
  # find markers for every cluster compared to all remaining cells, report only the positive ones
  seu_object.markers <- FindAllMarkers(object = seu_object, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
  
  #FC>2 -> logFC>1
  markers_use <- seu_object.markers[which(seu_object.markers$avg_logFC >1), ]
  write.csv(seu_object.markers, paste(outputname, "_pc", pc, "_res", res, "_marker.csv", sep = ""))
  write.csv(markers_use, paste(outputname, "m_pc", pc, "_res", res, "_marker_fc2.csv", sep = ""))
  
}
