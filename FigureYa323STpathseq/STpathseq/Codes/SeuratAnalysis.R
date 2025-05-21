suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(Matrix))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(hdf5r))

option_list = list(
  make_option(c("-n", "--samplename"), type = "character",  default = "CRC_16", 
              help="which sample in InputData/sample_sheet.csv to be process", dest = "samplename"),
  make_option(c("-w", "--workpath"), type = "character",  default= getwd(), 
              help="project path, default is current path", dest = "work.path")
)
opt = parse_args(OptionParser(option_list = option_list))
print(opt)

samplename <- opt$samplename
work.path <- opt$work.path

ST.seu <- Load10X_Spatial(file.path("Results", samplename, "outs"))
ST.seu$orig.ident <- samplename
p1 <- VlnPlot(ST.seu, features = "nCount_Spatial", pt.size = 0.1) + NoLegend()
p2 <- SpatialFeaturePlot(ST.seu, features = "nCount_Spatial") + theme(legend.position = "right")
# p1 + p2

ST.seu <- subset(ST.seu, subset = nCount_Spatial > 3 & nFeature_Spatial > 3)
# SpatialFeaturePlot(ST.seu, features = "nCount_Spatial") + theme(legend.position = "right")
ST.seu <- SCTransform(ST.seu, assay = "Spatial")

umi_table <- read.csv(sprintf("Results/%s/outs/%s.visium.raw_matrix.genus.csv", samplename, samplename),
                      sep=',', header=TRUE, row.names = 1)
umi_table$TotalMicrobe <- rowSums(umi_table, na.rm = T)
ST.seu <- AddMetaData(ST.seu, umi_table)
ST.seu$TotalMicrobe[is.na(ST.seu$TotalMicrobe)] = 0
SpatialFeaturePlot(ST.seu, features = "TotalMicrobe")
ggsave(sprintf("Results/%s/outs/%s.TotalMicrobe.pdf", samplename, samplename),
       width = 8, height = 8)
SpatialFeaturePlot(ST.seu, features = "Campylobacter")
ggsave(sprintf("Results/%s/outs/%s.Campylobacter.pdf", samplename, samplename),
       width = 8, height = 8)
saveRDS(ST.seu, file.path("Results", samplename, "outs", "ST.seu.rds"))

