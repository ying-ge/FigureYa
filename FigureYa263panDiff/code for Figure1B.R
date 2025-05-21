# 作者：大鱼海棠
# 单位：中国药科大学国家天然药物重点实验室，生物统计与计算药学研究中心

# 需求：分析铁死亡基因在肿瘤与正常组织中的差异表达情况并绘制拼接的点图与柱状图

# 设置工作路径
workdir <- "E:/IGBMC/externalproject/ExternalProject/FigurePancan1B"; setwd(workdir)

# 加载R包
library(ggplot2)
library(data.table)
library(randomcoloR)
library(ggpubr)
library(GSVA)
library(clusterProfiler)
library(impute)
library(ComplexHeatmap)
source("twoclasslimma.R")

# 设置颜色
blue <- "#4577FF"
red <- "#C2151A"
orange <- "#E45737"
green <- "#6F8B35"
darkblue <- "#303B7F"
darkred <- "#D51113"
yellow <- "#EECA1F"

# 获得同时有肿瘤和正常样本的肿瘤名
tumors <- c("BLCA","BRCA","CESC","CHOL","COAD",
            "ESCA","GBM","HNSC","KICH","KIRC",
            "KIRP","LIHC","LUAD","LUSC","PAAD",
            "PRAD","READ","STAD","THCA","UCEC")

# 获得感兴趣的基因集(TTC35是EMC2的同名)
frg <- c("CDKN1A","HSPA5","TTC35","SLC7A11","NFE2L2","MT1G","HSPB1","GPX4","FANCD2","CISD1","FDFT1","SLC1A5","SAT1",
         "TFRC","RPL8","NCOA4","LPCAT3","GLS2","DPP4","CS","CARS","ATP5G3","ALOX15","ACSL4","EMC2")

# 修正TCGA名称
# https://gdc.cancer.gov/about-data/publications/pancanatlas
rawAnno <- read.delim("merged_sample_quality_annotations.tsv",sep = "\t",row.names = NULL,check.names = F,stringsAsFactors = F,header = T) # 数据来自PanCanAtlas
rawAnno$simple_barcode <- substr(rawAnno$aliquot_barcode,1,15)
samAnno <- rawAnno[!duplicated(rawAnno$simple_barcode),c("cancer type", "simple_barcode")]
samAnno <- samAnno[which(samAnno$`cancer type` != ""),]
write.table(samAnno,"simple_sample_annotation.txt",sep = "\t",row.names = F,col.names = T,quote = F)

#-----------#
# Figure 1B #

# 快速读取表达谱
# https://gdc.cancer.gov/about-data/publications/pancanatlas
expr <- fread("EBPlusPlusAdjustPANCAN_IlluminaHiSeq_RNASeqV2.geneExp.tsv",sep = "\t",stringsAsFactors = F,check.names = F,header = T)
expr <- as.data.frame(expr); rownames(expr) <- expr[,1]; expr <- expr[,-1]
gene <- sapply(strsplit(rownames(expr),"|",fixed = T), "[",1) # 调整行名
expr$gene <- gene
expr <- expr[!duplicated(expr$gene),] # 移除重复样本
rownames(expr) <- expr$gene; expr <- expr[,-ncol(expr)]

comgene <- intersect(rownames(expr),frg) # 取部分表达谱（感兴趣的基因集）
expr_sub <- expr[comgene,]
colnames(expr_sub) <- substr(colnames(expr_sub),1,15)
expr_sub <- expr_sub[,!duplicated(colnames(expr_sub))]

exprTab <- ndegs <- NULL
log2fc.cutoff <- log2(1.5) # 设置差异表达的阈值 (FoldChange = 1.5)
fdr.cutoff <- 0.05 # 设置差异表达的阈值 (FDR = 0.05)
for (i in tumors) {
  message("--",i,"...")
  sam <- samAnno[which(samAnno$`cancer type` == i),"simple_barcode"]
  comsam <- intersect(colnames(expr_sub), sam)
  
  tumsam <- comsam[substr(comsam,14,14) == "0"] # 获得肿瘤样本
  norsam <- comsam[substr(comsam,14,14) == "1"] # 获得正常样本
  
  expr_subset <- expr_sub[,c(tumsam,norsam)]
  expr_subset[expr_subset < 0] <- 0 # 这份数据里存在负值，即便负值比较小，但也要矫正，如果使用其他泛癌表达谱根据情况而定
  expr_subset <- as.data.frame(impute.knn(as.matrix(expr_subset))$data)
  write.table(expr_subset, paste0("TCGA_",i,"_expr_subset.txt"),sep = "\t",row.names = T,col.names = NA,quote = F)
  
  subt <- data.frame(condition = rep(c("tumor","normal"),c(length(tumsam),length(norsam))),
                     row.names = colnames(expr_subset),
                     stringsAsFactors = F)
  twoclasslimma(subtype  = subt, # 亚型信息 (必须含有一列叫'condition')
                featmat  = expr_subset, # 表达谱 (会自动判断数据量级)
                treatVar = "tumor", # “治疗组”的名（就是要比较的组）
                ctrlVar  = "normal", # “对照组”的名（就是被比较的组）
                prefix   = paste0("TCGA_",i), # 差异表达的文件的前缀
                overwt   = T, # 是否覆盖已经存在的差异表达文件
                sort.p   = F, # 是否排序p值
                verbose  = TRUE, # 是否简化输出
                res.path = workdir) # 输出结果
  
  # 加载差异表达文件
  res <- read.table(paste0("TCGA_",i,"_limma_test_result.tumor_vs_normal.txt"),sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)
  upgene <- res[which(res$log2fc > log2fc.cutoff & res$padj < fdr.cutoff),] # 获取上调基因
  dngene <- res[which(res$log2fc < -log2fc.cutoff & res$padj < fdr.cutoff),] # 获取下调基因
  
  # 基因差异表达的数目以产生图片的上部注释
  if(nrow(upgene) > 0) {
    nup <- nrow(upgene)
  } else {nup <- 0}
  
  if(nrow(dngene) > 0) {
    ndn <- nrow(dngene)
  } else {ndn <- 0}
  
  exprTab <- rbind.data.frame(exprTab,
                              data.frame(gene = rownames(res),
                                         log2fc = res$log2fc,
                                         FDR = res$padj,
                                         tumor = i,
                                         stringsAsFactors = F),
                              stringsAsFactors = F)
  ndegs <- rbind.data.frame(ndegs,
                            data.frame(tumor = i,
                                       Group = c("UP","DOWN"),
                                       Number = c(nup,ndn),
                                       stringsAsFactors = F),
                            stringsAsFactors = F)
}
ndegs$Group <- factor(ndegs$Group, levels = c("UP","DOWN"))

# 产生上部注释
p_top <- ggplot(data = ndegs) +
  geom_bar(mapping = aes(x = tumor, y = Number, fill = Group), 
           stat = 'identity',position = 'stack') + 
  scale_fill_manual(values = c(orange,green)) +
  theme_classic() + 
  theme(axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        plot.margin = unit(c(1,0,0,1), "lines"))

# 产生下部泡泡图
exprTab$gene <- factor(exprTab$gene,
                       levels = rev(c("CDKN1A","HSPA5","TTC35","SLC7A11","NFE2L2","MT1G","HSPB1","GPX4","FANCD2","CISD1",
                                      "FDFT1","SLC1A5","SAT1","TFRC","RPL8","NCOA4","LPCAT3","GLS2","DPP4","CS","CARS","ATP5G3","ALOX15","ACSL4")))
my_palette <- colorRampPalette(c(green,"white",orange), alpha=TRUE)(n=128)
p_center <- ggplot(exprTab, aes(x=tumor,y=gene)) +
  geom_point(aes(size=-log10(FDR),color=log2fc)) +
  scale_color_gradientn('log2(FC)', 
                        colors=my_palette) + 
  theme_bw() +
  theme(panel.grid.minor = element_blank(), 
        panel.grid.major = element_blank(),
        axis.text.x = element_text(angle = 45, size = 12, hjust = 0.3, vjust = 0.5, color = "black"),
        axis.text.y = element_text(size = 12, color = rep(c(red,blue),c(14,10))),
        axis.title = element_blank(),
        panel.border = element_rect(size = 0.7, linetype = "solid", colour = "black"),
        plot.margin = unit(c(0,0,1,1), "lines")) 

# 排布图片
ggarrange(p_top,
          p_center, 
          nrow = 2, ncol = 1,
          align = "v",
          heights = c(2,6),
          common.legend = F)
ggsave("Figure 1B differential expression of interested genes in pancancer.pdf", width = 8,height = 8)
rm(expr); gc()
