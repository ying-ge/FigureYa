### 下载数据
##   https://xena-toil.hiplot.com.cn/download/tcga_RSEM_gene_tpm.gz
##   https://xena-toil.hiplot.com.cn/download/probeMap/gencode.v23.annotation.gene.probemap
##   https://xena-pancanatlas.hiplot.com.cn/download/Survival_SupplementalTable_S1_20171025_xena_sp.gz

rm(list = ls())
dd <- data.table::fread("tcga_RSEM_gene_tpm/tcga_RSEM_gene_tpm",data.table = F)
clin <- data.table::fread("Survival_SupplementalTable_S1_20171025_xena_sp/Survival_SupplementalTable_S1_20171025_xena_sp",data.table = F)
colnames(clin)[2:3] <- c("TCGA_id","type")
index <- clin$sample %in% colnames(dd)
clin <- clin[index,1:3]
dd <- dd[,c("sample",clin$sample)]
colnames(dd)[1] <- "gene_id"
#BiocManager::install("rtracklayer")
##注释信息
gtf1 <- rtracklayer::import('tcga_RSEM_gene_tpm/gencode.v23.annotation.gtf')
gtf_df <- as.data.frame(gtf1)
# save(gtf_df,file = "gtf_df.Rdata")
### 提取编码mRNA,基因注释,行列转置 
library(dplyr)
library(tibble)
mRNA_exprSet <- gtf_df %>% 
  #筛选gene,和编码指标
  filter(type=="gene",gene_type=="protein_coding") %>%
  #选取两列
  select(c(gene_name,gene_id)) %>% 
  #和表达量数据合并
  inner_join(dd,by ="gene_id") %>% 
  #去除一列
  select(-gene_id) %>% 
  #去重
  distinct(gene_name,.keep_all = T) %>% 
  #列名变成行名
  column_to_rownames("gene_name") %>% 
  #转置
  t() %>% 
  #转为数据框
  as.data.frame() %>% 
  # 行名变成列名
  rownames_to_column("sample")

mRNA_exprSet <- cbind(clin,subtype=substring(mRNA_exprSet$sample,14,15),mRNA_exprSet[,-1])

library(tidyverse)
mRNA_exprSet2= mRNA_exprSet %>%
  filter()
### 去除癌旁组织和重复组织
mRNA_exprSet <- mRNA_exprSet %>% 
  filter(subtype != "11") %>% 
  distinct(TCGA_id,.keep_all = T)

mRNA_exprSet[1:5,1:5]

### 转置一下
expr= mRNA_exprSet
rownames(expr)=expr$TCGA_id
expr=as.matrix(t(expr[,5:ncol(expr)]))
epic<-deconvo_tme(eset = expr, method = "epic", arrays = FALSE)
## save(epic,file = "epic.Rda")
## save(expr,file = "expr_pancancer.Rda")
## save(mRNA_exprSet,file = "pancancer_mRNA_exprSet.Rdata")
