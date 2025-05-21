#要求coxHR计算，单因素cox和多因素cox，输出为表格
#清除数据，安装数据包
rm(list = ls())
setwd('/Users/huaiy/Desktop/cox HR/')
if(!require("tidyverse")) install.packages("tidyverse")
if(!require('survival')) install.packages('survival')
library(tidyverse)
#临床样本信息，miRNA及mRNA表达谱数据来源于UCSC xena网站下载
#特别说明 mRNA为log2(reads+1) 后标准化的值(文章用什么标准化没陈述)
#miRNA为log2(RPM+1)后标准化的值(注:文章用的log2(RPM))
#临床信息已经在excel中进行过简单的整理，筛选后和文章陈述样本数量不一致，因为文章没有说明怎么筛选删除样本的。
#筛选样本处理原则，只要cox HR研究相关因素信息有缺失的，一律删除整个样本。
#文章补充材料提供了miRNA样本及表达值，也有点问题，给了276个样本，包括7个正常样本，去除重复肿瘤样本还有250多个样本。
#与文章表一240个样本不一致。疑问作者是怎么筛选样本数据，不是很清楚。大家如果看文献仔细的话，可以互相交流。
#另外文章验证数据集数据60样本没有上传，所以无法进行验证，只能从TCGA随机抽取60个样本作为验证。
clin<-read.csv('clin_coad.txt',sep = '\t')##已预处理后clin数据
surv<-readr::read_tsv('TCGA-COAD.survival.tsv')
mRNA<-readr::read_tsv('TCGA-COAD.htseq_counts.tsv')
miRNA<-readr::read_tsv('TCGA-COAD.mirna.tsv')
#首先临床病理信息和生存信息整合
surv<-surv[,c(1,2,6)]
colnames(surv) <- c("submitter_id.samples","EVENT","OS")
clin<-dplyr::inner_join(clin,surv,by='submitter_id.samples')
#按示例表格中的顺序排列
clin <- clin[,c(1,8,2:7,9:11)]
#mRNA miRNA共有的样本名；原则上研究应该是共有的样本。
sam<-intersect(colnames(mRNA),colnames(miRNA))
#mRNA miRNA共有样本和已筛选后的临床信息进行交集
clin_filter<-dplyr::filter(clin,submitter_id.samples%in%sam)
#其次筛选出肿瘤样本的临床信息
clin_filter_tumor<-dplyr::filter(clin_filter,
                                 substr(clin_filter$submitter_id.samples,14,14)==0)
#提取出miRNA hsa-mir-195表达数据
tumor_sam<-as.character(clin_filter_tumor$submitter_id.samples)
miRNA_tumor<-dplyr::select(miRNA,miRNA_ID,tumor_sam)%>%
  dplyr::filter(miRNA_ID=='hsa-mir-195')
miRNA_tumor<-as.data.frame(t(miRNA_tumor))
miRNA_tumor$submitter_id.samples<-rownames(miRNA_tumor)
miRNA_tumor$mir195<-miRNA_tumor$V1
miRNA_tumor<-miRNA_tumor[-1,]
rownames(miRNA_tumor)<-NULL
miRNA_tumor$V1<-NULL
#miRNA-195和临床信息合并
clin_miRNA<-dplyr::left_join(clin_filter_tumor,miRNA_tumor,by='submitter_id.samples')
#提取mRNA YAP1 == ENSG00000137693 ENSG00000137693.13
mRNA$Ensembl_ID<-substr(mRNA$Ensembl_ID,1,15)##可截取，也可以直接用.13
mRNA_tumor<-dplyr::select(mRNA,Ensembl_ID,tumor_sam)%>%
  dplyr::filter(Ensembl_ID=='ENSG00000137693')
mRNA_tumor<-as.data.frame(t(mRNA_tumor))
mRNA_tumor$submitter_id.samples<-rownames(mRNA_tumor)
mRNA_tumor$`YAP1`<-mRNA_tumor$V1
mRNA_tumor<-mRNA_tumor[-1,]
rownames(mRNA_tumor)<-NULL
mRNA_tumor$V1<-NULL
#YAP1和miRNA195 临床信息合并
clin_mRNA_miRNA<-dplyr::left_join(clin_miRNA,mRNA_tumor,by='submitter_id.samples')
#去除重复样本
clin_mRNA_miRNA$samp<-substr(clin_mRNA_miRNA$submitter_id.samples,1,12)
clin_mRNA_miRNA<-dplyr::filter(clin_mRNA_miRNA,!duplicated(clin_mRNA_miRNA$samp))
clin_mRNA_miRNA <- clin_mRNA_miRNA[,-samp]
#保存到文件
write.csv(clin_mRNA_miRNA,'easy_input.csv',quote = F, row.names = F)
