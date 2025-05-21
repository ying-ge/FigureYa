---
title: "FigureYa53 PPImodule"
author: "小丫画图出品"
date: "2018-11-11"
output: html_document
---
小丫微信: epigenomics  E-mail: figureya@126.com

作者：加大发

小丫编辑校验

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
```

## 需求描述

通过已知PPI数据库找到共表达蛋白中的module，输出到文件，用于下一步cytoscape展示。

![demo](FigureYa53PPImodule-0527更新.assets/demo.png)

出自<https://www.sciencedirect.com/science/article/pii/S1074761317300729?via%3Dihub>

## 应用场景

把自己的基因列表跟蛋白相互作用网络联系起来。

例如RNA-seq获得的差异表达基因，看他们在蛋白相互作用网络中，哪些基因处于同一module。


## 环境设置

安装需要的包

```r
#使用国内镜像安装包
options("repos"= c(CRAN="https://mirrors.tuna.tsinghua.edu.cn/CRAN/"))
install.packages("dynamicTreeCut")
install.packages("openxlsx")
install.packages("stringr")
install.packages("Matrix")
install.packages("WGCNA")
```

加载需要用到的包

```{r}
library(dynamicTreeCut)
library(openxlsx)
library(stringr)
library(Matrix)
library(WGCNA)

Sys.setenv(LANGUAGE = "en") #显示英文报错信息
options(stringsAsFactors = FALSE) #禁止chr转成factor
```

## 输入文件

需要3个文件：

- 共表达基因列表，此处用paper附件里的mm2.xlsx
- Uniprot ID转换文件：<ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/idmapping/by_organism/MOUSE_10090_idmapping_selected.tab.gz>
- 从STRING数据库下载PPI数据<https://stringdb-static.org/download/protein.actions.v10.5/10090.protein.actions.v10.5.txt.gz>。三个红框表示：数据库版本号，物种搜索框，相互作用数据（包含相互作用类型，如物理相互作用)。注意版本号，以及物种搜索框内为‘Mus musculus’，下载第三个框内文件。

![STRING](FigureYa53PPImodule-0527更新.assets/STRING.png)

## ！关于输入文件的更新

需要3个文件：

- ！共表达基因列表，此处用paper附件里的**mmc2.xlsx， Data S1. Whole Proteome and Phosphoproteome Analyses of WT and Rptor-Deficient T Cells (1st experiment)**.  https://ars.els-cdn.com/content/image/1-s2.0-S1074761317300729-mmc2.xlsx，下载之后重命名为**mmc2.xlsx**

- ！Uniprot ID转换文件：**链接更新为：https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/idmapping/by_organism/MOUSE_10090_idmapping_selected.tab.gz   把ftp改为https。**     《 **对应的人的为https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/idmapping/by_organism/HUMAN_9606_idmapping_selected.tab.gz** 》   《人和小鼠对应的版本号：HUMAN_9606， **MOUSE_10090**》  MOUSE_10090_idmapping_selected.tab.gz下载之后解压。

- **！从STRING数据库下载PPI数据：（此处更新较大！！）**

  在https://string-db.org/cgi/download网站看左上角版本号：**version 11.0**，输人或者小鼠，左下角没有找到小丫代码里面的文件。（如下图示：）

  ![STRINGv11-mm](FigureYa53PPImodule-0527更新.assets/STRINGv11-mm.png)

![STRINGv11-hg](FigureYa53PPImodule-0527更新.assets/STRINGv11-hg.png)

在https://stringdb-static.org/download/网站找到对应的**protein.actions.v11.0/** 路径，在此路径下找到 HUMAN_9606和**MOUSE_10090**分别对应的文件9606.protein.actions.v11.0.txt.gz https://stringdb-static.org/download/protein.actions.v11.0/9606.protein.actions.v11.0.txt.gz、**10090.protein.actions.v11.0.txt.gz https://stringdb-static.org/download/protein.actions.v11.0/10090.protein.actions.v11.0.txt.gz**

10090.protein.actions.v11.0.txt.gz下载之后解压，10090.protein.actions.v11.0.txt文件有318MB。

### 从paper附件获得WPC3基因及蛋白ID

```{r convert genename to ENSEMBLPROT}
wpc <- read.xlsx("mmc2.xlsx",cols=2:4,startRow = 3) #读入文件第一个sheet中的第2:4列，从第3行开始
wpc3 <- wpc[wpc$WPC==3&!is.na(wpc$WPC),]

uniprot_wpc3 <- str_split(wpc3$Protein.Accession,"[|]",3,simplify = T)[,2]
uniprot_wpc3 <- str_split(uniprot_wpc3,"[-]",2,simplify = T)[,1]
wpc3$uniprot <- uniprot_wpc3
```

### 根据Uniprot ID提取ENSEMBLPROT ID

```{r}
idmapping<-read.table("MOUSE_10090_idmapping_selected.tab",header = F,as.is=T,sep="\t")
rownames(idmapping)<-idmapping$V1
ids_wpc3 <- idmapping[intersect(uniprot_wpc3,rownames(idmapping)),][,c(1,21)]
dim(ids_wpc3)#由于uniprot版本问题会丢失1456-1415=41个蛋白
prots <- str_split(ids_wpc3$V21,"; ")
names(prots)<- rownames(ids_wpc3)
prots <- prots[prots!=""]
prots_tmp <- unlist(lapply(seq_along(prots), function(x) paste(names(prots)[x], prots[[x]],sep=";")))
prots_mat <- data.frame(str_split(prots_tmp,"[;]",2,simplify = T),stringsAsFactors = F)
colnames(prots_mat) <- c("uniprot","ensemblprot")
head(prots_mat)
gene_ensemblprot <- merge(wpc3,prots_mat,by="uniprot")
gene2ensemblprot <- gene_ensemblprot[,c(1,3,ncol(gene_ensemblprot))]
head(gene2ensemblprot)#后面会用到gene和ensemblprot id的对应关系
ensemblprot_wpc3 <- unique(prots_mat$ensemblprot)
```

### 从STRING数据库提取WPC3蛋白的PPI

从STRING数据库提取WPC3蛋白的PPI，把蛋白名转成基因名

#STRING数据文件10090.protein.actions.v10.5.txt有580M，我这里8G内存，读取它需要大概30s

**！更新后的STRING数据文件10090.protein.actions.v11.0.txt文件有318MB。**

```{r}
#ppi <- read.table("10090.protein.actions.v10.5.txt",header=T,sep = "\t")
ppi <- read.table("10090.protein.actions.v11.0.txt",header=T,sep = "\t")
head(ppi)
ppi$item_id_a <- str_replace(ppi$item_id_a,"10090.","")
ppi$item_id_b <- str_replace(ppi$item_id_b,"10090.","")
ppi_wpc3 <- ppi[(ppi$item_id_a %in% ensemblprot_wpc3)&(ppi$item_id_b %in% ensemblprot_wpc3),]
ppi_wpc3 <- ppi_wpc3[!duplicated(paste0(ppi_wpc3$item_id_a,ppi_wpc3$item_id_b)),]#去掉重复行

colnames(ppi_wpc3)[1]<-"ensemblprot"
ppi_gene1st <- merge(ppi_wpc3,gene2ensemblprot,by="ensemblprot")
colnames(gene2ensemblprot)[3] <- "item_id_b"
ppi_gene2nd <- merge(ppi_gene1st,gene2ensemblprot,by="item_id_b")
ppi_gene2nd <- ppi_gene2nd[!duplicated(paste0(ppi_gene2nd$Gene.Names.x,ppi_gene2nd$Gene.Names.y)),]
ppi_gene2nd <- na.omit(ppi_gene2nd)  #后面稀疏矩阵 NA's in (i,j) are not allowed
head(ppi_gene2nd)
genes <- sort(unique(union(ppi_gene2nd$Gene.Names.x,ppi_gene2nd$Gene.Names.y)))
#保存到文件
write.csv(data.frame(genes),file="ppi_gene.csv",row.names = F,quote = F)
```

以对称矩阵的形式表示网络的边，0表示两点之间无连边，1表示有边

```{r}
pos <- seq_along(genes)
names(pos) <- genes
ppi_gene2nd$posx <- pos[ppi_gene2nd$Gene.Names.x]
ppi_gene2nd$posy <- pos[ppi_gene2nd$Gene.Names.y]

#构建稀疏矩阵
ppi_mat <- sparseMatrix(ppi_gene2nd$posx, ppi_gene2nd$posy, x = 1, dims = c(length(genes), length(genes)))
#保存到文件
write.csv(as.matrix(ppi_mat),file="ppi_mat.csv",row.names = F,quote = F)
```

这时可以中场休息一下，下面继续

```{r}
NodeNames <- read.csv("ppi_gene.csv", header=T)
NodeNames <- as.vector(NodeNames[,1])
head(NodeNames)

ADJ <- read.csv("ppi_mat.csv", header=T)
ADJ <- as.matrix(ADJ)
rownames(ADJ) <- NodeNames
colnames(ADJ) <- NodeNames
ADJ[1:3,1:3]
```

## 计算TOM

用到函数TOMdist1，保存在TOMdist1.R文件中，位于当前文件夹。

来自：<https://horvath.genetics.ucla.edu/html/GeneralFramework/NetworkFunctions.txt>


```{r}
source("TOMdist1.R") #保存在当前文件夹

# We used the TOM matrix as a smoothed-out version of the adjacency matrix.
dissTOM <- TOMdist1(as.matrix(ADJ))
dissTOM[1:3, 1:3]
```

## 通过cut tree找module

用R包dynamicTreeCut实现。根据文章Supplemental Information描述：

Now we carry out hierarchical clustering with the TOM matrix. Branches of the resulting clustering tree will be used to define gene modules.


```{r}
library(dynamicTreeCut)

#基于TOM(topological overlap matrix)对基因进行层次聚类
#可以通过设置method参数，使用不同方法进行层次聚类，如median,complete, centroid等
hierTOM <- hclust(as.dist(dissTOM), method="average")
str(hierTOM)

#使用hybrid dynamic tree-cutting方法，将网络划分成模块
#可以通过设置minClusterSize来调整模块的最小大小
#minClusterSize的默认设置是20，因为通常太小的模块没有意义
#由于文章中的最小模块是2，所以在此我们设置为2
mods <- cutreeDynamic(hierTOM, cutHeight = NULL, minClusterSize = 2,method = "hybrid", distM = dissTOM)
table(mods)

# 将每个模块包含什么基因输出来，result_module.csv
list2<-c()
for (i in names(table(mods))){
    a <- paste(NodeNames[which(mods==i)], sep="", collapse=',')
    list2 <- c(list2,a)
}
list1 <- paste("Module_",names(table(mods)), sep="")
out_df <- data.frame(list1,list2)
write.csv(out_df,'result_module.csv',row.names=FALSE)

# 转成cytoscape，为fig2.F和fig2.G做准备
library(WGCNA)
#输出整个子网络的边和节点情况
exportNetworkToCytoscape(ADJ,
edgeFile = paste('./result_filtered_PPI_network_edge.txt'),
nodeFile = paste('./result_filtered_PPI_network_node.txt'),
weighted = TRUE,
threshold = 0.5,
nodeNames = NodeNames,
altNodeNames = NULL,
nodeAttr = mods,
includeColNames = TRUE)
```

这一步生成的两个文件，可以作为cytoscape的输入：

- result_filtered_PPI_network_edge.txt
- result_filtered_PPI_network_node.txt


## 参考文献

Dong J, Horvath S (2007) Understanding Network Concepts in Modules, BMC Systems Biology 2007, 1:24 PMID: 17547772 PMCID: PMC3238286

```{r}
sessionInfo()
```