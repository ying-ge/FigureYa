# 从GEO数据库下载表达矩阵
# 下载地址：ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE75nnn/GSE75041/matrix/GSE75041_series_matrix.txt.gz
# 130M
library(GEOquery)

eSet <- getGEO(filename = 'GSE75041_series_matrix.txt.gz',  destdir='./',getGPL = F)
pdata <- pData(eSet)
pdata[1:10,1:6]
##################################
#训练目标就在title上
pdata$title
group <- pdata[,c(1:2)];group
group$group <- grepl(group$title,pattern = "NR")
group$group <- ifelse(group$group=="TRUE","NR","R")
summary(as.factor(group$group))
group <- group[,-1]
colnames(group)[1] <- "ID"
#############################
#获取变量矩阵
eset <- exprs(eSet) 
eset[1:7,1:8]; dim(eset)
#########################
eset <- as.data.frame(t(eset));
eset <- rownames_to_column(eset,"ID")
summary(group$ID %in%eset$ID)
group <- group[order(group$ID),]
eset <- eset[order(eset$ID),]
data.frame(group$ID,eset$ID)
#######################################
#先进行变量的初筛-使用lapply批量进行统计检验-wilcoxon
#数据量有点大(48万个特征)，用lappy都跑了大概5分钟
result <- lapply(eset[,setdiff(colnames(eset),c("ID"))], 
                 function(x) wilcox.test(x ~ group$group, var.equal = TRUE))
#提取检验结果里面的Pvalue
pvalue <- data.frame(p.value = sapply(result, getElement, name = "p.value")); head(pvalue)
colnames(pvalue) <- "P_value"
pvalue <- rownames_to_column(pvalue,"features")
#选择显著有差异的位点 -- P < 0.005
fea_names <- pvalue[pvalue$P_value < 0.005, ]$features
fea_names#共获取到6000个特征
#构建训练矩阵-构建训练矩阵
train <- inner_join(group, eset[, c("ID",fea_names)], by="ID")
train[1:5,1:10]; dim(train)
rownames(train) <- train$ID;train - train[-1]
train$group <- as.factor(train$group)
#发现有NA值-使用均值替换
sum(is.na(train))
#循环使用均值来代替每列的缺失值
for(i in 2:ncol(train)){
  train[is.na(train[,i]), i] <- mean(train[!is.na(train[,i]),i])
}
sum(is.na(train))
write.csv(train,"easy_input.csv")
