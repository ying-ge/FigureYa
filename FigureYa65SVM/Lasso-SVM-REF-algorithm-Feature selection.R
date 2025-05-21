





#' @使用package-GEOquery读入从GEO下载的数据-130M
#' @下载链接：：ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE75nnn/GSE75041/matrix/
library(GEOquery)
eSet <- getGEO( filename = 'GSE75041_series_matrix.txt.gz',  destdir='./',getGPL = F)
pdata<-pData(eSet)
pdata[1:10,1:6]
##################################
#' @训练目标就在title上
pdata$title
group<-pdata[,c(1:2)];group
group$group<-grepl(group$title,pattern = "NR")
group$group<-ifelse(group$group=="TRUE","NR","R")
summary(as.factor(group$group))
group<-group[,-1]
colnames(group)[1]<-"ID"
#############################

#' @获取变量矩阵
eset <- exprs(eSet) 
eset[1:7,1:8];dim(eset)
#########################

library(tidyverse)
eset<-as.data.frame(t(eset));
eset<-rownames_to_column(eset,"ID")
summary(group$ID %in%eset$ID)
# train<-inner_join(eset,group,by="ID") #这一步运行很久还是不加进去了-直接按照排序来吧
group<-group[order(group$ID),]
eset<-eset[order(eset$ID),]
data.frame(group$ID,eset$ID)
#######################################


#' @先进行变量的初筛-使用lapply批量进行统计检验-wilcoxon
#' @数据量有点大(48万个特征)，用lappy都跑了大概5分钟

result<-lapply(eset[,setdiff(colnames(eset),c("ID"))], 
               function(x) wilcox.test(x ~ group$group, var.equal = TRUE))

#' @提取检验结果里面的Pvalue
pvalue<-data.frame(p.value = sapply(result, getElement, name = "p.value"));head(pvalue)
colnames(pvalue)<-"P_value"
pvalue<-rownames_to_column(pvalue,"features")


#' @选择显著有差异的位点--P<0.005
fea_names<-pvalue[pvalue$P_value<0.005,]$features
fea_names#共获取到6000个特征


#' @构建训练矩阵-构建训练矩阵
train<-inner_join(group,eset[,c("ID",fea_names)],by="ID")
train[1:5,1:10];dim(train)
rownames(train)<-train$ID;train<-train[-1]
train$group<-as.factor(train$group)

#' @发现有NA值-使用均值替换
sum(is.na(train))


#' @循环使用均值来代替每列的缺失值
for(i in 2:ncol(train)){
  train[is.na(train[,i]),i]<-mean(train[!is.na(train[,i]),i])
}
sum(is.na(train))
write.csv(train,"input.csv")
#################################

x<-as.matrix(train[,-1])
y<-ifelse(train$group == "NR", 0,1)
##################################

#' @进行lasso和SVM-REF的变量筛选
######################################

#' @第一部分
#' @LASSO-logitstic-Algorithm
library(glmnet)
fit = glmnet(x, y, family = "binomial",alpha = 1,lambda = NULL)


#' @如果想把所有标记都用探针号表示出来，请购买-小丫代码-31
plot(fit, xvar = "dev", label = TRUE)

cvfit = cv.glmnet(x, y, family = "binomial", type.measure = "class")
plot(cvfit)

cvfit$lambda.min
cvfit$lambda.1se

#' @LASSO选出来的特征如下
myCoefs <- coef(cvfit, s="lambda.min");
lasso_fea<-myCoefs@Dimnames[[1]][which(myCoefs != 0 )]
lasso_fea<-lasso_fea[-1];
lasso_fea

#' @可以使用选出来的特征进行模型的预测---因为没有验证数据-这里还是用了训练的数据来凑合一下
predict<-predict(cvfit, newx = x[1:nrow(x),], s = "lambda.min", type = "class")
table(predict,y)  #百分之百预测准确性
###################################


#' @第二部分
#' @使用SVM-REF-Algorithm-进行特征选择


#' @方法一：使用johncolby的代码
#' @https://github.com/johncolby/SVM-RFE
#' @http://www.colbyimaging.com/wiki/statistics/msvm-rfe
# Set up R environment
source('msvmRFE.R')  
library(e1071)
library(parallel)
input <- read.csv("input.csv",header = T,row.names = 1)  #这里默认第一列为目标-不需要分x,y
str(input)
#' @采用五折交叉验证
svmRFE(input, k=5, halve.above=100)
nfold = 5
nrows = nrow(input)
folds = rep(1:nfold, len=nrows)[sample(nrows)]
folds = lapply(1:nfold, function(x) which(folds == x))
results = lapply(folds, svmRFE.wrap, input, k=5, halve.above=100)
results
top.features = WriteFeatures(results, input, save=F)
head(top.features)

#' @画图--查看--选择多少个变量的时候预测错误率最低

# Each featsweep list element corresponds to using that many of the top features (i.e. featsweep[1] is using only the top feature, featsweep[2] is using the top 2 features, etc.). Within each, svm.list contains the generalization error estimates for each of the 10 folds in the external 5-fold CV. These accuracies are averaged as error.

#' @运行时间主要取决于选择变量的个数-下面选了前300个变量进行SVM模型构建-运行了大概9个小时(16G-台式机)
#' @（请会写并行运算的大神指点一下如何改进代码,是建模--节约时间的刚需）
#To run the code given bellow will take huge time. Thus the process
#data is given below. 

# featsweep = lapply(1:300, FeatSweep.wrap, results, input)
# featsweep
# save(featsweep,file = "featsweep.RData")
load("featsweep.RData")
######################################
no.info = min(prop.table(table(input[,1])))
errors = sapply(featsweep, function(x) ifelse(is.null(x), NA, x$error))
dev.new(width=4, height=4, bg='white')
# save(errors,file = "errors.RData")
load('errors.RData')
PlotErrors(errors, no.info=no.info)
Plotaccuracy(1-errors,no.info=no.info)

intersect(lasso_fea,top.features[1:212,"FeatureName"])
summary(lasso_fea%in%top.features[1:212,"FeatureName"])


# install.packages("VennDiagram")
library("VennDiagram")
venn.plot <- venn.diagram(list(LASSO=lasso_fea, 
                               SVM_RFE= as.character(top.features[1:212,"FeatureName"])), NULL, 
                          fill=c("#E31A1C","#E7B800"), 
                          alpha=c(0.5,0.5), cex = 4, cat.fontface=3, 
                          category.names=c("LASSO", "SVM_RFE"), 
                          main="Overlap")
# To plot the venn diagram we will use the grid.draw() function to plot the venn diagram
grid.draw(venn.plot)
#######################




#' @方法2：使用package:sigFeature
#' @https://github.com/pijush1285/sigFeature 
library(sigFeature);library(e1071)
#######################################
sigfeatureRankedList <- sigFeature(x, y)
sigFeature.model=svm(x[ ,sigfeatureRankedList[1:300]], y, 
                     type="C-classification", kernel="linear")
summary(sigFeature.model)
pred <- predict(sigFeature.model, x[ ,sigfeatureRankedList[1:300]])
table(pred,y)
#################################
featureRankedList = svmrfeFeatureRanking(x,y)
print(colnames(x)[featureRankedList[1:300]])
#################################
RFE.model=svm(x[ ,featureRankedList[1:500]], y, 
              type="C-classification", kernel="linear")
summary(RFE.model)
pred <- predict(RFE.model, x[ ,featureRankedList[1:300]])
table(pred,y)


#' @查看第一种方法和第二种方法获得的前300的变量是否大概一致
############################################
intersect(colnames(x)[featureRankedList[1:300]],top.features[1:300,"FeatureName"])


#' @以上两种方法挺相似，但是第二种方法已经开发成包了，后面的代码会报错，github上有人反映情况，但是作者还没有维护好，也说到了最后一步计算变量最佳组合数目需要花费大量时间的问题
#################################





#' @方法3：使用package:caret-RFE算法画原文的图

library(caret);library(randomForest)
set.seed(123)
control <- rfeControl(functions=rfFuncs, method="cv", number=5)


#' @使用第一个方法选择出来的前300个特征使用进行模型构建
rfe.train <- rfe(train[,as.character(top.features$FeatureName[1:300])], 
                 train[,1], sizes=1:300, rfeControl=control,method="svmlinear")
rfe.train
# Now, we’ll plot the result of RFE algorithm and obtain a variable importance chart.
plot(rfe.train, type=c("g", "o"), cex = 0.1, col = "firebrick3",lwd=3)
#Let’s extract the chosen features
predictors(rfe.train)

intersect(lasso_fea,predictors(rfe.train))
summary(lasso_fea%in%predictors(rfe.train))#有3/4的位点出现在REF算法过滤后的结果中
##################################


#' @使用第二个方法选择出来的前300个特征使用进行准确率评估
rfe.train.1 <- rfe(train[,as.character(colnames(x)[featureRankedList[1:300]])], 
                 train[,1], sizes=1:300, rfeControl=control,method="svmlinear")
rfe.train.1
# Now, we’ll plot the result of RFE algorithm and obtain a variable importance chart.
plot(rfe.train.1, type=c("g", "o"), cex = 0.1, col = "firebrick3",lwd=3)
#Let’s extract the chosen features
predictors(rfe.train.1)
intersect(lasso_fea,predictors(rfe.train.1))
summary(lasso_fea%in%predictors(rfe.train.1))
###############################



