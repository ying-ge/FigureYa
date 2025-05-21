#library
rm(list=ls())

library(GSA)
library(randomForest)
library(preprocessCore)
library(glmnet)

## input
x=read.csv('data.csv')
x <- na.omit(x)
#Only expression profiles of WT livers and tumours will be used to train model
y=c(rep(1,27),rep(2,11))
x=x[,c(1,14:16,26:28,32:34,38:40,44:46,59:61,65:67,71:73,89:91,56:58,86:88,96:100)]
z=read.csv('TFcluster.csv')

## significantly changed TF clusters
x_m=as.matrix(x[,-1])
rownames(x_m)=toupper(as.character(x$X))
x1=x_m
y1=y

gene_sets=list(1)
for (i in 1:dim(z)[2]){
  gene_sets[[i]]=toupper(unique(as.character(z[,i])))
}
genenames= toupper(rownames(x1))
GSA_temp=GSA(x1,y1, genesets=gene_sets, genenames,method="maxmean", 
             resp.type="Two class unpaired",minsize=40,maxsize=5000,nperms=1000)

GSA_result=cbind(GSA_temp$GSA.scores,GSA_temp$pvalues.lo,GSA_temp$pvalues.hi)
colnames(GSA_result)=c('score','pvalue_low','pvalue_high')
rownames(GSA_result)=toupper(colnames(z))
GSA_result1=GSA_result[!is.na(GSA_temp$GSA.scores),]
keep_set=(abs(GSA_result1[,1])>0.8 & (GSA_result1[,2]<0.05 | GSA_result1[,3]<0.05))
GSA_result2=GSA_result1[keep_set,]
##最后得出61个cluster
### heatmap
# code in TF.R

### weight matrix (genes*TF) using randomforest
#data normalization
x_f=normalize.quantiles(x_m)
rownames(x_f)=rownames(x_m)
colnames(x_f)=colnames(x_m)

#weight matrix
##算出TF对应的基因对结局变量的贡献度，获得的矩阵是TF*genes
M=matrix(0,dim(GSA_result2)[1],dim(x_f)[1])
rownames(M)=rownames(GSA_result2)
colnames(M)=rownames(x_f)
#random forest用的是normalization之后的数据
for (i in 1:dim(M)[1]){
  temp_TF=rownames(M)[i]
  temp_genes=unique(toupper(as.character(z[,temp_TF==colnames(z)])))
  temp_keep=(rownames(x_f)=='plplpl')
  for (j in 1:length(temp_genes)){
    temp_keep1=(temp_genes[j]==toupper(rownames(x_f)))
    temp_keep=(temp_keep | temp_keep1)
  }
  x_temp=x_f[temp_keep,]
  x_temp1=rbind(y,x_temp)
  x_forest=as.data.frame(t(x_temp1))
  
  sample.ind <- sample(2, nrow(x_forest),replace = T,prob = c(1,0))
  x_dev <- x_forest[sample.ind==1,]
  x_val <- x_forest[sample.ind==2,]
  
  varNames <- names(x_dev)
  varNames <- varNames[!varNames %in% c("y")]
  temp_var=varNames
  varNames1 <- paste('V',1:length(varNames), sep = "")
  colnames(x_dev)=c('y',varNames1)
  varNames2 <- paste(varNames1, collapse = "+")
  rf.form <- as.formula(paste("y", varNames2, sep = " ~ "))
  x_dev$y=as.factor(x_dev$y)
  cross.sell.rf <- randomForest(rf.form,x_dev,ntree=1000,importance=T)
  TF_M=cross.sell.rf$importance
  rownames(TF_M)=temp_var
  
  temp_listi=rownames(TF_M)
  temp_listj=colnames(M)
  for (temp_i in 1:dim(TF_M)[1]){
    for (temp_j in 1:dim(M)[2]){
      if (temp_listi[temp_i]==temp_listj[temp_j]){
        M[i,temp_j]=TF_M[temp_i,4]
      }
    }
  }
}
M_weight=t(M)

###deviation of samples: samples*TFs (normal liver and HCC)
D_normal=matrix(1000,dim(x_m)[2],dim(GSA_result2)[1])
rownames(D_normal)=colnames(x_m)
colnames(D_normal)=rownames(GSA_result2)
D_cancer=matrix(1000,dim(x_m)[2],dim(GSA_result2)[1])
rownames(D_cancer)=colnames(x_m)
colnames(D_cancer)=rownames(GSA_result2)

#deviatino of gene expression pattern
Ref_normal=x_m[,y==1]
Ref_cancer=x_m[,y==2]

for (temp_i in 1:dim(D_normal)[1]){
  temp_pro = x_m[,temp_i]
  normal_zscore = rep(0,length(temp_pro))
  cancer_zscore = rep(0,length(temp_pro))
  for (temp_genes in 1:length(temp_pro)){
    #normal
    temp_ave=mean(Ref_normal[temp_genes,])
    temp_sd=sd(Ref_normal[temp_genes,])
    if (temp_ave>0 & temp_sd>0){
      normal_zscore[temp_genes]=abs((temp_pro[temp_genes]-temp_ave)/temp_sd)
    }
      #cancer 
      temp_ave=mean(Ref_cancer[temp_genes,])
      temp_sd=sd(Ref_cancer[temp_genes,])
      if (temp_ave>0 & temp_sd>0){
        cancer_zscore[temp_genes]=abs((temp_pro[temp_genes]-temp_ave)/temp_sd)
      }
    }
  
  for (temp_j in 1:dim(D_normal)[2]){
    D_normal[temp_i,temp_j] = normal_zscore %*% M_weight[,temp_j]
    D_cancer[temp_i,temp_j] = cancer_zscore %*% M_weight[,temp_j]
  }
}

### adjusted activity of TF clusters (normal liver: -1, HCC:1)
A_TFactivity=matrix(0,dim(x_m)[2],dim(GSA_result2)[1])
rownames(A_TFactivity)=colnames(x_m)
colnames(A_TFactivity)=rownames(GSA_result2)

sign_TF=rep(0,dim(A_TFactivity)[2])
sign_TF[(GSA_result2[,1]>0)]=1
sign_TF[(GSA_result2[,1]<0)]=(-1)

for (ti in 1:dim(A_TFactivity)[1]){
  for (tj in 1:dim(A_TFactivity)[2]){
    A_TFactivity[ti,tj]=((D_normal[ti,tj]-D_cancer[ti,tj])/(D_normal[ti,tj]+D_cancer[ti,tj]))*sign_TF[tj]
  }
}
#step5:每个样本的tumor-promoting and -inhibiting strengths
## average activity of promoting and suppressing TF group
TFcluster=matrix(0,dim(A_TFactivity)[1],2)
colnames(TFcluster)=c('T_pro','T_supp')
rownames(TFcluster)=rownames(A_TFactivity)
temp_normalTF=A_TFactivity[,(sign_TF==(-1))]
temp_cancerTF=A_TFactivity[,(sign_TF==(1))]
TFcluster[,1]=rowMeans(temp_cancerTF)
TFcluster[,2]=rowMeans(temp_normalTF)

## weights of TF clusters in determining phenotypes
res_vector=y
res_vector[y==1]=(-1)
res_vector[y==2]=1
var_input=TFcluster
fit = glmnet(var_input, res_vector)
weig=coef(fit,s=0.1)


## index
Sample_index=predict(fit,newx=TFcluster,s=c(0.1,0.05))

#####
Readout=list(1)
Readout[[1]]=GSA_result2
Readout[[2]]=M_weight
Readout[[3]]=A_TFactivity
Readout[[4]]=TFcluster
Readout[[5]]=Sample_index
Readout[[6]]=fit

all_varlist=ls()
rm_var=(all_varlist[!(all_varlist=='Readout')])
rm(list=rm_var)
save(Readout, file = "Readout.RData")
