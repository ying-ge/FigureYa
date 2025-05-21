
#Input: training data, testing data
x=read.csv('data.csv')
x <- na.omit(x)
x=x[,c(1,14:16,26:28,32:34,38:40,44:46,59:61,65:67,71:73,89:91,56:58,86:88,96:100)]
y=c(rep(1,27),rep(2,11))
z=read.csv('data.csv')

###data normalization
x$X=toupper(as.character(x$X))
z=z[!is.na(z$X),]
z$X=toupper(as.character(z$X))

Data_all=merge(x,z,by='X',ALL=FALSE)
x_train=Data_all[,c(1:dim(x)[2])]
t_s=(dim(x)[2]+1)
t_e=dim(Data_all)[2]
z_test=Data_all[,c(1,t_s:t_e)]

x_train_m=as.matrix(x_train[,-1])
rownames(x_train_m)=toupper(as.character(x_train$X))
z_test_m=as.matrix(z_test[,-1])
rownames(z_test_m)=toupper(as.character(z_test$X))

x_ref=normalize.quantiles(x_train_m)
rownames(x_ref)=rownames(x_train_m)
colnames(x_ref)=colnames(x_train_m)
temp_sort=sort(x_ref[,1])
for (i in 1:dim(x_train_m)[2]){
  x_train_m[,i]=temp_sort[rank(x_train_m[,i])]
}
for (i in 1:dim(z_test_m)[2]){
  z_test_m[,i]=temp_sort[rank(z_test_m[,i])]
}

### deviation of TFs
### D_normal&D_camcer都是test，样本*TF
#step3:M_weight
GSA_result=Readout[[1]]
D_normal=matrix(1000,dim(z_test_m)[2],dim(GSA_result)[1])
rownames(D_normal)=colnames(z_test_m)
colnames(D_normal)=rownames(GSA_result)
D_cancer=matrix(1000,dim(z_test_m)[2],dim(GSA_result)[1])
rownames(D_cancer)=colnames(z_test_m)
colnames(D_cancer)=rownames(GSA_result)

#step4:TF activity.
#deviatino of gene expression pattern:基因*TF矩阵
Ref_normal=x_train_m[,y==1]
Ref_cancer=x_train_m[,y==2]

M_weight=Readout[[2]]
M_weight1=matrix(0,dim(z_test_m),dim(M_weight)[2])
rownames(M_weight1)=rownames(z_test_m)
colnames(M_weight1)=colnames(M_weight)

M1_genelist=rownames(M_weight1)
M_genelist=rownames(M_weight)
for (M1_i in 1:dim(M_weight1)[1]){
  for (M_i in 1:dim(M_weight)[1]){
    if (M1_genelist[M1_i]==M_genelist[M_i]){
      M_weight1[M1_i,]=M_weight[M_i,]
      break
    }
  }
}
# M_weight1其实就是训练模型的基因*TF，不过其他的基因顺序填充了而已。
for (temp_i in 1:dim(D_normal)[1]){
  temp_pro = z_test_m[,temp_i]#一个样本的基因表达值
  normal_zscore = rep(0,length(temp_pro))
  cancer_zscore = rep(0,length(temp_pro))
  for (temp_genes in 1:length(temp_pro)){
    #normal
    temp_ave=mean(Ref_normal[temp_genes,])#这个基因在normal的表达值
    temp_sd=sd(Ref_normal[temp_genes,])#这个基因在normal的sd
    if (temp_ave>0 & temp_sd>0){
      normal_zscore[temp_genes]=abs((temp_pro[temp_genes]-temp_ave)/temp_sd)#标准化
    }
    #cancer 
    temp_ave=mean(Ref_cancer[temp_genes,])
    temp_sd=sd(Ref_cancer[temp_genes,])
    if (temp_ave>0 & temp_sd>0){
      cancer_zscore[temp_genes]=abs((temp_pro[temp_genes]-temp_ave)/temp_sd)
    }
  }
  
  for (temp_j in 1:dim(D_normal)[2]){
    D_normal[temp_i,temp_j] = normal_zscore %*% M_weight1[,temp_j] #标准化*TF权重:样品*TF
    D_cancer[temp_i,temp_j] = cancer_zscore %*% M_weight1[,temp_j]
  }
}

### adjusted activity of TF clusters (normal liver: -1, HCC:1)
A_TFactivity=matrix(0,dim(z_test_m)[2],dim(GSA_result)[1])
rownames(A_TFactivity)=colnames(z_test_m)
colnames(A_TFactivity)=rownames(GSA_result)

sign_TF=rep(0,dim(A_TFactivity)[2])
sign_TF[(GSA_result[,1]>0)]=1
sign_TF[(GSA_result[,1]<0)]=(-1)

for (ti in 1:dim(A_TFactivity)[1]){
  for (tj in 1:dim(A_TFactivity)[2]){
    A_TFactivity[ti,tj]=((D_normal[ti,tj]-D_cancer[ti,tj])/(D_normal[ti,tj]+D_cancer[ti,tj]))*sign_TF[tj]
  }
}

# step5:计算每个样本的tumor-promoting and -inhibiting strengths. 
## average activity of promoting and suppressing TF group
TFcluster=matrix(0,dim(A_TFactivity)[1],2)
colnames(TFcluster)=c('T_pro','T_supp')
rownames(TFcluster)=rownames(A_TFactivity)

temp_normalTF=A_TFactivity[,(sign_TF==(-1))]
temp_cancerTF=A_TFactivity[,(sign_TF==(1))]
dim(temp_normalTF); dim(temp_cancerTF)
TFcluster[,1]=rowMeans(temp_cancerTF)
TFcluster[,2]=rowMeans(temp_normalTF)

## index
fit=Readout[[6]]
Sample_index=predict(fit,newx=TFcluster,s=c(0.1,0.05))
