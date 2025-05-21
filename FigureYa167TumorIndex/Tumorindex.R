Prepare_input <- function(x_train_m, z_test_m){
  Data_all <- merge(x_train_m,z_test_m,by='X',ALL=FALSE)
  x_train <- Data_all[,c(1:dim(x_train_m)[2])]; 
  x_train_1 <- as.matrix(x_train[,-1])
  rownames(x_train_1) <- toupper(as.character(x_train$X))
  t_s <- (dim(x_train_m)[2]+1)
  t_e <- dim(Data_all)[2]
  z_test <- Data_all[,c(1,t_s:t_e)]
  z_test_1 <- as.matrix(z_test[,-1])
  rownames(z_test_1) <- toupper(as.character(z_test$X))
  
  x_ref=normalize.quantiles(x_train_1)
  rownames(x_ref)=rownames(x_train_1)
  colnames(x_ref)=colnames(x_train_1)
  ##这一步的意义？？
  temp_sort=sort(x_ref[,1])
  for (i in 1:dim(x_train_1)[2]){
    x_train_1[,i]=temp_sort[rank(x_train_1[,i])]
  }
  for (i in 1:dim(z_test_1)[2]){
    z_test_1[,i]=temp_sort[rank(z_test_1[,i])]
  }
  return(list(x_train_m = x_train_1, z_test_m = z_test_1))
}
Prepare_M_weight <- function(z_test_m, Readout){
  M_weight <- Readout[[2]]
  M_weight1 <- matrix(0,dim(z_test_m),dim(M_weight)[2])
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
  return(M_weight1)
}

Calculate_TFactivity <- function(x_train_m, z_test_m, M_weight){
  GSA_result=Readout[[1]]
  D_normal=matrix(1000,dim(z_test_m)[2],dim(GSA_result)[1])
  rownames(D_normal)=colnames(z_test_m)
  colnames(D_normal)=rownames(GSA_result)
  D_cancer=matrix(1000,dim(z_test_m)[2],dim(GSA_result)[1])
  rownames(D_cancer)=colnames(z_test_m)
  colnames(D_cancer)=rownames(GSA_result)
  #deviatino of gene expression pattern
  Ref_normal=x_train_m[,y==1]
  Ref_cancer=x_train_m[,y==2]
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
      D_normal[temp_i,temp_j] = normal_zscore %*% M_weight[,temp_j] #标准化*TF权重:样品*TF
      D_cancer[temp_i,temp_j] = cancer_zscore %*% M_weight[,temp_j]
    }
  }
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
  return(list(A_TFactivity = A_TFactivity, sign_TF = sign_TF))
}

Calculate_TFcluster <- function(A_TFactivity = A_TFactivity, sign_TF = sign_TF){
  TFcluster=matrix(0,dim(A_TFactivity)[1],2)
  colnames(TFcluster)=c('T_pro','T_supp')
  rownames(TFcluster)=rownames(A_TFactivity)
  temp_normalTF=A_TFactivity[,(sign_TF==(-1))]
  temp_cancerTF=A_TFactivity[,(sign_TF==(1))]
  TFcluster[,1]=rowMeans(temp_cancerTF)
  TFcluster[,2]=rowMeans(temp_normalTF)
  return(TFcluster)}


