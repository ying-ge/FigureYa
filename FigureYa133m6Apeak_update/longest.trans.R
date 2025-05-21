args<-commandArgs(T)
library(dplyr)

top1_func<-function(data){
  data_sort <- arrange(data,desc(V7))
  data_sort_top1 <-data_sort[1,]
  return(data_sort_top1)
}

ref=read.table(args[1],head=F)

ref_list<-split(ref,ref$V6)
ref_list_top1<-lapply(ref_list,top1_func)
ref_ann<-do.call(rbind,lapply(ref_list_top1,data.frame))
ref_ann_bed<-arrange(ref_ann[,c(1,2,3,6)],V1,V2)

write.table(ref_ann_bed,args[2],sep="\t",quote=F,row.names=F,col.names=F)
