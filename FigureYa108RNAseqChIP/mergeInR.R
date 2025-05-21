aa <- read.table("Heatmap1sortedRegions.bed",head=F)
bb <- read.table("join.txt",head=F)
cc <- merge(aa,bb,by.x="V4",by.y="V2")
dd <- cc[,c(2,3,4,1,5:12,14)]
head(dd)
write.table(dd,"xx.bed",sep="\t",quote=F,row.names=F,col.name=F)
