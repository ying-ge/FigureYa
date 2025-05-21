### 画法二：ggplot2

# Remove upper value and prepare inupt file for ggplot2
est <- HyperGeoTest$estimate
p <- HyperGeoTest$pVal
est[est==Inf]= NA
p[p < 1.0e-150]=NA #or p[p == 0 ] = NA

est <- as.matrix(est) # estimate value
p <- as.matrix(p) # pvalue
r <- reorder(est,o)
p <- reorder(p,o)
r[lower.tri(r)]<- NA
p[lower.tri(p)]<- NA

# Reshape data from matrix to data.frame
library(reshape2)
est_melt <- melt(r)
p_melt <- melt(p)
colnames(est_melt)<- c("rowv","columnv","est")
est_melt$p<- p_melt$value

###########
# discrete value for estimate value (automatically)
###########
c<- est_melt$est
max_value<- range(c[!is.na(c)])
brks<- c(0,seq(1,20,l=10),seq(21,max_value[2],l=4))
labs<- c(as.character(round(c(1,seq(2,20,l=9),seq(21,max_value[2],l=4)))))

ranges<- cut(c, breaks=brks,include.lowest=TRUE, label=labs)

# Distribution of enrichment value
brks_num<- as.numeric(table(cut(c, breaks=c(0,seq(1,20,l=10),seq(21,max_value[2],l=4)))))
plot(1:length(levels(ranges)),y=brks_num,type="h")

# Setup up levels for enriched estimate value
est_melt$esitmate<- ranges
est_melt$esitmate<- factor(est_melt$esitmate, levels = rev(levels(ranges)))
levels(est_melt$esitmate)

###########
# discrete value for pvalue (automatically)
###########
c<- est_melt$p
c[is.na(c)]=1
ranges<- cut(c, breaks=c(1 %o% c(10^-(seq(1,max_value[2],l=10)),0)),include.lowest=TRUE, label=c(as.character(c(1 %o% c(10^-(seq(1,max_value[2],l=9)),0)))))

# log seq
# lseq <- function(from=1, to=100000, length.out=6) {
# logarithmic spaced sequence blatantly stolen from library("emdbook"), because need only this
# exp(seq(log(from), log(to), length.out = length.out))}

est_melt$pval<- ranges
est_melt$pval<- factor(est_melt$pval, levels = levels(ranges))
levels(est_melt$pval)

###########
# !!! discrete value for pvalue (manually)
###########
range(est_melt$p)
brks<- c(0, 10^-100, 10^-80, 10^-60,10^-40, 10^-20,10^-10,10^-5,10^-2, 1)
labs<- c("10^-100", "10^-80", "10^-60","10^-40", "10^-20", "10^-10","10^-5","10^-2", "1")
est_melt$pval <- cut(est_melt$p, breaks=brks, label=labs)  # Create column of significant labels

# Setup up order sequence
est_melt$rowv<- factor(est_melt$rowv, levels = levels(est_melt$rowv))
est_melt$columnv<- factor(est_melt$columnv, levels = levels(est_melt$rowv))

###########
# !!! Set up cluster information by manually
###########
est_melt$groups <- paste0(est_melt$rowv,"_",est_melt$columnv)
cluster <- rep("no_sig",length(est_melt$groups))
names(cluster) <- est_melt$groups
cluster_1 <- c("SRF_HNF4A","CREB1_HNF4A","RXRA_HNF4A","PPARG_HNF4A","FOXO1_HNF4A","SRF_FOXO1","CREB1_FOXO1","RXRA_FOXO1","PPARG_FOXO1","SRF_PPARG","CREB1_PPARG","RXRA_PPARG","SRF_PPARG","CREB1_PPARG","RXRA_PPARG","SRF_RXRA","CREB1_RXRA","SRF_CREB1")
cluster_2 <- c("JUN_PGR","CEBPA_PGR","CEBPB_PGR","NR3C1_PGR","AR_PGR", "JUN_AR","CEBPA_AR","CEBPB_AR","NR3C1_AR", "JUN_NR3C1","CEBPA_NR3C1","CEBPB_NR3C1", "JUN_CEBPB","CEBPA_CEBPB","JUN_CEBPA")

cluster[names(cluster) %in% cluster_1]<- "cluster_1"
cluster[names(cluster) %in% cluster_2]<- "cluster_2"

est_melt$clusters<- as.character(cluster)
est_melt$clusters<- factor(est_melt$clusters, levels = c("no_sig","cluster_1","cluster_2"))

## Final input for ggplot2
dat<- est_melt[!is.na(est_melt$est),]

p1<- ggplot(dat, aes(y = factor(rowv),x = factor(columnv))) +        
  geom_tile(aes(fill = pval),colour=dat$clusters,width=0.9, height=0.9,size=0.2) + 
  scale_color_manual(values=adjustcolor(c(colorRampPalette(brewer.pal(9,"RdBu"))(length(levels(dat$esitmate))),"white","darkgreen","darkred"), alpha.f = .95),limits = c(levels(dat$esitmate),levels(dat$clusters)))+
  scale_fill_manual("p-value",values=adjustcolor(rev(colorRampPalette(c("black","white"))(length(labs))), alpha.f = 0.7)) +
  labs(color = "Enrichment Value") +
  geom_point(aes(color=esitmate,size=esitmate)) +
  scale_size_manual("Enrichment Value by Size",values=(rev(seq(1,length(levels(dat$esitmate)),1))*0.2)) +
  xlab("TF genes") + 
  ylab("TF genes") +
  theme(legend.position="topright") + 
  ## remove backgroud
  theme_bw(base_size = 13) +
  theme(axis.text.x=element_text(angle=90,hjust=1)) +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))

pdf("Enrichment_ggplot2.v1.pdf", 14, 12)
p1
dev.off()

### 画法三：ggplot2

# Create column of significance labels
range(dat$p)
dat$pval <- cut(dat$p, breaks=c(0, 10^-50, 10^-35, 10^-25, 1), label=c("**", "*", "-", ""))  # Create column of significance labels

p2<- ggplot(dat, aes(y = factor(rowv),x = factor(columnv))) +        
  geom_tile(aes(fill = esitmate)) +  
  geom_tile(aes(fill = NA, colour=dat$clusters),width=0.85, height=0.85,size=0.3) + 
  scale_fill_manual("Enrichment Value",values=adjustcolor(colorRampPalette(brewer.pal(9,"RdBu"))(length(levels(dat$esitmate))), alpha.f = .95),limits = levels(dat$esitmate))+
  scale_color_manual("Cluster Groups",values=adjustcolor(c("white","darkgreen","darkred"), alpha.f = .95))+
  geom_text(aes(label=pval), color="black", size=5) + 
  labs(color = "p-value") +
  xlab("TF genes") + 
  ylab("TF genes") +
  ## remove backgroud
  theme_bw(base_size = 13) +
  theme(axis.text.x=element_text(angle=90,hjust=1)) +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
#geom_rangeframe() 
#theme_void()
#theme_tufte()

pdf("Enrichment_ggplot2.v2.pdf", 14, 12)
p2
dev.off()
