help = cat(
" 
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
Help text here 
Arguments in this order: 
    1) input file (chr\tstart\tend\treadsCount\tgeneID\tgene_windowNumber)
    2) IP file (gene_windowNumber\t readsCount)
    3) Totle reads in input
    4) Totle reads in ip
    5) PEAK bed file 
    6) all result
Author: Zhao Long zhaolong@genetics.ac.cn
Last Updata:2020 6 27 
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
\n\n") 


args<-commandArgs(T)

library(tidyverse)
cat(">>> Start...\n[==================================================]   0%", format(Sys.time(), "%Y-%m-%d %X"),"\n")

ipt<-read_tsv(args[1],
	col_names=c("chr","start","end","ipt","gene","transcript"),
	col_types=cols(
		chr=col_character(),start=col_integer(),end=col_integer(),ipt=col_integer(),gene=col_character(),transcript=col_character()
		)
	)
	
ip<-read_tsv(args[2],col_names=c("transcript","ip"),col_types=cols(transcript=col_character(),ip=col_integer()))

all<-inner_join(ipt,ip,by="transcript")

cat(">>> Read data Done...\n[>>>>>>>>>>========================================]  20%", format(Sys.time(), "%Y-%m-%d %X"),"\n")

all$ipt_totle<-args[3]
all$ip_totle<-args[4]

fisher_fun<-function(x){
	y<-c(as.numeric(x['ip']),as.numeric(x['ip_totle'])-as.numeric(x['ip']),
		as.numeric(x['ipt']),as.numeric(x['ipt_totle'])-as.numeric(x['ipt']))
	pvalue<-fisher.test(array(y,dim=c(2,2)))$p.value
	return(pvalue)
}

all$pvalue<-apply(all,1,fisher_fun)	
p.adjust(all$pvalue,method="fdr")->all$FDR

cat(">>> Fisher exact test Done...\n[>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>====================]  60%", format(Sys.time(), "%Y-%m-%d %X"),"\n")

all_median <- all %>% group_by(gene) %>% summarise(ip_median=median(ip),ipt_median=median(ipt),.groups = 'drop')
all <- all %>% left_join(all_median) %>% mutate(enrich_score=(ip*ipt_median)/(ipt*ip_median))

cat(">>> Enrich score Done...\n[>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>==========]  80%", format(Sys.time(), "%Y-%m-%d %X"),"\n")

peak_bed <- all %>% filter(FDR<=0.05 & log2(enrich_score) > 1) %>% select(chr,start,end,gene) %>% arrange(chr,start,end)

write.table(peak_bed,args[5],sep="\t",row.names=F,quote=F,col.names=F)
write.table(all,args[6],sep="\t",quote=F,row.names=F)
cat(">>> All jobs Done...\n[>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>] 100%", format(Sys.time(), "%Y-%m-%d %X"),"\n")

