#如果你下载的是TCGA数据，需要先将基因表达矩阵的行名基因ID转成gene symbol，便于检索和画图时显示基因名。
#not_easy_input_expr.txt，跟FigureYa71ssGSEA里的not_easy_input_expr.txt是同一个文件。

tcga_expr <- data.table::fread("not_easy_input_expr.txt")
names(tcga_expr)[1] <- "gene_id"
tcga_expr[1:3,1:3]


#读取gtf文件来转换gene_id。跟FigureYa71一样，下载v22版本gtf文件：<ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_22/gencode.v22.annotation.gtf.gz>
gtf_v22 <- rtracklayer::import("gencode.v22.annotation.gtf")
gtf_v22 <- as.data.frame(gtf_v22)

tcga_expr <- gtf_v22 %>% 
  filter(type == "gene",gene_type ==  "protein_coding") %>% 
  dplyr::select(c("gene_id", "gene_name")) %>% 
  separate(gene_id,into = c("gene_id"),sep = "\\.") %>% 
  distinct(gene_name,.keep_all = T) %>% 
  inner_join(tcga_expr,by="gene_id") %>% 
  select(-gene_id) 
tcga_expr[1:3, 1:3]

#选一个基因，例如FOXP3
write.table(tcga_expr["FOXP3",], file = "easy_input_expr.txt",sep = "\t", quote = F)