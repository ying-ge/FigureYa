Prepare_geodataset <- function(accession, destdir){
  library(GEOquery)
  gse <- getGEO(accession, destdir = destdir)
  gse <- gse[[1]]
  dat <- exprs(gse)
  pdata <- pData(gse)
  #anti-log if max < 50 in mixture file
  if(max(dat) < 50) {dat <- 2^dat}
  return(list(dat = dat, pdata = pdata))}

exprSet<-function(GEOexpr,GPL){ 
  library(GEOquery);library(stringr)
  ID<-rownames(GEOexpr);
  if (!all(is.numeric(GEOexpr[, 1]))){
    message("Column 1 in GEOexpr is not numeric, assume as probe name")
    colnames(GEOexpr)[1] <- "ID"}else{GEOexpr<-data.frame(ID,GEOexpr)}
  if (length(intersect(GEOexpr[, "ID"], GPL[, "ID"])) == 0){
    stop("Error: no intersection probeID between GEOexpr and the GPL")
  }
  GEOexpr<-merge(GPL,GEOexpr,by="ID");
  GEOexpr<-GEOexpr[,-1];
  formal<-as.formula(paste0(".~",colnames(GPL)[2]))
  GEOexpr<-aggregate(formal,GEOexpr,mean);
  rownames(GEOexpr)<-GEOexpr[,1];
  GEOexpr<-GEOexpr[,-1];
  GEOexpr<-GEOexpr[!str_detect(rownames(GEOexpr),"///"),]
  if (rownames(GEOexpr)[1] == "") {
    GEOexpr <- GEOexpr[-1, ]
  }
  return(GEOexpr)}
#GSE6764
gse6764 <-  Prepare_geodataset(accession = "GSE6764", destdir = getwd())
gse6764dat <- gse6764$dat
gse6764pdata <- gse6764$pdata
GPL570 <- getGEO("GPL570")
GPL570 <- Table(GPL570)[, c("ID", "Gene Symbol")]
colnames(GPL570) <- c("ID", "Gene_Symbol")
gse6764dat <- exprSet(GEOexpr = gse6764dat, GPL = GPL570)
#GSE56140
gse56140 <- Prepare_geodataset(accession = "GSE56140", destdir = getwd())
gse56140dat <- gse56140$dat 
GPL18461 <- getGEO("GPL18461")
GPL18461 <- Table(GPL18461)[, c("ID", "ILMN_GENE")]
gse56140dat <- exprSet(GEOexpr = gse56140dat, GPL = GPL18461)
gse56140pdata <- gse56140$pdata
datlist <- list(gse6764dat = gse6764dat, gse56140dat = gse56140dat)




