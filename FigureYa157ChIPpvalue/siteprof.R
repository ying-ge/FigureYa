#! /usr/bin/env Rscript

# ----- initiation -----
## choose by platform
## macOSX
#download.file("http://hgdownload.soe.ucsc.edu/admin/exe/macOSX.x86_64/bigWigSummary","bigWigSummary")
#download.file("http://hgdownload.soe.ucsc.edu/admin/exe/macOSX.x86_64/bigWigInfo","bigWigInfo")
## linux x86_64
#download.file("http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/bigWigSummary","bigWigSummary")
#download.file("http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/bigWigInfo","bigWigInfo")
## linux x86_64 v369
#download.file("http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64.v369/bigWigSummary","bigWigSummary")
#download.file("http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64.v369/bigWigInfo","bigWigInfo")

Sys.chmod(c("bigWigSummary","bigWigInfo"),755)


source("./corelib.R")

# ----- parameters -----
outputName <- "multiSitepro_example.pdf"
bigWigFile <- c("8cell.K4me3.bw")
bedFiles <- c("mm9_HCP_tss.bed","mm9_ICP_tss.bed","mm9_LCP_tss.bed")
labels <- c("HCP","ICP","LCP")
siteType <- "TSS"
resolution <- 10
span <- 2000
coreNumber <- detectCores()
normalization_constant <- 0 # 0 means use Fold for nromalization constant
colors <- c("darkorange","navy","darkgreen")
captureMethod="security" # security or speed, for bigwig generate from genomeCoverageBed with -bga option or MACS2, you can use speed method

# ----- capture signal -----
average_signal <- c()
for(i in 1:length(bedFiles)){
  
  # load bed file
  bed <- read.table(bedFiles[i])
  isNormalChrosome <- !grepl("_",bed$V1)
  bed <- bed[isNormalChrosome,]

  # caputre signal
  average_signal <- rbind(average_signal,signal_caputer_around_sites(bigWigFile,bed,resolution=resolution,span=span,captureMethod=captureMethod,cores=as.integer(coreNumber)))
}

if(normalization_constant == 0) {
  normalization_constant <- as.numeric(system(paste("./bigWigInfo", bigWigFile , "| grep mean | cut -d ' ' -f 2"), intern=T))
}
average_signal <- average_signal / normalization_constant


# ----- plot -----

minimum <- min(average_signal) - (max(average_signal) - min(average_signal)) * 0.1
maximum <- max(average_signal) + (max(average_signal) - min(average_signal)) * 0.1 * length(bedFiles)

datapoints <- span / resolution

pdf(outputName,width = 5,height = 5)
plot(1:(datapoints*2+1),average_signal[1,],type="l",col=colors[1],
     xaxt="n",xlab="",ylab="Normalized signal",ylim=c(minimum,maximum))
for(i in 2:length(bedFiles)){
  lines(1:(datapoints*2+1),average_signal[i,],col=colors[i])
}
legend("topleft",col=colors,legend=labels,lty=1)
axis(side=1,at=c(0,datapoints,datapoints*2)+1,labels=c(paste(-1*span/1000,"kb"),siteType,paste(span/1000,"kb")))
dev.off()
