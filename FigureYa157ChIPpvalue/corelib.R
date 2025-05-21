#! /usr/bin/env Rscript

regionCapture <- function(seqname,start,end,bw,datapoints,captureMethod="security",cores=4L) {
  # capture bigwig signal in region with given datapoints
  captureRegion <- data.frame(seqname=seqname,start=start,end=end)
  write.table(captureRegion,file="captureRegion.bed",quote=F,sep="\t",row.names=F,col.names=F)
  py_version <- system('python -V 2>&1', intern = T)
  py_version <- substr(py_version,8,nchar(py_version))
  py_version <- unlist(strsplit(py_version,'[.]'))
  if(py_version[1]=='3' && as.integer(py_version[2])>=6) {
    suppressMessages(system(paste("python getBigWigValue.py -b", "captureRegion.bed", "-w", bw, "-n regionCapture", "-p", cores,"-s", datapoints, "-m", captureMethod)))
  } else {
    if(py_version[1]=='2' && py_version[2]=='7') {
      suppressMessages(system(paste("python getBigWigValue_py2.py -b", "captureRegion.bed", "-w", bw, "-n regionCapture", "-p", cores,"-s", datapoints, "-m", captureMethod)))
    } else {
      stop('python version is either 2.7 nor 3.6+')
    }
  }
  signal <- read.table(gzfile("signal_regionCapture_siteprof1.gz"))
  delete_result <- file.remove("signal_regionCapture_siteprof1.gz")
  if(captureMethod=="security"){
    coverage <- read.table(gzfile("coverage_regionCapture_siteprof1.gz"))
    signal <- signal * coverage
    delete_result <- file.remove("coverage_regionCapture_siteprof1.gz")
  }
  
  return(signal)
}


signal_caputer_around_sites <- function(bigWigFile,bed,resolution=50,span=3000,captureMethod="security",cores=4L) {
  # options
  # bigWigFile: bigWigFile to be capture
  # bedFile: bedFile contain sites information
  # resolution: resolution for capture
  # span: span for upstream, downstream & gene body normalization length
  # cores: cores used for parallel
  
  datapoints <- span * 2 / resolution + 1
  
  seqname <- bed[,1]
  midpoints <- as.integer((bed[,2] + bed[,3] - 1) / 2)
  # get strand information
  if(ncol(bed)>5) { 
    strand <- bed[,6]
  } else {
    strand <- rep(".",nrow(bed))
  }
  
  # capture
  # IRange is 1-based
  signal <- regionCapture(seqname,midpoints-span-resolution/2,midpoints+span+resolution/2,bigWigFile,datapoints,captureMethod,cores)
  
  # correct orientation & get average signal
  write.table(signal,file="temp_signal.txt",quote=F,row.names=F,col.names=F,sep="\t")
  write.table(strand,file="temp_strand.txt",quote=F,row.names=F,col.names=F,sep="\t")
  system(paste0("paste temp_signal.txt temp_strand.txt | awk '{if($NF==",'"',"+",'"',") {for(i=1;i<NF-1;i++) printf $i",'"',"\t",'"',"; print $(NF-1)} else {for(i=NF-1;i>1;i--) printf $i",'"',"\t",'"',"; print $1}}' > temp_signal_rev.txt"))
  signal_rev <- read.table("temp_signal_rev.txt")
  delete_result <- file.remove("temp_signal.txt")
  delete_result <- file.remove("temp_strand.txt")
  delete_result <- file.remove("temp_signal_rev.txt")
  return(apply(signal_rev,2,mean))
  
}

average_signal_caputer_around_region <- function(bigWigFile,bed,captureMethod="security",cores=4L) {
  # options
  # bigWigFile: bigWigFile to be capture
  # bedFile: bedFile contain sites information
  # resolution: resolution for capture
  # span: span for upstream, downstream & gene body normalization length
  # cores: cores used for parallel
  
  datapoints <- 1
  
  seqname <- bed[,1]
  start <- bed[,2]
  end <- bed[,3]
  
  # capture
  signal <- regionCapture(seqname,start,end,bigWigFile,datapoints,captureMethod,cores)
  
  return(signal)
  
}
