#! /usr/bin/env Rscript
suppressPackageStartupMessages(require(rtracklayer))
suppressPackageStartupMessages(require(parallel))

regionCapture <- function(seqname,start,end,bw,datapoints) {
  # capture bigwig signal in region with given datapoints
  step <- (end-start) / datapoints
  regions <- GRanges(seqnames=rep(seqname,datapoints), ranges=IRanges(as.integer(seq(start,end,step))[1:datapoints],as.integer(seq(start,end,step))[2:(datapoints+1)]-1))
  return(as.numeric(mean(import(bw,which=regions,as="NumericList"))))
}

orientation_correct <- function(values,strand) {
  if(strand == "-"){
    return(rev(values))
  }
  else{
    return(values)
  }
}


signal_caputer_around_gene <- function(bigWigFile,refGene,resolution=50,span=3000,geneBodySpan=2000,cores=4L) {
  # options
  # bigWigFile: bigWigFile to be capture
  # refGene: genes for capture in refGene format
  # resolution: resolution for capture
  # span: span for upstream & downstream length
  # geneBodySpan: gene body normalization length
  # cores: cores used for parallel
  
  datapoints1 <- span / resolution
  datapoints2 <- geneBodySpan / resolution
  
  # capture
  # IRange is 1-based
  signalBody <- suppressWarnings(mcmapply(regionCapture,refGene[,3],refGene[,5]+1,refGene[,6],bigWigFile,datapoints2,mc.cores = getOption("mc.cores",cores)))
  signalUpstream <- suppressWarnings(mcmapply(regionCapture,refGene[,3],refGene[,5]-span+1,refGene[,5],bigWigFile,datapoints1,mc.cores = getOption("mc.cores",cores)))
  signalDownstream <- suppressWarnings(mcmapply(regionCapture,refGene[,3],refGene[,6]+1,refGene[,6]+span,bigWigFile,datapoints1,mc.cores = getOption("mc.cores",cores)))
  
  # correct orientation & get average signal
  signal <- data.frame(rbind(signalUpstream,signalBody,signalDownstream))
  
  signal_rev <-  matrix(mcmapply(orientation_correct,signal,refGene[,4],mc.cores = getOption("mc.cores",cores)),datapoints1*2+datapoints2,nrow(refGene),byrow=F)
  return(apply(signal_rev,1,mean))
  
}

signal_caputer_around_sites <- function(bigWigFile,bedFile,resolution=50,span=3000,cores=4L) {
  # options
  # bigWigFile: bigWigFile to be capture
  # bedFile: bedFile contain sites information
  # resolution: resolution for capture
  # span: span for upstream, downstream & gene body normalization length
  # cores: cores used for parallel
  
  datapoints <- span * 2 / resolution + 1
  
  # load bed file
  bed <- read.table(bedFile)
  isNormalChrosome <- !grepl("_",bed$V1)
  seqname <- bed[isNormalChrosome,1]
  midpoints <- as.integer(bed[isNormalChrosome,2]+bed[isNormalChrosome,3]-1) / 2
  if(ncol(bed)>5) { 
    strand <- bed[isNormalChrosome,6]
  } else {
    strand <- rep(".",sum(isNormalChrosome))
  }
  
  # capture
  # IRange is 1-based
  signal <- suppressWarnings(mcmapply(regionCapture,seqname,midpoints-span-resolution/2,midpoints+span+resolution/2,bigWigFile,datapoints,mc.cores = getOption("mc.cores",cores)))

  # correct orientation & get average signal
  signal_rev <- matrix(mcmapply(orientation_correct,signal,strand,mc.cores = getOption("mc.cores",cores)),datapoints,sum(isNormalChrosome),byrow=F)
  return(apply(signal_rev,1,mean))
  
}
