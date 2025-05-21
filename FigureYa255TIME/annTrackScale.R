annTrackScale <- function(indata=NULL, halfwidth=NULL, poolsd=FALSE) {
  if ( !(is.vector(indata) | is.matrix(indata) | is.data.frame(indata)) ) {
    stop("indata must be vector, matrix or data frame!")
  }
  
  if (is.vector(indata)) {
    if (sd(indata)!=0) {
      outdata = indata-median(indata, na.rm=T)  
      
      if (poolsd==FALSE) {      
        pos <- outdata[outdata>=0]
        neg <- outdata[outdata<0]
        sdpos <- sd(pos, na.rm = T)
        sdneg <- sd(neg, na.rm = T)
        pos <- pos/sdpos
        neg <- neg/sdneg
        outdata <- c(pos, neg)[names(outdata)]
      }else{    
        std <- sd(outdata, na.rm = T)
        outdata <- outdata/std
      }
      
      if (!is.null(halfwidth)) {
        outdata[outdata>halfwidth]=halfwidth
        outdata[outdata<(-halfwidth)]= -halfwidth
        outdata[outdata<0] <- halfwidth*outdata[outdata<0]/abs(min(outdata))
        outdata[outdata>0] <- halfwidth*outdata[outdata>0]/max(outdata)
      }
    }else{
      outdata <- rep(0, times=length(indata))
    }
    
  }
  
  if (is.matrix(indata) | is.data.frame(indata)) {    
    outdata = sweep(indata,1, apply(indata,1,median,na.rm=T))
    for (m in 1:nrow(outdata)) {
      tmp <- as.numeric(outdata[m, ]); names(tmp) <- colnames(outdata)
      
      if(sd(tmp)!=0) {
        if (poolsd==FALSE) {
          pos <- tmp[tmp>=0]
          neg <- tmp[tmp<0]
          sdpos <- sd(pos, na.rm = T)
          sdneg <- sd(neg, na.rm = T)
          pos <- pos/sdpos
          neg <- neg/sdneg
          outdata[m,] <- c(pos, neg)[names(tmp)]       
        }else{
          std <- sd(tmp, na.rm = T)
          outdata[m,] <- tmp/std
        }      
      }else{
        outdata[m,] <- rep(0, times=length(tmp))
      }
    }    
    
    if (!is.null(halfwidth)) {
      outdata[outdata>halfwidth]=halfwidth
      outdata[outdata<(-halfwidth)]= -halfwidth
      for (k in 1:nrow(outdata)) {
        rowk <- outdata[k, ]
        rowk[rowk<0] <- halfwidth*rowk[rowk<0]/abs(min(rowk))
        rowk[rowk>0] <- halfwidth*rowk[rowk>0]/max(rowk)
        outdata[k,] <- rowk
      }
    }
  }
  
  return(outdata)
}