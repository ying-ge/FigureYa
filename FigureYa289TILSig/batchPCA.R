batchPCA =function(indata, batch, fig.dir, PCA.fig.title, pos="bottomright", xy=c(1,2), cols=NULL, showID=FALSE, cex=1, showLegend=T) {
# indata is a data matrix with samples in columns and genes in rows.
# batch is a vector with the order matching the order in indata.
    library(ClassDiscovery)
    
    outfile = file.path(fig.dir, paste(PCA.fig.title, ".pdf",sep=""))
    N.batch = length(unique(batch))    
    if (is.null(cols)) { 
      cols <- rainbow(N.batch) 
    }else{
      if (length(cols) != N.batch) {stop("cols length not equal to batch length")}
    }           
    
    indata=na.omit(indata)
    pca<-SamplePCA(indata, usecor=F, center=T)
    pct1 <- round (pca@variances[xy[1]]/sum(pca@variances), digits=3)*100
    pct2 <- round (pca@variances[xy[2]]/sum(pca@variances), digits=3)*100
    xlab.text = paste("Comp ", xy[1], ": ", as.character(pct1), "% variance", sep="")
    ylab.text = paste("Comp ", xy[2], ": ", as.character(pct2), "% variance", sep="")    
    
    pdf(file=outfile)
    plot(pca@scores[,xy[1]], pca@scores[,xy[2]],  cex=0.7, xlab=xlab.text, ylab=ylab.text, col=cols[factor(batch)], pch=(1:N.batch)[factor(batch)],lwd=1.5, main=PCA.fig.title)
    abline(h=0, v=0, col="brown", lty=2)
    abline(h=0, v=0, col="brown", lty=2)
    center1<-tapply(pca@scores[,xy[1]], factor(batch), mean)
    center2<-tapply(pca@scores[,xy[2]], factor(batch), mean)
    for (ii in 1:length(center1)) {
        groupi<-pca@scores[as.numeric(factor(batch))==ii, xy]
        #  print(paste("Cluster", ii))
        if (inherits(groupi, "matrix")) {
            for (j in (1:nrow(groupi))) {
                segments( groupi[j,1], groupi[j,2], center1[ii], center2[ii], col=cols[ii] , lwd=0.3)
            }
        }else {
            segments( groupi[1], groupi[2], center1[ii], center2[ii], col=cols[ii] , lwd=0.3)
        }
    }
    points(center1, center2, pch=7, lwd=1.5,col=cols)
    if (showID) {
      text(pca@scores[,xy[1]], pca@scores[,xy[2]], colnames(indata), lwd=1, cex=cex)
    }
    if(showLegend){
      legend(pos,legend=names(table(factor(batch))), text.col=cols, pch=(1:N.batch), col=cols, lty=1)
    }
    invisible(dev.off())
}