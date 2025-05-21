library(BSgenome.Hsapiens.UCSC.hg19)
library(VariantAnnotation)
library(rgl)

get.complementary <- function(x) {
  if (x == "A") {
    y <- "T"
  } else if (x == "T") {
    y <- "A"
  } else if (x == "C") {
    y <- "G"
  } else if (x == "G") {
    y <- "C"
  }
  return(y)
}

get.mutationcontext <- function(x, genome, chrom, start, reverse) {
  chrom <- paste0("chr", as.character(x[chrom]))
  start <- as.numeric(x[start])
  reverse <- x[reverse]
  contextSeq <- getSeq(genome, chrom, start-1, start+1)
  contextSeq <- unlist(strsplit(as.character(contextSeq),""))
  if (reverse){
      return(paste0(get.complementary(contextSeq[3]), "_", get.complementary(contextSeq[1])))
    }else{
      return(paste0(contextSeq[1], "_", contextSeq[3]))
  }
}

get.mutationtype <- function(x, ref, alt) {
  ref <- x[ref]
  alt <- x[alt]
  if (ref %in% c("T", "C")) {
    return(list(paste0(ref, ">", alt), FALSE))
  } else {
    return(list(paste0(get.complementary(ref), ">", get.complementary(alt)), TRUE))
  }
}

stackplot.3d <- function(x, y, z, alpha=1, topcol="#078E53", sidecol="#aaaaaa"){
  
  ## These lines allow the active rgl device to be updated with multiple changes
  ## This is necessary to draw the sides and ends of the column separately  
  save <- par3d(skipRedraw=TRUE)
  on.exit(par3d(save))
  
  ## Determine the coordinates of each surface of the column and its edges
  x1 <- c(rep(c(x[1],x[2],x[2],x[1]),3),rep(x[1],4),rep(x[2],4))
  y1 <- c(y[1],y[1],y[2],y[2],rep(y[1],4),rep(y[2],4),rep(c(y[1],y[2],y[2],y[1]),2))
  z1 <- c(rep(0,4),rep(c(0,0,z,z),4))
  x2 <- c(rep(c(x[1],x[1],x[2],x[2]),2),rep(c(x[1],x[2],rep(x[1],3),rep(x[2],3)),2))
  y2 <- c(rep(y[1],4),rep(y[2],4),rep(c(rep(y[1],3),rep(y[2],3),y[1],y[2]),2) )
  z2 <- c(rep(c(0,z),4),rep(0,8),rep(z,8) )
  
  ## These lines create the sides of the column and its coloured top surface
  rgl.quads(x1,y1,z1,col=rep(sidecol,each=4),alpha=alpha,lit=FALSE)
  rgl.quads(c(x[1],x[2],x[2],x[1]),c(y[1],y[1],y[2],y[2]),rep(z,4),
            col=rep(topcol,each=4),alpha=1,lit=FALSE) 
  ## This line adds black edges to the column
  rgl.lines(x2, y2, z2, col="#000000",lit=FALSE)
}


legend.3d <- function(x, y, z, label_z, label, alpha=1, topcol="white", sidecol="white"){
  
  ## These lines allow the active rgl device to be updated with multiple changes
  ## This is necessary to draw the sides and ends of the column separately  
  save <- par3d(skipRedraw=TRUE)
  on.exit(par3d(save))
  
  ## Determine the coordinates of each surface of the column and its edges
  x1 <- c(rep(c(x[1],x[2],x[2],x[1]),3),rep(x[1],4),rep(x[2],4))
  y1 <- c(y[1],y[1],y[2],y[2],rep(y[1],4),rep(y[2],4),rep(c(y[1],y[2],y[2],y[1]),2))
  z1 <- c(rep(0,4),rep(c(0,0,z,z),4))
  x2 <- c(rep(c(x[1],x[1],x[2],x[2]),2),rep(c(x[1],x[2],rep(x[1],3),rep(x[2],3)),2))
  y2 <- c(rep(y[1],4),rep(y[2],4),rep(c(rep(y[1],3),rep(y[2],3),y[1],y[2]),2) )
  z2 <- c(rep(c(0,z),4),rep(0,8),rep(z,8) )
  
  ## These lines create the sides of the column and its coloured top surface
  rgl.quads(x1,y1,z1,col=rep(sidecol,each=4),alpha=alpha,lit=FALSE, shininess=0)
  rgl.quads(c(x[1],x[2],x[2],x[1]),c(y[1],y[1],y[2],y[2]),rep(z,4),
            col=rep(topcol,each=4),alpha=1,lit=FALSE, shininess=0) 
  ## This line adds black edges to the column
  rgl.lines(x2,y2,z2,col="#000000",lit=FALSE)
  rgl.texts((x[1]+x[2])/2, (y[1]+y[2])/2, label_z, text=label, color="black", lit=FALSE, adj=c(0.5, 0.8))
}


legoplot.3d <- function(z, alpha=1, scalexy=10, scalez=1, gap=0.1, zlab="Mutations", 
                        palettecolors=c("#805D3F", "#72549A", "#5EAFB2", 
                                        "#3F4F9D", "#F2EC3C", "#74B655")){
  ## These lines allow the active rgl device to be updated with multiple changes
  ## This is necessary to add each column sequentially
  save <- par3d(skipRedraw=TRUE)
  on.exit(par3d(save))

  ## Recreate Broad order
  types <- c("C>T", "C>A", "C>G", "T>G", "T>C", "T>A")
  contexts <- c("T_T", "C_T", "A_T", "G_T", "T_C", "C_C", "A_C", "G_C",
             "T_A", "C_A", "A_A", "G_A", "T_G", "C_G", "A_G", "G_G")
  typeorder <- c()
  for(type in types){
    typeorder <- c(typeorder, paste(type, contexts, sep="."))
  }
  z <- z[typeorder]
  
  ## Reorder data into 6 regions
  set1 <- c(20:17, 4:1, 24:21, 8:5, 28:25, 12:9, 32:29, 16:13)
  set2 <- set1 + 32
  set3 <- set1 + 64
  neworder <- c(set1, set2, set3)
  
  ## Scale column area and the gap between columns 
  x <- seq(1, 8) * scalexy
  y <- seq(1, 12) * scalexy
  gap <- gap * scalexy
  
  ## Scale z coordinate
  z <- z * scalez
  
  ## Set up colour palette
  legocolors <- as.vector(sapply(palettecolors, rep, 16))
  
  ## Plot each of the columns
  for(i in 1:12){
    for(j in 1:8){
      it=(i-1)*8+j # Variable to work out which column to plot; counts from 1:96
      stackplot.3d(c(gap + x[j], x[j] + scalexy),
                   c(-gap - y[i], -y[i] - scalexy),
                   z[neworder[it]],
                   alpha=alpha,
                   topcol=legocolors[neworder[it]],
                   sidecol=legocolors[neworder[it]])
    }
  }

  ## Plot lego legend
  legend_x <- seq(1, 4) * 1.2 * scalexy
  legend_y <- seq(1, 12) * 1.2 * scalexy
  legend_z <- 0.1 * scalexy
  legend_label_z <- 0.4 * scalexy
  a <- 0
  for(i in 5:8){
    for(j in 1:4){
      a = a + 1
      legend.3d(c(-gap - legend_x[j], -legend_x[j] - scalexy * 1.2),
                c(-gap - legend_y[i], -legend_y[i] - scalexy * 1.2),
                legend_z,
                label_z=legend_label_z,
                label=contexts[a],
                alpha=alpha,
                topcol="white",
                sidecol="white")
    }
  }

  # title3d(xlab='xlab', ylab='ylab', zlab='zlab')
  # axes3d(xlab="x", ylab="y",zlab="z")
  # add x,y line and label
  rgl.lines(c(-2.5 * scalexy,-5.5 * scalexy), c(-5.5 * scalexy, -5.5 * scalexy), c(0.1 * scalexy,0.1 * scalexy), col="#000000",lit=FALSE)
  rgl.lines(c(-6.8 * scalexy,-6.8 * scalexy), c(-7.0 * scalexy, -10.0 * scalexy), c(0.1 * scalexy,0.1 * scalexy), col="#000000",lit=FALSE)
  rgl.texts(-4.0 * scalexy, -5.0 * scalexy, 0.1 * scalexy, "X", color="black", adj=c(0.2, 0.8))
  rgl.texts(-7 * scalexy, -8.5 * scalexy, 0.1 * scalexy, "Y", color="black", adj=c(0.2, 0.8))

  ## Set the viewpoint and add axes and labels
  rgl.viewpoint(theta=50,phi=40,fov=0)
  axes3d("z++",labels=TRUE)
  mtext3d(zlab, edge="z++", line=2)
}
