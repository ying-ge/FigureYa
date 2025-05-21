#############################
## the main function
#############################
require(GenomicRanges);

## main ##
circos_plus <- function (mapping=mapping, xc=400, yc=400, R=400, W=W,  
                    cir="", type="n", col.v=3, B=FALSE, 
                    print.chr.lab=TRUE, col.bar=FALSE, col.bar.po = "topleft", cluster=FALSE, order=NULL,
                    scale=FALSE, cutoff = "n", zoom="", cex=1,
                    lwd=1, col=rainbow(10, alpha=0.8)[7], side=""
){
  
  ############################ main #############################
  # initial
  r        <- R;
  
  #### cir ####
  # In cir, the segments (chromosomes) should be sorted.
  cir.m <- 0;
  
  if (is.matrix(cir)==T|is.data.frame(cir)==T|class(cir)[1] == "GRanges") {
    # based on the input file, 
    if (class(cir)[1] == "GRanges"){
      chr.po    <- cbind(as.character(seqnames(cir)),
                         start(cir), names(cir), mcols(cir));   
    } else {
      chr.po    <- cir;
    }
    chr.po[,1]  <- gsub("chr","",chr.po[,1]);
    chr.num     <- nrow(chr.po);
    cir.m       <- 1;
  } else if (nchar(cir) >= 1){
    # by known genome
    pofile.s <- paste("UCSC.", cir, ".RData", sep="");
    pofile   <- system.file("data", pofile.s, package="OmicCircos");
    if (file.exists(pofile)){
    } else {
      stop ("cir name?");
    }
    chr.po     <- local(get(load(pofile)));
    chr.po[,1] <- gsub("chr","",chr.po[,1]);
    chr.num    <- nrow(chr.po);
  } else {
    stop("cir name ??")
  }
  # end cir
  
  ### zoom.in ###
  if (length(zoom)>1){
    chr1     <- which(chr.po[,1]==zoom[1]);
    chr2     <- which(chr.po[,1]==zoom[2]);
    chr.num  <- length(chr1:chr2);
    chr.po   <- zoom.in(cir.in=chr.po, zoom=zoom);
  }
  ### end zoom.in ###
  
  #### chr ###
  if (type == "chr"){
    chr.lw = W; 
    if (cir.m == 0){
      # for chromosomes have cytoband data
      chrfile.s <- paste("UCSC.", cir, ".chr.RData", sep="");
      chrfile   <- system.file("data", chrfile.s, package="OmicCircos");
      if (file.exists(chrfile)){
      } else {
        stop ("chrfile name?");
      }
      dat.c     <- local(get(load(chrfile)));
      
      dat.c[,1] <- gsub("chr", "", dat.c[,1]);
      
      for (chr.i in c(1:chr.num)){
        chr.s  <- chr.po[chr.i,1];
        
        v1 <- as.numeric(chr.po[chr.i,2]);
        v2 <- as.numeric(chr.po[chr.i,3]);
        v3 <- as.numeric(chr.po[chr.i,6]);
        v4 <- as.numeric(chr.po[chr.i,7]);
        
        dat.v <- subset(dat.c, dat.c[,1]==chr.s);
        dat.v <- dat.v[order(as.numeric(dat.v[,2])),];
        for (i in 1:nrow(dat.v)){  
          if (dat.v[i,5]=="gneg"){
            col <- colors()[351];
          } else if (dat.v[i,5]=="acen" | dat.v[i,5]=="gvar" | dat.v[i,5]=="stalk"){
            col <- colors()[26];
          } else {
            col.v <- gsub("gpos","",dat.v[i,5]);
            if (col.v >= 30){
              col <- colors()[300];
            } else {
              col <- colors()[351];
            }
          }
          w1 <- scale.v(dat.v[i,2], v1, v2, v3, v4);
          w2 <- scale.v(dat.v[i,3], v1, v2, v3, v4);
          draw.arc.s(xc, yc, r, w1, w2, col=col, lwd=chr.lw);
        }
        if (print.chr.lab){
          v1  <- as.numeric(chr.po[chr.i,2]);
          v2  <- as.numeric(chr.po[chr.i,3]);
          w.m <- (v1+v2)/2;
          r.j <- r/15;  #Ya#chromosome label, r/20 to r/15
          chr.t <- gsub("chr", "", chr.s);
          draw.text.rt(xc, yc, r+r.j, w.m, chr.t, cex=0.6);
        } 
        if (scale){
          total.num <- as.numeric(chr.po[nrow(chr.po),5]);
          do.scale.cir(xc=xc, yc=yc, the.r=r+3, total.num=total.num, #Ya#scale, r+10 to r+3
                       col="blue", lwd=0.001,
                       V1=v1, V2=v2, V3=v3, V4=v4);
        }
      }
    } else {
      if (!is.null(col)){
        col <- rep(col, nrow(chr.po))[c(1:nrow(chr.po))];    
      } else {
        col <- rep("blue", nrow(chr.po));
      }
      for (chr.i in c(1:chr.num)){
        w1 <- as.numeric(chr.po[chr.i,2]);
        w2 <- as.numeric(chr.po[chr.i,3]);
        draw.arc.s(xc, yc, r, w1, w2, col=col[chr.i], lwd=4);
        if (print.chr.lab){
          w.m <- (w1+w2)/2;
          r.j <- r/20;
          chr.t <- gsub("chr", "", chr.po[chr.i,1]);
          draw.text.rt(xc, yc, r+r.j, w.m, chr.t, cex=0.6);
        } 
        
        if (scale){
          v1 <- as.numeric(chr.po[chr.i,2]);
          v2 <- as.numeric(chr.po[chr.i,3]);
          v3 <- as.numeric(chr.po[chr.i,6]);
          v4 <- as.numeric(chr.po[chr.i,7]);
          
          total.num <- as.numeric(chr.po[nrow(chr.po),5]);
          do.scale.cir(xc=xc, yc=yc, the.r=r+10, total.num=total.num, 
                       col="blue", lwd=0.001,
                       V1=v1, V2=v2, V3=v3, V4=v4);
        }
      }
    } 
  }
  ### end chr ###
  
  ### type = chr2
  if (type == "chr2"){
    if (class(mapping)[1] == "GRanges"){
      po       <- as.integer((start(mapping)+end(mapping))/2+0.5);
      dat.in   <- cbind(as.character(seqnames(mapping)),
                        po, names(mapping), mcols(mapping));   
    } else {
      dat.in   <- mapping;
    }
    dat.in[,1] <- gsub("chr", "", dat.in[,1]);
    
    dat.in[dat.in[,1]==23,1] <- "X";
    dat.in[dat.in[,1]==24,1] <- "Y";
    
    chr.lw <- W;
    ##
    if (!is.null(col)){
      cols  <- rep(col, nrow(dat.in))[c(1:nrow(dat.in))];    
      col.i <- 0;
    }
    ## 
    if (!is.null(cex)){
      cexes <- rep(cex, nrow(dat.in))[c(1:nrow(dat.in))];    
      cex.i <- 0;
    }
    dat.c     <- mapping;
    dat.c[,1] <- gsub("chr", "", dat.c[,1]);
    chr.num   <- nrow(chr.po);
    for (chr.i in c(1:chr.num)){
      chr.s  <- chr.po[chr.i,1];
      
      v1 <- as.numeric(chr.po[chr.i,2]);
      v2 <- as.numeric(chr.po[chr.i,3]);
      v3 <- as.numeric(chr.po[chr.i,6]);
      v4 <- as.numeric(chr.po[chr.i,7]);
      
      dat.v <- subset(dat.c, dat.c[,1]==chr.s);
      dat.v <- dat.v[order(as.numeric(dat.v[,2])),];
      for (i in 1:nrow(dat.v)){  
        if (dat.v[i,5]=="gneg"){
          col <- colors()[351];
        } else if (dat.v[i,5]=="acen" | dat.v[i,5]=="gvar" | dat.v[i,5]=="stalk"){
          col <- colors()[26];
        } else {
          col.v <- gsub("gpos","",dat.v[i,5]);
          if (col.v >= 30){
            col <- colors()[300];
          } else {
            col <- colors()[351];
          }
        }
        w1 <- scale.v(dat.v[i,2], v1, v2, v3, v4);
        w2 <- scale.v(dat.v[i,3], v1, v2, v3, v4);
        draw.arc.s(xc, yc, r, w1, w2, col=col, lwd=chr.lw);
      }
      if (print.chr.lab){
        v1  <- as.numeric(chr.po[chr.i,2]);
        v2  <- as.numeric(chr.po[chr.i,3]);
        w.m <- (v1+v2)/2;
        r.j <- r/20;
        chr.t <- gsub("chr", "", chr.s);
        draw.text.rt(xc, yc, r+r.j, w.m, chr.t, cex=0.6);
      } 
      if (scale){
        total.num <- as.numeric(chr.po[nrow(chr.po),5]);
        do.scale.cir(xc=xc, yc=yc, the.r=r+10, total.num=total.num, 
                     col="blue", lwd=0.001,
                     V1=v1, V2=v2, V3=v3, V4=v4);
      }
    }
  }
  ### end chr2 ####
  
  ### type  ################
  ### l, h, ls, lh ###
  if (type == "l" | type == "h" | type == "ls" | type == "lh" ){
    if (class(mapping)[1] == "GRanges"){
      po       <- as.integer((start(mapping)+end(mapping))/2+0.5);
      dat.in   <- cbind(as.character(seqnames(mapping)),
                        po, names(mapping), mcols(mapping));   
    } else {
      dat.in   <- mapping;
    }
    dat.in[,1] <- gsub("chr", "", dat.in[,1]);
    
    dat.in[dat.in[,1]==23,1] <- "X";
    dat.in[dat.in[,1]==24,1] <- "Y";
    
    ##
    if (!is.null(col)){
      cols  <- rep(col, nrow(dat.in))[c(1:nrow(dat.in))];    
      col.i <- 0;
    }
    ## 
    if (!is.null(cex)){
      cexes <- rep(cex, nrow(dat.in))[c(1:nrow(dat.in))];    
      cex.i <- 0;
    }
    
    dat.min <- min(as.numeric(dat.in[,col.v]), na.rm=T);
    dat.max <- max(as.numeric(dat.in[,col.v]), na.rm=T);
    
    # data set for the chromosome
    for (chr.i in 1:chr.num){
      chr.s <- chr.po[chr.i,1];
      chr.s <- gsub("chr","",chr.s);      
      dat   <- subset(dat.in, dat.in[,1]==chr.s);
      dat   <- dat[order(as.numeric(dat[,2])),];
      
      v1 <- as.numeric(chr.po[chr.i,2]);
      v2 <- as.numeric(chr.po[chr.i,3]);
      v3 <- as.numeric(chr.po[chr.i,6]);
      v4 <- as.numeric(chr.po[chr.i,7]);
      
      # background line
      if (B){
        draw.arc.pg(xc, yc, v1, v2, r, r+W-5, col=colors()[245]);
      } else {
        draw.arc.s(xc, yc, r, v1, v2, col=colors()[245], lwd=lwd);
      }
      
      # if the chromosome has not dat.
      if (dim(dat)[1] == 0){
        next;
      }   
      
      my.R1 <- R + W/5;
      my.R2 <- R + W - W/5;
      my.v    <- as.numeric(dat[1, col.v]);      
      v.old   <- scale.v(my.v, my.R1, my.R2, dat.min, dat.max);
      po      <- as.numeric(dat[1,2]);
      w.from  <- scale.v(po, v1, v2, v3, v4);
      
      for (i in 2:nrow(dat)){
        col.i <- col.i + 1;
        cex.i <- cex.i + 1;
        col   <- cols[col.i];
        cex   <- cexes[cex.i];
        my.v  <- as.numeric(dat[i, col.v]);
        if (is.na(my.v)){
          next;
        } 
        my.R1 <- R + W/5;
        my.R2 <- R + W - W/5;
        
        v     <- scale.v(my.v, my.R1, my.R2, dat.min, dat.max);
        po    <- as.numeric(dat[i,2]);
        w.to  <- scale.v(po, v1, v2, v3, v4);
        
        if (type=="ls"){
          if (v.old > 0){
            draw.line(xc, yc, w.from, v, v.old, col=col, lwd=lwd, lend=2);
          }
          draw.arc.s(xc, yc, v, w.from, w.to, col=col, lwd=lwd);
          v.old  <- v;
        }
        
        if (type=="l"){
          if (v.old > 0){            
            draw.line3(xc, yc, w.from, w.to, v.old, v, col=col, lwd=lwd)
          } 
          v.old <- v;
        }
        
        if (type=="lh"){
          if (v.old > 0){            
            draw.line3(xc, yc, w.from, w.to, v, v, col=col, lwd=lwd)
          } 
          v.old <- v;
        }
        
        if (type=="h"){
          tmp.v <- v - r - 1;
          draw.arc.pg(xc, yc, w.from, w.to, r, r+tmp.v, col=col, border=col);
        }
        
        w.from <- w.to;
      }
    }
    if (scale){
      do.scale(xc, yc, dat.min, dat.max,  R+W/5, W-2*(W/5));
    }
    
  }
  ### end l, h, ls, lh ###
  
  ### type = b, s
  if (type == "b" | type == "s"){
    if (class(mapping)[1] == "GRanges"){
      po       <- as.integer((start(mapping)+end(mapping))/2+0.5);
      dat.in   <- cbind(as.character(seqnames(mapping)),
                        po, names(mapping), mcols(mapping));   
    } else {
      dat.in   <- mapping;
    }
    dat.in[,1] <- gsub("chr", "", dat.in[,1]);
    
    dat.in[dat.in[,1]==23,1] <- "X";
    dat.in[dat.in[,1]==24,1] <- "Y";
    
    ##
    if (!is.null(col)){
      cols  <- rep(col, nrow(dat.in))[c(1:nrow(dat.in))];    
      col.i <- 0;
    }
    ## 
    if (!is.null(cex)){
      cexes <- rep(cex, nrow(dat.in))[c(1:nrow(dat.in))];    
      cex.i <- 0;
    }
    if (!is.null(lwd)){
      lwds  <- rep(lwd, nrow(dat.in))[c(1:nrow(dat.in))];    
      lwd.i <- 0;
    }
    dat.min <- min(as.numeric(dat.in[,col.v]), na.rm=T);
    dat.max <- max(as.numeric(dat.in[,col.v]), na.rm=T);
    
    # data set for the chromosome
    for (chr.i in 1:chr.num){
      chr.s <- chr.po[chr.i,1];
      chr.s <- gsub("chr","",chr.s);      
      dat   <- subset(dat.in, dat.in[,1]==chr.s);
      dat   <- dat[order(as.numeric(dat[,2])),];
      
      v1 <- as.numeric(chr.po[chr.i,2]);
      v2 <- as.numeric(chr.po[chr.i,3]);
      v3 <- as.numeric(chr.po[chr.i,6]);
      v4 <- as.numeric(chr.po[chr.i,7]);
      
      # background line
      if (B){
        draw.arc.pg(xc, yc, v1, v2, r, r+W-5, col=colors()[245]);
      } else {
        draw.arc.s(xc, yc, r, v1, v2, col=colors()[245], lwd=lwd);
      }
      
      # if the chromosome has not dat.
      if (dim(dat)[1] == 0){
        next;
      }   
      
      my.R1 <- R + W/5;
      my.R2 <- R + W - W/5;
      
      for (i in 1:nrow(dat)){
        col.i <- col.i + 1;
        cex.i <- cex.i + 1;
        lwd.i <- lwd.i + 1;
        col   <- cols[col.i];
        cex   <- cexes[cex.i];
        lwd   <- lwds[lwd.i];
        my.v  <- as.numeric(dat[i, col.v]);
        if (is.na(my.v)){
          next;
        } 
        my.R1 <- R + W/5;
        my.R2 <- R + W - W/5;
        
        v     <- scale.v(my.v, my.R1, my.R2, dat.min, dat.max);
        po    <- as.numeric(dat[i,2]);
        w.to  <- scale.v(po, v1, v2, v3, v4);
        if (type=="b"){
          draw.line(xc, yc, w.to, r, v, col=col, lwd=lwd);
        }
        if (type=="s"){
          draw.point.w(xc, yc, v, w.to, col=col, cex=cex);
        }
      }
    }
    if (scale){
      do.scale(xc, yc, dat.min, dat.max, R+W/5, W-2*(W/5));
    }
    
  }
  ### end b, s
  
  ### b2
  if (type == "b2"){
    if (class(mapping)[1] == "GRanges"){
      po       <- as.integer((start(mapping)+end(mapping))/2+0.5);
      dat.in   <- cbind(as.character(seqnames(mapping)),
                        po, names(mapping), mcols(mapping));   
    } else {
      dat.in   <- mapping;
    }
    dat.in[,1] <- gsub("chr", "", dat.in[,1]);
    
    dat.in[dat.in[,1]==23,1] <- "X";
    dat.in[dat.in[,1]==24,1] <- "Y";
    
    ## colors
    m.col <- FALSE;
    if (!is.null(col)){
      col.l <- length(col);
      if (col.l==1){
        cols <- c(col, col);
      } else if (col.l==2){
        cols <- col;
      } else {
        cols  <- rep(col, nrow(dat.in))[c(1:nrow(dat.in))];
        m.col <- TRUE;
      }
    } else {
      cols <- c("red", "blue");
    }
    
    ## cutoff
    if (is.numeric(cutoff)){
      cutoff <- cutoff;
    } else {
      cutoff <- 0;
    }
    
    dat.min <- min(as.numeric(dat.in[,col.v]), na.rm=T);
    dat.max <- max(as.numeric(dat.in[,col.v]), na.rm=T);
    
    # data set for the chromosome
    for (chr.i in 1:chr.num){
      chr.s <- chr.po[chr.i,1];
      chr.s <- gsub("chr","",chr.s);      
      dat   <- subset(dat.in, dat.in[,1]==chr.s);
      dat   <- dat[order(as.numeric(dat[,2])),];
      
      v1 <- as.numeric(chr.po[chr.i,2]);
      v2 <- as.numeric(chr.po[chr.i,3]);
      v3 <- as.numeric(chr.po[chr.i,6]);
      v4 <- as.numeric(chr.po[chr.i,7]);
      
      # background line
      if (B){
        draw.arc.pg(xc, yc, v1, v2, r, r+W-5, col=colors()[245]);
      } else {
        draw.arc.s(xc, yc, r, v1, v2, col=colors()[245], lwd=lwd);
      }
      
      # if the chromosome has not dat.
      if (dim(dat)[1] < 1){
        next;
      }   
      
      my.R1 <- R + W/5;
      my.R2 <- R + W - W/5;
      my.Rm <- (my.R1 + my.R2)/2;
      
      for (i in 1:nrow(dat)){
        my.v  <- as.numeric(dat[i, col.v]);
        if (m.col){
          col.i <- col.i + 1;
          col   <- cols[col.i];
        } else {
          if (my.v > cutoff){
            col <- cols[1];
          } else {
            col <- cols[2];
          }
        }
        
        if (is.na(my.v)){
          next;
        } 
        
        if (my.v > cutoff){
          v     <- scale.v(my.v, my.Rm, my.R2, cutoff, dat.max);
          po    <- as.numeric(dat[i,2]);
          w.to  <- scale.v(po, v1, v2, v3, v4);
          draw.line(xc, yc, w.to, my.Rm, v, col=col, lwd=lwd);
        } else {
          v     <- scale.v(my.v, my.R1, my.Rm, cutoff, dat.min);
          po    <- as.numeric(dat[i,2]);
          w.to  <- scale.v(po, v1, v2, v3, v4);
          draw.line(xc, yc, w.to, v, my.Rm, col=col, lwd=lwd);
        }
        
      }
    }
    if (scale){
      do.scale(xc, yc, dat.min, dat.max, R+W/5, W-2*(W/5));
    }
  }
  ### end b2
  
  ### b3
  if (type == "b3"){
    if (class(mapping)[1] == "GRanges"){
      po       <- as.integer((start(mapping)+end(mapping))/2+0.5);
      dat.in   <- cbind(as.character(seqnames(mapping)),
                        po, names(mapping), mcols(mapping));   
    } else {
      dat.in   <- mapping;
    }
    dat.in[,1] <- gsub("chr", "", dat.in[,1]);
    
    dat.in[dat.in[,1]==23,1] <- "X";
    dat.in[dat.in[,1]==24,1] <- "Y";
    
    ## colors
    col.i <- 0;
    m.col <- FALSE;
    if (!is.null(col)){
      cols  <- rep(col, nrow(dat.in))[c(1:nrow(dat.in))];
    } 
    
    # data set for the chromosome
    for (chr.i in 1:chr.num){
      chr.s <- chr.po[chr.i,1];
      chr.s <- gsub("chr","",chr.s);      
      dat   <- subset(dat.in, dat.in[,1]==chr.s);
      dat   <- dat[order(as.numeric(dat[,2])),];
      
      v1 <- as.numeric(chr.po[chr.i,2]);
      v2 <- as.numeric(chr.po[chr.i,3]);
      v3 <- as.numeric(chr.po[chr.i,6]);
      v4 <- as.numeric(chr.po[chr.i,7]);
      
      # background line
      if (B){
        draw.arc.pg(xc, yc, v1, v2, r, r+W-5, col=colors()[245]);
      } else {
        draw.arc.s(xc, yc, r, v1, v2, col=colors()[245], lwd=lwd);
      }
      
      # if the chromosome has not dat.
      if (dim(dat)[1] == 0){
        next;
      }   
      
      my.R1 <- R + W/10;
      my.R2 <- R + W - W/5;
      
      for (i in 1:nrow(dat)){
        col.i <- col.i + 1;
        col   <- cols[col.i];
        po    <- as.numeric(dat[i, 2]);
        w.to  <- scale.v(po, v1, v2, v3, v4);
        draw.line(xc, yc, w.to, my.R1, my.R2, col=col, lwd=lwd);
      }
    }
    if (scale){
      #do.scale(xc, yc, my.R1, my.R2, R+W/5, W-2*(W/5));
    }
    
  }
  ### end b3
  
  ### label2
  if (type == "label2"){
    if (class(mapping)[1] == "GRanges"){
      po       <- as.integer((start(mapping)+end(mapping))/2+0.5);
      dat.in   <- cbind(as.character(seqnames(mapping)),
                        po, names(mapping), mcols(mapping));   
    } else {
      dat.in   <- mapping;
    }
    dat.in[,1] <- gsub("chr", "", dat.in[,1]);
    
    dat.in[dat.in[,1]==23, 1] <- "X";
    dat.in[dat.in[,1]==24, 1] <- "Y";
    
    ## colors
    col.i <- 0;
    m.col <- FALSE;
    if (!is.null(col)){
      cols  <- rep(col, nrow(dat.in))[c(1:nrow(dat.in))];
    } 
    
    my.R1 <- R + W/10;
    my.R2 <- R + W - W/5;
    
    # data set for the chromosome
    for (i in 1:nrow(dat.in)){
      chr.s   <- dat.in[i, 1];
      chr.i   <- which(chr.po[,1]==chr.s);
      
      v1 <- as.numeric(chr.po[chr.i,2]);
      v2 <- as.numeric(chr.po[chr.i,3]);
      v3 <- as.numeric(chr.po[chr.i,6]);
      v4 <- as.numeric(chr.po[chr.i,7]);
      
      col.i <- col.i + 1;
      col   <- cols[col.i];
      po    <- as.numeric(dat.in[i, 2]);
      label.n <- as.character(dat.in[i, col.v]);
      w.to  <- scale.v(po, v1, v2, v3, v4);
      
      if (side=="in"){
        draw.text.rt(xc, yc, my.R1, w.to, n=label.n, col=col, cex=cex, side="in");
      } else {
        draw.text.rt(xc, yc, my.R1, w.to, n=label.n, col=col, cex=cex, side="out");
      }
    }
  }
  ### end label2
  
  ### s2
  if (type == "s2"){
    if (class(mapping)[1] == "GRanges"){
      po       <- as.integer((start(mapping)+end(mapping))/2+0.5);
      dat.in   <- cbind(as.character(seqnames(mapping)),
                        po, names(mapping), mcols(mapping));   
    } else {
      dat.in   <- mapping;
    }
    dat.in[,1] <- gsub("chr", "", dat.in[,1]);
    
    dat.in[dat.in[,1]==23,1] <- "X";
    dat.in[dat.in[,1]==24,1] <- "Y";
    
    ## colors
    col.i <- 0;
    m.col <- FALSE;
    if (!is.null(col)){
      cols  <- rep(col, nrow(dat.in))[c(1:nrow(dat.in))];
    } 
    
    # data set for the chromosome
    for (chr.i in 1:chr.num){
      chr.s <- chr.po[chr.i,1];
      chr.s <- gsub("chr","",chr.s);      
      dat   <- subset(dat.in, dat.in[,1]==chr.s);
      dat   <- dat[order(as.numeric(dat[,2])),];
      
      v1 <- as.numeric(chr.po[chr.i,2]);
      v2 <- as.numeric(chr.po[chr.i,3]);
      v3 <- as.numeric(chr.po[chr.i,6]);
      v4 <- as.numeric(chr.po[chr.i,7]);
      
      # background line
      if (B){
        draw.arc.pg(xc, yc, v1, v2, r, r+W-5, col=colors()[245]);
      } else {
        draw.arc.s(xc, yc, r, v1, v2, col=colors()[245], lwd=lwd);
      }
      
      # if the chromosome has not dat.
      if (dim(dat)[1] == 0){
        next;
      }   
      
      my.R1 <- R + W/10;
      my.R2 <- R + W - W/5;
      my.Rm <- (my.R1 + my.R2)/2;
      
      for (i in 1:nrow(dat)){
        col.i <- col.i + 1;
        col   <- cols[col.i];
        po    <- as.numeric(dat[i, 2]);
        w.to  <- scale.v(po, v1, v2, v3, v4);
        draw.point.w(xc, yc, my.Rm, w.to, col=col, cex=cex);
      }
    }
    if (scale){
      do.scale(xc, yc, dat.min, dat.max, R+W/5, W-2*(W/5));
    }
    
  }
  ### end s2
  
  ### label
  if (type == "label"){
    if (class(mapping)[1] == "GRanges"){
      po       <- as.integer((start(mapping)+end(mapping))/2+0.5);
      dat.in   <- cbind(as.character(seqnames(mapping)),
                        po, names(mapping), mcols(mapping));   
    } else {
      dat.in   <- mapping;
    }
    r <- R;
    w <- W;
    
    if (!is.null(col)){
      cols <- rep(col, nrow(dat.in))[c(1:nrow(dat.in))];    
      col.i <- 0;
    }
    if (!is.null(col)){
      cexes <- rep(cex, nrow(dat.in))[c(1:nrow(dat.in))];    
      cex.i <- 0;
    }
    
    dat.in[,1] <- gsub("chr", "", dat.in[,1]);
    
    # data set for the chromosome
    for (chr.i in 1:chr.num){
      chr.s <- chr.po[chr.i,1];
      chr.s <- gsub("chr", "", chr.s);
      dat   <- subset(dat.in, dat.in[,1]==chr.s);
      
      # if the chromosome has not dat.
      if (dim(dat)[1] == 0){
        next;
      }  
      
      v1 <- as.numeric(chr.po[chr.i,2]);
      v2 <- as.numeric(chr.po[chr.i,3]);
      v3 <- as.numeric(chr.po[chr.i,6]);
      v4 <- as.numeric(chr.po[chr.i,7]);      
      
      # sort by position
      dat <- dat[order(dat[,2]),];
      label.num <- nrow(dat);
      
      v.angle <- v2-v1;
      # if the label number > maxlab number then label the fist maxlab
      maxlab  <- as.integer((v2-v1)*3);
      if (label.num > maxlab){
        label.num == maxlab;
      }
      
      v.wide  <- v.angle/label.num;
      gap.len <- 2;
      if (v.wide > gap.len){
        v.wide <- gap.len;
        w.po  <- v1+v.angle/2-(label.num/2)*gap.len;
      } else {
        w.po  <- v1;
      }
      
      for (i in 1:label.num){
        col.i <- col.i + 1;
        col   <- cols[col.i];
        cex.i <- cex.i + 1;
        cex   <- cexes[col.i];
        label.n <- dat[i,3];
        po      <- as.numeric(dat[i,2]);
        w.to    <- scale.v(po, v1, v2, v3, v4);
        if (side == "in"){             
          draw.line(xc,  yc, w.to,   r,   r-w/3,          col=col, lwd=lwd);
          draw.line(xc,  yc, w.po,   r-w+w/3, r-w,        col=col, lwd=lwd);
          draw.line3(xc, yc, w.to, w.po, r-w/3, r-w+w/3,  col=col, lwd=lwd);                      
          draw.text.rt(xc, yc, w.po, r=r-w-40, n=label.n, col=col, cex=cex, side="in");
        } else if (side == "out"){
          draw.line(xc,  yc, w.to,   r,  r+w/3,           col=col, lwd=lwd);
          draw.line(xc,  yc, w.po,   r+w-w/3, r+w,        col=col, lwd=lwd);
          draw.line3(xc, yc, w.to, w.po, r+w/3, r+w-w/3,  col=col, lwd=lwd);                      
          draw.text.rt(xc, yc, w.po, r=r+w+10, n=label.n, col=col, cex=cex, side="out");
        }
        w.po <- w.po + v.wide;
      }   
    }
  }
  # end label
  
  ### link
  if (type == "link"){
    chr.po[,4] <- gsub("chr", "", chr.po[,4]);
    dat.in <- mapping;
    dat.in[,1] <- gsub("chr", "", dat.in[,1]);
    dat.in[,4] <- gsub("chr", "", dat.in[,4]);
    
    #dat.in[dat.in[,1]==23,1] <- "X";
    #dat.in[dat.in[,1]==24,1] <- "Y";
    #dat.in[dat.in[,4]==23,1] <- "X";
    #dat.in[dat.in[,4]==24,1] <- "Y";
    
    dat    <- dat.in;
    
    if (!is.null(col)){
      col <- rep(col, nrow(dat.in))[c(1:nrow(dat.in))];    
    }
    if (!is.null(lwd)){
      lwd <- rep(lwd, nrow(dat.in))[c(1:nrow(dat.in))];    
    } else {
      lwd <- rep(1, nrow(dat.in));
    }
    
    r <- R;
    for (i in 1:nrow(dat)){
      chr1.s   <- dat[i,1];
      chr2.s   <- dat[i,4];
      po1      <- dat[i,2];
      po2      <- dat[i,5];
      
      chr1     <- which(chr.po[,1]==chr1.s);
      chr2     <- which(chr.po[,1]==chr2.s);
      
      v1 <- as.numeric(chr.po[chr1,2]);
      v2 <- as.numeric(chr.po[chr1,3]);
      v3 <- as.numeric(chr.po[chr1,6]);
      v4 <- as.numeric(chr.po[chr1,7]);
      
      w1 <- scale.v(as.numeric(po1), v1, v2, v3, v4);
      
      v1 <- as.numeric(chr.po[chr2,2]);
      v2 <- as.numeric(chr.po[chr2,3]);
      v3 <- as.numeric(chr.po[chr2,6]);
      v4 <- as.numeric(chr.po[chr2,7]);
      
      w2 <- scale.v(as.numeric(po2), v1, v2, v3, v4);
      
      if (chr1 == chr2){
        draw.link(xc, yc, r, w1, w2, col=col[i], lwd=lwd);
      } else {
        draw.link(xc, yc, r, w1, w2, col=col[i], lwd=lwd);
      }
    }
  }
  ### end link
  
  ### link2
  if (type == "link2"){
    chr.po[,4] <- gsub("chr", "", chr.po[,4]);
    dat.in <- mapping;
    dat.in[,1] <- gsub("chr", "", dat.in[,1]);
    dat.in[,4] <- gsub("chr", "", dat.in[,4]);
    
    #dat.in[dat.in[,1]==23,1] <- "X";
    #dat.in[dat.in[,1]==24,1] <- "Y";
    #dat.in[dat.in[,4]==23,1] <- "X";
    #dat.in[dat.in[,4]==24,1] <- "Y";
    
    dat    <- dat.in;
    
    if (!is.null(col)){
      col <- rep(col, nrow(dat.in))[c(1:nrow(dat.in))];    
    }
    if (!is.null(lwd)){
      lwd <- rep(lwd, nrow(dat.in))[c(1:nrow(dat.in))];    
    } else {
      lwd <- rep(1, nrow(dat.in));
    }
    
    r <- R;
    for (i in 1:nrow(dat)){
      chr1.s   <- dat[i,1];
      chr2.s   <- dat[i,4];
      po1      <- dat[i,2];
      po2      <- dat[i,5];
      
      chr1     <- which(chr.po[,1]==chr1.s);
      chr2     <- which(chr.po[,1]==chr2.s);
      
      v1 <- as.numeric(chr.po[chr1,2]);
      v2 <- as.numeric(chr.po[chr1,3]);
      v3 <- as.numeric(chr.po[chr1,6]);
      v4 <- as.numeric(chr.po[chr1,7]);
      
      w1 <- scale.v(as.numeric(po1), v1, v2, v3, v4);
      
      v1 <- as.numeric(chr.po[chr2,2]);
      v2 <- as.numeric(chr.po[chr2,3]);
      v3 <- as.numeric(chr.po[chr2,6]);
      v4 <- as.numeric(chr.po[chr2,7]);
      
      w2 <- scale.v(as.numeric(po2), v1, v2, v3, v4);
      
      if (chr1 == chr2){
        draw.link2(xc, yc, r, w1, w2, col=col[i], lwd=lwd);
      } else {
        draw.link(xc, yc, r, w1, w2, col=col[i], lwd=lwd);
      }
    }
  }
  ### end link2
  
  ### links of sub-segments
  
  if (type=="link.pg"){
    thefile  <- system.file("data", "UCSC.chr.colors.RData", package="OmicCircos");
    if (file.exists(thefile)){
    } else {
      stop ("UCSC.chr.colors ?");
    }
    dat.col     <- local(get(load(thefile)))
    dat.col[,1] <- gsub("chr", "", dat.col[,1]);
    if (class(mapping)[1] == "GRanges"){
      po       <- as.integer((start(mapping)+end(mapping))/2+0.5);
      dat.in   <- cbind(as.character(seqnames(mapping)),
                        po, names(mapping), mcols(mapping));   
    } else {
      dat.in   <- mapping;
    }
    
    dat.in[,1] <- gsub("chr", "", dat.in[,1]);
    dat.in[,4] <- gsub("chr", "", dat.in[,4]);
    
    #dat.in[dat.in[,1]==23,1] <- "X";
    #dat.in[dat.in[,1]==24,1] <- "Y";
    
    if (!is.null(col)){
      col <- rep(col, nrow(dat.in))[c(1:nrow(dat.in))];    
    }
    
    if (!is.null(lwd)){
      lwd <- rep(lwd, nrow(dat.in))[c(1:nrow(dat.in))];    
    } else {
      lwd <- rep(1, nrow(dat.in));
    }
    
    for (i in 1:nrow(dat.in)){
      the.d1    <- as.numeric(dat.in[i,3])-as.numeric(dat.in[i,2]);
      the.d2    <- as.numeric(dat.in[i,6])-as.numeric(dat.in[i,5]);
      
      #if (the.d < 10000000){
      #next;
      #}
      
      chr1.s    <- dat.in[i,1];
      chr2.s    <- dat.in[i,4];
      
      #cir.i      <- sum(grep("mm", cir));
      #if (cir.i == 1){
      # chr2.s    <- paste("mm.", dat.in[i,4], sep="");
      #} 
      #cir.i      <- sum(grep("mac", cir));
      #if (cir.i == 1){
      # chr2.s    <- paste("mac.", dat.in[i,4], sep="");
      #}
      #the.color <- as.character(dat.col[dat.col[,1]==chr1.s,2]);
      #the.color <- paste(the.color, "50", sep="");
      if (is.null(col)){
        if (dat.col[i,1]==chr1.s){
          the.color <- rainbow(10, alpha=0.8)[7];
        } else {
          the.color <- rainbow(10, alpha=0.8)[1];
        }
      } else {
        the.color <- col[i];
      }
      po1.1    <- as.numeric(dat.in[i,2]);
      po1.2    <- as.numeric(dat.in[i,3]);
      
      po2.1    <- as.numeric(dat.in[i,5]);
      po2.2    <- as.numeric(dat.in[i,6]);
      
      chr1     <- which(chr.po[,1]==chr1.s);
      chr2     <- which(chr.po[,1]==chr2.s);
      
      v1 <- as.numeric(chr.po[chr1,2]);
      v2 <- as.numeric(chr.po[chr1,3]);
      v3 <- as.numeric(chr.po[chr1,6]);
      v4 <- as.numeric(chr.po[chr1,7]);
      
      w1.1 <- scale.v(po1.1, v1, v2, v3, v4);
      w1.2 <- scale.v(po1.2, v1, v2, v3, v4);
      
      v1 <- as.numeric(chr.po[chr2,2]);
      v2 <- as.numeric(chr.po[chr2,3]);
      v3 <- as.numeric(chr.po[chr2,6]);
      v4 <- as.numeric(chr.po[chr2,7]);
      
      w2.1 <- scale.v(po2.1, v1, v2, v3, v4);
      w2.2 <- scale.v(po2.2, v1, v2, v3, v4);
      
      if (po1.1 > po1.2 & po2.2 < po2.1){
        draw.link.pg(xc, yc, r, w1.2, w1.1, w2.2, w2.1, col=the.color, lwd=lwd[i]);
      } else if (po1.1 < po1.2 & po2.2 > po2.1){
        draw.link.pg(xc, yc, r, w1.1, w1.2, w2.1, w2.2, col=the.color, lwd=lwd[i]);
      } else {
        draw.link.pg(xc, yc, r, w1.1, w1.2, w2.1, w2.2, col=the.color, lwd=lwd[i]);
      }
    }
  }
  ### end links of sub-segments
  
  if (type == "highlight.link"){
    xc   <- mapping[1];
    yc   <- mapping[2];
    r    <- mapping[3];
    w1.1 <- mapping[4];
    w1.2 <- mapping[5];
    w2.1 <- mapping[6];
    w2.2 <- mapping[7];
    draw.link.pg(xc, yc, r, w1.1, w1.2, w2.1, w2.2, col=col, lwd=lwd);
  }
  
  ### heatmap
  if (type == "heatmap"){
    if (class(mapping)[1] == "GRanges"){
      po       <- as.integer((start(mapping)+end(mapping))/2+0.5);
      dat.in   <- cbind(as.character(seqnames(mapping)),
                        po, names(mapping), mcols(mapping));   
    } else {
      dat.in   <- mapping;
    }
    dat.in[,1] <- gsub("chr", "", dat.in[,1]);
    
    dat.in[dat.in[,1]==23,1] <- "X";
    dat.in[dat.in[,1]==24,1] <- "Y";
    
    dat.i      <- c(col.v:ncol(dat.in));
    dat.m      <- dat.in[,dat.i];
    dat.m      <- as.matrix(dat.m);
    
    if (cluster){
      dat.d      <- dist(t(dat.m));
      dat.h      <- hclust(dat.d);
      
      dat.m      <- dat.m[,dat.h$order];
    }
    
    ## for the matrix colors
    num.col <- ncol(dat.m);
    num.row <- nrow(dat.m);
    
    dat.c <- as.numeric(unlist(dat.m));
    
    dat.min <- min(dat.c, na.rm=T);
    dat.max <- max(dat.c, na.rm=T);
    
    dat.n <- matrix(dat.c, ncol=num.col, nrow=num.row);
    colnames(dat.n) <- colnames(dat.n);
    
    rbPal <- colorRampPalette(c("blue","white","red"))
    col.i <- rbPal(20)[cut(dat.c, breaks = 20)];
    col.m <- matrix(col.i, ncol=num.col, nrow=num.row);
    
    w          <- W-W/25;
    
    w.i        <- w/length(dat.i);
    for (i in 1:nrow(dat.in)){
      dat.i <- as.numeric(dat.m[i,]);
      chr.i <- which(chr.po[,1]==dat.in[i,1]);
      
      v1 <- as.numeric(chr.po[chr.i,2]);
      v2 <- as.numeric(chr.po[chr.i,3]);
      v3 <- as.numeric(chr.po[chr.i,6]);
      v4 <- as.numeric(chr.po[chr.i,7]);
      w.to <- scale.v(as.numeric(dat.in[i,2]), v1, v2, v3, v4);
      
      r.from <- r + w - w/25;            
      for (j in 1:length(dat.i)){
        w.to <- scale.v(as.numeric(dat.in[i,2]), v1, v2, v3, v4);
        r.to <- r.from - w.i;
        draw.line(xc,  yc, w.to, r.from, r.to, col=col.m[i,j], lwd=lwd);
        r.from <- r.to;
      } 
    }
    
    if (col.bar){
      if (col.bar.po == "topleft"){
        color.bar(40,700,50,780, v.min=dat.min, v.max=dat.max);
      }
      if (col.bar.po == "topright"){
        color.bar(770,700,780,780, v.min=dat.min, v.max=dat.max);
      }
      if (col.bar.po == "bottomleft"){
        color.bar(40, 10, 50, 90, v.min=dat.min, v.max=dat.max);
      }
      if (col.bar.po == "bottomright"){
        color.bar(770, 10, 780, 90, v.min=dat.min, v.max=dat.max);
      }
    }
    if (cluster){
      heatmap.cluster(0, 20, 150, 70, dat.m);
    }
  }
  ### end heatmap
  
  ### heatmap2 ## 
  if (type == "heatmap2"){
    if (class(mapping)[1] == "GRanges"){
      po       <- as.integer((start(mapping)+end(mapping))/2+0.5);
      dat.in   <- cbind(as.character(seqnames(mapping)),
                        po, names(mapping), mcols(mapping));   
    } else {
      dat.in   <- mapping;
    }
    dat.in[,1] <- gsub("chr", "", dat.in[,1]);
    
    dat.in[dat.in[,1]==23,1] <- "X";
    dat.in[dat.in[,1]==24,1] <- "Y";
    
    dat.i      <- c(col.v:ncol(dat.in));
    dat.m      <- dat.in[,dat.i];
    dat.m      <- as.matrix(dat.m);
    
    if (cluster){
      dat.d      <- dist(t(dat.m));
      dat.h      <- hclust(dat.d);
      
      dat.m      <- dat.m[,dat.h$order];
    }
    
    ## for the matrix colors
    num.col <- ncol(dat.m);
    num.row <- nrow(dat.m);
    
    dat.c <- as.numeric(unlist(dat.m));
    
    dat.min <- min(dat.c, na.rm=T);
    dat.max <- max(dat.c, na.rm=T);
    
    dat.n <- matrix(dat.c, ncol=num.col, nrow=num.row);
    colnames(dat.n) <- colnames(dat.n);
    
    rbPal  <- colorRampPalette(c("blue","white","red"))
    col.i  <- rbPal(20)[cut(dat.c, breaks = 20)];
    col.m  <- matrix(col.i, ncol=num.col, nrow=num.row);
    col.df <- cbind(dat.in[,c(1:2)], col.m);
    
    w          <- W-W/25;
    w          <- w-w/5;
    w5         <- w/5;
    
    lwd1       <- 360/(num.row+24)*1.5;
    
    w.i        <- w/length(dat.i);
    
    for (chr.i in 1:chr.num){
      chr.s <- chr.po[chr.i,1];
      chr.s <- gsub("chr","",chr.s);
      
      dat   <- subset(dat.in, dat.in[,1]==chr.s);
      dat   <- dat[order(as.numeric(dat[,2])),];
      col.c <- subset(col.df, col.df[,1]==chr.s);
      col.c <- col.c[order(as.numeric(col.c[,2])),];   
      col.c <- col.c[,-c(1:2)];   
      
      if (dim(dat)[1]==0){
        next;
      }
      v1 <- as.numeric(chr.po[chr.i,2]);
      v2 <- as.numeric(chr.po[chr.i,3]);
      v3 <- as.numeric(chr.po[chr.i,6]);
      v4 <- as.numeric(chr.po[chr.i,7]);
      
      ## for heatmap position
      v.angle <- v2-v1;
      v.wide  <- v.angle/nrow(dat);
      w.po    <- v1;
      
      for (i in 1:nrow(dat)){   
        r.from <- r + w - w/25; 
        
        ## chr position
        w.to   <- scale.v(as.numeric(dat[i,2]), v1, v2, v3, v4);
        
        for (j in 1:length(dat.i)){
          r.to <- r.from - w.i;
          if (v.wide > 0.5){
            draw.arc.pg(xc, yc, w.po, w.po+v.wide, r.from, r.to, col=col.c[i,j], border=col.c[i,j]);
          } else {
            draw.line(xc,  yc, w.po, r.from, r.to, col=col.c[i,j], lwd=lwd1);
          }
          r.from <- r.to;
        }
        
        ## link heatmap position and chromosome position
        draw.line(xc,  yc, w.to,   r+w+w5*2/3, r+w+w5,       col=col, lwd=lwd);
        draw.line(xc,  yc, w.po,   r+w,        r+w+w5/3,     col=col, lwd=lwd);
        draw.line3(xc, yc, w.to, w.po, r+w+w5*2/3, r+w+w5/3, col=col, lwd=lwd);                      
        
        w.po <- w.po + v.wide;
      } 
    }
    
    if (col.bar){
      if (col.bar.po == "topleft"){
        color.bar(40,700,50,780, v.min=dat.min, v.max=dat.max);
      }
      if (col.bar.po == "topright"){
        color.bar(770,700,780,780, v.min=dat.min, v.max=dat.max);
      }
      if (col.bar.po == "bottomleft"){
        color.bar(40, 10, 50, 90, v.min=dat.min, v.max=dat.max);
      }
      if (col.bar.po == "bottomright"){
        color.bar(770, 10, 780, 90, v.min=dat.min, v.max=dat.max);
      }
    }
    if (cluster){
      heatmap.cluster(0, 20, 150, 70, dat.m);
    }
  }
  ### end heatmap2
  
  ### boxplot
  if (type == "box"){
    
    if (class(mapping)[1] == "GRanges"){
      po       <- as.integer((start(mapping)+end(mapping))/2+0.5);
      dat.in   <- cbind(as.character(seqnames(mapping)),
                        po, names(mapping), mcols(mapping));   
    } else {
      dat.in   <- mapping;
    }
    dat.in  <- na.omit(dat.in);
    
    dat.in[,1] <- gsub("chr", "", dat.in[,1]);
    
    dat.in[dat.in[,1]==23,1] <- "X";
    dat.in[dat.in[,1]==24,1] <- "Y";
    
    # data set for the chromosome
    for (chr.i in 1:chr.num){
      chr.s <- chr.po[chr.i,1];
      chr.s <- gsub("chr","",chr.s);      
      dat   <- subset(dat.in, dat.in[,1]==chr.s);
      dat   <- dat[order(as.numeric(dat[,2])),];
      v1 <- as.numeric(chr.po[chr.i,2]);
      v2 <- as.numeric(chr.po[chr.i,3]);
      v3 <- as.numeric(chr.po[chr.i,6]);
      v4 <- as.numeric(chr.po[chr.i,7]);
      
      # background line
      if (B){
        draw.arc.pg(xc, yc, v1, v2, r, r+W-5, col=colors()[245]);
      } else {
        draw.arc.s(xc, yc, r, v1, v2, col=colors()[245], lwd=lwd);
      }
    }
    ## 
    
    dat.i      <- c(col.v:ncol(dat.in));
    dat.m      <- dat.in[,dat.i];
    
    dat.m      <- data.matrix(dat.m);
    
    dat.min    <- min(as.numeric(dat.m), na.rm=T);
    dat.max    <- max(as.numeric(dat.m), na.rm=T);
    
    mat.quant  <- apply(dat.m, 1, quantile, 
                        probs=c(0.05, 0.25, 0.5, 0.75, 0.95));
    
    #lwd <- 360/(nrow(dat.in)+length(unique(dat.in[,1])));
    lwd <- 360/nrow(dat.in)
    if (lwd < 3){
      lwd <- 3;
    }
    
    ## for color
    color1 <- rainbow(10, alpha=0.5)[7];
    color2 <- rainbow(10)[7];
    
    for (i in 1:nrow(dat.in)){
      
      chr.i <- which(dat.in[i,1]==chr.po[,1]);
      
      v1 <- as.numeric(chr.po[chr.i,2]);
      v2 <- as.numeric(chr.po[chr.i,3]);
      v3 <- as.numeric(chr.po[chr.i,6]);
      v4 <- as.numeric(chr.po[chr.i,7]);
      w.to <- scale.v(as.numeric(dat.in[i,2]), v1, v2, v3, v4);
      
      my.R1 <- r;
      my.R2 <- r + W;
      
      mat.num <- length(mat.quant[,i]);
      if (mat.num < 5){
        next;
      }
      
      q1     <- scale.v(mat.quant[1,i], my.R1, my.R2, dat.min, dat.max);
      q2     <- scale.v(mat.quant[2,i], my.R1, my.R2, dat.min, dat.max);
      qm     <- scale.v(mat.quant[3,i], my.R1, my.R2, dat.min, dat.max);
      q3     <- scale.v(mat.quant[4,i], my.R1, my.R2, dat.min, dat.max);
      q4     <- scale.v(mat.quant[5,i], my.R1, my.R2, dat.min, dat.max);
      
      draw.line(xc,  yc, w.to, q1, q4, col=color1, lwd=0.2);
      draw.line(xc,  yc, w.to, q2, q3, col=color1, lwd=lwd);
      draw.line(xc,  yc, w.to, qm-0.1, qm+0.1, col=color2, lwd=lwd);
      
    }
    
    if (scale){
      do.scale(xc, yc, dat.min, dat.max, R+W/5, W-2*(W/5));
    }
    
  }
  ### end boxplot
  
  ### quantile plot
  if (type == "quant90"){
    if (class(mapping)[1] == "GRanges"){
      po       <- as.integer((start(mapping)+end(mapping))/2+0.5);
      dat.in   <- cbind(as.character(seqnames(mapping)),
                        po, names(mapping), mcols(mapping));   
    } else {
      dat.in   <- mapping;
    }
    dat.in  <- na.omit(dat.in);
    
    dat.in[,1] <- gsub("chr", "", dat.in[,1]);
    
    dat.in[dat.in[,1]==23,1] <- "X";
    dat.in[dat.in[,1]==24,1] <- "Y";
    
    # data set for the chromosome
    for (chr.i in 1:chr.num){
      chr.s <- chr.po[chr.i,1];
      chr.s <- gsub("chr","",chr.s);      
      dat   <- subset(dat.in, dat.in[,1]==chr.s);
      dat   <- dat[order(as.numeric(dat[,2])),];
      v1 <- as.numeric(chr.po[chr.i,2]);
      v2 <- as.numeric(chr.po[chr.i,3]);
      v3 <- as.numeric(chr.po[chr.i,6]);
      v4 <- as.numeric(chr.po[chr.i,7]);
      
      # background line
      if (B){
        draw.arc.pg(xc, yc, v1, v2, r, r+W-5, col=colors()[245]);
      } else {
        draw.arc.s(xc, yc, r, v1, v2, col=colors()[245], lwd=lwd);
      }
    }
    ## 
    
    dat.i      <- c(col.v:ncol(dat.in));
    dat.m      <- dat.in[,dat.i];
    
    dat.m      <- data.matrix(dat.m);
    
    dat.min    <- min(as.numeric(dat.m), na.rm=T);
    dat.max    <- max(as.numeric(dat.m), na.rm=T);
    
    mat.quant  <- apply(dat.m, 1, quantile, 
                        probs=c(0.05, 0.25, 0.5, 0.75, 0.95));
    
    dat.qu     <- cbind(dat.in[,c(1,2)], t(mat.quant));
    
    ## for color
    colors <- rainbow(10, alpha=0.5);
    
    for (chr.i in 1:chr.num){
      chr.s <- chr.po[chr.i,1];
      chr.s <- gsub("chr","",chr.s);
      
      dat   <- subset(dat.qu, dat.qu[,1]==chr.s);
      dat   <- dat[order(as.numeric(dat[,2])),];
      
      v1 <- as.numeric(chr.po[chr.i,2]);
      v2 <- as.numeric(chr.po[chr.i,3]);
      v3 <- as.numeric(chr.po[chr.i,6]);
      v4 <- as.numeric(chr.po[chr.i,7]);
      
      my.R1 <- r;
      my.R2 <- r + W;
      
      i <- 1; 
      q1.from <- scale.v(as.numeric(dat[i,3]), my.R1, my.R2, dat.min, dat.max);
      q2.from <- scale.v(as.numeric(dat[i,4]), my.R1, my.R2, dat.min, dat.max);
      qm.from <- scale.v(as.numeric(dat[i,5]), my.R1, my.R2, dat.min, dat.max);
      q3.from <- scale.v(as.numeric(dat[i,6]), my.R1, my.R2, dat.min, dat.max);
      q4.from <- scale.v(as.numeric(dat[i,7]), my.R1, my.R2, dat.min, dat.max);
      w.from  <- scale.v(as.numeric(dat[i,2]), v1, v2, v3, v4);
      
      for (i in 2:nrow(dat)){
        
        q1.to     <- scale.v(as.numeric(dat[i,3]), my.R1, my.R2, dat.min, dat.max);
        q2.to     <- scale.v(as.numeric(dat[i,4]), my.R1, my.R2, dat.min, dat.max);
        qm.to     <- scale.v(as.numeric(dat[i,5]), my.R1, my.R2, dat.min, dat.max);
        q3.to     <- scale.v(as.numeric(dat[i,6]), my.R1, my.R2, dat.min, dat.max);
        q4.to     <- scale.v(as.numeric(dat[i,7]), my.R1, my.R2, dat.min, dat.max);
        w.to      <- scale.v(as.numeric(dat[i,2]), v1, v2, v3, v4);
        
        draw.line3(xc,  yc, w.from, w.to, q1.from, q1.to, col=colors[1], lwd=lwd);
        draw.line3(xc,  yc, w.from, w.to, qm.from, qm.to, col=colors[4], lwd=lwd);
        draw.line3(xc,  yc, w.from, w.to, q4.from, q4.to, col=colors[7], lwd=lwd);
        
        w.from  <- w.to;
        q1.from <- q1.to;
        q2.from <- q2.to;
        qm.from <- qm.to;
        q3.from <- q3.to;
        q4.from <- q4.to;
        
      }
    }
    
    if (scale){
      do.scale(xc, yc, dat.min, dat.max, R+W/5, W-2*(W/5));
    }
    
  }
  ### end quantile plot
  
  ### quantile plot
  if (type == "quant75"){
    if (class(mapping)[1] == "GRanges"){
      po       <- as.integer((start(mapping)+end(mapping))/2+0.5);
      dat.in   <- cbind(as.character(seqnames(mapping)),
                        po, names(mapping), mcols(mapping));   
    } else {
      dat.in   <- mapping;
    }
    dat.in  <- na.omit(dat.in);
    
    dat.in[,1] <- gsub("chr", "", dat.in[,1]);
    
    dat.in[dat.in[,1]==23,1] <- "X";
    dat.in[dat.in[,1]==24,1] <- "Y";
    
    # data set for the chromosome
    for (chr.i in 1:chr.num){
      chr.s <- chr.po[chr.i,1];
      chr.s <- gsub("chr","",chr.s);      
      dat   <- subset(dat.in, dat.in[,1]==chr.s);
      dat   <- dat[order(as.numeric(dat[,2])),];
      v1 <- as.numeric(chr.po[chr.i,2]);
      v2 <- as.numeric(chr.po[chr.i,3]);
      v3 <- as.numeric(chr.po[chr.i,6]);
      v4 <- as.numeric(chr.po[chr.i,7]);
      
      # background line
      if (B){
        draw.arc.pg(xc, yc, v1, v2, r, r+W-5, col=colors()[245]);
      } else {
        draw.arc.s(xc, yc, r, v1, v2, col=colors()[245], lwd=lwd);
      }
    }
    ## 
    
    dat.i      <- c(col.v:ncol(dat.in));
    dat.m      <- dat.in[,dat.i];
    
    dat.m      <- data.matrix(dat.m);
    
    dat.min    <- min(as.numeric(dat.m), na.rm=T);
    dat.max    <- max(as.numeric(dat.m), na.rm=T);
    
    mat.quant  <- apply(dat.m, 1, quantile, 
                        probs=c(0.05, 0.25, 0.5, 0.75, 0.95));
    
    dat.qu     <- cbind(dat.in[,c(1,2)], t(mat.quant));
    
    ## for color
    colors <- rainbow(10, alpha=0.5);
    
    for (chr.i in 1:chr.num){
      chr.s <- chr.po[chr.i,1];
      chr.s <- gsub("chr","",chr.s);
      
      dat   <- subset(dat.qu, dat.qu[,1]==chr.s);
      dat   <- dat[order(as.numeric(dat[,2])),];
      
      v1 <- as.numeric(chr.po[chr.i,2]);
      v2 <- as.numeric(chr.po[chr.i,3]);
      v3 <- as.numeric(chr.po[chr.i,6]);
      v4 <- as.numeric(chr.po[chr.i,7]);
      
      my.R1 <- r;
      my.R2 <- r + W;
      
      i <- 1; 
      q1.from <- scale.v(as.numeric(dat[i,3]), my.R1, my.R2, dat.min, dat.max);
      q2.from <- scale.v(as.numeric(dat[i,4]), my.R1, my.R2, dat.min, dat.max);
      qm.from <- scale.v(as.numeric(dat[i,5]), my.R1, my.R2, dat.min, dat.max);
      q3.from <- scale.v(as.numeric(dat[i,6]), my.R1, my.R2, dat.min, dat.max);
      q4.from <- scale.v(as.numeric(dat[i,7]), my.R1, my.R2, dat.min, dat.max);
      w.from  <- scale.v(as.numeric(dat[i,2]), v1, v2, v3, v4);
      
      for (i in 2:nrow(dat)){
        
        q1.to     <- scale.v(as.numeric(dat[i,3]), my.R1, my.R2, dat.min, dat.max);
        q2.to     <- scale.v(as.numeric(dat[i,4]), my.R1, my.R2, dat.min, dat.max);
        qm.to     <- scale.v(as.numeric(dat[i,5]), my.R1, my.R2, dat.min, dat.max);
        q3.to     <- scale.v(as.numeric(dat[i,6]), my.R1, my.R2, dat.min, dat.max);
        q4.to     <- scale.v(as.numeric(dat[i,7]), my.R1, my.R2, dat.min, dat.max);
        w.to      <- scale.v(as.numeric(dat[i,2]), v1, v2, v3, v4);
        
        draw.line3(xc,  yc, w.from, w.to, q2.from, q2.to, col=colors[1], lwd=lwd);
        draw.line3(xc,  yc, w.from, w.to, qm.from, qm.to, col=colors[4], lwd=lwd);
        draw.line3(xc,  yc, w.from, w.to, q3.from, q3.to, col=colors[7], lwd=lwd);
        
        w.from  <- w.to;
        q1.from <- q1.to;
        q2.from <- q2.to;
        qm.from <- qm.to;
        q3.from <- q3.to;
        q4.from <- q4.to;
        
      }
    }
    
    if (scale){
      do.scale(xc, yc, dat.min, dat.max, R+W/5, W-2*(W/5));
    }
    
  }
  ### end quantile 75% plot
  
  ### sample 95% confidence intervals plot
  if (type == "ci95"){
    if (class(mapping)[1] == "GRanges"){
      po       <- as.integer((start(mapping)+end(mapping))/2+0.5);
      dat.in   <- cbind(as.character(seqnames(mapping)),
                        po, names(mapping), mcols(mapping));   
    } else {
      dat.in   <- mapping;
    }
    dat.in  <- na.omit(dat.in);
    
    dat.in[,1] <- gsub("chr", "", dat.in[,1]);
    
    dat.in[dat.in[,1]==23,1] <- "X";
    dat.in[dat.in[,1]==24,1] <- "Y";
    
    # data set for the chromosome
    for (chr.i in 1:chr.num){
      chr.s <- chr.po[chr.i,1];
      chr.s <- gsub("chr","",chr.s);      
      dat   <- subset(dat.in, dat.in[,1]==chr.s);
      dat   <- dat[order(as.numeric(dat[,2])),];
      v1 <- as.numeric(chr.po[chr.i,2]);
      v2 <- as.numeric(chr.po[chr.i,3]);
      v3 <- as.numeric(chr.po[chr.i,6]);
      v4 <- as.numeric(chr.po[chr.i,7]);
      
      # background line
      if (B){
        draw.arc.pg(xc, yc, v1, v2, r, r+W-5, col=colors()[245]);
      } else {
        draw.arc.s(xc, yc, r, v1, v2, col=colors()[245], lwd=lwd);
      }
    }
    ## 
    
    dat.i      <- c(col.v:ncol(dat.in));
    dat.m      <- dat.in[,dat.i];
    
    dat.m      <- data.matrix(dat.m);
    
    dat.min    <- min(as.numeric(dat.m), na.rm=T);
    dat.max    <- max(as.numeric(dat.m), na.rm=T);
    
    get.conf.int <- function(x) t.test(x)$conf.int
    mat.ci     <- apply(dat.m, 1, get.conf.int);
    mat.me     <- apply(dat.m, 1, mean, rm.na=T);
    
    dat.ci     <- cbind(dat.in[,c(1,2)], t(mat.ci), mat.me);
    
    ## for color
    colors <- rainbow(10, alpha=0.5);
    
    for (chr.i in 1:chr.num){
      chr.s <- chr.po[chr.i,1];
      chr.s <- gsub("chr","",chr.s);
      
      dat   <- subset(dat.ci, dat.ci[,1]==chr.s);
      dat   <- dat[order(as.numeric(dat[,2])),];
      
      v1 <- as.numeric(chr.po[chr.i,2]);
      v2 <- as.numeric(chr.po[chr.i,3]);
      v3 <- as.numeric(chr.po[chr.i,6]);
      v4 <- as.numeric(chr.po[chr.i,7]);
      
      my.R1 <- r;
      my.R2 <- r + W;
      
      i <- 1; 
      ci1.from <- scale.v(as.numeric(dat[i,3]), my.R1, my.R2, dat.min, dat.max);
      ci2.from <- scale.v(as.numeric(dat[i,4]), my.R1, my.R2, dat.min, dat.max);
      m.from   <- scale.v(as.numeric(dat[i,5]), my.R1, my.R2, dat.min, dat.max);
      w.from   <- scale.v(as.numeric(dat[i,2]), v1, v2, v3, v4);
      
      for (i in 2:nrow(dat)){
        
        ci1.to     <- scale.v(as.numeric(dat[i,3]), my.R1, my.R2, dat.min, dat.max);
        ci2.to     <- scale.v(as.numeric(dat[i,4]), my.R1, my.R2, dat.min, dat.max);
        m.to       <- scale.v(as.numeric(dat[i,5]), my.R1, my.R2, dat.min, dat.max);
        w.to       <- scale.v(as.numeric(dat[i,2]), v1, v2, v3, v4);
        
        draw.line3(xc,  yc, w.from, w.to, ci1.from, ci1.to, col=colors[1], lwd=lwd);
        draw.line3(xc,  yc, w.from, w.to, m.from,   m.to, col=colors[4], lwd=lwd);
        draw.line3(xc,  yc, w.from, w.to, ci2.from, ci2.to, col=colors[7], lwd=lwd);
        
        w.from   <- w.to;
        ci1.from <- ci1.to;
        ci2.from <- ci2.to;
        m.from   <- m.to;
      }
    }
    
    if (scale){
      do.scale(xc, yc, dat.min, dat.max, R+W/5, W-2*(W/5));
    }
    
  }
  ### end quantile 75% plot
  
  ### sv
  if (type == "sv"){
    if (class(mapping)[1] == "GRanges"){
      po       <- as.integer((start(mapping)+end(mapping))/2+0.5);
      dat.in   <- cbind(as.character(seqnames(mapping)),
                        po, names(mapping), mcols(mapping));   
    } else {
      dat.in   <- mapping;
    }
    dat.in[,1] <- gsub("chr", "", dat.in[,1]);
    
    #
    dat.in[dat.in[,1]==23,1] <- "X";
    dat.in[dat.in[,1]==24,1] <- "Y";
    
    # data set for the chromosome
    for (chr.i in 1:chr.num){
      chr.s <- chr.po[chr.i,1];
      chr.s <- gsub("chr","",chr.s);      
      dat   <- subset(dat.in, dat.in[,1]==chr.s);
      dat   <- dat[order(as.numeric(dat[,2])),];
      v1 <- as.numeric(chr.po[chr.i,2]);
      v2 <- as.numeric(chr.po[chr.i,3]);
      v3 <- as.numeric(chr.po[chr.i,6]);
      v4 <- as.numeric(chr.po[chr.i,7]);
      
      # background line
      if (B){
        draw.arc.pg(xc, yc, v1, v2, r, r+W-5, col=colors()[245]);
      } else {
        draw.arc.s(xc, yc, r, v1, v2, col=colors()[245], lwd=lwd);
      }
    }
    
    my.R1 <- R + W/5;
    my.R2 <- R + W - W/5;
    ## 
    dat.i   <- c(col.v:ncol(dat.in));
    dat.m   <- dat.in[,dat.i];
    dat.m   <- as.matrix(dat.m);
    dat.min <- min(as.numeric(dat.m), na.rm=T);
    dat.max <- max(as.numeric(dat.m), na.rm=T);
    dat.me  <- mean(as.numeric(dat.m), na.rm=T);
    
    var1    <- function(x) var(x,na.rm=T)
    all.var <- apply(dat.m, 1, var1)
    var.min <- min(all.var, na.rm=T);
    var.max <- max(all.var, na.rm=T);
    num.col <- ncol(dat.in);
    
    the.cex <- 360/nrow(dat.in);
    if (the.cex < 3){
      the.cex <- 3;
    }
    
    for (chr.i in 1:chr.num){
      chr.s <- chr.po[chr.i,1];
      chr.s <- gsub("chr","",chr.s);
      
      dat   <- subset(dat.in, dat.in[,1]==chr.s);
      dat   <- dat[order(as.numeric(dat[,2])),];
      
      v1 <- as.numeric(chr.po[chr.i,2]);
      v2 <- as.numeric(chr.po[chr.i,3]);
      v3 <- as.numeric(chr.po[chr.i,6]);
      v4 <- as.numeric(chr.po[chr.i,7]);
      
      
      for (i in 1:nrow(dat)){           
        my.v   <- as.numeric(dat[i,dat.i]); 
        d.v    <- var(my.v, na.rm = T);
        d.m    <- mean(my.v, na.rm = T);
        
        if (d.m > dat.me){
          col <- rainbow(10, alpha=0.5)[1];
        } else {
          col <- rainbow(10, alpha=0.5)[7];
        }
        
        po      <- as.numeric(dat[i,2]);
        w.to    <- scale.v(po, v1, v2, v3, v4);
        
        the.v  <- scale.v(d.m, my.R1, my.R2, dat.min, dat.max);
        the.c  <- scale.v(d.v, 0.1, the.cex, var.min, var.max);
        
        draw.point.w(xc, yc, the.v, w.to, col=col, cex=the.c);
        
      }   
    }    
    
    if (scale){
      do.scale(xc, yc, dat.min, dat.max, R+W/5, W-2*(W/5));
    }
    
  }
  ### end sv
  
  ### sd
  if (type == "s.sd"){
    if (class(mapping)[1] == "GRanges"){
      po       <- as.integer((start(mapping)+end(mapping))/2+0.5);
      dat.in   <- cbind(as.character(seqnames(mapping)),
                        po, names(mapping), mcols(mapping));   
    } else {
      dat.in   <- mapping;
    }
    dat.in[,1] <- gsub("chr", "", dat.in[,1]);
    
    #
    dat.in[dat.in[,1]==23,1] <- "X";
    dat.in[dat.in[,1]==24,1] <- "Y";
    
    # data set for the chromosome
    for (chr.i in 1:chr.num){
      chr.s <- chr.po[chr.i,1];
      chr.s <- gsub("chr","",chr.s);      
      dat   <- subset(dat.in, dat.in[,1]==chr.s);
      dat   <- dat[order(as.numeric(dat[,2])),];
      v1 <- as.numeric(chr.po[chr.i,2]);
      v2 <- as.numeric(chr.po[chr.i,3]);
      v3 <- as.numeric(chr.po[chr.i,6]);
      v4 <- as.numeric(chr.po[chr.i,7]);
      
      # background line
      if (B){
        draw.arc.pg(xc, yc, v1, v2, r, r+W-5, col=colors()[245]);
      } else {
        draw.arc.s(xc, yc, r, v1, v2, col=colors()[245], lwd=lwd);
      }
    }
    
    my.R1 <- R + W/5;
    my.R2 <- R + W - W/5;
    ## 
    dat.i   <- c(col.v:ncol(dat.in));
    dat.m   <- dat.in[,dat.i];
    dat.m   <- as.matrix(dat.m);
    dat.min <- min(as.numeric(dat.m), na.rm=T);
    dat.max <- max(as.numeric(dat.m), na.rm=T);
    dat.me  <- mean(as.numeric(dat.m), na.rm=T);
    
    sd1    <- function(x) sd(x,na.rm=T)
    all.sd <- apply(dat.m, 1, sd1)
    sd.min <- min(all.sd,na.rm=T);
    sd.max <- max(all.sd,na.rm=T);
    
    num.col <- ncol(dat.in);
    
    the.cex <- 360/nrow(dat.in);
    if (the.cex < 3){
      the.cex <- 3;
    }
    
    for (chr.i in 1:chr.num){
      chr.s <- chr.po[chr.i,1];
      chr.s <- gsub("chr","",chr.s);
      
      dat   <- subset(dat.in, dat.in[,1]==chr.s);
      dat   <- dat[order(as.numeric(dat[,2])),];
      
      v1 <- as.numeric(chr.po[chr.i,2]);
      v2 <- as.numeric(chr.po[chr.i,3]);
      v3 <- as.numeric(chr.po[chr.i,6]);
      v4 <- as.numeric(chr.po[chr.i,7]);
      
      
      for (i in 1:nrow(dat)){           
        my.v   <- as.numeric(dat[i,dat.i]); 
        sd.v    <- sd(my.v, na.rm = T);
        d.m     <- mean(my.v, na.rm = T);
        
        if (d.m > dat.me){
          col <- rainbow(10, alpha=0.5)[1];
        } else {
          col <- rainbow(10, alpha=0.5)[7];
        }
        
        po      <- as.numeric(dat[i,2]);
        w.to    <- scale.v(po, v1, v2, v3, v4);
        
        the.v  <- scale.v(d.m, my.R1, my.R2, dat.min, dat.max);
        the.c  <- scale.v(sd.v, 0.1, the.cex, sd.min, sd.max);
        
        draw.point.w(xc, yc, the.v, w.to, col=col, cex=the.c);
        
      }   
    }    
    
    if (scale){
      do.scale(xc, yc, dat.min, dat.max, R+W/5, W-2*(W/5));
    }
    
  }
  ### end s sd
  
  ### ss
  if (type == "ss"){
    if (class(mapping)[1] == "GRanges"){
      po       <- as.integer((start(mapping)+end(mapping))/2+0.5);
      dat.in   <- cbind(as.character(seqnames(mapping)),
                        po, names(mapping), mcols(mapping));   
    } else {
      dat.in   <- mapping;
    }
    dat.in[,1] <- gsub("chr", "", dat.in[,1]);
    
    #
    dat.in[dat.in[,1]==23,1] <- "X";
    dat.in[dat.in[,1]==24,1] <- "Y";
    
    # data set for the chromosome
    for (chr.i in 1:chr.num){
      chr.s <- chr.po[chr.i,1];
      chr.s <- gsub("chr","",chr.s);      
      dat   <- subset(dat.in, dat.in[,1]==chr.s);
      dat   <- dat[order(as.numeric(dat[,2])),];
      v1 <- as.numeric(chr.po[chr.i,2]);
      v2 <- as.numeric(chr.po[chr.i,3]);
      v3 <- as.numeric(chr.po[chr.i,6]);
      v4 <- as.numeric(chr.po[chr.i,7]);
      
      # background line
      if (B){
        draw.arc.pg(xc, yc, v1, v2, r, r+W-5, col=colors()[245]);
      } else {
        draw.arc.s(xc, yc, r, v1, v2, col=colors()[245], lwd=lwd);
      }
    }
    
    my.R1 <- R + W/5;
    my.R2 <- R + W - W/5;
    ## 
    dat.i   <- c(col.v:ncol(dat.in));
    dat.m   <- dat.in[,dat.i];
    dat.m   <- as.matrix(dat.m);
    dat.min <- min(as.numeric(dat.m), na.rm=T);
    dat.max <- max(as.numeric(dat.m), na.rm=T);
    dat.me  <- mean(as.numeric(dat.m), na.rm=T);
    
    num.col <- ncol(dat.in);
    
    the.cex <- 360/nrow(dat.in);
    if (the.cex < 3){
      the.cex <- 3;
    }
    
    for (chr.i in 1:chr.num){
      chr.s <- chr.po[chr.i,1];
      chr.s <- gsub("chr","",chr.s);
      
      dat   <- subset(dat.in, dat.in[,1]==chr.s);
      dat   <- dat[order(as.numeric(dat[,2])),];
      
      v1 <- as.numeric(chr.po[chr.i,2]);
      v2 <- as.numeric(chr.po[chr.i,3]);
      v3 <- as.numeric(chr.po[chr.i,6]);
      v4 <- as.numeric(chr.po[chr.i,7]);
      
      
      for (i in 1:nrow(dat)){           
        my.v   <- as.numeric(dat[i,dat.i]); 
        d.m    <- mean(my.v, na.rm = T);
        
        if (d.m > dat.me){
          col <- rainbow(10, alpha=0.5)[1];
        } else {
          col <- rainbow(10, alpha=0.5)[7];
        }
        
        po      <- as.numeric(dat[i,2]);
        w.to    <- scale.v(po, v1, v2, v3, v4);
        
        the.v  <- scale.v(d.m, my.R1, my.R2, dat.min, dat.max);
        the.c  <- scale.v(d.m, 0.1, the.cex, dat.min, dat.max);
        
        draw.point.w(xc, yc, the.v, w.to, col=col, cex=the.c);
        
      }   
    }    
    
    if (scale){
      do.scale(xc, yc, dat.min, dat.max, R+W/5, W-2*(W/5));
    }
    
  }
  ### end sv
  
  ### ml
  if (type == "ml"){
    if (class(mapping)[1] == "GRanges"){
      po       <- as.integer((start(mapping)+end(mapping))/2+0.5);
      dat.in   <- cbind(as.character(seqnames(mapping)),
                        po, names(mapping), mcols(mapping));   
    } else {
      dat.in   <- mapping;
    }
    dat.in[,1] <- gsub("chr", "", dat.in[,1]);
    
    #
    dat.in[dat.in[,1]==23,1] <- "X";
    dat.in[dat.in[,1]==24,1] <- "Y";
    
    # data set for the chromosome
    for (chr.i in 1:chr.num){
      chr.s <- chr.po[chr.i,1];
      chr.s <- gsub("chr","",chr.s);      
      dat   <- subset(dat.in, dat.in[,1]==chr.s);
      dat   <- dat[order(as.numeric(dat[,2])),];
      v1 <- as.numeric(chr.po[chr.i,2]);
      v2 <- as.numeric(chr.po[chr.i,3]);
      v3 <- as.numeric(chr.po[chr.i,6]);
      v4 <- as.numeric(chr.po[chr.i,7]);
      
      # background line
      if (B){
        draw.arc.pg(xc, yc, v1, v2, r, r+W-5, col=colors()[245]);
      } else {
        draw.arc.s(xc, yc, r, v1, v2, col=colors()[245], lwd=lwd);
      }
    }
    ## 
    dat.i   <- c(col.v:ncol(dat.in));
    dat.m   <- dat.in[,dat.i];
    dat.m   <- as.matrix(dat.m);
    dat.min <- min(as.numeric(dat.m), na.rm=T);
    dat.max <- max(as.numeric(dat.m), na.rm=T);
    
    my.R1 <- R + W/5;
    my.R2 <- R + W - W/5;
    
    ## for the matrix colors
    num.col <- ncol(dat.m);
    num.row <- nrow(dat.m);
    
    if (length(col) == num.col){
      colors <- col;
    } else {
      colors <- rainbow(num.col, alpha=0.5);
    }
    
    for (chr.i in 1:chr.num){
      chr.s <- chr.po[chr.i,1];
      chr.s <- gsub("chr","",chr.s);
      
      dat   <- subset(dat.in, dat.in[,1]==chr.s);
      dat   <- dat[order(as.numeric(dat[,2])),];
      
      v1 <- as.numeric(chr.po[chr.i,2]);
      v2 <- as.numeric(chr.po[chr.i,3]);
      v3 <- as.numeric(chr.po[chr.i,6]);
      v4 <- as.numeric(chr.po[chr.i,7]);
      
      # colors 
      col.i <- 0;
      for (j in col.v:ncol(dat)){
        col.i <- col.i + 1;
        col   <- colors[col.i];
        
        my.v      <- as.numeric(dat[1,j]); 
        dat.i.old <- my.v;
        v.old   <- scale.v(my.v, my.R1, my.R2, dat.min, dat.max);
        
        po      <- as.numeric(dat[1,2]);
        w.from  <- scale.v(po, v1, v2, v3, v4);
        
        for (k in 1:nrow(dat)){
          
          dat.i <- as.numeric(dat[k,j]);
          
          if (is.na(dat.i)){
            next;
          }
          
          v    <- scale.v(dat.i, my.R1, my.R2, dat.min, dat.max);
          w.to <- scale.v(as.numeric(dat[k,2]), v1, v2, v3, v4);
          
          if (w.from > 0){            
            draw.line3(xc, yc, w.from, w.to, v.old, v, col=col, lwd=lwd)
          } 
          
          dat.i.old <- dat.i;
          w.from    <- w.to;
          v.old     <- v;
        } # end the row
      }   # end the col
    }     # end the chr/segment
    
    if (scale){
      do.scale(xc, yc, dat.min, dat.max, R+W/5, W-2*(W/5));
    }
    
  }
  ### end ml
  
  ### arc
  if (type == "arc"){
    if (class(mapping)[1] == "GRanges"){
      po       <- as.integer((start(mapping)+end(mapping))/2+0.5);
      dat.in   <- cbind(as.character(seqnames(mapping)),
                        po, names(mapping), mcols(mapping));   
    } else {
      dat.in   <- mapping;
    }
    dat.in[,1] <- gsub("chr", "", dat.in[,1]);
    
    if (!is.null(col)){
      colors <- rep(col, nrow(dat.in))[c(1:nrow(dat.in))];    
      col.i  <- 0; 
    }
    #
    dat.in[dat.in[,1]==23,1] <- "X";
    dat.in[dat.in[,1]==24,1] <- "Y";
    
    # data set for the chromosome
    for (chr.i in 1:chr.num){
      chr.s <- chr.po[chr.i,1];
      chr.s <- gsub("chr","",chr.s);      
      dat   <- subset(dat.in, dat.in[,1]==chr.s);
      dat   <- dat[order(as.numeric(dat[,2])),];
      v1 <- as.numeric(chr.po[chr.i,2]);
      v2 <- as.numeric(chr.po[chr.i,3]);
      v3 <- as.numeric(chr.po[chr.i,6]);
      v4 <- as.numeric(chr.po[chr.i,7]);
      
      # background line
      if (B){
        draw.arc.pg(xc, yc, v1, v2, r, r+W-5, col=colors()[245]);
      } else {
        draw.arc.s(xc, yc, r, v1, v2, col=colors()[245], lwd=lwd);
      }
    }
    ## 
    my.R1 <- R + W/10;
    my.R2 <- R + W - W/6;
    my.Rm <- (my.R1 + my.R2)/2;
    
    dat.i   <- c(col.v:ncol(dat.in));
    dat.m   <- dat.in[,dat.i];
    dat.m   <- as.matrix(dat.m);
    dat.min <- min(as.numeric(dat.m), na.rm=T);
    dat.max <- max(as.numeric(dat.m), na.rm=T);
    
    for (i in 1:nrow(dat.in)){
      col.i <- col.i + 1;
      col   <- colors[col.i];
      chr.i <- which(chr.po[,1]==dat.in[i,1]);
      v1 <- as.numeric(chr.po[chr.i,2]);
      v2 <- as.numeric(chr.po[chr.i,3]);
      v3 <- as.numeric(chr.po[chr.i,6]);
      v4 <- as.numeric(chr.po[chr.i,7]);
      
      the.r   <- scale.v(as.numeric(dat.in[i,col.v]), my.R1, my.R2, dat.min, dat.max);
      
      w1 <- scale.v(as.numeric(dat.in[i,2]), v1, v2, v3, v4);
      w2 <- scale.v(as.numeric(dat.in[i,3]), v1, v2, v3, v4);
      
      draw.arc.s(xc, yc, the.r, w1, w2, col=col, lwd=lwd)
    }
    
    if (scale){
      do.scale(xc, yc, dat.min, dat.max, R+W/5, W-2*(W/5));
    }
  }
  ### end arc
  
  ### arc2
  if (type == "arc2"){
    if (class(mapping)[1] == "GRanges"){
      po       <- as.integer((start(mapping)+end(mapping))/2+0.5);
      dat.in   <- cbind(as.character(seqnames(mapping)),
                        po, names(mapping), mcols(mapping));   
    } else {
      dat.in   <- mapping;
    }
    dat.in[,1] <- gsub("chr", "", dat.in[,1]);
    
    if (!is.null(col)){
      colors <- rep(col, nrow(dat.in))[c(1:nrow(dat.in))];    
      col.i  <- 0; 
    }
    #
    dat.in[dat.in[,1]==23,1] <- "X";
    dat.in[dat.in[,1]==24,1] <- "Y";
    
    # data set for the chromosome
    for (chr.i in 1:chr.num){
      chr.s <- chr.po[chr.i,1];
      chr.s <- gsub("chr","",chr.s);      
      dat   <- subset(dat.in, dat.in[,1]==chr.s);
      dat   <- dat[order(as.numeric(dat[,2])),];
      v1 <- as.numeric(chr.po[chr.i,2]);
      v2 <- as.numeric(chr.po[chr.i,3]);
      v3 <- as.numeric(chr.po[chr.i,6]);
      v4 <- as.numeric(chr.po[chr.i,7]);
      
      # background line
      if (B){
        draw.arc.pg(xc, yc, v1, v2, r, r+W-5, col=colors()[245]);
      } else {
        draw.arc.s(xc, yc, r, v1, v2, col=colors()[245], lwd=1);
      }
    }
    ## 
    my.R1 <- R + W/10;
    my.R2 <- R + W - W/6;
    my.Rm <- (my.R1 + my.R2)/2;
    
    for (i in 1:nrow(dat.in)){
      col.i <- col.i + 1;
      col   <- colors[col.i];
      chr.i <- which(chr.po[,1]==dat.in[i,1]);
      v1 <- as.numeric(chr.po[chr.i,2]);
      v2 <- as.numeric(chr.po[chr.i,3]);
      v3 <- as.numeric(chr.po[chr.i,6]);
      v4 <- as.numeric(chr.po[chr.i,7]);
      
      w1 <- scale.v(as.numeric(dat.in[i,2]), v1, v2, v3, v4);
      w2 <- scale.v(as.numeric(dat.in[i,3]), v1, v2, v3, v4);
      
      draw.arc.s(xc, yc, my.Rm, w1, w2, col=col, lwd=lwd)
    }
    
    if (scale){
      do.scale(xc, yc, dat.min, dat.max, R+W/5, W-2*(W/5));
    }
  }
  ### end arc2
  
  ### ml2
  if (type == "ml2"){
    if (class(mapping)[1] == "GRanges"){
      po       <- as.integer((start(mapping)+end(mapping))/2+0.5);
      dat.in   <- cbind(as.character(seqnames(mapping)),
                        po, names(mapping), mcols(mapping));   
    } else {
      dat.in   <- mapping;
    }
    dat.in[,1] <- gsub("chr", "", dat.in[,1]);
    
    dat.in[dat.in[,1]==23,1] <- "X";
    dat.in[dat.in[,1]==24,1] <- "Y";
    
    # data set for the chromosome
    for (chr.i in 1:chr.num){
      chr.s <- chr.po[chr.i,1];
      chr.s <- gsub("chr","",chr.s);      
      dat   <- subset(dat.in, dat.in[,1]==chr.s);
      dat   <- dat[order(as.numeric(dat[,2])),];
      v1 <- as.numeric(chr.po[chr.i,2]);
      v2 <- as.numeric(chr.po[chr.i,3]);
      v3 <- as.numeric(chr.po[chr.i,6]);
      v4 <- as.numeric(chr.po[chr.i,7]);
      
      # background line
      if (B){
        draw.arc.pg(xc, yc, v1, v2, r, r+W-5, col=colors()[245]);
      } else {
        draw.arc.s(xc, yc, r, v1, v2, col=colors()[245], lwd=lwd);
      }
    }
    ## 
    dat.in[,1] <- gsub("chr", "", dat.in[,1]);
    dat.i      <- c(col.v:ncol(dat.in));
    dat.m      <- dat.in[,dat.i];
    dat.m      <- as.matrix(dat.m);
    dat.min    <- min(as.numeric(dat.m), na.rm=T);
    dat.max    <- max(as.numeric(dat.m), na.rm=T);
    
    ## for the matrix colors
    num.col <- ncol(dat.m);
    num.row <- nrow(dat.m);
    
    if (cutoff == "n"){
      if (length(col) == num.col){
        colors <- col;
      } else {
        colors <- rainbow(num.col, alpha=0.5);
      }
    } else {
      if (length(col) == 2){
        colors <- col;
      } else {
        colors <- rainbow(10, alpha=0.5)[c(1,7)];
      }
    }
    
    my.R1 <- R + W/5;
    my.R2 <- R + W - W/5;
    
    for (chr.i in 1:chr.num){
      chr.s <- chr.po[chr.i,1];
      chr.s <- gsub("chr","",chr.s);
      
      dat   <- subset(dat.in, dat.in[,1]==chr.s);
      dat   <- dat[order(as.numeric(dat[,2])),];
      
      v1 <- as.numeric(chr.po[chr.i,2]);
      v2 <- as.numeric(chr.po[chr.i,3]);
      v3 <- as.numeric(chr.po[chr.i,6]);
      v4 <- as.numeric(chr.po[chr.i,7]);
      
      col.i <- 0;
      for (j in col.v:ncol(dat)){           
        my.v      <- as.numeric(dat[1,j]);
        if (cutoff == "n"){
          col.i <- col.i + 1;
          col   <- colors[col.i];
        } else {
          if (my.v > cutoff){
            col <- colors[1];
          } else {
            col <- colors[2];
          }
        }
        
        dat.i.old <- my.v;
        v.old   <- scale.v(my.v, my.R1, my.R2, dat.min, dat.max);
        po      <- as.numeric(dat[1,2]);
        w.from  <- scale.v(po, v1, v2, v3, v4);
        
        for (k in 1:nrow(dat)){
          
          dat.i <- as.numeric(dat[k,j]);
          
          if (is.na(dat.i)){
            next;
          }
          
          v    <- scale.v(dat.i, my.R1, my.R2, dat.min, dat.max);
          w.to <- scale.v(as.numeric(dat[k,2]), v1, v2, v3, v4);
          
          draw.arc.s(xc, yc, v, w.from, w.to, col=col, lwd=lwd);
          
          dat.i.old <- dat.i;
          w.from    <- w.to;
          v.old     <- v;
        } # end the row
      }   # end the col
    }     # end the chr/segment
    
    
    if (scale){
      do.scale(xc, yc, dat.min, dat.max, R+W/5, W-2*(W/5));
    }
    
  }
  ### end ml2
  
  ### ml3
  if (type == "ml3"){
    if (class(mapping)[1] == "GRanges"){
      po       <- as.integer((start(mapping)+end(mapping))/2+0.5);
      dat.in   <- cbind(as.character(seqnames(mapping)),
                        po, names(mapping), mcols(mapping));   
    } else {
      dat.in   <- mapping;
    }
    dat.in[,1] <- gsub("chr", "", dat.in[,1]);
    
    dat.in[dat.in[,1]==23,1] <- "X";
    dat.in[dat.in[,1]==24,1] <- "Y";
    
    # data set for the chromosome
    for (chr.i in 1:chr.num){
      
      v1 <- as.numeric(chr.po[chr.i,2]);
      v2 <- as.numeric(chr.po[chr.i,3]);
      v3 <- as.numeric(chr.po[chr.i,6]);
      v4 <- as.numeric(chr.po[chr.i,7]);
      
      # background line
      if (B){
        draw.arc.pg(xc, yc, v1, v2, r, r+W-5, col=colors()[245]);
      } else {
        draw.arc.s(xc, yc, r, v1, v2, col=colors()[245], lwd=lwd);
      }
    }
    ## 
    dat.in[,1] <- gsub("chr", "", dat.in[,1]);
    dat.i      <- c(col.v:ncol(dat.in));
    dat.m      <- dat.in[,dat.i];
    dat.m      <- as.matrix(dat.m);
    dat.min    <- min(as.numeric(dat.m), na.rm=T);
    dat.max    <- max(as.numeric(dat.m), na.rm=T);
    
    ## for the matrix colors
    num.col <- ncol(dat.m);
    num.row <- nrow(dat.m);
    
    if (length(col) == 2){
      colors <- col;
    } else {
      colors <- rainbow(10, alpha=0.5)[c(1,7)];
    }
    
    my.R1 <- R + W/5;
    my.R2 <- R + W - W/5;
    
    for (chr.i in 1:chr.num){
      chr.s <- chr.po[chr.i,1];
      chr.s <- gsub("chr","",chr.s);
      
      dat   <- subset(dat.in, dat.in[,1]==chr.s);
      dat   <- dat[order(as.numeric(dat[,2])),];
      
      v1 <- as.numeric(chr.po[chr.i,2]);
      v2 <- as.numeric(chr.po[chr.i,3]);
      v3 <- as.numeric(chr.po[chr.i,6]);
      v4 <- as.numeric(chr.po[chr.i,7]);
      
      for (j in col.v:ncol(dat)){     
        
        my.v      <- as.numeric(dat[1,j]); 
        dat.i.old <- my.v;
        if (is.na(my.v)==T){
          next;
        }
        v.old   <- scale.v(my.v, my.R1, my.R2, dat.min, dat.max);
        po      <- as.numeric(dat[1,2]);
        w.from  <- scale.v(po, v1, v2, v3, v4);
        
        for (k in 2:nrow(dat)){
          
          dat.i <- as.numeric(dat[k,j]);
          
          if (is.na(dat.i)){
            next;
          }
          
          if (dat.i > cutoff){
            col <- colors[1];
          } else {
            col <- colors[2];
          }
          
          v    <- scale.v(dat.i, my.R1, my.R2, dat.min, dat.max);
          w.to <- scale.v(as.numeric(dat[k,2]), v1, v2, v3, v4);
          
          if (v.old > 0){
            tmp.v <- dat.i * dat.i.old;
            if (tmp.v < 0){
              the0 <- scale.v(0, my.R1, my.R2, dat.min, dat.max);
              if (dat.i > cutoff){
                col <- colors[1];
                draw.line(xc, yc, w.from, v, the0, col=col, lwd=lwd, lend=2);
                col <- colors[2];
                draw.line(xc, yc, w.from, the0, v.old, col=col, lwd=lwd, lend=2);
              } else {
                col <- colors[2];
                draw.line(xc, yc, w.from, v, the0, col=col, lwd=lwd, lend=2);
                col <- colors[1];
                draw.line(xc, yc, w.from, the0, v.old, col=col, lwd=lwd, lend=2);
              }
            } else {
              draw.line(xc, yc, w.from, v, v.old, col=col, lwd=lwd, lend=2);
            }
          }
          if (dat.i > cutoff){
            col <- colors[1];
          } else {
            col <- colors[2];
          }
          draw.arc.s(xc, yc, v, w.from, w.to, col=col, lwd=lwd);
          
          dat.i.old <- dat.i;
          w.from    <- w.to;
          v.old     <- v;
        } # end the row
      }   # end the col
    }     # end the chr/segment
    
    if (scale){
      do.scale(xc, yc, dat.min, dat.max, R+W/5, W-2*(W/5));
    }
    
  }       # end the function
  ### end ml3
  
  ### ms
  if (type == "ms"){
    if (class(mapping)[1] == "GRanges"){
      po       <- as.integer((start(mapping)+end(mapping))/2+0.5);
      dat.in   <- cbind(as.character(seqnames(mapping)),
                        po, names(mapping), mcols(mapping));   
    } else {
      dat.in   <- mapping;
    }
    dat.in[,1] <- gsub("chr", "", dat.in[,1]);
    
    dat.in[dat.in[,1]==23,1] <- "X";
    dat.in[dat.in[,1]==24,1] <- "Y";
    
    # data set for the chromosome
    for (chr.i in 1:chr.num){
      chr.s <- chr.po[chr.i,1];
      chr.s <- gsub("chr","",chr.s);      
      dat   <- subset(dat.in, dat.in[,1]==chr.s);
      dat   <- dat[order(as.numeric(dat[,2])),];
      v1 <- as.numeric(chr.po[chr.i,2]);
      v2 <- as.numeric(chr.po[chr.i,3]);
      v3 <- as.numeric(chr.po[chr.i,6]);
      v4 <- as.numeric(chr.po[chr.i,7]);
      
      # background line
      if (B){
        draw.arc.pg(xc, yc, v1, v2, r, r+W-5, col=colors()[245]);
      } else {
        draw.arc.s(xc, yc, r, v1, v2, col=colors()[245], lwd=lwd);
      }
    }
    ## 
    dat.in[,1] <- gsub("chr", "", dat.in[,1]);
    dat.i      <- c(col.v:ncol(dat.in));
    dat.m      <- dat.in[,dat.i];
    dat.m      <- as.matrix(dat.m);
    
    ## for the matrix colors
    num.col <- ncol(dat.m);
    num.row <- nrow(dat.m);
    
    my.R1 <- R + W/5;
    my.R2 <- R + W - W/5;
    dat.min <- min(as.numeric(dat.m), na.rm=T);
    dat.max <- max(as.numeric(dat.m), na.rm=T);
    dat.med <- median(as.numeric(dat.m), na.rm=T);
    
    colors <- rainbow(10, alpha=0.3)[c(1,7)];
    
    for (j in 1:num.col){
      
      for (i in 1:num.row){
        dat.i <- as.numeric(dat.m[i,j]);
        if (is.na(dat.i)){
          next;
        }
        chr.i <- which(chr.po[,1]==dat.in[i,1]);
        
        if (dat.i > dat.med){
          col <- colors[1];
        } else {
          col <- colors[2];
        }
        
        v     <- scale.v(dat.i, my.R1, my.R2, dat.min, dat.max);
        
        v1 <- as.numeric(chr.po[chr.i,2]);
        v2 <- as.numeric(chr.po[chr.i,3]);
        v3 <- as.numeric(chr.po[chr.i,6]);
        v4 <- as.numeric(chr.po[chr.i,7]);
        
        w.to <- scale.v(as.numeric(dat.in[i,2]), v1, v2, v3, v4);
        
        draw.point.w(xc, yc, v, w.to, col=col, cex=0.1);
        
      }
    }
    
    if (scale){
      do.scale(xc, yc, dat.min, dat.max, R+W/5, W-2*(W/5));
    }
    
  }
  ### end ms
  
  ### hist
  if (type == "hist"){
    
    if (class(mapping)[1] == "GRanges"){
      po       <- as.integer((start(mapping)+end(mapping))/2+0.5);
      dat.in   <- cbind(as.character(seqnames(mapping)),
                        po, names(mapping), mcols(mapping));   
    } else {
      dat.in   <- mapping;
    }
    dat.in[,1] <- gsub("chr", "", dat.in[,1]);
    
    dat.in[dat.in[,1]==23,1] <- "X";
    dat.in[dat.in[,1]==24,1] <- "Y";
    
    dat.i      <- c(col.v:ncol(dat.in));
    dat.m      <- dat.in[,dat.i];
    dat.m      <- data.matrix(dat.m, );
    
    dat.min    <- min(as.numeric(dat.m), na.rm=T);
    dat.max    <- max(as.numeric(dat.m), na.rm=T);
    
    dat.df     <- as.data.frame(t(dat.m));
    
    hist.max   <- ncol(dat.m)*0.75;    
    
    cut.num    <- 16;
    region.num <- cut.num-1;
    ## breaks
    breaks     <- seq(dat.min, dat.max, 
                      length.out=cut.num);
    ## get freq
    out        <- lapply(dat.df, cut, br=breaks, labels=c(1:region.num));
    
    ## for color
    color1 <- rainbow(10, alpha=0.5)[7];
    
    w.i    <- (W-W/10)/region.num;
    
    for (i in 1:nrow(dat.in)){
      out.i <- table(out[[i]]);     
      
      chr.i <- which(dat.in[i,1]==chr.po[,1]);
      
      v1 <- as.numeric(chr.po[chr.i,2]);
      v2 <- as.numeric(chr.po[chr.i,3]);
      v3 <- as.numeric(chr.po[chr.i,6]);
      v4 <- as.numeric(chr.po[chr.i,7]);
      w.to <- scale.v(as.numeric(dat.in[i,2]), v1, v2, v3, v4);
      draw.line(xc,  yc, w.to, r, r+W-5, col=color1, lwd=0.01);
      
      my.R <- r;
      for (i in 1:region.num){
        my.R <- my.R + w.i;
        the.n <- out.i[i];
        if (the.n == 0){
          next;
        }
        w.from <- scale.v(the.n, w.to, w.to+360/nrow(dat.m), 0, hist.max-1);
        draw.arc.s(xc, yc, my.R, w.from, w.to, col="lightblue", 
                   lwd=w.i, lend=1);
        
      }
    }
    
    if (scale){
      do.scale(xc, yc, dat.min, dat.max, R+W/5, W-2*(W/5));
    }
    
  }
  ### end hist
  
  ### highlight
  if (type == "hl"){
    if (class(mapping)[1] == "GRanges"){
      po       <- as.integer((start(mapping)+end(mapping))/2+0.5);
      dat.in   <- cbind(as.character(seqnames(mapping)),
                        po, names(mapping), mcols(mapping));   
    } else {
      dat.in   <- mapping;
    }
    
    r1       <- as.numeric(dat.in[1]);
    r2       <- as.numeric(dat.in[2]);
    seg1     <- dat.in[3];
    seg2     <- dat.in[5];
    po1      <- as.numeric(dat.in[4]);
    po2      <- as.numeric(dat.in[6]);
    the.col1 <- dat.in[7];
    the.col2 <- dat.in[8];
    
    seg.i   <- which(seg1==chr.po[,1]);
    
    v1 <- as.numeric(chr.po[seg.i,2]);
    v2 <- as.numeric(chr.po[seg.i,3]);
    v3 <- as.numeric(chr.po[seg.i,6]);
    v4 <- as.numeric(chr.po[seg.i,7]);
    
    w1 <- scale.v(as.numeric(po1), v1, v2, v3, v4);
    
    seg.i   <- which(seg2==chr.po[,1]);
    
    v1 <- as.numeric(chr.po[seg.i,2]);
    v2 <- as.numeric(chr.po[seg.i,3]);
    v3 <- as.numeric(chr.po[seg.i,6]);
    v4 <- as.numeric(chr.po[seg.i,7]);
    
    w2 <- scale.v(as.numeric(po2), v1, v2, v3, v4);
    
    draw.arc.pg (xc, yc, w1, w2, r1, r2, col=the.col1, border=the.col2, lwd=lwd);
  }
  ### end highlight
  ### end type  ################
}

##############################
### functions
############################## 
bezierCurve <- function(x, y, n=10)	{	
  outx <- NULL	
  outy <- NULL 	
  i <- 1	
  for (t in seq(0, 1, length.out=n))		{		
    b <- bez(x, y, t)		
    outx[i] <- b$x		
    outy[i] <- b$y 		
    i <- i+1		
  } 	
  return (list(x=outx, y=outy))	
} 

###
bez <- function(x, y, t)	{	
  outx <- 0	
  outy <- 0	
  n <- length(x)-1	
  for (i in 0:n)		{		
    outx <- outx + choose(n, i)*((1-t)^(n-i))*t^i*x[i+1]		
    outy <- outy + choose(n, i)*((1-t)^(n-i))*t^i*y[i+1]		
  } 	
  return (list(x=outx, y=outy))	
}

###########################################
# one value : from a to b 
scale.v <- function(v, a, b, min.v, max.v) {
  v <- v-min.v; 
  v <- v/(max.v-min.v); 
  v <- v*(b-a);  
  v+a
}

## permutation index
perm_list <- function (n, r, v = 1:n){ 
  if (r == 1) 
    X <- matrix(v, n, 1) 
  else if (n == 1) 
    X <- matrix(v, 1, r) 
  else { 
    X <- NULL 
    for (i in 1:n){ 
      X <- rbind(X, cbind(v[i], perm_list(n-1 , r-1 , v[-i]))) 
    } 
  } 
  return(X);
} 

### draw.link
draw.link <- function(xc, yc, r, w1, w2, col=col, lwd=lwd) {
  # for translocation
  w3  <- (w1+w2)/2;
  w1  <- w1/360*2*pi;
  w2  <- w2/360*2*pi;
  w3  <- w3/360*2*pi;
  x0  <- xc+r*cos(w1);
  y0  <- yc-r*sin(w1);
  x1  <- xc+r*cos(w2);
  y1  <- yc-r*sin(w2);
  x <- c(x0,xc,xc,x1);
  y <- c(y0,yc,yc,y1);
  points(bezierCurve(x,y,60), type="l", col=col, lwd=lwd, lend="butt")
}

### draw.link2
draw.link2 <- function(xc, yc, r, w1, w2, col=col, lwd=lwd) {
  # for translocation
  w3  <- (w1+w2)/2;
  w1  <- w1/360*2*pi;
  w2  <- w2/360*2*pi;
  w3  <- w3/360*2*pi;
  x2  <- xc+r/2*cos(w3);
  y2  <- yc-r/2*sin(w3);
  x0  <- xc+r*cos(w1);
  y0  <- yc-r*sin(w1);
  x1  <- xc+r*cos(w2);
  y1  <- yc-r*sin(w2);
  x <- c(x0, x2, x2, x1);
  y <- c(y0, y2, y2, y1);
  points(bezierCurve(x,y,60), type="l", col=col, lwd=lwd, lend="butt")
}
### 

###
draw.link.pg <- function(xc, yc, r, w1.1, w1.2, w2.1, w2.2, col=col, lwd=lwd) {
  #####################################################
  w1 <- w1.1;
  w2 <- w2.2;
  w3  <- (w1+w2)/2;
  w1  <- w1/360*2*pi;
  w2  <- w2/360*2*pi;
  w3  <- w3/360*2*pi;
  x0  <- xc+r*cos(w1);
  y0  <- yc-r*sin(w1);
  x1  <- xc+r*cos(w2);
  y1  <- yc-r*sin(w2);
  x <- c(x0,xc,xc,x1);
  y <- c(y0,yc,yc,y1);
  bc1 <- bezierCurve(x,y,60);
  
  ang.d <- abs(w1.1-w1.2);
  pix.n <- ang.d * 10;
  if (pix.n < 10){
    pix.n <- 10;
  }
  
  ang.seq <- rev(seq(w1.1,w1.2,length.out=pix.n));
  ang.seq <- ang.seq/360*2*pi;
  
  fan.1.x <- xc + cos(ang.seq) * r;
  fan.1.y <- yc - sin(ang.seq) * r;
  
  ######################################################
  w1 <- w1.2;
  w2 <- w2.1;
  w3  <- (w1+w2)/2;
  w1  <- w1/360*2*pi;
  w2  <- w2/360*2*pi;
  w3  <- w3/360*2*pi;
  x0  <- xc+r*cos(w1);
  y0  <- yc-r*sin(w1);
  x1  <- xc+r*cos(w2);
  y1  <- yc-r*sin(w2);
  x <- c(x0,xc,xc,x1);
  y <- c(y0,yc,yc,y1);
  bc2 <- bezierCurve(x,y,60);
  
  ang.d <- abs(w2.1-w2.2);
  pix.n <- ang.d * 10;
  if (pix.n < 10){
    pix.n <- 10;
  }
  
  ang.seq <- rev(seq(w2.1,w2.2,length.out=pix.n));
  ang.seq <- ang.seq/360*2*pi;
  
  fan.2.x <- xc + cos(ang.seq) * r;
  fan.2.y <- yc - sin(ang.seq) * r;
  
  polygon(c(bc1$x, fan.2.x, rev(bc2$x), rev(fan.1.x)), 
          c(bc1$y, fan.2.y, rev(bc2$y), rev(fan.1.y)), 
          fillOddEven=TRUE, border=col, col=col, lwd=lwd); 
}

###
draw.point.w <- function(xc, yc, r, w, col=col, cex=cex){
  w <- w/360*2*pi;
  x <- xc+r*cos(w);
  y <- yc-r*sin(w);
  points(x, y, pch=20, col=col, cex=cex);
}

###
draw.text.w <- function(xc, yc, r, w, n, col="black", cex=1){
  w <- w%%360;
  w <- w/360*2*pi;
  x <- xc+r*cos(w);
  y <- yc-r*sin(w);
  text(x,y,labels=n, col=col, cex=cex);
}

###
draw.text.rt <- function(xc, yc, r, w, n, col="black", cex=1, side="out"){
  
  w     <- w%%360;
  the.o <- w;
  
  if (w <= 90){
    #w=0.0219780219780059+0.978021978021979*w;
  } else if (w > 90 & w <= 180){
    #w=3.91304347826081+0.978260869565218*w;
  } else if (w > 180 & w <= 270){
    #w=3.91304347826091+0.978260869565217*w;
  } else if (w > 270 & w <= 360){
    #w=7.82608695652158+0.978260869565218*w;
  }
  
  the.w <- 360-w;
  w     <- w/360*2*pi;
  x     <- xc+r*cos(w);
  y     <- yc-r*sin(w);
  
  
  num2  <- 26;
  
  if (side=="out"){
    if (the.w <= 90 ){
      the.pos <- 4;
    } else if (the.w > 90 & the.w <= 180) {
      the.w <- the.w + 180;
      the.pos <- 2;
    } else if (the.w > 180 & the.w <= 270){
      the.w <- the.w%%180;
      the.pos <- 2;
    } else if (the.w > 270 & the.w <= 360){
      the.pos <- 4;
    }
    
    if (the.pos==2){
      x <- x+num2;
    }
    if (the.pos==4){
      x <- x-num2;
    }
  } 
  
  if (side=="in"){
    if (the.w <= 90 ){
      the.pos <- 4;
    } else if (the.w > 90 & the.w <= 180) {
      the.w <- the.w + 180;
      the.pos <- 2;
    } else if (the.w > 180 & the.w <= 270){
      the.w <- the.w%%180;
      the.pos <- 2;
    } else if (the.w > 270 & the.w <= 360){
      the.pos <- 4;
    }
    
    if (the.pos==2){
      x <- x+num2;
    }
    if (the.pos==4){
      x <- x-num2;
    }
  }
  
  text(x, y, adj=0, offset=1, labels=n, srt=the.w, 
       pos=the.pos, col=col, cex=cex);
}

###strokeLine2
draw.line <- function (xc, yc, w, l1, l2, col=col, lwd=lwd, lend=1) {
  w  <- (w/360)*2*pi;
  x1 <- xc+l1*cos(w);
  y1 <- yc-l1*sin(w);
  x2 <- xc+l2*cos(w);
  y2 <- yc-l2*sin(w);
  segments(x1, y1, x2, y2, col=col, lwd=lwd, lend=lend);
}

###strokeLine3
draw.line2 <- function (xc, yc, w, r, l, col=col, lwd=lwd){
  line_w   <- l;
  theangle <- w;
  l1       <- r;
  theangle <- (theangle/360)*2*pi;
  x0       <- xc+l1*cos(theangle);
  y0       <- yc+l1*sin(theangle);
  w1       <- 45/360*2*pi;
  x1 = xc + sin(w1) * (x0);
  y1 = yc + cos(w1) * (y0);
  x2 = xc - sin(w1) * (x0);
  y2 = yc - cos(w1) * (y0);
  segments(x1, y1, x2, y2, col=col, lwd=lwd, lend="butt");
}

###strokeLine by two angles
draw.line3 <- function (xc, yc, w1, w2, r1, r2, col=col, lwd=lwd){
  theangle1 <- w1;
  theangle2 <- w2;
  l1        <- r1;
  l2        <- r2;
  
  theangle1 <- (theangle1/360)*2*pi;
  x1        <- xc+l1*cos(theangle1);
  y1        <- yc-l1*sin(theangle1);
  
  theangle2 <- (theangle2/360)*2*pi;
  x2        <- xc+l2*cos(theangle2);
  y2        <- yc-l2*sin(theangle2);
  
  segments(x1, y1, x2, y2, col=col, lwd=lwd, lend="butt");
}

### plot fan or sector that likes a piece of doughnut (plotFan)
draw.arc.pg <- function (xc, yc, 
                         w1, w2, r1, r2, col="lightblue", border="lightblue", lwd=0.01
){
  
  ang.d <- abs(w1-w2);
  pix.n <- ang.d * 10;
  if (pix.n < 10){
    pix.n <- 10;
  }
  
  ang.seq <- rev(seq(w1,w2,length.out=pix.n));
  ang.seq <- ang.seq/360*2*pi;
  
  fan.i.x <- xc + cos(ang.seq) * r1;
  fan.i.y <- yc - sin(ang.seq) * r1;
  
  
  fan.o.x <- xc + cos(ang.seq) * r2;
  fan.o.y <- yc - sin(ang.seq) * r2;
  
  polygon(c(rev(fan.i.x), fan.o.x ), c(rev(fan.i.y), fan.o.y), 
          fillOddEven=TRUE, border=border, col=col, lwd=lwd, lend=1)
  
}

###
draw.arc.s <- function (xc, yc, r, w1, w2, col="lightblue", lwd=1, lend=1){
  # s = simple
  # r = radius
  
  ang.d <- abs(w1-w2);
  pix.n <- ang.d * 5;
  if (pix.n < 2){
    pix.n <- 2;
  }
  
  ang.seq <- rev(seq(w1,w2,length.out=pix.n));
  ang.seq <- ang.seq/360*2*pi;
  
  fan.i.x <- xc + cos(ang.seq) * r;
  fan.i.y <- yc - sin(ang.seq) * r;
  ## lend=0(round); lend=1(butt); lend=2(square)
  lines(fan.i.x, fan.i.y, col=col, lwd=lwd, type="l", lend=lend);
  #points(fan.i.x, fan.i.y, col=col, lwd=lwd, type="l", lend=lend);
}

#########################################
## segment to angle and position
## segAnglePo
#########################################

# get angle if given seg and position
# seg=chromosome in genome, shuch as: 1:22, 23="X", 24="Y";

setGeneric("segAnglePo", function(seg.dat=seg.dat, seg=seg, angle.start=angle.start, 
                                  angle.end=angle.end){standardGeneric("segAnglePo")})

setMethod("segAnglePo", "data.frame", function(seg.dat=seg.dat, seg=seg, 
                                               angle.start=angle.start, angle.end=angle.end){
  ## 
  seg.dat  <- seg.dat;
  .segAnglePo(seg.dat=seg.dat, seg=seg, angle.start=angle.start, angle.end=angle.end);
})

setMethod("segAnglePo", "GRanges", function(seg.dat=seg.dat){
  ## 
  seg.dat  <- cbind(as.character(seqnames(seg.dat)), start(seg.dat), end(seg.dat), 
                    names(seg.dat), as.character(strand(seg.dat)));
  colnames(seg.dat) <- c("seg.name","seg.Start","seg.End","name","gieStain");
  .segAnglePo(seg.dat=seg.dat, seg=seg, angle.start=angle.start, angle.end=angle.end);
})

## seg should be ordered by user
.segAnglePo <- function (seg.dat=seg.dat, seg=seg, angle.start=angle.start, 
                         angle.end=angle.end){
  
  if (missing(angle.start)){
    angle.start <- 0;
  }
  if (missing(angle.end)){
    angle.end <- 360;
  }
  ## check data.frame?
  colnames(seg.dat) <- c("seg.name","seg.Start","seg.End","name","gieStain");
  
  ## get length of the segomosomes
  seg.l   <- c();
  seg.min <- c();
  seg.sum <- c();
  seg.s   <- 0;
  seg.num   <- length(seg);
  seg.names <- seg;
  
  ########################################################
  ########################################################
  for (i in 1:seg.num){
    seg.n <- seg.names[i];
    dat.m <- subset(seg.dat, seg.dat[,1]==seg.n);
    seg.full.l <- max(as.numeric(dat.m[,"seg.End"]));
    seg.full.m <- min(as.numeric(dat.m[,"seg.Start"]));
    seg.l      <- cbind(seg.l, seg.full.l);
    seg.min    <- cbind(seg.min, seg.full.m);
    seg.s      <- seg.s + seg.full.l;
    seg.sum    <- cbind(seg.sum, seg.s);
  }
  
  ## initial parameters
  
  gap.angle.size <- 2;
  seg.angle.from <- angle.start + 270;
  
  seg.full.l  <- sum(as.numeric(seg.l));
  angle.range <- angle.end - angle.start;
  cir.angle.r <- (angle.range - seg.num * gap.angle.size)/seg.full.l;
  
  out.s     <- c();
  l.old     <- 0;
  gap.angle <- 0;
  for (i in 1:seg.num){
    seg.n <- seg.names[i];
    dat.m <- subset(seg.dat, seg.dat[,1]==seg.n);
    len   <- seg.sum[i];
    w1    <- cir.angle.r*l.old + gap.angle;
    w2    <- cir.angle.r*len   + gap.angle;
    out.s     <- rbind(out.s, c(seg.n, w1+seg.angle.from, w2+seg.angle.from, l.old, len, seg.min[i], seg.l[i]));
    gap.angle <- gap.angle + gap.angle.size;
    l.old     <- len;
  }
  
  colnames(out.s) <- c("seg.name","angle.start", "angle.end", "seg.sum.start", "seg.sum.end","seg.start", "seg.end");
  #write.table(out.s, outfile, quote=F, row.names=F, sep="\t");
  return(out.s);
}

sim.circos <- function (seg=10, po=c(20,50), ind=10, link=10, 
                        link.pg=10){
  seg.num    <- seg;
  seg.po     <- po;
  ind.num    <- ind;
  link.num   <- link;
  link.w.num <- link.pg;
  
  ################################################
  ## segment, poiter, individual 
  ################################################
  seg.out <- c();
  seg.v   <- c();
  
  for (i in 1:seg.num){
    po.num <- sample(seg.po, 1);
    i.s    <- paste("chr", i, sep="");
    for (j in 1:po.num){
      seg.out <- rbind(seg.out, c(i.s, j-1, j, "NA", "NA"));
      c.v <- c();
      for (k in 1:ind.num){
        if (k > ind.num/2){
          the.v   <- round(rnorm(1) + i + i/ind.num , 3);
        } else {
          the.v   <- round(rnorm(1) + i,3);
        }
        c.v     <- c(c.v, the.v);
      }
      seg.v   <- rbind(seg.v,   c(i.s, j, c.v));
    }
  }
  
  colnames(seg.out) <- c("seg.name", "seg.Start", "seg.End", "the.v", "NO")
  seg.out <- as.data.frame(seg.out);
  names   <- paste("name", 1:ind.num, sep="");
  colnames(seg.v) <- c("seg.name", "seg.po", names)
  seg.v   <- as.data.frame(seg.v);
  
  ################################################
  ## link data
  ################################################
  link.out <- c();
  for (i in 1:link.num){
    name <- paste("n", i, sep="");
    c1 <- sample(c(1:nrow(seg.out)), 1);
    c2 <- sample(c(1:nrow(seg.out)), 1);
    link.out <- rbind(link.out, c(seg.out[c1,1], seg.out[c1,2], name, 
                                  seg.out[c2,1], seg.out[c2,2], name, name));
  }
  colnames(link.out) <- c("seg1", "po1", "name1", "seg2", "po2", "name2", "name3")
  link.out <- as.data.frame(link.out);
  
  ################################################
  ## link data with wide
  ################################################
  link.w.out <- c();
  for (i in 1:link.w.num){
    n.i <- sample(1:seg.num, 2);
    i1    <- paste("chr", n.i[1], sep="");
    i2    <- paste("chr", n.i[2], sep="");
    n1    <- subset(seg.out, seg.out[,1]==i1);
    n2    <- subset(seg.out, seg.out[,1]==i2);
    re1   <- sample(n1[,2], 2);
    re2   <- sample(n2[,2], 2);
    link.w.out <- rbind(link.w.out, c(i1, re1, i2, re2));
  }
  colnames(link.w.out) <- c("seg1", "start1", "end1", "seg2", "start2", "end2");
  link.w.out <- as.data.frame(link.w.out)
  
  sim.out <- list(seg.frame=seg.out, seg.mapping=seg.v, 
                  seg.link = link.out, seg.link.pg = link.w.out);
  
  return(sim.out);
}


## color bar
## xl=xleft, yb=ybottom, xr=xright, yt=ytop
color.bar <- function(xl, yb, xr, yt, v.min, v.max) {
  
  nticks <- 11; 
  ticks  <- seq(v.min, v.max, len=nticks);
  
  lut   <- colorRampPalette(c("blue", "white", "red"))(100);
  scale <- (length(lut)-1)/(v.max-v.min);     
  ys    <- (yt-yb)/(length(lut)-1);
  
  yb.old <- yb;
  for (i in 1:(length(lut)-1)) {
    rect(xl, yb.old, xr, yb.old+ys, col=lut[i], border=NA)
    yb.old <- yb.old + ys;
  }
  v.med <- (v.max+v.min)/2;
  
  ym  <- (yt+yb)/2;
  
  v.max <- round(v.max,2);
  v.min <- round(v.min,2);
  v.med <- round(v.med,2);
  text(xl, yt+30, "heatmap", cex=1);
  text(xl-30, yt, v.max, cex=0.6);
  text(xl-30, yb, v.min, cex=0.6);
  text(xl-30, ym, v.med, cex=0.6);
  segments(xl-5, yb, xl-5,  yt, col="black", lwd=0.8);
  segments(xl-5, yb, xl-15, yb, col="black", lwd=0.8);
  segments(xl-5, yt, xl-15, yt, col="black", lwd=0.8);
  segments(xl-5, ym, xl-15, ym, col="black", lwd=0.8);
  
}

### heatmap.cluster
heatmap.cluster <- function(x1, y1, x2, y2, dat=dat){
  dat.d  <- dist(t(dat));
  dat.h  <- hclust(dat.d);
  
  ####################
  c.x1   <- x1;
  c.x2   <- x2;
  c.y1   <- y1;
  c.y2   <- y2;
  max.h  <- max(dat.h$height);
  max.n  <- length(dat.h$labels);
  
  lab   = (c.x2-c.x1)/max.n;
  ratio = (c.y2-c.y1)/max.h;
  
  ####################  
  sub2h  <- c();
  sub2po <- c();
  y0 <- c.y1 - 2;
  for (i in 1:nrow(dat.h$merge)){
    le <- dat.h$merge[i,1];
    ri <- dat.h$merge[i,2];
    if (le < 0 && ri < 0){
      po1 <- which(dat.h$order==abs(le));
      po2 <- which(dat.h$order==abs(ri));
      y1  <- y0;
      y2  <- y0;
    } else if (le < 0){ 
      po1 <- which(dat.h$order==abs(le));
      po2 <- sub2po[ri,2];
      y1  <- y0;
      y2  <- c.y1 + ratio * sub2h[ri,2];
    } else if (ri < 0){
      po1 <- sub2po[le, 2];
      po2 <- which(dat.h$order==abs(ri));
      y1  <- c.y1 + ratio * sub2h[le,2];
      y2  <- y0;
    } else {
      po1 = sub2po[le, 2];
      po2 = sub2po[ri, 2];
      y1  = c.y1 + ratio * sub2h[le,2];
      y2  = c.y1 + ratio * sub2h[ri,2];
    }
    sub2po <- rbind(sub2po, c(i, (po1+po2)/2));
    sub2h  <- rbind(sub2h,  c(i, dat.h$height[i]));
    x1 = c.x1 + lab*po1;
    x2 = c.x1 + lab*po2;
    y3 = c.y1 + ratio*dat.h$height[i];
    segments(x1,y3,x2,y3, lwd=0.1);
    # le
    segments(x1,y3,x1,y1, lwd=0.1);
    # ri
    segments(x2,y3,x2,y2, lwd=0.1);
  }
  
  y0 <- c.y1 - 8;
  for (i in 1:length(dat.h$order)){
    n  <- dat.h$label[dat.h$order[i]];
    x0 <- c.x1 + lab * i - 24;
    text(x0, y0, n, srt=270, col="blue", cex=0.2, pos=4, adj=0, offset=1);
  }
  ####################
}
### end heatmap.cluster

### start do.scale
do.scale <- function(xc=xc, yc=yc, dat.min=dat.min, dat.max=dat.max, 
                     R=R, W=W, s.n=1, col="blue"){
  dat.m   <- round((dat.min+dat.max)/2, s.n);
  dat.min <- round(dat.min, s.n);
  dat.max <- round(dat.max, s.n);
  y1      <- yc + R ;
  y2      <- yc + R + W/2;
  y3      <- yc + R + W;
  x1      <- xc - W/20;
  x2      <- x1 - (W/20)*1.2;
  x3      <- x1 - (W/20)*3;
  segments(x1, y1, x1, y3, lwd=0.01, col=col);
  segments(x1, y1, x2, y1, lwd=0.01, col=col);
  segments(x1, y2, x2, y2, lwd=0.01, col=col);
  segments(x1, y3, x2, y3, lwd=0.01, col=col);
  text(x3, y1, dat.min, cex=0.2, col=col);
  text(x3, y2, dat.m,   cex=0.2, col=col);
  text(x3, y3, dat.max, cex=0.2, col=col);
}
### end do.scale

### start do.scale.cir
do.scale.cir <- function (xc=xc, yc=yc, the.r=the.r, total.num=total.num, 
                          col="blue", lwd=0.001,
                          V1=V1, V2=V2, V3=V3, V4=V4
){
  
  ### scale initial start ##################
  ## scale label number = 150 in 360 degree circumference
  ## the.r  = R
  ## sum.po = total point number 
  ## should return: 1) scale number 2) scale unit
  sum.po  <- as.numeric(total.num);         # total number
  scale.w <- sum.po/150;                         # one scale size
  scale.l <- nchar(as.integer(scale.w))-1;       # length of the number
  scale.i <- as.integer(scale.w/(10^scale.l));   # the first digital
  scale.d <- scale.i * 10^scale.l;               # the first digital, then 0
  scale.m <- 1 * 10^scale.l;                     # the unit of the scale   
  ### scale initial end   ##################
  
  ### scale start ######################
  draw.arc.s(xc, yc, w1=V1, w2=V2, r=the.r, col=col, lwd=lwd);
  
  start.p  <- 0;
  w        <- scale.v(as.numeric(start.p), V1, V2, V3, V4);  
  scale.n  <- as.integer(as.numeric(V4)/scale.d+1);
  
  for (j in 1:scale.n){
    po.s <- as.integer(start.p/scale.m);
    
    draw.line(xc,    yc,   w, the.r, the.r+4, col=col, lwd=lwd);
    draw.text.rt(xc, yc,  the.r+6, w, po.s, col=col, cex=0.2);
    
    start.p <- start.p + scale.d/2;
    w        <- scale.v(as.numeric(start.p), V1, V2, V3, V4);
    
    if (w <= V2){  
      draw.line(xc, yc, w, the.r, the.r+2, col=col, lwd=lwd);
    }
    
    start.p <- start.p + scale.d/2;
    w        <- scale.v(as.numeric(start.p), V1, V2, V3, V4);
  }
  ### scale end ######################
  
}
### end do.scale.cir

### zoom ###
zoom.in <- function(cir.in=cir.in, zoom=zoom){
  
  out.cir  <- c();
  out.dat  <- c();
  
  ## for all data sets
  chr1     <- which(cir.in[,1]==zoom[1]);
  chr2     <- which(cir.in[,1]==zoom[2]);
  
  v1 <- as.numeric(cir.in[chr1,2]);
  v2 <- as.numeric(cir.in[chr1,3]);
  v3 <- as.numeric(cir.in[chr1,6]);
  v4 <- as.numeric(cir.in[chr1,7]);
  w3 <- scale.v(zoom[3], v1, v2, v3, v4);
  
  v1 <- as.numeric(cir.in[chr2,2]);
  v2 <- as.numeric(cir.in[chr2,3]);
  v3 <- as.numeric(cir.in[chr2,6]);
  v4 <- as.numeric(cir.in[chr2,7]);
  w4 <- scale.v(zoom[4], v1, v2, v3, v4);
  
  for (i in chr1:chr2){
    w1 <- as.numeric(zoom[5])+270;
    w2 <- as.numeric(zoom[6])+270;
    
    v1 <- scale.v(as.numeric(cir.in[i,2]), w1, w2, w3, w4);
    v2 <- scale.v(as.numeric(cir.in[i,3]), w1, w2, w3, w4);
    out.cir <- rbind(out.cir, c(as.character(cir.in[i,1]), v1, v2, cir.in[i,4:ncol(cir.in)]));
  }
  return(out.cir);
} 
### end zoom ###

###
# chrOrder <- c(paste("chr",1:22,sep=""),"chrX","chrY","chrM")
# df$chr   <- factor(df$chr, levels=chrOrder)
# df[order(df$chr),]
###

###################################################################
## from gtools package
mixedorder <- function (x) 
{
  delim = "\\$\\@\\$"
  numeric <- function(x) {
    optwarn = options("warn")
    on.exit(options(optwarn))
    options(warn = -1)
    as.numeric(x)
  }
  nonnumeric <- function(x) {
    optwarn = options("warn")
    on.exit(options(optwarn))
    options(warn = -1)
    ifelse(is.na(as.numeric(x)), toupper(x), NA)
  }
  x <- as.character(x)
  which.nas <- which(is.na(x))
  which.blanks <- which(x == "")
  if (length(which.blanks) > 0) 
    x[which.blanks] <- -Inf
  if (length(which.nas) > 0) 
    x[which.nas] <- Inf
  delimited <- gsub("([+-]{0,1}[0-9]+([eE][\\+\\-]{0,1}[0-9]+){0,1})", 
                    paste(delim, "\\1", delim, sep = ""), x)
  step1 <- strsplit(delimited, delim)
  step1 <- lapply(step1, function(x) x[x > ""])
  step1.numeric <- lapply(step1, numeric)
  step1.character <- lapply(step1, nonnumeric)
  maxelem <- max(sapply(step1, length))
  step1.numeric.t <- lapply(1:maxelem, function(i) sapply(step1.numeric, 
                                                          function(x) x[i]))
  step1.character.t <- lapply(1:maxelem, function(i) sapply(step1.character, 
                                                            function(x) x[i]))
  rank.numeric <- sapply(step1.numeric.t, rank)
  rank.character <- sapply(step1.character.t, function(x) as.numeric(factor(x)))
  rank.numeric[!is.na(rank.character)] <- 0
  rank.character <- t(t(rank.character) + apply(matrix(rank.numeric), 
                                                2, max, na.rm = TRUE))
  rank.overall <- ifelse(is.na(rank.character), rank.numeric, 
                         rank.character)
  order.frame <- as.data.frame(rank.overall)
  if (length(which.nas) > 0) 
    order.frame[which.nas, ] <- Inf
  retval <- do.call("order", order.frame)
  return(retval)
}

