#The function TOMdist1 computes a dissimilarity 
# based on the topological overlap matrix (Ravasz et al)
# Input: an Adjacency matrix with entries in [0,1] for PPI network.
# The following function is the source code (NetworkFunctions.txt) provided by Horvath Group
#(https://horvath.genetics.ucla.edu/html/ModuleConformity/ModuleNetworks/).
if(exists("TOMdist1")) rm(TOMdist1);
TOMdist1 <- function(adjmat1, maxADJ=FALSE) {
  diag(adjmat1)=0;
  adjmat1[is.na(adjmat1)]=0;
  maxh1=max(as.dist(adjmat1) ); minh1=min(as.dist(adjmat1) ); 
  if (maxh1>1 | minh1 < 0 ) {print(paste("ERROR: the adjacency matrix contains entries that are larger than 1 or smaller than 0!!!, max=",maxh1,", min=",minh1)) } else { 
    if (  max(c(as.dist(abs(adjmat1-t(adjmat1)))))>10^(-12)   ) {print("ERROR: non-symmetric adjacency matrix!!!") } else { 
      adjmat1= (adjmat1+ t(adjmat1) )/2
      kk=apply(adjmat1,2,sum)
      maxADJconst=1
      if (maxADJ==TRUE) maxADJconst=max(c(as.dist(adjmat1 ))) 
      Dhelp1=matrix(kk,ncol=length(kk),nrow=length(kk))
      denomTOM= pmin(as.dist(Dhelp1),as.dist(t(Dhelp1)))   +as.dist(maxADJconst-adjmat1); 
      gc();gc();
      numTOM=as.dist(adjmat1 %*% adjmat1 +adjmat1);
      #TOMmatrix=numTOM/denomTOM
      # this turns the TOM matrix into a dissimilarity 
      out1=1-as.matrix(numTOM/denomTOM) 
      diag(out1)=1 
      # setting the diagonal to 1 is unconventional (it should be 0)
      # but it leads to nicer looking TOM plots... 
      out1
    }}}
