
## Enrichment
## By Haitao Wang
## Dr. Haitao Wang: JNU,Guagnzhou => GIBH-CAS,Guagnzhou => BNU,Beijing => IMCB,Singapore => UMAC,Macao => NCCS,Singapore


##########################################################
### genelist filtering based on given background genes ###
##########################################################
# function to match and filter geneset from your background genelist，and remove duplicates
# useage filtered.GeneSets <- GeneSetsFilterByBackground(GeneSets,bglist$GeneSym,100,2000)

GeneSetsFilterByBackground <- function(GSet, ## A list of gene vectors.
                                       total.genelist, ## A vector with all genes, including the genes of interest.
                                       minSize=2, ## minimum gene number.
                                       maxSize=1000){ ## maximum gene number.
  finalGSet<-list(); ## A list of gene vectors from GSet which filtered by background genelist (total.genelist).
  total.genelist<-toupper(total.genelist);
  for(i in 1:length(GSet)){
    print(paste("Processing", i));
    probes<-GSet[[i]] 
    
    allProbes<-intersect(toupper(probes),total.genelist);
    if(length(allProbes)>=minSize && length(allProbes)<=maxSize){
      finalGSet[[names(GSet)[i]]]<-allProbes
    }
  }
  return(finalGSet);
}

# based on your clustering to reorder matrix
reorder <- function(x, o){
  v <- x;
  v[lower.tri(v)] <- t(x)[lower.tri(x)]
  v <- v[o,o];
  y <- x;
  y[upper.tri(v)] <- t(x)[upper.tri(x)]
  y <- y[o,o];
  Z <- v;
  Z[lower.tri(Z)] <- y[lower.tri(Z)]
  return(Z)
}

