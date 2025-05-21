# Copyright (C) 2011  John Colby
# http://github.com/johncolby/SVM-RFE

svmRFE.wrap <- function(test.fold, X, ...) {
# Wrapper to run svmRFE function while omitting a given test fold
    train.data = X[-test.fold, ]
    test.data  = X[test.fold, ]
    
    # Rank the features
    features.ranked = svmRFE(train.data, ...)
    
    return(list(feature.ids=features.ranked, train.data.ids=row.names(train.data), test.data.ids=row.names(test.data)))
}

svmRFE <- function(X, k=1, halve.above=5000) {
# Feature selection with Multiple SVM Recursive Feature Elimination (RFE) algorithm
    n = ncol(X) - 1
    
    # Scale data up front so it doesn't have to be redone each pass
    cat('Scaling data...')
    X[, -1] = scale(X[, -1])
    cat('Done!\n')
    flush.console()
    
    pb = txtProgressBar(1, n, 1, style=3)

    i.surviving = 1:n
    i.ranked    = n
    ranked.list = vector(length=n)

    # Recurse through all the features
    while(length(i.surviving) > 0) {
        if(k > 1) {
            # Subsample to obtain multiple weights vectors (i.e. mSVM-RFE)            
            folds = rep(1:k, len=nrow(X))[sample(nrow(X))]
            folds = lapply(1:k, function(x) which(folds == x))
            
            # Obtain weights for each training set
            w = lapply(folds, getWeights, X[, c(1, 1+i.surviving)])
            w = do.call(rbind, w)
            
            # Normalize each weights vector
            w = t(apply(w, 1, function(x) x / sqrt(sum(x^2))))
            
            # Compute ranking criteria
            v    = w * w
            vbar = apply(v, 2, mean)
            vsd  = apply(v, 2, sd)
            c    = vbar / vsd
        } else {
            # Only do 1 pass (i.e. regular SVM-RFE)
            w = getWeights(NULL, X[, c(1, 1+i.surviving)])
            c = w * w
        }

        # Rank the features
        ranking = sort(c, index.return=T)$ix
        if(length(i.surviving) == 1) {
            ranking = 1
        }
        
        if(length(i.surviving) > halve.above) {
            # Cut features in half until less than halve.above
            nfeat = length(i.surviving)
            ncut  = round(nfeat / 2)
            n     = nfeat - ncut
            
            cat('Features halved from', nfeat, 'to', n, '\n')
            flush.console()
            
            pb = txtProgressBar(1, n, 1, style=3)
            
        } else ncut = 1
        
        # Update feature list
        ranked.list[i.ranked:(i.ranked-ncut+1)] = i.surviving[ranking[1:ncut]]
        i.ranked    = i.ranked - ncut
        i.surviving = i.surviving[-ranking[1:ncut]]

        setTxtProgressBar(pb, n-length(i.surviving))
        flush.console()
    }
    
    close(pb)

    return (ranked.list)
}

getWeights <- function(test.fold, X) {
# Fit a linear SVM model and obtain feature weights
    train.data = X
    if(!is.null(test.fold)) train.data = X[-test.fold, ]
    
    svmModel = svm(train.data[, -1], train.data[, 1], cost=10, cachesize=500,
    scale=F, type="C-classification", kernel="linear")
    
    t(svmModel$coefs) %*% svmModel$SV
}

WriteFeatures <- function(results, input, save=T, file='features_ranked.txt') {
# Compile feature rankings across multiple folds
    featureID = sort(apply(sapply(results, function(x) sort(x$feature, index.return=T)$ix), 1, mean), index=T)$ix
    avg.rank  = sort(apply(sapply(results, function(x) sort(x$feature, index.return=T)$ix), 1, mean), index=T)$x
    feature.name = colnames(input[, -1])[featureID]
    features.ranked = data.frame(FeatureName=feature.name, FeatureID=featureID, AvgRank=avg.rank)
    if(save==T) {
        write.table(features.ranked, file=file, quote=F, row.names=F)
    } else {
        features.ranked
    }
}

FeatSweep.wrap <- function(i, results, input) {
# Wrapper to estimate generalization error across all hold-out folds, for a given number of top features
    svm.list = lapply(results, function(x) tune(svm,
                      train.x      = input[x$train.data.ids, 1+x$feature.ids[1:i]],
                      train.y      = input[x$train.data.ids, 1],
                      validation.x = input[x$test.data.ids, 1+x$feature.ids[1:i]],
                      validation.y = input[x$test.data.ids, 1],
                      # Optimize SVM hyperparamters
                      ranges       = tune(svm,
                                          train.x = input[x$train.data.ids, 1+x$feature.ids[1:i]],
                                          train.y = input[x$train.data.ids, 1],
                                          ranges  = list(gamma=2^(-12:0), cost=2^(-6:6)))$best.par,
                      tunecontrol  = tune.control(sampling='fix'))$perf)
        
    error = mean(sapply(svm.list, function(x) x$error))
    return(list(svm.list=svm.list, error=error))
}


PlotErrors <- function(errors, errors2=NULL, no.info=0.5, 
                       ylim=range(c(errors, errors2), na.rm=T), 
                       xlab='Number of Features',  ylab='5 x CV Error') {
# Makes a plot of average generalization error vs. number of top features
    AddLine <- function(x, col='dodgerblue') {
        lines(which(!is.na(errors)), na.omit(x), col=col,lwd=3)
        points(which.min(x), min(x, na.rm=T), col='firebrick3')
        text(which.min(x), min(x, na.rm=T), paste(which.min(x), '-', 
                                                  format(min(x, na.rm=T), dig=3)), pos=2, col='red', cex=1.15)
    }
    
    plot(errors, type='n', ylim=ylim, xlab=xlab, ylab=ylab)
    AddLine(errors)
    if(!is.null(errors2)) AddLine(errors2, 'gray30')
    abline(h=no.info, lty=2)
}


Plotaccuracy <- function(errors, errors2=NULL, no.info=0.5, 
                       ylim=range(c(errors, errors2), na.rm=T), 
                       xlab='Number of Features',  ylab='5 x CV Accuracy') {
  # Makes a plot of average generalization error vs. number of top features
  AddLine <- function(x, col='dodgerblue') {
    lines(which(!is.na(errors)), na.omit(x), col=col,lwd=3)
    points(which.max(x), max(x, na.rm=T), col='firebrick3')
    text(which.max(x), max(x, na.rm=T), paste(which.max(x), '-', 
                                              format(max(x, na.rm=T), dig=3)), pos=2, col='red', cex=1.15)
  }
  
  plot(errors, type='n', ylim=ylim, xlab=xlab, ylab=ylab)
  AddLine(errors)
  if(!is.null(errors2)) AddLine(errors2, 'gray30')
  abline(h=no.info, lty=2)
}



