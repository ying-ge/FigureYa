#' This function performs variable selection on gene expression matrices.
#' It can, for instance, remove genes with low variation.
#' @param exprMat A matrix of gene expression levels. rownames() are genes, and colnames() are samples.
#' @param removeLowVaryingGenes The proportion of low varying genes to be removed.The default is .2
#' @return A vector of row/genes to keep.
#' @export
doVariableSelection <- function(exprMat, removeLowVaryingGenes=.2)
{
  vars <- apply(exprMat, 1, var)
  return(order(vars, decreasing=TRUE)[seq(1:as.integer(nrow(exprMat)*(1-removeLowVaryingGenes)))])
}
#'This function takes two gene expression matrices (like trainExprMat and testExprMat) and returns homogenized versions of the matrices by employing the homogenization method specified.
#'By default, the Combat method from the sva library is used.
#'In both matrices, genes are row names and samples are column names.
#'It will deal with duplicated gene names, as it subsets and orders the matrices correctly.
#'@param testExprMat A gene expression matrix for samples on which we wish to predict a phenotype.Genes are rows, samples are columns.
#'@param trainExprMat A gene expression matrix for samples for which the phenotype is already known.Genes are rows, samples are columns.
#'@param batchCorrect The type of batch correction to be used. Options are 'eb' for Combat, 'none', or 'qn' for quantile normalization.
#'#The default is 'eb'.
#'@param selection This parameter can be used to specify how duplicates are handled. The default value of -1 means to ask the user.
#'#Other options include '1' to summarize duplicates by their mean, and '2'to discard all duplicated genes.
#'@param printOutput To suppress output, set to false. Default is TRUE.
#'@import sva
#'@import preprocessCore
#'@keywords Homogenize gene expression data.
#'@return A list containing two entries $train and $test, which are the homogenized input matrices.
#'@export
homogenizeData<-function (testExprMat, trainExprMat, batchCorrect = "eb", selection = -1, printOutput = TRUE)
{
  #Check the batchCorrect parameter
  if (!(batchCorrect %in% c("eb", "qn", "none",
                            "rank", "rank_then_eb", "standardize")))
    stop("\"batchCorrect\" must be one of \"eb\", \"qn\", \"rank\", \"rank_then_eb\", \"standardize\" or \"none\"")
  
  #Check if both row and column names have been specified.
  if (is.null(rownames(trainExprMat)) || is.null(rownames(testExprMat))) {
    stop("ERROR: Gene identifiers must be specified as \"rownames()\" on both training and test expression matrices. Both matices must have the same type of gene identifiers.")
  }
  
  #Check that some of the row names overlap between both datasets (print an error if none overlap)
  if (sum(rownames(trainExprMat) %in% rownames(testExprMat)) == 0) {
    stop("ERROR: The rownames() of the supplied expression matrices do not match. Note that these are case-sensitive.")
  }
  else {
    if (printOutput)
      message(paste("\n", sum(rownames(trainExprMat) %in% rownames(testExprMat)), " gene identifiers overlap between the supplied expression matrices... \n", paste = ""))
  }
  
  #If there are duplicate gene names, give the option of removing them or summarizing them by their mean.
  if ((sum(duplicated(rownames(trainExprMat))) > 0) || sum(sum(duplicated(rownames(testExprMat))) > 0)) {
    if (selection == -1) {
      message("\nExpression matrix contain duplicated gene identifiers (i.e. duplicate rownames()), how would you like to proceed:")
      message("\n1. Summarize duplicated gene ids by their mean value (acceptable in most cases)")
      message("\n2. Disguard all duplicated genes (recommended if unsure)")
      message("\n3. Abort (if you want to deal with duplicate genes ids manually)\n")
    }
    while (is.na(selection) | selection <= 0 | selection > 3) {
      selection <- readline("Selection: ")
      selection <- ifelse(grepl("[^1-3.]", selection), -1, as.numeric(selection))
    }
    
    message("\n")
    
    if (selection == 1) #Summarize duplicates by their mean.
    {
      if ((sum(duplicated(rownames(trainExprMat))) > 0)) {
        trainExprMat <- summarizeGenesByMean(trainExprMat)
      }
      if ((sum(duplicated(rownames(testExprMat))) > 0)) {
        testExprMat <- summarizeGenesByMean(testExprMat)
      }
    }
    else if (selection == 2) #Disguard all duplicated genes.
    {
      if ((sum(duplicated(rownames(trainExprMat))) > 0)) {
        keepGenes <- names(which(table(rownames(trainExprMat)) == 1))
        trainExprMat <- trainExprMat[keepGenes, ]
      }
      if ((sum(duplicated(rownames(testExprMat))) > 0)) {
        keepGenes <- names(which(table(rownames(testExprMat)) == 1))
        testExprMat <- testExprMat[keepGenes, ]
      }
    }
    else {
      stop("Exectution Aborted!")
    }
  }
  
  #Subset and order gene ids on the expression matrices.
  commonGenesIds <- rownames(trainExprMat)[rownames(trainExprMat) %in%
                                             rownames(testExprMat)]
  trainExprMat <- trainExprMat[commonGenesIds, ]
  testExprMat <- testExprMat[commonGenesIds, ]
  
  #Subset and order the 2 expression matrices.
  if (batchCorrect == "eb") {
    #Subset to common genes and batch correct using ComBat.
    dataMat <- cbind(trainExprMat, testExprMat)
    mod <- data.frame(`(Intercept)` = rep(1, ncol(dataMat)))
    rownames(mod) <- colnames(dataMat)
    whichbatch <- as.factor(c(rep("train", ncol(trainExprMat)),
                              rep("test", ncol(testExprMat))))
    
    # Added
    # Filter out genes with low variances to make sure comBat run correctly
    dataMat <- cbind(trainExprMat, testExprMat)
    gene_vars = apply(dataMat, 1, var)
    genes<-as.vector(gene_vars)
    
    if (length(which(genes <= 1e-3) != 0)){ #If some genes have low variances (if the variance is not 0), remove them.
      dataMat = dataMat[-(which(genes <= 1e-3)),]
    }
    
    # End added
    
    combatout <- ComBat(dataMat, whichbatch, mod = mod)
    return(list(train = combatout[, whichbatch == "train"],
                test = combatout[, whichbatch == "test"], selection = selection))
  }
  else if (batchCorrect == "standardize") #Standardize to mean 0 and variance 1 in each dataset using a non EB based approach.
  {
    for (i in 1:nrow(trainExprMat)) {
      row <- trainExprMat[i, ]
      trainExprMat[i, ] <- ((row - mean(row))/sd(row))
    }
    for (i in 1:nrow(testExprMat)) {
      row <- testExprMat[i, ]
      testExprMat[i, ] <- ((row - mean(row))/sd(row))
    }
    return(list(train = trainExprMat, test = testExprMat,
                selection = selection))
  }
  else if (batchCorrect == "rank") #The random-rank transform approach, that may be better when applying models to RNA-seq data.
  {
    for (i in 1:nrow(trainExprMat)) {
      trainExprMat[i, ] <- rank(trainExprMat[i, ], ties.method = "random")
    }
    for (i in 1:nrow(testExprMat)) {
      testExprMat[i, ] <- rank(testExprMat[i, ], ties.method = "random")
    }
    return(list(train = trainExprMat, test = testExprMat,
                selection = selection))
  }
  else if (batchCorrect == "rank_then_eb") #Rank-transform the RNAseq data, then apply ComBat
  {
    #First, rank transform the RNA-seq data.
    for (i in 1:nrow(testExprMat)) {
      testExprMat[i, ] <- rank(testExprMat[i, ], ties.method = "random")
    }
    #Subset to common genes and batch correct using ComBat.
    dataMat <- cbind(trainExprMat, testExprMat)
    mod <- data.frame(`(Intercept)` = rep(1, ncol(dataMat)))
    rownames(mod) <- colnames(dataMat)
    whichbatch <- as.factor(c(rep("train", ncol(trainExprMat)),
                              rep("test", ncol(testExprMat))))
    combatout <- ComBat(dataMat, whichbatch, mod = mod)
    return(list(train = combatout[, whichbatch == "train"],
                test = combatout[, whichbatch == "test"], selection = selection))
  }
  else if (batchCorrect == "qn")
  {
    dataMat <- cbind(trainExprMat, testExprMat)
    dataMatNorm <- normalize.quantiles(dataMat)
    whichbatch <- as.factor(c(rep("train", ncol(trainExprMat)),
                              rep("test", ncol(testExprMat))))
    return(list(train = dataMatNorm[, whichbatch == "train"],
                test = dataMatNorm[, whichbatch == "test"],
                selection = selection))
  }
  else {
    return(list(train = trainExprMat, test = testExprMat,
                selection = selection))
  }
}
#'This function takes a gene expression matrix and if duplicate genes are measured, summarizes them by their means.
#'@param exprMat A gene expression matrix with genes as rownames() and samples as colnames().
#'@return A gene expression matrix that does not contain duplicate genes.
#'@keywords Summarize duplicate genes by their mean.
#'@export
summarizeGenesByMean <- function(exprMat)
{
  geneIds <- rownames(exprMat)
  t <- table(geneIds) #How many times is each gene name duplicated.
  allNumDups <- unique(t)
  allNumDups <- allNumDups[-which(allNumDups == 1)]
  
  #Create a *new* gene expression matrix with everything in the correct order....
  #Starrt by just adding stuff that isn't duplicated
  exprMatUnique <- exprMat[which(geneIds %in% names(t[t == 1])), ]
  gnamesUnique <- geneIds[which(geneIds %in% names(t[t == 1]))]
  
  #Add all the duplicated genes to the bottom of "exprMatUniqueHuman", summarizing as you go
  for(numDups in allNumDups)
  {
    geneList <- names(which(t == numDups))
    
    for(i in 1:length(geneList))
    {
      exprMatUnique <- rbind(exprMatUnique, colMeans(exprMat[which(geneIds == geneList[i]), ]))
      gnamesUnique <- c(gnamesUnique, geneList[i])
      # print(i)
    }
  }
  
  if(class(exprMatUnique) == "numeric")
  {
    exprMatUnique <- matrix(exprMatUnique, ncol=1)
  }
  
  rownames(exprMatUnique) <- gnamesUnique
  return(exprMatUnique)
}
#'This function predicts a phenotype (drug sensitivity score) when provided with microarray or bulk RNAseq gene expression data of different platforms.
#'The imputations are performed using ridge regression, training on a gene expression matrix where phenotype is already known.
#'This function integrates training and testing datasets via a user-defined procedure, and power transforming the known phenotype.
#'@param trainingExprData The training data. A matrix of expression levels. rownames() are genes, colnames() are samples (cell line names or cosmic ides, etc.). rownames() must be specified and must contain the same type of gene ids as "testExprData"
#'@param trainingPtype The known phenotype for "trainingExprData". This data must be a matrix of drugs/rows x cell lines/columns or cosmic ids/columns. This matrix can contain NA values, that is ok (they are removed in the calcPhenotype() function).
#'@param testExprData The test data where the phenotype will be estimated. It is a matrix of expression levels, rows contain genes and columns contain samples, "rownames()" must be specified and must contain the same type of gene ids as "trainingExprData".
#'@param batchCorrect How should training and test data matrices be homogenized. Choices are "eb" (default) for ComBat, "qn" for quantiles normalization or "none" for no homogenization.
#'@param powerTransformPhenotype Should the phenotype be power transformed before we fit the regression model? Default to TRUE, set to FALSE if the phenotype is already known to be highly normal.
#'@param removeLowVaryingGenes What proportion of low varying genes should be removed? 20 percent be default
#'@param minNumSamples How many training and test samples are required. Print an error if below this threshold
#'@param selection How should duplicate gene ids be handled. Default is -1 which asks the user. 1 to summarize by their or 2 to disguard all duplicates.
#'@param printOutput Set to FALSE to supress output.
#'@param pcr Indicates whether or not you'd like to use pcr for feature (gene) reduction. Options are 'TRUE' and 'FALSE'. If you indicate 'report_pc=TRUE' you need to also indicate 'pcr=TRUE'
#'@param removeLowVaringGenesFrom Determine method to remove low varying genes. Options are 'homogenizeData' and 'rawData'.
#'@param report_pc Indicates whether you want to output the training principal components. Options are 'TRUE' and 'FALSE'.
#'@param cc Indicate if you want correlation coefficients for biomarker discovery.
#'@param percent Indicate percent variability (of the training data) you'd like principal components to reflect if pcr=TRUE. Default is 80 for 80%
#'These are the correlations between a given gene of interest across all samples vs. a given drug response across samples.
#'These correlations can be ranked to obtain a ranked correlation to determine highly correlated drug-gene associations.
#'@param rsq Indicate whether or not you want to output the R^2 values for the data you train on from true and predicted values.
#'These values represent the percentage in which the optimal model accounts for the variance in the training data.
#'Options are 'TRUE' and 'FALSE'.
#'@return .txt files will be saved into your working directory. Depending on the parameter specified, the .txt file outputs of this function can include the estimated phenotype/sensitivity predictions,
#'the R^2 data, and the correlation coefficients. Principal components are stored as .RData files for each drug in your drug dataset.
#'@import sva
#'@import ridge
#'@import car
#'@import ridge
#'@import glmnet
#'@import tidyverse
#'@import utils
#'@import stats
#'@importFrom pls explvar explvar
#'@keywords predict drug sensitivity and phenotype
#'@export
calcPhenotype<-function (trainingExprData,
                         trainingPtype,
                         testExprData,
                         batchCorrect,
                         powerTransformPhenotype=TRUE,
                         removeLowVaryingGenes=0.2,
                         minNumSamples,
                         selection=1,
                         printOutput,
                         pcr=FALSE,
                         removeLowVaringGenesFrom,
                         report_pc=FALSE,
                         cc=FALSE,
                         percent=80,
                         rsq=FALSE)
{
  
  #Initiate empty lists for each data type you'd like to collect.
  #_______________________________________________________________
  DrugPredictions<-list() #Collects drug predictions.
  rsqs<-list() #Collects R^2 values.
  cors<-list() #Collects correlation coefficient for each gene across all samples vs. each drug across all samples.
  pvalues<-list() #Collects p-values for the correlation coefficients. 
  
  #vs=c()
  
  drugs<-colnames(trainingPtype) #Store all the possible drugs in a vector.
  
  #Check the supplied data and parameters.
  #_______________________________________________________________
  if (class(testExprData)[1] != "matrix")
    stop("\nERROR: \"testExprData\" must be a matrix.")
  if (class(trainingExprData)[1] != "matrix")
    stop("\nERROR: \"trainingExprData\" must be a matrix.")
  if (class(trainingPtype)[1] != "matrix")
    stop("\nERROR: \"trainingPtype\" must be a matrix.")
  
  if (report_pc)
    if (pcr == FALSE)
      stop("\nERROR: pcr must be TRUE if report_pc is TRUE")
  
  if (pcr)
    if (cc)
      stop("\nERROR: pcr must be FALSE if cc is TRUE")
  
  #Make sure training samples are equivalent in both matrices.
  if (!any(colnames(trainingExprData) %in% rownames(trainingPtype)))
    stop("\nERROR: No Cell Lines Found in Common: Sample names must be consistent in training matrices")
  
  #Subset and order the training Expr and trainingPtype to the cell lines in common (and order them)
  
  commonCellLines<-colnames(trainingExprData)[colnames(trainingExprData) %in% rownames(trainingPtype)]
  
  trainingExprData <- trainingExprData[,commonCellLines]
  trainingPtype <- trainingPtype[commonCellLines,]
  
  #Check if an adequate number of training and test samples have been supplied.
  #_______________________________________________________________
  if ((nrow(trainingExprData) < minNumSamples) || (nrow(testExprData) < minNumSamples)) {
    stop(paste("\nThere are less than", minNumSamples, "samples in your test or training set. It is strongly recommended that you use larger numbers of samples in order to (a) correct for batch effects and (b) fit a reliable model. To supress this message, change the \"minNumSamples\" parameter to this function."))
  }
  
  #Get the homogenized data.
  #_______________________________________________________________
  homData <- homogenizeData(testExprMat=testExprData, trainExprMat=trainingExprData, batchCorrect, selection, printOutput)
  
  #Remove low varying genes.
  #_______________________________________________________________
  #Do variable selection if specified. By default, we remove 20% of least varying genes from the homogenized dataset.
  #We can also remove the intersection of the lowest 20% from both training and test sets (for the gene ids remaining in the homogenized data).
  #Otherwise, keep all genes.
  
  #Check batchCorrect parameter.
  if (!(removeLowVaringGenesFrom %in% c("homogenizeData", "rawData"))) {
    stop("\nremoveLowVaringGenesFrom\" must be one of \"homogenizeData\", \"rawData\"")
  }
  
  keepRows <- seq(1:nrow(homData$train)) #By default we will keep all the genes.
  if (removeLowVaryingGenes > 0 && removeLowVaryingGenes < 1) { #If the proportion of variability to keep is between 0 and 1.
    if (removeLowVaringGenesFrom == "homogenizeData") { #If you're filtering based on homogenized data.
      keepRows <- doVariableSelection(cbind(homData$test, homData$train), removeLowVaryingGenes = removeLowVaryingGenes)
      
      numberGenesRemoved <- nrow(homData$test) - length(keepRows)
      if (printOutput) message(paste("\n", numberGenesRemoved, "low variabilty genes filtered."));
    }
    else if (removeLowVaringGenesFrom == "rawData") { #If we are filtering based on the raw data i.e. the intersection of the things filtered from both datasets.
      evaluabeGenes <- rownames(homData$test)
      keepRowsTrain <- doVariableSelection(trainingExprData[evaluabeGenes,], removeLowVaryingGenes = removeLowVaryingGenes)
      keepRowsTest <- doVariableSelection(testExprData[evaluabeGenes,], removeLowVaryingGenes = removeLowVaryingGenes)
      keepRows <- intersect(keepRowsTrain, keepRowsTest)
      numberGenesRemoved <- nrow(homData$test) - length(keepRows)
      if (printOutput)
        message(paste("\n", numberGenesRemoved, "low variabilty genes filtered."));
    }
  }
  
  #Predict for each drug.
  #_______________________________________________________________
  for(a in 1:length(drugs)){ #For each drug...
    
    #Modify trainingPtype and trainingExprData so that you only use cell lines for which you have expression and response data for.
    #_______________________________________________________________
    trainingPtype2<-trainingPtype[,a] #Obtain the response data for the compound of interest.
    NonNAindex <- which(!is.na(trainingPtype2)) #Get the indices of the non NAs. You only want the cell lines/cosmic ids that you have drug info for.
    
    samps<-rownames(trainingPtype)[NonNAindex] #Obtain cell lines you have expression and response data for.
    
    if (length(samps) == 1){ #Make sure training data has more than just 1 sample. If the drug has one sample, it will be skipped.
      drugs = drugs[-a]
      message(paste("\n", drugs[a], "is skipped due to insufficient cell lines to fit the model."))
      next
    } else {
      
      trainingPtype4<-as.numeric(trainingPtype2[NonNAindex]) #This makes sure you use all the response data for the drug without NA values.
      
      #PowerTransform the phenotype if specified.
      #_______________________________________________________________
      offset = 0
      if (powerTransformPhenotype){
        if (min(trainingPtype4) < 0){ #All numbers must be positive for a powertransform to work, so make them positive.
          offset <- -min(trainingPtype4) + 1
          trainingPtype4 <- trainingPtype4 + offset
        }
        transForm <- powerTransform(trainingPtype4)[[6]]
        trainingPtype4 <- trainingPtype4^transForm
      }
      
      #Create the ridge regression model on the training data using pcr (keeping the number of components required for 80% variance).
      #_______________________________________________________________
      if (pcr){
        
        #There are many ways for pcr in R...here I will use pcr() (principal component regression).
        
        #Use pcr to predict for testing data.
        #_______________________________________________________________
        train_x<-(t(homData$train)[samps,keepRows]) #samps represent the cell lines that have been filtered, keepRows represents the genes.
        train_y<-trainingPtype4
        
        test_x<-(t(homData$test)[,keepRows])
        
        #Remove genes that you have no transcriptome data for (aka columns are filled with 0's).
        #If you don't, it causes problems in linearRidge().
        #This is weird, but there were 2 genes in CTRPv2 data that had 0's for one drug (I suppose this occurs depending on sample filtration).
        #Make sure you remove those same genes from the test data...want to make sure the genes are the same in both train_x and test_x
        x<-as.vector(colSums(train_x)) #Sum each column.
        bad<-which(x == 0) #Column/genes indices that are filled with only 0's.
        if(length(bad) != 0){
          train_x<-train_x[,-bad]
          test_x<-data.frame(test_x[,-bad])
        }
        #Remove genes with 0 variance. If you don't, this will also cause problems in linearRidge().
        variance<-c()
        for(i in 1:ncol(train_x)){
          variance[i]<-var(as.vector(train_x[,i]))
        }
        bi<-which(variance %in% 0) #Bad index...gene has 0 variance.
        if(length(bi) != 0){
          train_x<-train_x[,-bi]
          test_x<-data.frame(test_x[,-bi])
        }
        
        #Check to make sure you have more than 1 training sample for the drug's model.
        #_______________________________________________________________
        trainFrame<-try(data.frame(Resp=train_y, train_x), silent = TRUE)
        if (dim(trainFrame)[1] == 1){ #Make sure you have more than 1 sample.
          drugs = drugs[-a]
          message(paste("\n", drugs[a], "is skipped due to insufficient cell lines to fit the model."))
          next
        } else {
          
          pcr_model<-pcr(Resp~., data=trainFrame, validation='CV')
          
          v=cumsum(explvar(pcr_model)) #A vector of all the pcs and their percent of variance.
          ncomp=min(which(v > percent)) #Identify which pcs will represent 80% variance.
          
          #vs[a]<-ncomp
          
          if(printOutput) message("\nCalculating predicted phenotype using pcr...")
          preds<-predict(pcr_model, newdata=test_x, ncomp=ncomp)
          
          #You can compute an R^2 value for the data you train on from true and predicted values.
          #The rsq value represents the percentage in which the optimal model accounts for the variance in the training data.
          #_______________________________________________________________
          if (rsq){
            
            if (dim(train_x)[1] < 4){ #The code will result in an error if you have 3 samples (which is enough for the model fitting but not when you do a 70/30% split)...
              message(paste("\n", drugs[a], 'is skipped for R^2 analysis'))
            }else{
              
              #trainFrame<-data.frame(Resp=train_y, train_x)
              
              data<-(cbind(train_x, train_y)) #Rows are samples, columns are genes.
              
              dt<-sort(sample(nrow(data), nrow(data)*.7)) #sample() randomly picks 70% of rows/samples from the dataset. It samples without replacement.
              
              #Prepare the training data (70% of original training data)
              train_x<-data[dt,]
              ncol<-dim(train_x)[2]
              train_y<-train_x[,ncol]
              train_x<-train_x[,-ncol]
              
              #Prepare the testing data (30% of original training data)
              test_x<-data[-dt,]
              ncol<-dim(test_x)[2]
              test_y<-test_x[,ncol]
              test_x<-test_x[,-ncol]
              
              #Remove genes that you have no transcriptome data for aka columns are filled with 0's.
              #If you don't, it causes problems in linearRidge()
              #Make sure you remove these same genes from the test data...want to make sure the genes are the same in both.
              x<-as.vector(colSums(train_x))
              bad<-which(x == 0)
              if (length(bad) != 0){
                train_x<-train_x[,-bad]
                test_x<-data.frame(test_x[,-bad])
              }
              #Remove genes with 0 variance. If you don't, this will also cause problems in linearRidge()
              variance<-c()
              for(i in 1:ncol(train_x)){
                variance[i]<-var(as.vector(train_x[,i]))
              }
              bi<-which(variance %in% 0) #Bad index...gene has 0 variance.
              if (length(bi) != 0){
                train_x<-train_x[,-bi]
              }
              
              data<-data.frame(Resp=train_y, train_x)
              
              pcr_model<-pcr(Resp~., data=data, validation='CV') #Set validation argument to CV.
              
              v=cumsum(explvar(pcr_model)) #A vector of all the pcs and their percent of variance.
              ncomp=min(which(v > percent)) #Identify which pcs will represent 80% variance.
              
              pcr_pred<-predict(pcr_model, test_x, ncomp=ncomp)
              
              if (printOutput) message("\nCalculating R^2...")
              sst<-sum((test_y - mean(test_y))^2) #Compute the sum of squares total.
              sse<-sum((pcr_pred - test_y)^2) #Compute the sum of squares error.
              rsq_value<-1 - sse/sst #Compute the rsq value.
            }
          }
        }
        
        if (report_pc){
          if (printOutput) message("\nObtaining principal components...")
          pcs<-coef(pcr_model, comps = ncomp) #comps: numeric, which components to return.
          dir.create("./calcPhenotype_Output")
          path<-paste('./calcPhenotype_Output/', drugs[a], '.RData', sep="")
          save(pcs, file=path)
        }
        
      } else {
        
        #Create the ridge regression model on our training data to predict for our actual testing data without pcr.
        #_______________________________________________________________
        if(printOutput) message("\nFitting Ridge Regression model...");
        
        expression<-(t(homData$train)[samps,keepRows]) #samps represent the cell lines that have been filtered, keepRows represents the genes.
        
        test_x<-(t(homData$test)[,keepRows])
        
        #Remove genes that you have no transcriptome data for (aka columns are filled with 0's).
        #If you don't, it causes problems in linearRidge().
        #This is weird, but there were 2 genes in CTRPv2 data that had 0's for one drug (I suppose this occurs depending on sample filtration).
        #Make sure you remove those same genes from the test data...want to make sure the genes are the same in both train_x and test_x
        x<-as.vector(colSums(expression))
        bad<-which(x == 0) #Column/genes indices that are filled with only 0's.
        if (length(bad) != 0){
          expression<-expression[,-bad]
          test_x<-data.frame(test_x[,-bad])
        }
        
        #Remove genes with 0 variance. If you don't, this will also cause problems in linearRidge().
        variance<-c()
        for (i in 1:ncol(expression)){
          variance[i]<-var(as.vector(expression[,i]))
        }
        bi<-which(variance %in% 0) #Bad index...gene has 0 variance.
        if (length(bi) != 0){ #If there are actually bad indices/genes that had 0 variance across all samples, remove them.
          expression<-expression[,-bi]
          test_x<-data.frame(test_x[,-bi])
        }
        
        #Check to make sure you have more than 1 training sample for the drug's model.
        #_______________________________________________________________
        trainFrame<-try(data.frame(Resp=trainingPtype4, expression), silent = TRUE)
        if (dim(trainFrame)[1] == 1){ #Make sure you have more than 1 sample.
          drugs = drugs[-a]
          message(paste("\n", drugs[a], "is skipped due to insufficient cell lines to fit the model."))
          next
        } else {
          
          if(printOutput) message("\nCalculating predicted phenotype...")
          
          rrModel<-linearRidge(Resp ~., data=trainFrame)
          
          preds<-predict(rrModel, newdata=data.frame(test_x))
          
          #You can compute an R^2 value for the data you train on from true and predicted values.
          #The R^2 value represents the percentage in which the optimal model accounts for the variance in the training data.
          #_______________________________________________________________
          if(rsq){
            
            if (dim(expression)[1] < 4){ #The code will result in an error if you have 3 samples...this makes sure you have more than 3 samples.
              #It results in an error because there is a 70/30% split of training data.
              message(paste("\n", drugs[a], 'is skipped for R^2 analysis'))
            } else {
              expression<-(cbind(expression, trainingPtype4))
              dt<-sort(sample(nrow(expression), nrow(expression)*.7)) #sample() randomly picks 70% of rows/samples from the dataset. It samples without replacement.
              
              #Prepare the training data (70% of original training data)
              train_x<-expression[dt,]
              ncol<-dim(train_x)[2]
              train_y<-train_x[,ncol]
              train_x<-train_x[,-ncol]
              
              #Prepare the testing data (30% of original training data)
              test_x<-expression[-dt,]
              ncol<-dim(test_x)[2]
              test_y<-test_x[,ncol]
              test_x<-test_x[,-ncol]
              
              #Remove genes that you have no transcriptome data for (aka columns are filled with 0's).
              #If you don't, it causes problems in linearRidge().
              #Make sure you remove those same genes from the test data...want to make sure the genes are the same in both train_x and test_x
              x<-as.vector(colSums(train_x))
              bad<-which(x == 0) #Column/genes indices that are filled with only 0's.
              if (length(bad) != 0){
                train_x<-train_x[,-bad]
                test_x<-data.frame(test_x[,-bad])
              }
              
              #Remove genes with 0 variance. If you don't, this will also cause problems in linearRidge().
              variance<-c()
              for (i in 1:ncol(train_x)){
                variance[i]<-var(as.vector(train_x[,i]))
              }
              bi<-which(variance %in% 0) #Bad index...gene has 0 variance.
              if (length(bi) != 0){ #If there are actually bad indices/genes that had 0 variance across all samples, remove them.
                train_x<-train_x[,-bi]
                test_x<-data.frame(test_x[,-bi])
              }
              
              trainFrame<-data.frame(Resp=train_y, train_x)
              rrModel<-linearRidge(Resp ~., data=trainFrame)
              
              testFrame<-data.frame(test_x)
              
              pred<-predict(rrModel, newdata=testFrame)
              
              if(printOutput) message("\nCalculating R^2...")
              
              sst<-sum((test_y - mean(test_y))^2) #Compute the sum of squares total.
              sse<-sum((pred - test_y)^2) #Compute the sum of squares error.
              rsq_value<-1 - sse/sst #Compute the rsq value.
              
            }
          }
        }
      }
      
      #If the response variable was transformed (aka powerTransformPhenotype=TRUE), untransform it here.
      #_______________________________________________________________
      if(powerTransformPhenotype) {
        preds <- preds^(1/transForm)
        preds <- preds - offset
      }
      
      #Find correlation between imputed response and expression of a gene.
      #Each gene is corrected differently; therefore, it may not be ideal to determine weight or percentage or weight in which each gene contributes to the prediction.
      #_______________________________________________________________
      if(cc){ #You can only collect correlations if pcr=FALSE!
        if(pcr){
          stop('ERROR: pcr must equal FALSE in order to compute correlations') #It doesn't make sense to compute correlations when the features have changed from genes to pcs.
        }
        
        if(printOutput) message("\nCalculating correlation coefficients...") #This is only relevant if you aren't using pcr.
        
        cors_vec<-c() #This vector will store the correlation coefficients. 
        cors_vec2<-c() #This vector will store the p values.
        matrix<-homData$test[keepRows,] #Matrix of genes x cell lines/cosmic ids.
        for(d in 1:nrow(matrix)){ #For each gene...
          cors_vec[d]<-cor.test(as.vector(matrix[d,]), as.vector(preds))$estimate #Compute correlation coefficient for expression of a given gene across patients vs.imputed values for a given drug across patients.
          cors_vec2[d]<-cor.test(as.vector(matrix[d,]), as.vector(preds))$p.value
        }
        #indices<-order(cor, decreasing=TRUE)
        #ordered_cor<-cor[indices] #Order the correlation coefficients from big to small.
        #genes<-rownames(matrix)
        #ordered_genes<-genes[indices]
        #names(ordered_cor)<-ordered_genes
      }
      
      if(printOutput) message(paste("\nDone making prediction for drug", a, "of", ncol(trainingPtype)))
      
      #Store the data in your lists.
      #_______________________________________________________________
      DrugPredictions[[a]]<-preds
      
      if(rsq){
        rsqs[[a]]<-rsq_value
      }
      
      if(cc){
        cors[[a]]<-cors_vec
        pvalues[[a]]<-cors_vec2
      }
    }
  }
  
  #Time to save the data!
  #_______________________________________________________________
  #Save drug prediction data to your home directory as a .txt file.
  names(DrugPredictions)<-drugs
  DrugPredictions_mat<-do.call(cbind, DrugPredictions)
  colnames(DrugPredictions_mat)<-drugs
  rownames(DrugPredictions_mat)<-colnames(testExprData)
  dir.create("./calcPhenotype_Output")
  write.csv(DrugPredictions_mat, file="./calcPhenotype_Output/DrugPredictions.csv", row.names = TRUE, col.names = TRUE)
  
  #If rsq=TRUE, save R^2 data.
  if(rsq){
    names(rsqs)<-drugs
    rsqs_mat<-do.call(cbind, rsqs)
    dir.create("./calcPhenotype_Output")
    write.table(rsqs_mat, file="./calcPhenotype_Output/R^2.txt")
  }
  
  #If CC=TRUE, save correlation coefficient data.
  if(cc){
    names(cors)<-drugs
    cor_mat<-do.call(cbind, cors)
    rownames(cor_mat)<-rownames(homData$train[keepRows,NonNAindex])
    colnames(cor_mat)<-drugs
    dir.create("./calcPhenotype_Output")
    write.table(cor_mat, file="./calcPhenotype_Output/cors.txt")
    
    names(pvalues)<-drugs
    p_mat<-do.call(cbind, pvalues)
    rownames(p_mat)<-rownames(homData$train[keepRows, NonNAindex])
    colnames(p_mat)<-drugs
    dir.create("./calcPhenotype_Output")
    write.table(p_mat, file="./calcPhenotype_Output/pvalues.txt")
  }
  
  #print(vs)
}