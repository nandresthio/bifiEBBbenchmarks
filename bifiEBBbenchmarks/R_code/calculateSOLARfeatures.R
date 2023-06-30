library(flacco)
library(lhs)
library(plyr)
library(numDeriv)
library(e1071)
library(Rcpp)
library(RANN)
library(mlr)
library(mda)
library(class)
library(pracma)
library(ParamHelpers)
library(stringr)
source("customFeatureCalculation.R")
tic()

args = commandArgs(trailingOnly=TRUE)
indexNum = as.numeric(args[[1]])

sourceCpp(paste0('../cpp_code/testSave.cpp'), rebuild = TRUE)

Y <- read.table("../solar/points/evaluatedPoints.txt", header = TRUE, sep = " ")




Yreplaced <- Y
Yreplaced[Yreplaced == 1e+20] <- 2500

Yremoved <- Y
for(fid in c("0.10", '0.20', '0.30', '0.40', '0.50', '0.60', '0.70', '0.80', '0.90', "1.00")){
  Yremoved <- Yremoved[Yremoved[paste0("fid", fid)] != 1e+20, ]
}

# X <- Yreplaced[, paste0("x", 1:5)]
# X$x1 <- (X$x1 - 793) / (995 - 793)
# X$x2 <- (X$x2 - 2) / (50 - 2)
# X$x3 <- (X$x3 - 2) / (30 - 2)
# X$x4 <- (X$x4 - 0.01) / (5 - 0.01)
# X$x5 <- (X$x5 - 0.01) / (5 - 0.01)
# write.table(X, 'solar/points/evaluatedPointsReplaced.txt', quote = FALSE, row.names = FALSE, col.names = FALSE)
# Y <- Yreplaced$fid1.00
# write.table(Y, 'solar/points/evaluatedPointsReplacedValues.txt', quote = FALSE, row.names = FALSE, col.names = FALSE)
# 
# X <- Yremoved[, paste0("x", 1:5)]
# X$x1 <- (X$x1 - 793) / (995 - 793)
# X$x2 <- (X$x2 - 2) / (50 - 2)
# X$x3 <- (X$x3 - 2) / (30 - 2)
# X$x4 <- (X$x4 - 0.01) / (5 - 0.01)
# X$x5 <- (X$x5 - 0.01) / (5 - 0.01)
# write.table(X, 'solar/points/evaluatedPointsRemoved.txt', quote = FALSE, row.names = FALSE, col.names = FALSE)
# Y <- Yremoved$fid1.00
# write.table(Y, 'solar/points/evaluatedPointsRemovedValues.txt', quote = FALSE, row.names = FALSE, col.names = FALSE)

index <- 0
for(type in c("replaced", "removed")){
  for(fid in c("0.10", '0.20', '0.30', '0.40', '0.50', '0.60', '0.70', '0.80', '0.90')){
    index <- index + 1
    if(type == "replaced"){instanceName <- paste0("SOLAR", fid)}
    else{instanceName <- paste0("SOLAR", fid, type)}
    if(type == "replaced"){
      Yhigh <- Yreplaced$fid1.00
      Ylow <- Yreplaced[, paste0("fid", fid)]
      Ymid <- Yreplaced[, paste0("fid", fid, "Diff")]
      X <- Yreplaced[, paste0("x", 1:5)]
    }else{
      Yhigh <- Yremoved$fid1.00
      Ylow <- Yremoved[, paste0("fid", fid)]
      Ymid <- Yremoved[, paste0("fid", fid, "Diff")]
      X <- Yremoved[, paste0("x", 1:5)]
    }
    
    if(index != indexNum){next}
  
    X$x1 <- (X$x1 - 793) / (995 - 793)
    X$x2 <- (X$x2 - 2) / (50 - 2)
    X$x3 <- (X$x3 - 2) / (30 - 2)
    X$x4 <- (X$x4 - 0.01) / (5 - 0.01)
    X$x5 <- (X$x5 - 0.01) / (5 - 0.01)
    
    X <- data.frame( t(X))
    disp("Basic features:")
    basicFeatures <- functionBasicFeaturesWithSampleAndVals(paste0("SOLAR", fid), X, Yhigh, Ylow)
    disp(basicFeatures)
    X <- data.frame( t(X))
    
    suffixes <- c("_0_1", "_0_2", "_0_3", "_0_4", "_0_5", "_0_6", "_0_7", "_0_8", "_0_9", "_0_95", "_0_975", "_mean", "_sd", "_coeff")
    basicFeaturesNames <- c("instances",
                            "feature_dimension",
                            "feature_CC",
                            "feature_C",
                            "feature_RRMSE", 
                            paste0("feature_LCC", suffixes),
                            paste0("feature_LC", suffixes))
    
    
    dim <- functionDimension(paste0("SOLAR", fid))
    
    initialVals <- c(instanceName, dim, basicFeatures)
    
    allFeatures <- initialVals
    names(allFeatures) <- basicFeaturesNames
    write.table(allFeatures, paste0("../data/clusterResults/featureRun", index, "_SOLAR.txt"), quote = FALSE, row.names = FALSE)
    toc()
    
    testFunctionHigh <- function(x){
      return(sampleFunction(paste0("SOLAR", fid, index), x, 0))
    }
    
    testFunctionLow <- function(x){
      return(sampleFunction(paste0("SOLAR", fid, index), x, 1))
    }
    
    testFunctionMid <- function(x){
      return(sampleFunction(paste0("SOLAR", fid, index), x, 2))
    }
    
    
    featobjectHigh = createFeatureObject(X = X, y = Yhigh, fun = testFunctionHigh)
    featobjectLow = createFeatureObject(X = X, y = Ylow, fun = testFunctionLow)
    featobjectMid = createFeatureObject(X = X, y = Ymid, fun = testFunctionMid)
    
    featureSets <- c("ela_conv", "ela_distr", "ela_level", "ela_meta", "ic", "basic", "disp", "pca", "nbc")
    
    for(set in featureSets){
      disp(paste0("Working on feature set ", set))
      out <- tryCatch({
        localFeatures <- c()
        # Special cases
        if(set == "pca" & max(featobjectHigh$env$init$y) == min(featobjectHigh$env$init$y)){
          disp("Performing special pca")
          localFeatures <- calculatePrincipalComponentFeaturesWithConstantY(featobjectHigh)
          
        }else if(set == "ela_level" & max(featobjectHigh$env$init$y) == min(featobjectHigh$env$init$y)){
          disp("Performing special ela_level")
          lcalFeatures <- calculateLevelsetFeaturesWithConstantY(featobjectHigh)
          
        }else if(set == "ela_distr"){
          disp("Performing special ela_distr")
          localFeatures <- calculateDistributionFeaturesCustom(featobjectHigh)
          
        }else if(set == "ela_meta"){
          disp("Performing special ela_meta")
          localFeatures <- calculateMetaModelFeaturesCustom(featobjectHigh)
          
        }else if(set == "disp"){
          disp("Performing special disp")
          localFeatures <- calculateDispersionFeaturesCustom(featobjectHigh)
          
        }else if(set == "nbc"){
          disp("Performing special nbc")
          localFeatures <- calculateNearestBetterFeaturesCustom(featobjectHigh)
          
        }else if(set == "ic"){
          disp("Performing special ic")
          localFeatures <- calculateInformationContentFeaturesCustom(featobjectHigh)
          
        }else{
          localFeatures <- calculateFeatureSet(feat.object = featobjectHigh, set = set)
        }
        
        localFeaturesNames <- names(localFeatures)
        localFeaturesNames <- gsub("\\.", "_", localFeaturesNames)
        localFeaturesNames <- gsub("-", "_", localFeaturesNames)
        localFeaturesNames <- paste0('feature_high_', localFeaturesNames)
        
        
        names(localFeatures) <- localFeaturesNames
        allFeatures <- c(allFeatures, localFeatures)
        
        write.table(allFeatures, paste0("../data/clusterResults/featureRun", index, "_SOLAR.txt"), quote = FALSE, row.names = FALSE)
        toc()
      },
      error=function(cond) {
        disp(paste0("Problem with high fi function and set ", set))
        print(cond)
        disp(featobjectHigh$env$init$y)
        return(NA)
      },
      finally={}
      )
      out <- tryCatch({
        localFeatures <- c()
        # Special cases
        if(set == "pca" & max(featobjectLow$env$init$y) == min(featobjectLow$env$init$y)){
          disp("Performing special pca")
          localFeatures <- calculatePrincipalComponentFeaturesWithConstantY(featobjectLow)
          
        }else if(set == "ela_level" & max(featobjectLow$env$init$y) == min(featobjectLow$env$init$y)){
          disp("Performing special ela_level")
          localFeatures <- calculateLevelsetFeaturesWithConstantY(featobjectLow)
          
        }else if(set == "ela_distr"){
          disp("Performing special ela_distr")
          localFeatures <- calculateDistributionFeaturesCustom(featobjectLow)
          
        }else if(set == "ela_meta"){
          disp("Performing special ela_meta")
          localFeatures <- calculateMetaModelFeaturesCustom(featobjectLow)
          
        }else if(set == "disp"){
          disp("Performing special disp")
          localFeatures <- calculateDispersionFeaturesCustom(featobjectLow)
          
        }else if(set == "nbc"){
          disp("Performing special nbc")
          localFeatures <- calculateNearestBetterFeaturesCustom(featobjectLow)
          
        }else if(set == "ic"){
          disp("Performing special ic")
          localFeatures <- calculateInformationContentFeaturesCustom(featobjectLow)
          
        }else{
          localFeatures <- calculateFeatureSet(feat.object = featobjectLow, set = set)
        }
        
        localFeaturesNames <- names(localFeatures)
        localFeaturesNames <- gsub("\\.", "_", localFeaturesNames)
        localFeaturesNames <- gsub("-", "_", localFeaturesNames)
        localFeaturesNames <- paste0('feature_low_', localFeaturesNames)
        
        names(localFeatures) <- localFeaturesNames
        allFeatures <- c(allFeatures, localFeatures)
        
        
        write.table(allFeatures, paste0("../data/clusterResults/featureRun", index, "_SOLAR.txt"), quote = FALSE, row.names = FALSE)
        toc()
      },
      error=function(cond) {
        disp(paste0("Problem with low fi function and set ", set))
        disp(featobjectLow$env$init$y)
        print(cond)
        return(NA)
      },
      finally={}
      )
      out <- tryCatch({
        localFeatures <- c()
        # Special cases
        if(set == "pca" & max(featobjectMid$env$init$y) == min(featobjectMid$env$init$y)){
          disp("Performing special pca")
          localFeatures <- calculatePrincipalComponentFeaturesWithConstantY(featobjectMid)
          
        }else if(set == "ela_level" & max(featobjectMid$env$init$y) == min(featobjectMid$env$init$y)){
          disp("Performing special ela_level")
          localFeatures <- calculateLevelsetFeaturesWithConstantY(featobjectMid)
          
        }else if(set == "ela_distr"){
          disp("Performing special ela_distr")
          localFeatures <- calculateDistributionFeaturesCustom(featobjectMid)
          
        }else if(set == "ela_meta"){
          disp("Performing special ela_meta")
          localFeatures <- calculateMetaModelFeaturesCustom(featobjectMid)
          
        }else if(set == "disp"){
          disp("Performing special disp")
          localFeatures <- calculateDispersionFeaturesCustom(featobjectMid)
          
        }else if(set == "nbc"){
          disp("Performing special nbc")
          localFeatures <- calculateNearestBetterFeaturesCustom(featobjectMid)
          
        }else if(set == "ic"){
          disp("Performing special ic")
          localFeatures <- calculateInformationContentFeaturesCustom(featobjectMid)
          
        }else{
          localFeatures <- calculateFeatureSet(feat.object = featobjectMid, set = set)
        }
        localFeaturesNames <- names(localFeatures)
        localFeaturesNames <- gsub("\\.", "_", localFeaturesNames)
        localFeaturesNames <- gsub("-", "_", localFeaturesNames)
        localFeaturesNames <- paste0('feature_mid_', localFeaturesNames)
        
        names(localFeatures) <- localFeaturesNames
        allFeatures <- c(allFeatures, localFeatures)
        
        write.table(allFeatures, paste0("../data/clusterResults/featureRun", index, "_SOLAR.txt"), quote = FALSE, row.names = FALSE)
        toc()
      },
      error=function(cond) {
        disp(paste0("Problem with mid function and set ", set))
        print(cond)
        disp(featobjectMid$env$init$y)
        return(NA)
      },
      finally={}
      )
    }
    
    # At this point should have all the features, if this is the first index, turn matrix into a data frame
    allFeatures <- as.data.frame(allFeatures)
    write.table(allFeatures, paste0("../data/clusterResults/featureRun", index, "_SOLAR.txt"), quote = FALSE, row.names = FALSE)
  }
  toc()
}

disp("Success!")

toc()
