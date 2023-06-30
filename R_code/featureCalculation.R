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

args = commandArgs(trailingOnly=TRUE)

indexNum = as.numeric(args[[1]])
indexMult = as.numeric(args[[2]])
indexAdd = as.numeric(args[[3]])

arrayNumStart = indexAdd + (indexNum - 1) * indexMult + 1
arrayNumEnd = indexAdd + indexNum * indexMult

tic()

litFunctions <- read.table("data/availableFunctions/literatureBiSourceAll.txt", header = FALSE, sep = " ", fill = TRUE)[[1]]
newFunctions <- read.table("data/availableFunctions/disturbanceFunctions.txt", header = FALSE, sep = " ", fill = TRUE)[[1]]

litFunctions <- litFunctions[str_which(litFunctions, "SOLAR", negate = TRUE)]

functions <- c(litFunctions, newFunctions)

functionNames  <- functions[arrayNumStart:arrayNumEnd]
functionNames <- functionNames[!is.na(functionNames)]

sourceCpp('cpp_code/rFeatureAnalysis.cpp')



# Want to do something here on repeat, so get a for loop going
for(index in 1:length(functionNames)){
  disp("")
  # Get the function name
  functionName <- functionNames[[index]]
  disp(functionName)
  # Get the range of the function, save time when initialising COCO and Disturbance Based functions
  # Make sure we extract the seed to keep things consistent
  bits <- str_split(functionName, "-")[[1]]
  seedChar <- str_which(bits, "seed", negate = FALSE)
  if(length(seedChar) == 0){
    seed <- 1
  }else{
    seed <- as.numeric(substring(bits[seedChar], 5))
  }
  
  disp("Get range:")
  range <- functionMinMax(functionName, seed)
  disp(range)
  
  fmin <- range[[1]]
  fmax <- range[[2]]
  
  
  testFunctionHigh <- function(x){
    return(sampleFunction(functionName, x, 0, TRUE, fmin, fmax))
  }
  
  testFunctionLow <- function(x){
    return(sampleFunction(functionName, x, 1, TRUE, fmin, fmax))
  }
  
  testFunctionMid <- function(x){
    return(sampleFunction(functionName, x, 2, TRUE, fmin, fmax))
  }
  
  dim <- functionDimension(functionName, TRUE, fmin, fmax)
  lowerBound <- functionLowerBound(functionName, TRUE, fmin, fmax)
  upperBound <- functionUpperBound(functionName, TRUE, fmin, fmax)
  
  # disp(paste0("Function ", functionName, " of dimension ", dim))
  disp("Dimension: ")
  disp(dim)
  disp("Lower bound:")
  disp(lowerBound)
  disp("Upper bound:")
  disp(upperBound)
  
  disp("Get sample:")
  X <- read.table(paste0("data/samplePlans/LHS-dim", dim, "-n", 1000*dim,".txt"))
  # Need to unscale it (i.e. take it away from the [0,1]^d hypercube)
  for(i in 1:dim){
    X[i] <- lowerBound[[i]] + X[i] * (upperBound[[i]] - lowerBound[[1]])
  }
  toc()
  
  disp("Basic features:")
  X <- data.frame(t(X))
  basicFeatures <- functionBasicFeaturesWithSample(functionName, X, TRUE, fmin, fmax)
  X <- data.frame(t(X))
  # basicFeatures <- functionBasicFeatures(functionName, seed, TRUE, fmin, fmax)
  disp(basicFeatures)

  suffixes <- c("_0_1", "_0_2", "_0_3", "_0_4", "_0_5", "_0_6", "_0_7", "_0_8", "_0_9", "_0_95", "_0_975", "_mean", "_sd", "_coeff")
  # basicFeaturesNames <- c("instances",
  #                         "feature_CC",
  #                         "feature_C",
  #                         "feature_RRMSE",
  #                         paste0("feature_LCC", suffixes),
  #                         paste0("feature_LC", suffixes))
  basicFeaturesNames <- c("instances",
                          "feature_dimension",
                          "feature_CC",
                          "feature_RRMSE",
                          paste0("feature_LCC", suffixes))

  initialVals <- c(functionName, dim, basicFeatures)

  if(index == 1){
    allFeatures <- initialVals
    names(allFeatures) <- basicFeaturesNames
  }else{
    allFeatures[index, basicFeaturesNames] <- initialVals
  }
  write.table(allFeatures, paste0("data/clusterResults/featureRun_arrayJob", arrayNumStart, "-", arrayNumEnd, ".txt"), quote = FALSE, row.names = FALSE)
  toc()
  
  # disp("Create feature sample")
  # ctrl = list(init_sample.type = "lhs",
  #             init_sample.lower = lowerBound,
  #             init_sample.upper = upperBound)
  # 
  # X = createInitialSample(n.obs = 1000 * dim, dim = dim, control = ctrl)
  # toc()
  
  featobjectHigh = createFeatureObject(X = X, fun = testFunctionHigh)
  featobjectLow = createFeatureObject(X = X, fun = testFunctionLow)
  featobjectMid = createFeatureObject(X = X, fun = testFunctionMid)
  
  featureSets <- c("ela_meta", "ic", "basic", "disp", "pca", "nbc")
  
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

        if(index == 1){
          names(localFeatures) <- localFeaturesNames
          allFeatures <- c(allFeatures, localFeatures)
        }else{
          allFeatures[index, localFeaturesNames] <- localFeatures
        }
        write.table(allFeatures, paste0("data/clusterResults/featureRun_arrayJob", arrayNumStart, "-", arrayNumEnd, ".txt"), quote = FALSE, row.names = FALSE)
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

        if(index == 1){
          names(localFeatures) <- localFeaturesNames
          allFeatures <- c(allFeatures, localFeatures)
        }else{
          allFeatures[index, localFeaturesNames] <- localFeatures
        }

        write.table(allFeatures, paste0("data/clusterResults/featureRun_arrayJob", arrayNumStart, "-", arrayNumEnd, ".txt"), quote = FALSE, row.names = FALSE)
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

        if(index == 1){
          names(localFeatures) <- localFeaturesNames
          allFeatures <- c(allFeatures, localFeatures)
        }else{
          allFeatures[index, localFeaturesNames] <- localFeatures
        }
        write.table(allFeatures, paste0("data/clusterResults/featureRun_arrayJob", arrayNumStart, "-", arrayNumEnd, ".txt"), quote = FALSE, row.names = FALSE)
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
}

disp("Success!")

toc()
