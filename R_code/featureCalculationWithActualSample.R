options(warn = -1)
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
source("R_code/customFeatureCalculation.R")


# RUN THIS FROM TERMINAL AS FOLLOWS:
# Rscript R_code/featureCalculationWithActualSample.R 1 10 0

# Select functions for which sample features will be calculated
functions <- read.table("data/availableFunctions/chosenTestSuiteN207.txt", header = FALSE, sep = " ", fill = TRUE)[[1]] 
# functions <- read.table("data/availableFunctions/literatureBiSourceDim5.txt", header = FALSE, sep = " ", fill = TRUE)[[1]] 
# functions <- functions[str_which(functions, "SOLAR")]
# Select sizes for low and high fidelities
highFiSizes <- seq(2, 20, by = 2)
lowFiSizes <- seq(4, 20, by = 4)
relativeSize <- TRUE
seeds <- 1:40
# NOTHING SHOULD BE CHANGED PAST THIS POINT



args = commandArgs(trailingOnly=TRUE)

indexNum = as.numeric(args[[1]])
indexMult = as.numeric(args[[2]])
indexAdd = as.numeric(args[[3]])

arrayNumStart = indexAdd + (indexNum - 1) * indexMult + 1
arrayNumEnd = indexAdd + indexNum * indexMult

disp(paste0("Running from index ", arrayNumStart, " to index ", arrayNumEnd))

tic()


sourceCpp('cpp_code/rFeatureAnalysis.cpp')

functionNames  <- functions[arrayNumStart:arrayNumEnd]
functionNames <- functionNames[!is.na(functionNames)]

# Want to do something here on repeat, so get a for loop going
currIndex <- 0
for(index in 1:length(functionNames)){
  disp("")
  functionName <- functionNames[[index]]
  bits <- str_split(functionName, "-")[[1]]
  seedChar <- str_which(bits, "seed", negate = FALSE)
  if(length(seedChar) == 0){
    seed <- 1
  }else{
    seed <- as.numeric(substring(bits[seedChar], 5))
  }
  disp(functionName)
  disp("Get range:")
  # If dealing with SOLAR instance, finding these will take forever and
  # never be actually used, so skip
  if(substr(functionName, 1, 5) == "SOLAR"){
    disp("Actually skipping SOLAR range (takes too long and is not used)")
    fmin = 0
    fmax = 0
  }else{
    range <- functionMinMax(functionName, seed)
    fmin <- range[[1]]
    fmax <- range[[2]]
    disp(range)
  }
  
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
  
  disp("Dimension: ")
  disp(dim)
  disp("Lower bound:")
  disp(lowerBound)
  disp("Upper bound:")
  disp(upperBound)
  
  for(highFiSize in highFiSizes){
    for(lowFiSize in lowFiSizes){
      if(lowFiSize < highFiSize){next}
      for(seed in seeds){
        if(relativeSize){
          highFiBudget <- dim * highFiSize
          lowFiBudget <- dim * lowFiSize
        }else{
          highFiBudget <- highFiSize
          lowFiBudget <- lowFiSize
        }
        currIndex <- currIndex + 1
        disp(paste0("Working on seed ", seed, " with high fi budget ", highFiBudget, " and low fi budget ", lowFiBudget))
        # At this stage now need the sample locations
        samples <- functionSample(functionName, seed, highFiBudget, lowFiBudget, TRUE, fmin, fmax)
        lowFiSamples <- samples[[1]]
        highFiSamples <- samples[[2]]
        
        disp("Have obtained samples. Start with basic features: ")
        basicFeatures <- functionBasicFeaturesWithSample(functionName, highFiSamples, TRUE, fmin, fmax)
        disp(basicFeatures)
        
        suffixes <- c("_0_1", "_0_2", "_0_3", "_0_4", "_0_5", "_0_6", "_0_7", "_0_8", "_0_9", "_0_95", "_0_975", "_mean", "_sd", "_coeff")
        basicFeaturesNames <- c("instances",
                                "feature_sample_highFiBudget",
                                "feature_sample_lowFiBudget",
                                "feature_sample_highFiBudgetRatio",
                                "feature_sample_lowFiBudgetRatio",
                                "feature_sample_budgetRatio",
                                "feature_sample_CC",
                                "feature_sample_RRMSE", 
                                paste0("feature_sample_LCC", suffixes))
        
        initialVals <- c(paste0("(", functionName, ",", highFiBudget, ",", lowFiBudget, ",", seed, ')'), 
                         as.numeric(highFiBudget),
                         as.numeric(lowFiBudget),
                         as.numeric(highFiBudget) / as.numeric(dim),
                         as.numeric(lowFiBudget) / as.numeric(dim),
                         as.numeric(highFiBudget) /as.numeric(lowFiBudget),
                         basicFeatures)
        
        if(currIndex == 1){
          allFeatures <- initialVals
          names(allFeatures) <- basicFeaturesNames
        }else{
          allFeatures[currIndex, basicFeaturesNames] <- initialVals
        }
        if(arrayNumStart == arrayNumEnd){write.table(allFeatures, paste0("data/clusterResults/sampleFeatureRun_arrayJob", arrayNumStart, ".txt"), quote = FALSE, row.names = FALSE)}
        else{write.table(allFeatures, paste0("data/clusterResults/sampleFeatureRun_arrayJob", arrayNumStart, "-", arrayNumEnd, ".txt"), quote = FALSE, row.names = FALSE)}
        toc()
        
        disp("Create feature sample")
        highFiSamples <- do.call(rbind, highFiSamples)
        lowFiSamples <- do.call(rbind, lowFiSamples)
        
        featobjectHigh = createFeatureObject(X = highFiSamples, fun = testFunctionHigh)
        featobjectLow = createFeatureObject(X = lowFiSamples, fun = testFunctionLow)
        featobjectMid = createFeatureObject(X = highFiSamples, fun = testFunctionMid)
        
        # Reducing the number of groups for CV in ela_level for a higher probability that smaller samples will not be NA
        control <- list(ela.level.resample_iterations = 5)
        featureSets <- c("ela_distr", "ela_level", "ela_meta", "ic", "basic", "disp", "nbc", "pca")
        
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
              localFeatures <- calculateFeatureSet(feat.object = featobjectHigh, set = set, control = control)
            }
            
            localFeaturesNames <- names(localFeatures)
            localFeaturesNames <- gsub("\\.", "_", localFeaturesNames)
            localFeaturesNames <- gsub("-", "_", localFeaturesNames)
            localFeaturesNames <- paste0('feature_sample_high_', localFeaturesNames)
            
            if(currIndex == 1){
              names(localFeatures) <- localFeaturesNames
              allFeatures <- c(allFeatures, localFeatures)
            }else{
              allFeatures[currIndex, localFeaturesNames] <- localFeatures
            }
            
            if(arrayNumStart == arrayNumEnd){write.table(allFeatures, paste0("data/clusterResults/sampleFeatureRun_arrayJob", arrayNumStart, ".txt"), quote = FALSE, row.names = FALSE)}
            else{write.table(allFeatures, paste0("data/clusterResults/sampleFeatureRun_arrayJob", arrayNumStart, "-", arrayNumEnd, ".txt"), quote = FALSE, row.names = FALSE)}
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
              localFeatures <- calculateFeatureSet(feat.object = featobjectLow, set = set, control = control)
            }
            localFeaturesNames <- names(localFeatures)
            localFeaturesNames <- gsub("\\.", "_", localFeaturesNames)
            localFeaturesNames <- gsub("-", "_", localFeaturesNames)
            localFeaturesNames <- paste0('feature_sample_low_', localFeaturesNames)
            
            if(currIndex == 1){
              names(localFeatures) <- localFeaturesNames
              allFeatures <- c(allFeatures, localFeatures)
            }else{
              allFeatures[currIndex, localFeaturesNames] <- localFeatures
            }
            if(arrayNumStart == arrayNumEnd){write.table(allFeatures, paste0("data/clusterResults/sampleFeatureRun_arrayJob", arrayNumStart, ".txt"), quote = FALSE, row.names = FALSE)}
            else{write.table(allFeatures, paste0("data/clusterResults/sampleFeatureRun_arrayJob", arrayNumStart, "-", arrayNumEnd, ".txt"), quote = FALSE, row.names = FALSE)}
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
              localFeatures <- calculateFeatureSet(feat.object = featobjectMid, set = set, control = control)
            }
            localFeaturesNames <- names(localFeatures)
            localFeaturesNames <- gsub("\\.", "_", localFeaturesNames)
            localFeaturesNames <- gsub("-", "_", localFeaturesNames)
            localFeaturesNames <- paste0('feature_sample_mid_', localFeaturesNames)
            
            if(currIndex == 1){
              names(localFeatures) <- localFeaturesNames
              allFeatures <- c(allFeatures, localFeatures)
            }else{
              allFeatures[currIndex, localFeaturesNames] <- localFeatures
            }
            if(arrayNumStart == arrayNumEnd){write.table(allFeatures, paste0("data/clusterResults/sampleFeatureRun_arrayJob", arrayNumStart, ".txt"), quote = FALSE, row.names = FALSE)}
            else{write.table(allFeatures, paste0("data/clusterResults/sampleFeatureRun_arrayJob", arrayNumStart, "-", arrayNumEnd, ".txt"), quote = FALSE, row.names = FALSE)}
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
        if(currIndex == 1){allFeatures <- as.data.frame(allFeatures)}
      }
    }
  }
}

disp("Success!")

toc()
