######################################################################################

# Copyright 2023, Nicolau Andrés-Thió

# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:

# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.

# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.

######################################################################################

# Code which calculates both the real and the sample features
# of the SOLAR simulation engine.  New users should not run this code unless 
# new features need to be calculated. This code relies on the objective function
# values stored in the data/misc folder, to be taken from datashare.
# For information on the simulator, consult https://github.com/bbopt/solar


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
tic()

# RUN THIS FROM TERMINAL AS FOLLOWS:
# Rscript R_code/calculateSOLARfeatures.R 1 9 0

args = commandArgs(trailingOnly=TRUE)
indexNum = as.numeric(args[[1]])
indexMult = as.numeric(args[[2]])

indexStart = (indexNum - 1) * indexMult + 1
indexEnd = indexNum * indexMult

print(indexStart)
print(indexEnd)

if(indexStart == indexEnd){
  outputFilename <- paste0("data/clusterResults/featureRun", indexStart, "_SOLAR.txt")
}else{
  outputFilename <- paste0("data/clusterResults/featureRun", indexStart, "-", indexEnd, "_SOLAR.txt")
}

print(outputFilename)

sourceCpp('cpp_code/rFeatureAnalysis.cpp')

Y <- read.table("data/misc/solarEvaluatedPoints.txt", header = TRUE, sep = " ")

Yreplaced <- Y
Yreplaced[Yreplaced == 1e+20] <- 2500

Yremoved <- Y
for(fid in c("0.10", '0.20', '0.30', '0.40', '0.50', '0.60', '0.70', '0.80', '0.90', "1.00")){
  Yremoved <- Yremoved[Yremoved[paste0("fid", fid)] != 1e+20, ]
}

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
    
    if(!(index %in% indexStart:indexEnd)){next}
    # Scale input as internally the interface with the executable takes the domain to be in [0,1]^5
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
                          "feature_RRMSE",
                          paste0("feature_LCC", suffixes),
                          paste0("feature_LCCrel", suffixes))
    
    
    dim <- functionDimension(paste0("SOLAR", fid))
    
    initialVals <- c(instanceName, dim, basicFeatures)
    
    allFeatures <- initialVals
    names(allFeatures) <- basicFeaturesNames
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
    if(index != indexStart){
      existingFeatures <- read.table(outputFilename, header = TRUE, sep = " ")
      allFeatures <- rbind.fill(existingFeatures, allFeatures)
    }
    write.table(allFeatures, outputFilename, quote = FALSE, row.names = FALSE)
  }
  toc()
}
# Skip the next 2 so that the first 20 are real feautres, and the rest are sample features
index <- index + 2


# Move onto sample features
highFiSizes <- seq(2, 20, by = 2)
lowFiSizes <- seq(4, 20, by = 4)
seeds <- 1:40
for(seed in seeds){
  for(fid in c("0.90", '0.80', '0.70', '0.60', '0.50', '0.40', '0.30', '0.20', '0.10')){
    for(lowFiSize in lowFiSizes){
      for(highFiSize in highFiSizes){
        if(lowFiSize < highFiSize){next}
        highFiBudget <- 5 * highFiSize
        lowFiBudget <- 5 * lowFiSize
        index <- index + 1
        if(!(index %in% indexStart:indexEnd)){next}
        # Need an instance name to pass to c++ so that the point file does not get
        # written by more than one cluster run
        functionName <- paste0("SOLAR", fid, index)
        instanceRealName <- paste0("(SOLAR", fid, ",", highFiBudget, ",", lowFiBudget, ",", seed,  ")")
        
        disp(paste0("Working on index ", index, " with seed ", seed, " for fidelity ", fid, "with high fi budget ", highFiBudget, " and low fi budget ", lowFiBudget))
        
        fmin = 0
        fmax = 0
        
        # testFunctionHigh <- function(x){
        #   return(sampleFunction(functionName, x, 0, TRUE, fmin, fmax))
        # }
        # 
        # testFunctionLow <- function(x){
        #   return(sampleFunction(functionName, x, 1, TRUE, fmin, fmax))
        # }
        # 
        # testFunctionMid <- function(x){
        #   return(sampleFunction(functionName, x, 2, TRUE, fmin, fmax))
        # }
        
        dim <- functionDimension(functionName, TRUE, fmin, fmax)
        lowerBound <- functionLowerBound(functionName, TRUE, fmin, fmax)
        upperBound <- functionUpperBound(functionName, TRUE, fmin, fmax)
        disp("Dimension: ")
        disp(dim)
        disp("Lower bound:")
        disp(lowerBound)
        disp("Upper bound:")
        disp(upperBound)
        
        # At this stage now need the sample locations
        samples <- functionSample(functionName, seed, highFiBudget, lowFiBudget, TRUE, fmin, fmax)
        lowFiSamples <- samples[[1]]
        highFiSamples <- samples[[2]]
        
        # Ok, going to add something here to read in the objective function values instead
        values <- read.table(paste0("data/misc/SOLARn",lowFiBudget, "s", seed, ".txt"), header = TRUE, sep = " ")
        # The low fidelity values should be all of them
        Ylow <- values[, paste0("fid", fid)]
        # The high fi values need to be found
        Yhigh <- c()
        Ymid <- c()
        YlowLimited <- c()
        for(sample in highFiSamples){
          # Find the right point
          point <- values[abs(values$x1 -sample[[1]]) < 0.000001 &
                            abs(values$x2 - sample[[2]]) < 0.000001 &
                            abs(values$x3 - sample[[3]]) < 0.000001 &
                            abs(values$x4 - sample[[4]]) < 0.000001 &
                            abs(values$x5 - sample[[5]]) < 0.000001,]
          if(nrow(point) != 1){
            stop("Weird, did not find exactly one corresponding point in the file!")
          }
          Yhigh <- c(Yhigh, point$fid1.00)
          Ymid <- c(Ymid, point$fid1.00 - point[, paste0("fid", fid)])
          YlowLimited <- c(YlowLimited, point[, paste0("fid", fid)])
        }
        
        disp("Have obtained samples. Start with basic features: ")
        # basicFeatures <- functionBasicFeaturesWithSample(functionName, highFiSamples, TRUE, fmin, fmax)
        basicFeatures <- functionBasicFeaturesWithSampleAndVals(paste0("SOLAR", fid), highFiSamples, Yhigh, YlowLimited)
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
                                paste0("feature_sample_LCC", suffixes),
                                paste0("feature_sample_LCCrel", suffixes))
        
        initialVals <- c(instanceRealName, 
                         as.numeric(highFiBudget),
                         as.numeric(lowFiBudget),
                         as.numeric(highFiBudget) / as.numeric(dim),
                         as.numeric(lowFiBudget) / as.numeric(dim),
                         as.numeric(highFiBudget) /as.numeric(lowFiBudget),
                         basicFeatures)
        
        allFeatures <- initialVals
        names(allFeatures) <- basicFeaturesNames
        toc()
        
        disp("Create feature sample")
        # Do not specify the function as it should not be needed
        highFiSamples <- do.call(rbind, highFiSamples)
        lowFiSamples <- do.call(rbind, lowFiSamples)
        featobjectHigh = createFeatureObject(X = highFiSamples, y = Yhigh)
        featobjectLow = createFeatureObject(X = lowFiSamples, y = Ylow)
        featobjectMid = createFeatureObject(X = highFiSamples, y = Ymid)
        
        # disp("Create feature sample")
        # highFiSamples <- do.call(rbind, highFiSamples)
        # lowFiSamples <- do.call(rbind, lowFiSamples)
        # featobjectHigh = createFeatureObject(X = highFiSamples, fun = testFunctionHigh)
        # featobjectLow = createFeatureObject(X = lowFiSamples, fun = testFunctionLow)
        # featobjectMid = createFeatureObject(X = highFiSamples, fun = testFunctionMid)
        
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
            
            names(localFeatures) <- localFeaturesNames
            allFeatures <- c(allFeatures, localFeatures)
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
            
            names(localFeatures) <- localFeaturesNames
            allFeatures <- c(allFeatures, localFeatures)
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
            
            names(localFeatures) <- localFeaturesNames
            allFeatures <- c(allFeatures, localFeatures)
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
        allFeatures <- as.data.frame(allFeatures)
        if(index != indexStart){
          existingFeatures <- read.table(outputFilename, header = TRUE, sep = " ")
          allFeatures <- rbind.fill(existingFeatures, allFeatures)
        }
        write.table(allFeatures, outputFilename, quote = FALSE, row.names = FALSE)
      }
    }
  }
}








disp("Success!")

toc()
