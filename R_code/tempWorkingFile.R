data <- data.frame(matrix(ncol = 3, nrow = 0))
row <- 0
for(i in 1:nrow(Xfeat)){
  for(j in (i+1):nrow(Xfeat)){
    row <- row + 1
    print(paste0(i, "-", j, ":", sqrt(sum((Xfeat[i,] - Xfeat[j,])^2))))
    data[row, ] <- c(i, j, sqrt(sum((Xfeat[i,] - Xfeat[j,])^2)))
  }
}

temp <- allFeatures[str_which(allFeatures$instances, "DisturbanceBasedFunction1-"),]


Xfeat <- temp[, str_which(colnames(temp), "feature_")]
Xfeat <- data.matrix(Xfeat)




val <- purifyInstances(Xfeat[1:220, ], 16)

temp <- allFeaturesNormalised[1:220, ]
temp <- temp[!(val$preclude), ]
temp

testDim <- allFeaturesNormalised[[2]]
testDim[testDim > testDim + 5*IQR(testDim)] <- testDim + 5*IQR(testDim)
testDim[testDim < testDim - 5*IQR(testDim)] <- testDim - 5*IQR(testDim)

data <- allFeatures
chosenData <- standarisedData[standarisedData$instances %in% finalInstances$instances, ]
for(i in 2:ncol(data)){
  feat <- colnames(data)[[i]]
  print(i)
  png(paste0(i, "-", feat, ".png"))
  plot(data$feature_dimension, data[feat][[1]])
  dev.off()
  
  png(paste0(i, "-", feat, "processed.png"))
  plot(data$feature_dimension, standarisedData[feat][[1]], ylim = c(-4, 4))
  dev.off()
  
  # int <- (max(standarisedData[feat]) - min(standarisedData[feat])) / 20
  # xVals <- seq(min(standarisedData[feat]), max(standarisedData[feat]), int)
  # yVals <- c()
  # for(val in xVals){
  #   yVals <- c(yVals, sum(standarisedData[feat] <= val & standarisedData[feat] > (val - int)) / nrow(standarisedData))
  # }
  # png(paste0(i, "-", feat, "processedDist.png"))
  # plot(xVals, yVals, type = 'l', ylim = c(0,1))
  # dev.off()
  # 
  # int <- (max(data[feat]) - min(data[feat])) / 20
  # xVals <- seq(min(data[feat]), max(data[feat]), int)
  # yVals <- c()
  # for(val in xVals){
  #   yVals <- c(yVals, sum(data[feat] <= val & data[feat] > (val - int)) / nrow(data))
  # }
  # png(paste0(i, "-", feat, "Dist.png"))
  # plot(xVals, yVals, type = 'l', ylim = c(0,1))
  # dev.off()
  # 
  # 
  # 
  # png(paste0(i, "-", feat, "processedFinal.png"))
  # plot(chosenData$feature_dimension, chosenData[feat][[1]])
  # dev.off()
  # 
  # int <- (max(standarisedData[feat]) - min(standarisedData[feat])) / 20
  # xVals <- seq(min(standarisedData[feat]), max(standarisedData[feat]), int)
  # yVals <- c()
  # for(val in xVals){
  #   yVals <- c(yVals, sum(chosenData[feat] <= val & chosenData[feat] > (val - int)) / nrow(chosenData))
  # }
  # png(paste0(i, "-", feat, "processedFinalDist.png"))
  # plot(xVals, yVals, type = 'l', ylim = c(0,1))
  # dev.off()
  
  # png(paste0(i, "-", feat, "chosenSubset.png"))
  # plot(chosenData$feature_dimension, chosenData[feat][[1]])
  # dev.off()
  # 
  # 
  # png(paste0(i, "-", colnames(allFeatures)[[i]], "Norm.png"))
  # plot(allFeaturesNormalised$feature_dimension, allFeaturesNormalised[i][[1]])
  # dev.off()
  # test <- allFeaturesNormalised[[i]]
  # test[test > test + 5*IQR(test)] <- test + 5*IQR(test)
  # test[test < test - 5*IQR(test)] <- test - 5*IQR(test)
  # png(paste0(i, "-", colnames(allFeatures)[[i]], "NormUnif.png"))
  # plot(testDim, test)
  # dev.off()
  
  # print(paste0(colnames(allFeaturesNormalised)[[i]], ": ", min(allFeaturesNormalised[i]), " ", max(allFeaturesNormalised[i])))

}
library(ggplot2)
test <- allFeaturesNormalised
test$feature_dimension <- allFeatures$feature_dimension

ggplot(allFeatures, aes(x=as.factor(feature_dimension), y=feature_CC)) + 
  geom_boxplot()



  
  
plot(allFeaturesNormalised$feature_dimension, allFeaturesNormalised)

test <- allFeaturesNormalised[c("instances", "feature_mid_disp_ratio_median_25")]
test <- allFeatures[c("instances", "feature_mid_disp_ratio_median_25")]

test <- finalInstancesOriginal[str_which(finalInstancesOriginal$instances, "COCO", negate = TRUE), ]

for(dim in 1:20){
  if(sum(finalInstancesOriginal$feature_dimension == dim) == 0){next}
  print(paste0(dim, " ", sum(finalInstancesOriginal$feature_dimension == dim)))
}



test <- finalInstances[str_which(finalInstances$instances, "amp0.5"), ]






ctrl = list(init_sample.type = "lhs",
            init_sample.lower = rep(0,15),
            init_sample.upper = rep(0,15))

X = createInitialSample(n.obs = 1000 * 15, dim = 15, control = ctrl)
write.table()
paste0("../data/samplePlans/LHS-dim", dim, "-n", 1000*dim,".txt")




names <- colnames(allFeatures)[!(colnames(allFeatures) %in% c(featuresBound01, featuresBound11, features_unbound, notSure))]
temp <- allFeatures[names]




corrs <- c()
for(feat in colnames(standarisedData[-1])){
  corrs <- c(corrs, cor(allFeatures$feature_dimension, standarisedData[feat]))
  print(paste0(feat, ": ", cor(allFeatures$feature_dimension, standarisedData[feat])))
}

lowCorrFeats <- colnames(standarisedData[-1])[abs(corrs) < 0.3]




temp <- read.table("data/clusterResults/sampleFeatureRun_arrayJob1-10.txt", header = TRUE, sep = " ")


temp <- allFeatures[c(1, str_which(colnames(allFeatures), "level"))]

allFeatures$feature_dimension <- as.factor(allFeatures$feature_dimension)
ggscatter(allFeatures[allFeatures$feature_dimension == 4, ], x = "feature_LCC_0_5", y = "feature_CC", shape = "feature_dimension")
plot()


temp <- features[features$feature_dimension == 10, 1:20]
plot(temp$feature_LCC_0_7, temp$feature_CC)


allFeatures <- features


c("instances", featuresBound01, featuresBound11, features_norm, features_scale)[!(c("instances", featuresBound01, featuresBound11, features_norm, features_scale) %in% colnames(allFeatures))]


temp <- features[c("instances", "feature_low_ela_level_lda_qda_50")]




