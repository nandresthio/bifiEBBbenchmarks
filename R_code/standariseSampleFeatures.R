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

# The code that follows is divided into three sections, similar to the file
# standariseFeaturesAndFilterInstances.R, done for sample features, each of 
# which can be run for reproducibility purposes. Note no instance filtering
# is performed as this requires algorithm performance.

# - The first combines the files from the cluster into a single file.
#   This file is given with the repo, so new users should not require to run this

# - The second cleans the feature values by removing undesired
#   features and removing features with None and Inf.

# - The third standarises the feature values.
#   Bounded features are linearly scaled to lie in [-2,2].
#   Unbounded features are normalised, given a standard mean and variance,
#   and then bound within [-4,4]. This is done so that the feature
#   values are comparable when running the instance filtering. "Difference" 
#   features i.e. the difference between the processed "real" and "sample" 
#   features is also calculated.

source("R_code/dataProcessor.R")

# # COMBINE ALL THE FEATURE FILES INTO A SINGLE FILE
# features <- combineArrayResults("sampleFeatureRun", 1, 312)
# solarFeatures <- read.table("data/clusterResults/featureRun1-2700_SOLAR.txt", header = TRUE, sep = " ", fill = TRUE)
# solarFeatures <- rbind(solarFeatures, read.table("data/clusterResults/featureRun2701-5400_SOLAR.txt", header = TRUE, sep = " ", fill = TRUE))
# solarFeatures <- rbind(solarFeatures, read.table("data/clusterResults/featureRun5401-8100_SOLAR.txt", header = TRUE, sep = " ", fill = TRUE))
# solarFeatures <- rbind(solarFeatures, read.table("data/clusterResults/featureRun8101-10800_SOLAR.txt", header = TRUE, sep = " ", fill = TRUE))
# allFeatures <- rbind.fill(solarFeatures, features)
# write.table(allFeatures, "data/features/sampleFeatures.txt", quote = FALSE, row.names = FALSE)
# # COMPLETED COMBINING FEATURE FILES



# # REMOVE FEATURES WITH NAS, INFINITES, AND UNDESIRED FEATURES
# allFeatures <- read.table("data/features/sampleFeatures.txt", header = TRUE, sep = " ")
# # First of all, set default level_mmce_lda behaviour to 1 (i.e. error is 1)
# # when not enough data exists to train a linear model
# for(feat in colnames(allFeatures)[str_which(colnames(allFeatures), "ela_level_mmce_lda")]){
#   allFeatures[is.na(allFeatures[feat]), feat] <- 1
# }
# # So here want to remove any features with nas, infinites, or non changing values
# print(colnames(allFeatures[, colSums(is.na(allFeatures)) > 0]))
# allFeatures <- allFeatures[, colSums(is.na(allFeatures)) == 0]
# print(is.finite(colSums(allFeatures[, -1]))[is.finite(colSums(allFeatures[, -1])) == FALSE])
# allFeatures <- allFeatures[, c(TRUE, is.finite(colSums(allFeatures[, c(-1, -ncol(allFeatures))])), TRUE)]
# # Do not use runtime or function evaluations, or the basic set of features
# allFeatures <- allFeatures[str_which(colnames(allFeatures), "runtime", negate = TRUE)]
# allFeatures <- allFeatures[str_which(colnames(allFeatures), "fun_evals", negate = TRUE)]
# allFeatures <- allFeatures[str_which(colnames(allFeatures), "basic", negate = TRUE)]
# # Remove the C features as they are almost identical to CC features
# allFeatures <- allFeatures[str_which(colnames(allFeatures), "feature_LC_", negate = TRUE)]
# allFeatures <- allFeatures[colnames(allFeatures) != "feature_C"]
# # Remove features which always have the same value
# for(name in colnames(allFeatures[, str_which(colnames(allFeatures), "feature_", negate = FALSE)])){
#   if(min(allFeatures[, name]) == max(allFeatures[, name])){
#     print(name)
#     allFeatures <- allFeatures[, colnames(allFeatures) != name]
#   }
# }
# # Remove other features which will mess with the instance selection
# # Convex / linearity features which are heavily impacted by function range
# allFeatures <- allFeatures[str_which(colnames(allFeatures), "lin_dev_orig", negate = TRUE)]
# allFeatures <- allFeatures[str_which(colnames(allFeatures), "lin_dev_abs", negate = TRUE)]
# # Also heavily affected by function range: linear intercept and coefficients of a linear / quadratic model
# allFeatures <- allFeatures[str_which(colnames(allFeatures), "lin_simple_intercept", negate = TRUE)]
# allFeatures <- allFeatures[str_which(colnames(allFeatures), "lin_simple_coef_min", negate = TRUE)]
# allFeatures <- allFeatures[-str_which(colnames(allFeatures), "lin_simple_coef_max", negate = FALSE)[!(str_which(colnames(allFeatures), "lin_simple_coef_max", negate = FALSE) %in% str_which(colnames(allFeatures), "lin_simple_coef_max_by_min", negate = FALSE))]]
# # Take out disp diff features, they are heavily affected by the domain of the functions
# allFeatures <- allFeatures[str_which(colnames(allFeatures), "disp_diff", negate = TRUE)]
# # Take out pca _x features, which only use the sample locations for PCA
# # as this is dependent on dimension and the samples are "large"
# allFeatures <- allFeatures[str_which(colnames(allFeatures), "_x", negate = TRUE)]
# # Also take out PC1 features as this is heavily correlated with dimension,
# # keep var_cov and var_cor as these are the relative number of components needed
# allFeatures <- allFeatures[str_which(colnames(allFeatures), "_PC1", negate = TRUE)]
# # Remove conv probability as it needs extra sampling
# allFeatures <- allFeatures[str_which(colnames(allFeatures), "conv", negate = TRUE)]
# # Save this "cleaned" set of features
# write.table(allFeatures, "data/features/sampleFeaturesClean.txt", quote = FALSE, row.names = FALSE)
# 
# # Add "real" feature values
# realFeatures <- read.table("data/features/featuresClean.txt", header = TRUE, sep = " ")
# # Add column to read in features with the name of the instance
# allFeatures$functionName <- gsub('[(]', '', sapply(strsplit(allFeatures$instance, ","), "[[", 1))
# order <- match(allFeatures$functionName, realFeatures$instances)
# realFeatureNames <- colnames(realFeatures[-1])
# allFeatures[realFeatureNames] <- realFeatures[order, realFeatureNames]
# allFeatures <- allFeatures[str_which(colnames(allFeatures), "functionName", negate = TRUE)]
# write.table(allFeatures, "data/features/sampleAndRealFeaturesClean.txt", quote = FALSE, row.names = FALSE)
# # COMPLETED CLEANING OF FEATURE VALUES





# # STANDARISE FEATURE VALUES FOR INSTANCE FILTERING
# allFeatures <- read.table("data/features/sampleFeaturesClean.txt", header = TRUE, sep = " ")
# # Do this using matlab automatic boxox
# library(R.matlab)
# Matlab$startServer()
# matlab <- Matlab()
# isOpen <- open(matlab)
# temp <- allFeatures[str_which(colnames(allFeatures), "feature_")]
# Xtemp <- data.matrix(temp)
# ncols <- ncol(Xtemp)
# setVariable(matlab, X = Xtemp)
# setVariable(matlab, ncols = ncols)
# # Normalise data using Box-Cox and Z transformations
# evaluate(matlab, "out.minX = min(X,[],1);")
# evaluate(matlab, "X = bsxfun(@minus,X,out.minX)+1;")
# for(i in 1:ncols){
#   print(i)
#   setVariable(matlab, i = i)
#   evaluate(matlab, "aux = X(:,i);")
#   evaluate(matlab, "idx = isnan(aux);")
#   evaluate(matlab, "[aux, temp] = boxcox(aux(~idx));")
#   evaluate(matlab, "[aux, temp, temp2] = zscore(aux);")
#   evaluate(matlab, "X(~idx,i) = aux;")
# }
# Xout <- getVariable(matlab, "X")[[1]]
# allFeaturesNormalised <- allFeatures
# allFeaturesNormalised[str_which(colnames(allFeaturesNormalised), "feature_")] <- Xout
# close(matlab)
# 
# # List the feature values based on their ranges
# # Features with the range [0,1]
# featuresBound01 <- c("feature_sample_CC",
#                      "feature_sample_LCC_0_1",
#                      "feature_sample_LCC_0_2",
#                      "feature_sample_LCC_0_3",
#                      "feature_sample_LCC_0_4",
#                      "feature_sample_LCC_0_5",
#                      "feature_sample_LCC_0_6",
#                      "feature_sample_LCC_0_7",
#                      "feature_sample_LCC_0_8",
#                      "feature_sample_LCC_0_9",
#                      "feature_sample_LCC_0_95",
#                      "feature_sample_LCC_0_975",
#                      "feature_sample_LCC_sd",
#                      "feature_sample_LCC_mean",
#                      "feature_sample_LCCrel_0_1",
#                      "feature_sample_LCCrel_0_2",
#                      "feature_sample_LCCrel_0_3",
#                      "feature_sample_LCCrel_0_4",
#                      "feature_sample_LCCrel_0_5",
#                      "feature_sample_LCCrel_0_6",
#                      "feature_sample_LCCrel_0_7",
#                      "feature_sample_LCCrel_0_8",
#                      "feature_sample_LCCrel_0_9",
#                      "feature_sample_LCCrel_0_95",
#                      "feature_sample_LCCrel_0_975",
#                      "feature_sample_LCCrel_sd",
#                      "feature_sample_LCCrel_mean",
#                      "feature_sample_high_ela_level_mmce_lda_10",
#                      "feature_sample_high_ela_level_mmce_qda_10",
#                      "feature_sample_high_ela_level_mmce_lda_25",
#                      "feature_sample_high_ela_level_mmce_qda_25",
#                      "feature_sample_high_ela_level_mmce_lda_50",
#                      "feature_sample_high_ela_level_mmce_qda_50",
#                      "feature_sample_low_ela_level_mmce_lda_10",
#                      "feature_sample_low_ela_level_mmce_qda_10",
#                      "feature_sample_low_ela_level_mmce_lda_25",
#                      "feature_sample_low_ela_level_mmce_qda_25",
#                      "feature_sample_low_ela_level_mmce_lda_50",
#                      "feature_sample_low_ela_level_mmce_qda_50",
#                      "feature_sample_mid_ela_level_mmce_lda_10",
#                      "feature_sample_mid_ela_level_mmce_qda_10",
#                      "feature_sample_mid_ela_level_mmce_lda_25",
#                      "feature_sample_mid_ela_level_mmce_qda_25",
#                      "feature_sample_mid_ela_level_mmce_lda_50",
#                      "feature_sample_mid_ela_level_mmce_qda_50",
#                      "feature_sample_high_ela_meta_lin_simple_adj_r2",
#                      "feature_sample_high_ela_meta_lin_w_interact_adj_r2",
#                      "feature_sample_high_ela_meta_quad_simple_adj_r2",
#                      "feature_sample_high_ela_meta_quad_w_interact_adj_r2",
#                      "feature_sample_low_ela_meta_lin_simple_adj_r2",
#                      "feature_sample_low_ela_meta_lin_w_interact_adj_r2",
#                      "feature_sample_low_ela_meta_quad_simple_adj_r2",
#                      "feature_sample_low_ela_meta_quad_w_interact_adj_r2",
#                      "feature_sample_mid_ela_meta_lin_simple_adj_r2",
#                      "feature_sample_mid_ela_meta_lin_w_interact_adj_r2",
#                      "feature_sample_mid_ela_meta_quad_simple_adj_r2",
#                      "feature_sample_mid_ela_meta_quad_w_interact_adj_r2",
#                      "feature_sample_high_pca_expl_var_cov_init",
#                      "feature_sample_high_pca_expl_var_cor_init",
#                      "feature_sample_low_pca_expl_var_cov_init",
#                      "feature_sample_low_pca_expl_var_cor_init",
#                      "feature_sample_mid_pca_expl_var_cov_init",
#                      "feature_sample_mid_pca_expl_var_cor_init",
#                      "feature_sample_high_nbc_nn_nb_sd_ratio",
#                      "feature_sample_high_nbc_nn_nb_mean_ratio",
#                      "feature_sample_low_nbc_nn_nb_sd_ratio",
#                      "feature_sample_low_nbc_nn_nb_mean_ratio",
#                      "feature_sample_mid_nbc_nn_nb_sd_ratio",
#                      "feature_sample_mid_nbc_nn_nb_mean_ratio",
#                      "feature_sample_high_ic_h_max",
#                      "feature_sample_high_ic_m0",
#                      "feature_sample_low_ic_h_max",
#                      "feature_sample_low_ic_m0",
#                      "feature_sample_mid_ic_h_max",
#                      "feature_sample_mid_ic_m0")
# 
# for(feat in featuresBound01){
#   if(!(feat %in% colnames(allFeatures))){
#     print(paste0("Note feat ", feat, " is not being used as it conatined either inf or na values"))
#     featuresBound01 <- featuresBound01[featuresBound01 != feat]
#   }
# }
# 
# # Features with the range [-1,1]
# featuresBound11 <- c("feature_sample_high_nbc_nn_nb_cor",
#                      "feature_sample_high_nbc_nb_fitness_cor",
#                      "feature_sample_low_nbc_nn_nb_cor",
#                      "feature_sample_low_nbc_nb_fitness_cor",
#                      "feature_sample_mid_nbc_nn_nb_cor",
#                      "feature_sample_mid_nbc_nb_fitness_cor")
# 
# for(feat in featuresBound11){
#   if(!(feat %in% colnames(allFeatures))){
#     print(paste0("Note feat ", feat, " is not being used as it conatined either inf or na values"))
#     featuresBound11 <- featuresBound11[featuresBound11 != feat]
#   }
# }
# 
# # Features with an unbound range which need to be normalised
# features_norm <- c("feature_sample_RRMSE",
#                    "feature_sample_LCC_coeff",
#                    "feature_sample_LCCrel_coeff",
#                    "feature_sample_high_ela_distr_skewness",
#                    "feature_sample_high_ela_distr_kurtosis",
#                    "feature_sample_high_ela_distr_number_of_peaks",
#                    "feature_sample_low_ela_distr_skewness",
#                    "feature_sample_low_ela_distr_kurtosis",
#                    "feature_sample_low_ela_distr_number_of_peaks",
#                    "feature_sample_mid_ela_distr_skewness",
#                    "feature_sample_mid_ela_distr_kurtosis",
#                    "feature_sample_mid_ela_distr_number_of_peaks",
#                    "feature_sample_high_ela_meta_lin_simple_coef_max_by_min",
#                    "feature_sample_high_ela_meta_quad_simple_cond",
#                    "feature_sample_low_ela_meta_lin_simple_coef_max_by_min",
#                    "feature_sample_low_ela_meta_quad_simple_cond",
#                    "feature_sample_mid_ela_meta_lin_simple_coef_max_by_min",
#                    "feature_sample_mid_ela_meta_quad_simple_cond",
#                    "feature_sample_high_nbc_dist_ratio_coeff_var",
#                    "feature_sample_low_nbc_dist_ratio_coeff_var",
#                    "feature_sample_mid_nbc_dist_ratio_coeff_var",
#                    "feature_sample_high_ic_eps_s",
#                    "feature_sample_high_ic_eps_max",
#                    "feature_sample_high_ic_eps_ratio",
#                    "feature_sample_low_ic_eps_s",
#                    "feature_sample_low_ic_eps_max",
#                    "feature_sample_low_ic_eps_ratio",
#                    "feature_sample_mid_ic_eps_s",
#                    "feature_sample_mid_ic_eps_max",
#                    "feature_sample_mid_ic_eps_ratio",
#                    "feature_sample_high_ela_level_lda_qda_10",
#                    "feature_sample_high_ela_level_lda_qda_25",
#                    "feature_sample_high_ela_level_lda_qda_50",
#                    "feature_sample_low_ela_level_lda_qda_10",
#                    "feature_sample_low_ela_level_lda_qda_25",
#                    "feature_sample_low_ela_level_lda_qda_50")
# for(feat in features_norm){
#   if(!(feat %in% colnames(allFeatures))){
#     print(paste0("Note feat ", feat, " is not being used as it conatined either inf or na values"))
#     features_norm <- features_norm[features_norm != feat]
#   }
# }
# 
# # Features with an unbound value which only need to be scaled
# features_scale <- c("feature_sample_lowFiBudget",
#                     "feature_sample_highFiBudget",
#                     "feature_sample_lowFiBudgetRatio",
#                     "feature_sample_highFiBudgetRatio",
#                     "feature_sample_budgetRatio",
#                     "feature_sample_high_disp_ratio_mean_02",
#                     "feature_sample_high_disp_ratio_mean_05",
#                     "feature_sample_high_disp_ratio_mean_10",
#                     "feature_sample_high_disp_ratio_mean_25",
#                     "feature_sample_high_disp_ratio_median_02",
#                     "feature_sample_high_disp_ratio_median_05",
#                     "feature_sample_high_disp_ratio_median_10",
#                     "feature_sample_high_disp_ratio_median_25",
#                     "feature_sample_low_disp_ratio_mean_02",
#                     "feature_sample_low_disp_ratio_mean_05",
#                     "feature_sample_low_disp_ratio_mean_10",
#                     "feature_sample_low_disp_ratio_mean_25",
#                     "feature_sample_low_disp_ratio_median_02",
#                     "feature_sample_low_disp_ratio_median_05",
#                     "feature_sample_low_disp_ratio_median_10",
#                     "feature_sample_low_disp_ratio_median_25",
#                     "feature_sample_mid_disp_ratio_mean_02",
#                     "feature_sample_mid_disp_ratio_mean_05",
#                     "feature_sample_mid_disp_ratio_mean_10",
#                     "feature_sample_mid_disp_ratio_mean_25",
#                     "feature_sample_mid_disp_ratio_median_02",
#                     "feature_sample_mid_disp_ratio_median_05",
#                     "feature_sample_mid_disp_ratio_median_10",
#                     "feature_sample_mid_disp_ratio_median_25")
# 
# for(feat in features_scale){
#   if(!(feat %in% colnames(allFeatures))){
#     print(paste0("Note feat ", feat, " is not being used as it conatined either inf or na values"))
#     features_scale <- features_scale[features_scale != feat]
#   }
# }
# 
# 
# standarisedData <- allFeatures[c("instances", featuresBound01, featuresBound11, features_norm, features_scale)]
# 
# # First bound r2 features to lie in [0,1]
# for(feat in colnames(standarisedData[str_which(colnames(standarisedData), "r2")])){
#   print(nrow(standarisedData[standarisedData[feat] < 0, ]))
#   standarisedData[standarisedData[feat] < 0, feat] <- 0
# }
# 
# 
# for(feat in featuresBound01){
#   standarisedData[feat] <- standarisedData[feat] * 4 - 2
# }
# for(feat in featuresBound11){
#   standarisedData[feat] <- standarisedData[feat] * 2
# }
# for(feat in features_norm){
#   standarisedData[feat] <- allFeaturesNormalised[feat]
#   standarisedData[standarisedData[feat] > 4, feat] <- 4
#   standarisedData[standarisedData[feat] < -4, feat] <- -4
# 
# }
# for(feat in features_scale){
#   standarisedData[feat] <- (standarisedData[feat] - min(standarisedData[feat])) / (max(standarisedData[feat]) - min(standarisedData[feat]))
#   standarisedData[feat] <- standarisedData[feat] * 4 - 2
# }
# 
# if(ncol(allFeatures[!(colnames(allFeatures) %in% c("instances", featuresBound01, featuresBound11, features_norm, features_scale))]) > 0){
#   print(paste0("NOTE: Not including features ",
#                colnames(allFeatures)[!(colnames(allFeatures) %in% c("instances", featuresBound01, featuresBound11, features_norm, features_scale))],
#   " as not specified what transformation to apply to them!"))
# }
# write.table(standarisedData, "data/features/sampleFeaturesCleanStandarised.txt", quote = FALSE, row.names = FALSE)
# 
# 
# # Now want diff features i.e. the difference between real features and sample features
# # when these were available
# # Start by reading in the actual feature data and adding it
# realFeaturesStandarised <- read.table("data/features/featuresCleanStandarised.txt", header = TRUE, sep = " ")
# # Add column to read in features with the name of the instance
# standarisedData$functionName <- gsub('[(]', '', sapply(strsplit(standarisedData$instance, ","), "[[", 1))
# order <- match(standarisedData$functionName, realFeaturesStandarised$instances)
# realFeatureNames <- colnames(realFeaturesStandarised[-1])
# standarisedData[realFeatureNames] <- realFeaturesStandarised[order, realFeatureNames]
# standarisedData <- standarisedData[str_which(colnames(standarisedData), "functionName", negate = TRUE)]
# 
# nonSampleFeat <- c()
# for(feat in colnames(standarisedData[str_which(colnames(standarisedData), "feature")])){
#   if(substr(feat, 1, 14) == "feature_sample"){next}
#   nonSampleFeat <- c(nonSampleFeat, substr(feat, 8, nchar(feat)))
#   # Get actual feature name
#   actualFeat <- substr(feat, 8, nchar(feat))
#   # Check an equivalent sample feature exists
#   if(length(str_which(colnames(standarisedData), paste0("feature_sample", actualFeat))) == 0){
#     print(paste0("Skipping ", actualFeat, " as do not have sample feature value"))
#     next
#   }else if(length(str_which(colnames(standarisedData), paste0("feature_", actualFeat))) == 1){
#     print(paste0("Skipping ", actualFeat, " as do not have real feature value"))
#     next
#   }
#   
#   standarisedData[paste0("feature_diff", actualFeat)] <- standarisedData[paste0("feature", actualFeat)] - standarisedData[paste0("feature_sample", actualFeat)]
#   print(paste0("Max and min diff of ", actualFeat, " ", min(standarisedData[paste0("feature_diff", actualFeat)]), " ", max(standarisedData[paste0("feature_diff", actualFeat)])))
# }
# 
# # Add "real" keyword to features calculated with large sample
# colnames(standarisedData)[colnames(standarisedData) %in% paste0("feature", nonSampleFeat)] <- paste0("feature_real", nonSampleFeat)
# 
# write.table(standarisedData, "data/features/sampleAndRealFeaturesCleanStandarised.txt", quote = FALSE, row.names = FALSE)
# # COMPLETED STANDARISE FEATURE VALUES



