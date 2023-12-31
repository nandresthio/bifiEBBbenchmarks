######################################################################################

# Copyright 2023, Nicolau Andres-Thio

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

# The code that follows is divided into four sections, each of which can be run
# for reproducibility purposes:

# - The first combines the files from the cluster into a single file.
#   This file is given with the repo, so new users should not require to run this

# - The second cleans the feature values by removing undesired
#   features and removing features with None and Inf.

# - The third standarises the feature values.
#   Bounded features are linearly scaled to lie in [-2,2].
#   Unbounded features are normalised, given a standard mean and variance,
#   and then bound within [-4,4]. This is done so that the feature
#   values are comparable when running the instance filtering

# - The fourth runs an instance filtering technique to create as varied a subset
#   as possible to conduct further testing.

source("R_code/dataProcessor.R")

# # COMBINE ALL THE FEATURE FILES INTO A SINGLE FILE
# features <- combineArrayResults("featureRun", 1, 8120, 10)
# solarFeatures <- read.table("data/clusterResults/featureRun1_SOLAR.txt", header = TRUE, sep = " ", fill = TRUE)
# for(i in 2:9){
#   solarFeatures <- rbind(solarFeatures, read.table(paste0("data/clusterResults/featureRun", i, "_SOLAR.txt"), header = TRUE, sep = " ", fill = TRUE))
# }
# solarFeatures <- solarFeatures[1:9, ]
# allFeatures <- rbind.fill(solarFeatures, features)
# write.table(allFeatures, "data/features/features.txt", quote = FALSE, row.names = FALSE)
# # COMPLETED COMBINING FEATURE FILES



# # REMOVE FEATURES WITH NAS, INFINITES, AND UNDESIRED FEATURES
# allFeatures <- read.table("data/features/features.txt", header = TRUE, sep = " ")
# # First of all, set default level_mmce_lda behaviour to 1 (i.e. error is 1)
# # when not enough data exists to train a linear or a quadratic model
# for(feat in colnames(allFeatures)[str_which(colnames(allFeatures), "ela_level_mmce_lda")]){
#   allFeatures[is.na(allFeatures[feat]), feat] <- 1
# }
# for(feat in colnames(allFeatures)[str_which(colnames(allFeatures), "ela_level_mmce_qda")]){
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
# write.table(allFeatures, "data/features/featuresClean.txt", quote = FALSE, row.names = FALSE)
# # COMPLETED CLEANING OF FEATURE VALUES




# # STANDARISE FEATURE VALUES FOR INSTANCE FILTERING
# allFeatures <- read.table("data/features/featuresClean.txt", header = TRUE, sep = " ")
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
# featuresBound01 <- c("feature_CC",
#                      "feature_LCC_0_1",
#                      "feature_LCC_0_2",
#                      "feature_LCC_0_3",
#                      "feature_LCC_0_4",
#                      "feature_LCC_0_5",
#                      "feature_LCC_0_6",
#                      "feature_LCC_0_7",
#                      "feature_LCC_0_8",
#                      "feature_LCC_0_9",
#                      "feature_LCC_0_95",
#                      "feature_LCC_0_975",
#                      "feature_LCC_sd",
#                      "feature_LCC_mean",
#                      "feature_LCCrel_0_1",
#                      "feature_LCCrel_0_2",
#                      "feature_LCCrel_0_3",
#                      "feature_LCCrel_0_4",
#                      "feature_LCCrel_0_5",
#                      "feature_LCCrel_0_6",
#                      "feature_LCCrel_0_7",
#                      "feature_LCCrel_0_8",
#                      "feature_LCCrel_0_9",
#                      "feature_LCCrel_0_95",
#                      "feature_LCCrel_0_975",
#                      "feature_LCCrel_sd",
#                      "feature_LCCrel_mean",
#                      "feature_high_ela_level_mmce_lda_10",
#                      "feature_high_ela_level_mmce_qda_10",
#                      "feature_high_ela_level_mmce_lda_25",
#                      "feature_high_ela_level_mmce_qda_25",
#                      "feature_high_ela_level_mmce_lda_50",
#                      "feature_high_ela_level_mmce_qda_50",
#                      "feature_low_ela_level_mmce_lda_10",
#                      "feature_low_ela_level_mmce_qda_10",
#                      "feature_low_ela_level_mmce_lda_25",
#                      "feature_low_ela_level_mmce_qda_25",
#                      "feature_low_ela_level_mmce_lda_50",
#                      "feature_low_ela_level_mmce_qda_50",
#                      "feature_mid_ela_level_mmce_lda_10",
#                      "feature_mid_ela_level_mmce_qda_10",
#                      "feature_mid_ela_level_mmce_lda_25",
#                      "feature_mid_ela_level_mmce_qda_25",
#                      "feature_mid_ela_level_mmce_lda_50",
#                      "feature_mid_ela_level_mmce_qda_50",
#                      "feature_high_ela_meta_lin_simple_adj_r2",
#                      "feature_high_ela_meta_lin_w_interact_adj_r2",
#                      "feature_high_ela_meta_quad_simple_adj_r2",
#                      "feature_high_ela_meta_quad_w_interact_adj_r2",
#                      "feature_low_ela_meta_lin_simple_adj_r2",
#                      "feature_low_ela_meta_lin_w_interact_adj_r2",
#                      "feature_low_ela_meta_quad_simple_adj_r2",
#                      "feature_low_ela_meta_quad_w_interact_adj_r2",
#                      "feature_mid_ela_meta_lin_simple_adj_r2",
#                      "feature_mid_ela_meta_lin_w_interact_adj_r2",
#                      "feature_mid_ela_meta_quad_simple_adj_r2",
#                      "feature_mid_ela_meta_quad_w_interact_adj_r2",
#                      "feature_high_pca_expl_var_cov_init",
#                      "feature_high_pca_expl_var_cor_init",
#                      "feature_low_pca_expl_var_cov_init",
#                      "feature_low_pca_expl_var_cor_init",
#                      "feature_mid_pca_expl_var_cov_init",
#                      "feature_mid_pca_expl_var_cor_init",
#                      "feature_high_nbc_nn_nb_sd_ratio",
#                      "feature_high_nbc_nn_nb_mean_ratio",
#                      "feature_low_nbc_nn_nb_sd_ratio",
#                      "feature_low_nbc_nn_nb_mean_ratio",
#                      "feature_mid_nbc_nn_nb_sd_ratio",
#                      "feature_mid_nbc_nn_nb_mean_ratio",
#                      "feature_high_ic_h_max",
#                      "feature_high_ic_m0",
#                      "feature_low_ic_h_max",
#                      "feature_low_ic_m0",
#                      "feature_mid_ic_h_max",
#                      "feature_mid_ic_m0")
# 
# # Features with the range [-1,1]
# featuresBound11 <- c("feature_high_nbc_nn_nb_cor",
#                      "feature_high_nbc_nb_fitness_cor",
#                      "feature_low_nbc_nn_nb_cor",
#                      "feature_low_nbc_nb_fitness_cor",
#                      "feature_mid_nbc_nn_nb_cor",
#                      "feature_mid_nbc_nb_fitness_cor")
# 
# # Features with an unbound range which need to be normalised
# features_norm <- c("feature_dimension",
#                    "feature_RRMSE",
#                    "feature_LCC_coeff",
#                    "feature_LCCrel_coeff",
#                    "feature_high_ela_distr_skewness",
#                    "feature_high_ela_distr_kurtosis",
#                    "feature_high_ela_distr_number_of_peaks",
#                    "feature_low_ela_distr_skewness",
#                    "feature_low_ela_distr_kurtosis",
#                    "feature_low_ela_distr_number_of_peaks",
#                    "feature_mid_ela_distr_skewness",
#                    "feature_mid_ela_distr_kurtosis",
#                    "feature_mid_ela_distr_number_of_peaks",
#                    "feature_high_ela_meta_lin_simple_coef_max_by_min",
#                    "feature_high_ela_meta_quad_simple_cond",
#                    "feature_low_ela_meta_lin_simple_coef_max_by_min",
#                    "feature_low_ela_meta_quad_simple_cond",
#                    "feature_mid_ela_meta_lin_simple_coef_max_by_min",
#                    "feature_mid_ela_meta_quad_simple_cond",
#                    "feature_high_nbc_dist_ratio_coeff_var",
#                    "feature_low_nbc_dist_ratio_coeff_var",
#                    "feature_mid_nbc_dist_ratio_coeff_var",
#                    "feature_high_ic_eps_s",
#                    "feature_high_ic_eps_max",
#                    "feature_high_ic_eps_ratio",
#                    "feature_low_ic_eps_s",
#                    "feature_low_ic_eps_max",
#                    "feature_low_ic_eps_ratio",
#                    "feature_mid_ic_eps_s",
#                    "feature_mid_ic_eps_max",
#                    "feature_mid_ic_eps_ratio",
#                    "feature_high_ela_level_lda_qda_10",
#                    "feature_high_ela_level_lda_qda_25",
#                    "feature_high_ela_level_lda_qda_50",
#                    "feature_low_ela_level_lda_qda_10",
#                    "feature_low_ela_level_lda_qda_25",
#                    "feature_low_ela_level_lda_qda_50")
# 
# # Features with an unbound value which only need to be linearly scaled
# features_scale <- c("feature_high_disp_ratio_mean_02",
#                     "feature_high_disp_ratio_mean_05",
#                      "feature_high_disp_ratio_mean_10",
#                     "feature_high_disp_ratio_mean_25",
#                      "feature_high_disp_ratio_median_02",
#                     "feature_high_disp_ratio_median_05",
#                      "feature_high_disp_ratio_median_10",
#                     "feature_high_disp_ratio_median_25",
#                      "feature_low_disp_ratio_mean_02",
#                     "feature_low_disp_ratio_mean_05",
#                      "feature_low_disp_ratio_mean_10",
#                     "feature_low_disp_ratio_mean_25",
#                      "feature_low_disp_ratio_median_02",
#                     "feature_low_disp_ratio_median_05",
#                      "feature_low_disp_ratio_median_10",
#                     "feature_low_disp_ratio_median_25",
#                      "feature_mid_disp_ratio_mean_02",
#                     "feature_mid_disp_ratio_mean_05",
#                      "feature_mid_disp_ratio_mean_10",
#                     "feature_mid_disp_ratio_mean_25",
#                      "feature_mid_disp_ratio_median_02",
#                     "feature_mid_disp_ratio_median_05",
#                      "feature_mid_disp_ratio_median_10",
#                     "feature_mid_disp_ratio_median_25")
# 
# 
# names <- c("instances", featuresBound01, featuresBound11, features_norm, features_scale)
# for(name in names){
#   if(!(name %in% colnames(allFeatures))){
#     print(paste0("Note, not including feature ", name, " as it does not exist in allFeatures!"))
#     names <- names[names != name]
#     featuresBound01 <- featuresBound01[featuresBound01 != name]
#     featuresBound11 <- featuresBound11[featuresBound11 != name]
#     features_norm <- features_norm[features_norm != name]
#     features_scale <- features_scale[features_scale != name]
#   }
# }
# 
# standarisedData <- allFeatures[names]
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
# 
# write.table(standarisedData, "data/features/featuresCleanStandarised.txt", quote = FALSE, row.names = FALSE)
# # COMPLETED STANDARISE FEATURE VALUES FOR INSTANCE FILTERING




# # BEGIN INSTANCE FILTERING
# # What follows is the implementation of the instance filtering proposed by
# # Hossein Alipour et al. in "Enhanced instance space analysis for the maximum flow problem"
# # Define function for instance filtering based on feature values only
# purifyInstances <- function(featureData, eps){
#   # First create a dataframe in which to store all the information
#   labels <- data.frame(matrix(ncol = 0, nrow = nrow(featureData)))
#   labels$preclude <- FALSE
#   labels$dissimlar <- TRUE
#   for(i in 1:(nrow(labels)-1)){
#     if(labels[i, "preclude"]){next}
#     cat(paste0("\rWorking on row ", i, "/", nrow(labels)))
#     for(j in (i+1):nrow(labels)){
#       if(labels[j, "preclude"]){next}
#       # cat(paste0("\rWorking on row ", i, "/", nrow(labels), " and second row ", j, "/", nrow(labels)))
#       dist <- sqrt(sum((featureData[i,] - featureData[j,])^2))
#       if(dist > eps){next}
#       labels[j, "dissimlar"] <- FALSE
#       labels[j, "preclude"] <- TRUE
#       labels[j, "removed"] <- i
#     }
#   }
#   cat(paste0("\rWorking on row ", i, "/", nrow(labels), " - done\n"))
#   return(labels)
# }
# 
# # Define function which calculates the Coefficient of variation of the nearest neighbor distances
# # Essentially tells you how spread out the instances are based on feature values
# calculateCVNND <- function(featureData){
#   labels <- data.frame(matrix(ncol = 0, nrow = nrow(featureData)))
#   labels$minDist <- 0
#   if(nrow(labels) == 1){return(c(0,1))}
#   for(i in 1:nrow(labels)){
#     cat(paste0("\rWorking on row ", i, "/", nrow(labels)))
#     temp <- featureData[-i, ]
#     temp <- sweep(x = temp, MARGIN = 2, STATS = featureData[i, ], FUN = "-")
#     temp <- temp^2
#     temp <- rowSums(temp)
#     temp <- sqrt(temp)
#     labels[i, "minDist"] <- min(temp)
#   }
#   CV = sd(labels$minDist)/mean(labels$minDist)
#   Uniformity = 1 - CV
#   cat(paste0("\rWorking on row ", i, "/", nrow(labels), " - done\n"))
#   return(c(CV, Uniformity))
# }
# 
# 
# # Next prepare the standarised feature data to run through the instance filtering
# standarisedData <- read.table("data/features/featuresCleanStandarised.txt", header = TRUE, sep = " ")
# specialAllFeatures <- standarisedData
# # Separate into types of instances, namely SOLAR (which are not run through the filtering
# # as they are a "real" black-box and we want to keep all of them), literature
# # instances (found in the literature), disturbance instances (created by adding
# # a disturbance to a high fidelity function from the literature), and COCO disturbance
# # instances (created by adding a disturbance to a COCO function).
# instancesSolar <- specialAllFeatures[str_which(specialAllFeatures$instances, "SOLAR"), ]
# 
# instancesDist <- specialAllFeatures[str_which(specialAllFeatures$instances, "Disturbance"), ]
# 
# instancesCoco <- specialAllFeatures[str_which(specialAllFeatures$instances, "COCO"), ]
# 
# instancesLit <- specialAllFeatures[-c(str_which(specialAllFeatures$instances, "SOLAR"),
#                                           str_which(specialAllFeatures$instances, "Disturbance"),
#                                           str_which(specialAllFeatures$instances, "COCO")), ]
# # Scramble each of the set of instances so that order does not impact
# # the creation of the test suite. For reproducibility purposes, change seed if
# # want different randomisation
# set.seed(1)
# instancesDist <- instancesDist[sample(1:nrow(instancesDist)), ]
# instancesCoco <- instancesCoco[sample(1:nrow(instancesCoco)), ]
# instancesLit <- instancesLit[sample(1:nrow(instancesLit)), ]
# # Put them together again
# instancesOrdered <- rbind(instancesLit, instancesDist, instancesCoco)
# # To avoid overrepresentation of a type of instance, select a subset for this
# # pass of instance selection
# useFeatures <- c("feature_CC",
#                  "feature_RRMSE",
#                  "feature_LCCrel_0_5",
#                  "feature_LCCrel_0_95",
#                  "feature_LCCrel_sd",
#                  "feature_high_ela_distr_skewness",
#                  "feature_low_ela_distr_skewness",
#                  "feature_mid_ela_distr_skewness",
#                  "feature_high_ela_distr_kurtosis",
#                  "feature_low_ela_distr_kurtosis",
#                  "feature_mid_ela_distr_kurtosis",
#                  "feature_high_ela_distr_number_of_peaks",
#                  "feature_low_ela_distr_number_of_peaks",
#                  "feature_mid_ela_distr_number_of_peaks",
#                  "feature_high_ela_level_mmce_lda_25",
#                  "feature_low_ela_level_mmce_lda_25",
#                  "feature_mid_ela_level_mmce_lda_25",
#                  "feature_high_ela_meta_quad_simple_adj_r2",
#                  "feature_low_ela_meta_quad_simple_adj_r2",
#                  "feature_mid_ela_meta_quad_simple_adj_r2",
#                  "feature_high_ela_meta_quad_simple_cond",
#                  "feature_low_ela_meta_quad_simple_cond",
#                  "feature_mid_ela_meta_quad_simple_cond",
#                  "feature_high_ic_h_max",
#                  "feature_low_ic_h_max",
#                  "feature_mid_ic_h_max",
#                  "feature_high_ic_m0",
#                  "feature_low_ic_m0",
#                  "feature_mid_ic_m0",
#                  "feature_high_ic_eps_max",
#                  "feature_low_ic_eps_max",
#                  "feature_mid_ic_eps_max",
#                  "feature_high_ic_eps_ratio",
#                  "feature_low_ic_eps_ratio",
#                  "feature_mid_ic_eps_ratio",
#                  "feature_high_disp_ratio_mean_10",
#                  "feature_low_disp_ratio_mean_10",
#                  "feature_mid_disp_ratio_mean_10",
#                  "feature_high_nbc_nn_nb_sd_ratio",
#                  "feature_low_nbc_nn_nb_sd_ratio",
#                  "feature_mid_nbc_nn_nb_sd_ratio",
#                  "feature_high_nbc_nn_nb_mean_ratio",
#                  "feature_low_nbc_nn_nb_mean_ratio",
#                  "feature_mid_nbc_nn_nb_mean_ratio",
#                  "feature_high_pca_expl_var_cov_init",
#                  "feature_low_pca_expl_var_cov_init",
#                  "feature_mid_pca_expl_var_cov_init")
# 
# instancesOrdered <- instancesOrdered[c("instances", useFeatures)]
# results <- data.frame(matrix(nrow = 0, ncol = 0))
# numResults <- 0
# # Note this range should be modified if the features used is different;
# # if using more features, use larger epsilon values, if using less features
# # use smaller epsilon values
# epsRange <- c(10, 9.5, 9, 8.5, 8, 7.5, 7, 6.5, 6, 5.5, 5, 4.5, 4, 3.5, 3)
# for(eps in epsRange){
#   print(eps)
#   numResults <- numResults + 1
#   results[numResults, "epsilon"] <- eps
#   Xfeat <- instancesOrdered[, str_which(colnames(instancesOrdered), "feature_")]
#   Xfeat <- data.matrix(Xfeat)
#   labels <- purifyInstances(Xfeat, eps)
#   results[numResults, "InstancesRatio"] <- sum(!labels$preclude) / nrow(labels)
#   results[numResults, "InstancesNumber"] <- sum(!labels$preclude)
#   Xfeat <- Xfeat[!labels$preclude, ]
#   # Now calculate uniformity
#   vals <- calculateCVNND(Xfeat)
#   results[numResults, c("CV", "UniformVal")] <- vals
#   print(results)
#   # Save results
#   write.table(labels, paste0("data/features/filtering/filteringLabelsEps", eps, ".txt"), quote = FALSE, row.names = FALSE)
#   # Also save the data
#   write.table(results, "data/features/instancePurificationResults.txt", quote = FALSE, row.names = FALSE)
# }
# # What follows is code that just read in all of the processed data, so
# # a new user can see the different datasets and choose their own
# # data subset
# for(i in 7:12){
#   eps <- results[i, "epsilon"]
#   labels <- read.table(paste0("data/features/filtering/filteringLabelsEps", eps, ".txt"), header = TRUE, sep = " ")
#   finalInstances <- instancesOrdered[!labels$preclude,]
#   # Write final instances to file!
#   finalInstancesOriginal <- allFeatures[allFeatures$instances %in% finalInstances$instances, ]
#   finalInstancesLit <- allFeatures[allFeatures$instances %in% instancesLit$instances, ]
#   finalInstancesLit <- finalInstancesLit[finalInstancesLit$instances %in% finalInstancesOriginal$instances, ]
#   print(eps)
#   for(dim in 1:20){
#     if(sum(finalInstancesOriginal$feature_dimension == dim) == 0){next}
#     print(paste0(dim, " ", sum(finalInstancesOriginal$feature_dimension == dim), " ", sum(finalInstancesLit$feature_dimension == dim)))
#   }
#   write.table(finalInstances$instances, paste0("data/availableFunctions/chosenTestSuiteN", nrow(finalInstances), ".txt"), quote = FALSE, col.names = FALSE, row.names = FALSE)
# }
# # INSTANCE FILTERING COMPLETE

