# Modification to code in flacco repo, 
# https://github.com/kerschke/flacco/blob/master/R/feature_misc_principal_component.R
# Idea is to modify it so that it can be used when the objective function values are 0,
# to get a value instead of a crash
calculatePrincipalComponentFeaturesWithConstantY = function(feat.object) {
  # assertClass(feat.object, "FeatureObject")
  # assertList(control)
  measureTime(expression({
    # prop.cov_x = control_parameter(control, "pca.cov_x", 0.9)
    # prop.cor_x = control_parameter(control, "pca.cor_x", 0.9)
    # prop.cov_init = control_parameter(control, "pca.cov_init", 0.9)
    # prop.cor_init = control_parameter(control, "pca.cor_init", 0.9)
    # init = feat.object$env$init
    X = as.matrix(subset(feat.object$env$init, select = feat.object$feature.names))
    d = feat.object$dim
    
    explainVariance = function(data, cov = TRUE) {
      if (cov) {
        ev = eigen(cov(data))$values
      } else {
        ev = eigen(cor(data))$values
      }
      cumsum(ev) / sum(ev)
    }
    
    cov_x = explainVariance(X, cov = TRUE)
    cor_x = explainVariance(X, cov = FALSE)
    # cov_init = explainVariance(init, cov = TRUE)
    # cor_init = explainVariance(init, cov = FALSE)
    
    return(list(
      # pca.expl_var.cov_x = min(which(cov_x >= prop.cov_x)) / d,
      # pca.expl_var.cor_x = min(which(cor_x >= prop.cor_x)) / d,
      # pca.expl_var.cov_init = min(which(cov_init >= prop.cov_init)) / (d + 1),
      # pca.expl_var.cor_init = min(which(cor_init >= prop.cor_init)) / (d + 1),
      pca.expl_var.cov_x = min(which(cov_x >= 0.9)) / d,
      pca.expl_var.cor_x = min(which(cor_x >= 0.9)) / d,
      pca.expl_var.cov_init = 1 / (d + 1),
      pca.expl_var.cor_init = 1 / (d + 1),
      pca.expl_var_PC1.cov_x = cov_x[1],
      pca.expl_var_PC1.cor_x = cor_x[1],
      # pca.expl_var_PC1.cov_init = cov_init[1],
      # pca.expl_var_PC1.cor_init = cor_init[1]
      pca.expl_var_PC1.cov_init = 1,
      pca.expl_var_PC1.cor_init = 1
      ))
  }), "pca")
}



calculateDistributionFeaturesCustom = function(feat.object) {
  measureTime(expression({
    y = feat.object$env$init[, feat.object$objective.name]
    # NICO: Added the case when the objective function is flat
    if(min(y) == max(y)){
      n <- length(y)
      skewness = (n - 1)^1.5/n
      kurtosis = (n - 1)^2/n - 3
      number_of_peaks = 1
    }else{
      skewness = e1071::skewness(y, type = 3L)
      kurtosis = e1071::kurtosis(y, type = 3L)
      number_of_peaks = number_of_peaksCustom(y)
    }
    list(ela_distr.skewness = skewness,
        ela_distr.kurtosis = kurtosis,
        ela_distr.number_of_peaks = number_of_peaks)
  }), "ela_distr")
}


number_of_peaksCustom = function(x, smoothing.bandwidth = "SJ", 
                           modemass.threshold = 0.01) {
  intdens = function(a, b) {
    mean(y[a:b]) * diff(d$x[c(a, b)])
  }
  # NICO: Using this setting instead of "SJ" seems to work
  d = density(x, bw = "nrd0")
  y = d$y
  n = length(y)
  index = 2L : (n - 1L)
  # NICO: Changed from 
  # min.index = c(1L, which((y[index] < y[index - 1L]) & y[index] < y[index + 1L]), n + 1L)
  min.index = c(1L, 1 + which((y[index] < y[index - 1L]) & y[index] < y[index + 1L]), n + 1L)
  modemass = vapply(1L : (length(min.index) - 1L), function(i) 
    intdens(min.index[i], min.index[i + 1L] - 1L), double(1))
  sum(modemass > modemass.threshold)
}



calculateLevelsetFeaturesWithConstantY = function(feat.object) {
  X = as.matrix(subset(feat.object$env$init, select = feat.object$feature.names))
  y = feat.object$env$init[, feat.object$objective.name]
  measureTime(expression({
    
    result <- c()
    names <- c()
    for(method in c("lda", "qda", "mda")){
      for(quantile in c("10", "25", "50")){
        result <- c(result, 0)
        names <- c(names, paste0("ela_level.mmce_", method, "_", quantile))
      }
    }
    for(method1 in c("lda", "qda", "mda")){
      for(method2 in c("lda", "qda", "mda")){
        if(method1 == method2){next}
        for(quantile in c("10", "25", "50")){
          result <- c(result, 1)
          names <- c(names, paste0("ela_level.", method1, "_", method2, "_", quantile))
        }
      }
    }
    result = as.vector(result, mode = "list")
    names(result) = names
    return(result)
  }), "ela_level")
}





calculateMetaModelFeaturesCustom = function(feat.object) {
  measureTime(expression({
    X = as.matrix(subset(feat.object$env$init, select = feat.object$feature.names))
    y = feat.object$env$init[, feat.object$objective.name]
    df = as.data.frame(X)
    # NICO: If values are constant, all adj_r2 = 1 as any of the models should predict exactly, set all the values
    if(min(y) == max(y)){
      res = list(ela_meta.lin_simple.adj_r2 = 1,
                 ela_meta.lin_simple.intercept = min(y),
                 ela_meta.lin_simple.coef.min = 0,
                 ela_meta.lin_simple.coef.max = 0,
                 ela_meta.lin_simple.coef.max_by_min = 1,
                 ela_meta.lin_w_interact.adj_r2 = 1,
                 ela_meta.quad_simple.adj_r2 = 1,
                 ela_meta.quad_simple.cond = 1,
                 ela_meta.quad_w_interact.adj_r2 = 1)
    }else{
      
      ## Simple linear model:
      lin.model = lm(y ~ ., data = df)
      model.coeff = coef(lin.model)
      # NICO: Slight modification which ignores NA coefficients in case of overfitting
      res = list(ela_meta.lin_simple.adj_r2 = calculateAdjustedR2Custom(lin.model),
                 ela_meta.lin_simple.intercept = as.numeric(model.coeff[1L]),
                 ela_meta.lin_simple.coef.min = min(abs(model.coeff[!is.na(model.coeff)][-1L])),
                 ela_meta.lin_simple.coef.max = max(abs(model.coeff[!is.na(model.coeff)][-1L])),
                 ela_meta.lin_simple.coef.max_by_min = max(abs(model.coeff[!is.na(model.coeff)][-1L])) / min(abs(model.coeff[!is.na(model.coeff)][-1L]))
      )
      ## Linear interactions:
      lin.model = lm(y ~ .^2, data = df) ## Include interactions
      res$ela_meta.lin_w_interact.adj_r2 = calculateAdjustedR2Custom(lin.model)
      
      ## Simple quadratic model:
      cns = names(df)
      df = cbind(df, df^2)
      cns.squared = sprintf("%s_squared", cns)
      names(df) = c(cns, cns.squared)
      quad.model = lm(y ~ ., data = df)
      res$ela_meta.quad_simple.adj_r2 = calculateAdjustedR2Custom(quad.model)
      # NICO: If all quadratic coefficients are NA (i.e. overfitting), set this to 1 (as all coefficients are 0)
      quadraticCoeff <- coef(quad.model)[cns.squared]
      quadraticCoeff <- quadraticCoeff[!is.na(quadraticCoeff)]
      quadraticCoeff <- quadraticCoeff[quadraticCoeff != 0]
      
      if(length(quadraticCoeff) == 0){
        res$ela_meta.quad_simple.cond <- 1
      }else{
        quad.model_cond = range(abs(quadraticCoeff))
        res$ela_meta.quad_simple.cond = quad.model_cond[2L] / quad.model_cond[1L]
      }
      
      ## Quadratic interactions:
      quad.model_matrix = model.matrix(~ .^2, data = df)
      quad.model = lm.fit(quad.model_matrix, y)
      res$ela_meta.quad_w_interact.adj_r2 = calculateAdjustedR2Custom(quad.model)
    }
    res
  }), "ela_meta")
}

## Calculate adjusted R^2
calculateAdjustedR2Custom = function(mod) {
  pred = fitted(mod)
  resi = residuals(mod)
  SS_reg = crossprod(pred - mean(pred))
  SS_res = crossprod(resi)
  SS_total = SS_reg + SS_res
  n = length(resi)
  print("Info")
  print(SS_reg)
  print(SS_res)
  print(SS_total)
  print(n)
  # NICO: Removing code below as it should only happen with a perfect fit,
  # which should be simplified to adjusted R2 of 1
  # NICO: Some of the coefficients might be 0, should not count them in this
  # coeffs <- mod$coefficients
  # coeffs <- coeffs[-1]
  # coeffs <- coeffs[coeffs != 0]
  # p = length(coeffs) - 1L
  
  if(sum(resi) == 0){return(1)}
  p = length(mod$coefficients) - 1L
  print(p)
  print(drop(1 - (SS_res / SS_total) / ((n - p - 1) / (n - 1))))
  print((SS_res / SS_total) / ((n - p - 1) / (n - 1)))
  print("Done")
  
  drop(1 - (SS_res / SS_total) / ((n - p - 1) / (n - 1)))
}





calculateDispersionFeaturesCustom = function(feat.object) {
  measureTime(expression({
    X = as.matrix(subset(feat.object$env$init, select = feat.object$feature.names))
    y = feat.object$env$init[, feat.object$objective.name]
    quantiles = c(0.02, 0.05, 0.1, 0.25)
    index = lapply(quantile(y, quantiles), 
                   function(quant) which(y <= quant))
    
    dists = lapply(seq_along(index), 
                   function(i) as.numeric(dist(X[index[[i]], ], method = "euclidean")))  
    
    dists.full_sample = as.numeric(dist(X, method = "euclidean"))
    means = vapply(dists, mean, double(1))
    medians = vapply(dists, median, double(1))
    # NICO: If there is only a single point in the quantile, the dists will return
    # and empty vector; really the mean and median distance is 0 (as there is only 1)
    means[is.na(means)] <- 0
    medians[is.na(medians)] <- 0
    
    res = c(means / mean(dists.full_sample), medians / median(dists.full_sample),
            means - mean(dists.full_sample), medians - median(dists.full_sample))
    res = as.vector(res, mode = "list")
    names(res) = c(sprintf("disp.ratio_mean_%02i", quantiles * 100), 
                   sprintf("disp.ratio_median_%02i", quantiles * 100),
                   sprintf("disp.diff_mean_%02i", quantiles * 100), 
                   sprintf("disp.diff_median_%02i", quantiles * 100))
    return(res)
  }), "disp")
}






# NICO: Just forcing the feature to calculate using the "slow" version, as
# the "fast" version gives an error.
calculateNearestBetterFeaturesCustom = function(feat.object) {
  measureTime(expression({
    meth = "euclidean"
    cor_na = "pairwise.complete.obs"
    init = feat.object$env$init
    X = init[, feat.object$feature.names]
    y = init[, feat.object$objective.name]
    distmat = as.matrix(dist(X, method = meth, diag = TRUE, upper = TRUE))  
    nb.stats = computeNearestBetterStatsCustom(distmat = distmat, objectives = y)
    nn.dists = nb.stats$nearDist
    nb.dists = nb.stats$nbDist
    # cure global optima
    nb.dists[is.na(nb.stats$nbID)] = nn.dists[is.na(nb.stats$nbID)]
    dist_ratio = nn.dists / nb.dists
    
    # NICO: Corr will not be calculated if values are constant, set to 1
    if(min(y) == max(y)){
      corrVal = 1
    }else{
      corrVal = cor(nb.stats$toMe_count, y, use = cor_na)
    }
    
    return(list(
      nbc.nn_nb.sd_ratio = sd(nn.dists, na.rm = TRUE) / sd(nb.dists, na.rm = TRUE),
      nbc.nn_nb.mean_ratio = mean(nn.dists, na.rm = TRUE) / mean(nb.dists, na.rm = TRUE),
      nbc.nn_nb.cor = cor(nn.dists, nb.dists, use = cor_na),
      nbc.dist_ratio.coeff_var = 
        sd(dist_ratio, na.rm = TRUE) / mean(dist_ratio, na.rm = TRUE),
      nbc.nb_fitness.cor = corrVal
    ))
  }), "nbc")

}

# compute various distance measures and ratios wrt the nearest and nearest
# better elements:
computeNearestBetterStatsCustom = function(distmat, objectives) {
  result = data.frame(ownID = BBmisc::seq_row(distmat))
  result = cbind(result, t(vapply(result$ownID, function(row) {
    rowDists = as.numeric(distmat[row, ])
    # first look for elements with better fitness
    better = which(objectives < objectives[row])
    # if no better elements are available check for elements with equal fitness
    if (length(better) == 0L) {
      better = which(objectives == objectives[row])
      better = setdiff(better, row)
    }
    # select the nearest among the (equal-or-)better elements
    if (length(better) > 0L) {
      nb = better[selectMin(rowDists[better])]
      return(c(nbID = nb, nbDist = rowDists[nb], nearDist = min(rowDists[-row])))
    } else {
      return(c(nbID = NA_real_, nbDist = NA_real_, nearDist = min(rowDists[-row])))
    }
  }, double(3))))
  # compute ratio of nearestBetter and nearest
  result$nb_near_ratio = result$nbDist / result$nearDist
  result$fitness = objectives
  result = cbind(result, t(vapply(result$ownID, function(row) {
    x = which(result$nbID == result$ownID[row])
    # number of elements, which have ownID[row] as nearest better
    count = length(x)
    # median distance to the elements, to which ownID[row] is nearest better
    toMe_dist = median(result$nbDist[x])
    if (count > 0L) {
      return(c(toMe_count = count, toMe_dist_median = toMe_dist,
               nb_median_toMe_ratio = result$nbDist[row] / toMe_dist))
    } else {
      return(c(toMe_count = 0, toMe_dist_median = NA_real_, 
               nb_median_toMe_ratio = NA_real_))
    }
  }, double(3))))
  return(result)
}

selectMin = function(x, tie.breaker = "sample") {
  i = which(x == min(x))
  if (length(i) > 1L) {
    if (tie.breaker == "first") {
      return(i[1L])
    } else if (tie.breaker == "last") {
      return(i[length(i)])
    } else if (tie.breaker == "sample") {
      return(sample(i, 1L))
    }
  }
  return(i)
}





calculateInformationContentFeaturesCustom = function(feat.object) {
  measureTime(expression({
    res = computeInfoContentStatisticsCustom(feat.object)
    # h.max = "maximum information content" - cf. equation (5)
    # eps.s = "settling sensitivity" - cf. equation (6)
    # eps.max (created due to a comment from Mario Andres Munoz), this
    # is the epsilon which holds H(epsilon_max) = h.max
    # eps.ratio = "half partial information sensitivity" - cf. equation (8)
    # M0 = "initial partial information" - cf. equation (7)
    return(list(ic.h.max = res$Hmax,
                ic.eps.s = res$eps.S,
                ic.eps.max = res$eps.max,
                ic.eps.ratio = res$eps05,
                ic.m0 = res$M0
    ))
  }), "ic")
}

computeInfoContentStatisticsCustom = function(feat.object) {
  epsilon <- c(0, 10^(seq(-5, 15, length.out = 1000)))
  epsilon = unique(epsilon)
  X = as.matrix(subset(feat.object$env$init, select = feat.object$feature.names))
  y = feat.object$env$init[, feat.object$objective.name]
  n = feat.object$n.obs
  
  # NICO: Have default values when the 
  
  seed = 1
  set.seed(seed)

  start = sample.int(n, 1)
  hood = 20L
  res = constructSequence(X = X, start = start, hood = hood)
  permutation = res$permutation
  d = res$distance
  
  psi.eps = vapply(epsilon, function(eps) {
    computePsi(permutation = permutation, xdists = d, y = y, eps = eps)
  }, integer(length(permutation) - 1L))
  
  H = apply(psi.eps, 2, computeH)
  M = apply(psi.eps, 2, computeM)
  
  # calculate eps.S, cf. equation (6) ("settling sensitivity")
  settl.sens = 0.05
  if(settl.sens < 0 | settl.sens > .Machine$double.xmax){
    stop("settl.send is larger than the machine's double max!")
  }
  eps.S = epsilon[which(H < settl.sens)]
  
  # NICO: Replacing -Inf (if min(eps.S) = 0) with the second smallest entry
  if(length(eps.S) > 0){
    if(min(eps.S) == 0){
      eps.S <- log10(min(eps.S[eps.S != min(eps.S)]))
    }else{
      eps.S <- log10(min(eps.S))
    }
  }else{
    eps.S <- NA_real_
  }
  
  # calculate M0, cf. equation (7) ("initial partial information")
  M0 = M[epsilon == 0]
  
  # calculate epsilon05 [Eq. (8)] ("half partial information sensitivity")
  inf.sens = 0.5

  eps05 = which(M > inf.sens * M0)
  # NICO: Also replacing eps05 = NaN (if max(eps05) is empty) with the second lowest entry
  if(length(eps05) > 0 & max(epsilon[eps05]) != 0){
    eps05 <- log10(max(epsilon[eps05]))
  }else{
    eps05 <- log10(min(epsilon[epsilon != min(epsilon)]))
  }
  # eps05 = ifelse(length(eps05) > 0, log10(max(epsilon[eps05])), NA_real_)
  
  return(list(H = H, M = M, eps = epsilon, eps05 = eps05, eps.S = eps.S,
              M0 = M0, Hmax = max(H), eps.max = median(epsilon[H == max(H)])))
}

## construct a path through the landscape - starting with an initial
## observation and walking (greedily) from an observation to its nearest
## (not-yet-visited) neighbour;
## the output is a matrix with two columns: the first returns the index
## of the elements (in which they've been visited) and the second one
## returns the distance from neighbour to neighbour
constructSequence = function(X, start, hood) {
  nn.list = RANN::nn2(X, k = min(hood, nrow(X)))
  n = nrow(X)
  # add first candidate (random) and initialise permutation vector (avoids
  # continuous allocation of space)
  if (missing(start))
    current = sample.int(n, 1L)
  else
    current = as.integer(start)
  candidates = seq_len(n)[-current]
  permutation = c(current, rep(NA_integer_, n - 1L))
  dists = rep(NA_real_, n)
  
  # successively add next candidates
  for (i in seq_len(n)[-1L]) {
    currents = nn.list$nn.idx[permutation[i - 1L], ]
    current = intersect(currents, candidates)
    if (length(current) > 0L) {
      current = current[1L]
      permutation[i] = current
      candidates = candidates[-which(candidates == current)]
      dists[i] = nn.list$nn.dists[permutation[i - 1], currents == current]
    } else {
      # list of nearest (yet unvisited) neighbor
      nn.list2 = RANN::nn2(X[candidates, , drop = FALSE],
                           query = X[permutation[i - 1L], , drop = FALSE],
                           k = min(nrow(X), 1L))
      current = as.integer(candidates[nn.list2$nn.idx])
      permutation[i] = current
      candidates = candidates[-which(candidates == current)]
      dists[i] = as.numeric(nn.list2$nn.dists)
    }
  }
  return(list(permutation = permutation, distance = dists[-1L]))
}

## cf. equation (9)
computePsi = function(permutation, xdists, y, eps) {
  y = y[permutation]
  ratio = diff(y) / xdists
  ifelse(abs(ratio) < eps, 0L, as.integer(sign(ratio)))
}

## cf. equation(2)
computeH = function(psi) {
  a = psi[-length(psi)]
  b = psi[-1]
  probs = c(
    #neg_neg = mean((a == -1) & (b == -1)), 
    neg_neu = mean((a == -1) & (b == 0)), 
    neg_pos = mean((a == -1) & (b == 1)),
    neu_neg = mean((a == 0) & (b == -1)),
    #neu_neu = mean((a == 0) & (b == 0)),
    neu_pos = mean((a == 0) & (b == 1)),
    pos_neg = mean((a == 1) & (b == -1)),
    pos_neu = mean((a == 1) & (b == 0))#,
    #pos_pos = mean((a == 1) & (b == 1))
  )
  -sum(ifelse(probs == 0, 0, probs * log(probs, base = 6)))
}

## cf. equation(3)
computeM = function(psi) {
  n = length(psi)
  psi = psi[psi != 0]
  psi = psi[c(FALSE, diff(psi) != 0)]
  length(psi) / (n - 1)
}







