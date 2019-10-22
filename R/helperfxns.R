type1ErrTest <- function(sigma_sq_, sigma_b_sq_, nsubjects_){
  ICC <- createMissingVarianceParam(sigma_sq = sigma_sq_,
    sigma_b_sq = sigma_b_sq_)
  nobs <- nsubjects_
  clusters <- unlist(lapply(nsubjects_, length))
  for (i in 1:length(nobs)){
    for (j in 1:length(nobs[[i]])){
      nobstemp <- nobs[[i]][[j]]
      ICCtemp <- ICC[i]
      clusterstemp <- clusters[i]  
      if ((nobstemp < 5 & ICCtemp < 0.1 & clusterstemp < 20) | 
      (nobstemp < 10  & ICCtemp < 0.05 & clusterstemp < 20) |
      (nobstemp < 20  & ICCtemp < 0.02 & clusterstemp < 20) |
      (nobstemp < 50  & ICCtemp < 0.01 & clusterstemp < 10)){
        warning("WARNING: Type 1 error rate increases when ICC and nclusters are small. True power may be lower than estimates.")
      }
    }
  }
}

