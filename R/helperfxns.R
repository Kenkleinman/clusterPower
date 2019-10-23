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
        warning("ClusterPower uses maximum likelihood inference, which depends on asymptotic results.  These results hove known limitations in finite settings, such as when there are small numbers of clusters or small numbers of observation per cluster.  In the parameters you have chosen, the power estimate may be conservative.")
      }
    }
  }
}

#calculate the confidence intervals using binomial exact test
confint.calc <- function(nsim = nsim, alpha = alpha,
  p.val = p.val, names.power = names.power){
  sig.val <-  ifelse(p.val < alpha, 1, 0)
  pval.power <- apply (sig.val, 2, FUN = sum)
  power.parms <- list()
  for (q in 1:length(pval.power)){
    power.parms[[q]] <- binom.test(p = 0.05, n = nsim, x = pval.power[q], 
               alternative = "greater")
  }
  Power <- vector(length = length(pval.power))
  Lower.95.CI <- vector(length = length(pval.power))
  Upper.95.CI <- vector(length = length(pval.power))
  for (o in 1:length(pval.power)){
    Power[o] = power.parms[[o]]$estimate
    Lower.95.CI[o] = power.parms[[o]]$conf.int[1]
    Upper.95.CI[o] = power.parms[[o]]$conf.int[2]
  }
  power.parms <- data.frame(Power, Lower.95.CI, Upper.95.CI)
  rownames(power.parms) <- names.power
  return(power.parms)
}

   
# Create wholenumber function
is.wholenumber = function(x, tol = .Machine$double.eps^0.5)  abs(x - round(x)) < tol

# create proportion of F-test rejections fxn
prop_H0_rejection <- function (alpha = alpha, nsim = nsim, 
                               LRT.holder.abbrev = LRT.holder.abbrev){
f.test <- binom.test(p = 0.05, n = nsim, x = LRT.holder.abbrev, 
                               alternative = "greater")
  Power = f.test$estimate
  Lower.95.CI = f.test$conf.int[1]
  Upper.95.CI = f.test$conf.int[2]
  Ftest <- data.frame(Power, Lower.95.CI, Upper.95.CI)
  return(Ftest)
}
