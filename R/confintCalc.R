#' Power simulations for cluster-randomized trials: get exact confidence intervals.
#'
#' This function uses binom.test() to get exact
#' confidence intervals for the power in Monte Carlo power estimation.  
#' It is called by some of the clusterPower functions. 
#' 
#' @author Alexandria C. Sakrejda (\email{acbro0@@umass.edu} and Ken Kleinman (\email{ken.kleinman@@gmail.com})
#' 
#' @param alpha A numeric; the user-selected alpha cutoff.
#' @param p.val A vector of p-values.
#' 
#' @return A dataframe
#' \describe{
#'   \item{power.parms}{Exact confidence intervals produced using \code{binom.test()}}
#' }
#' 
#' @export confintCalc
confintCalc <- function(alpha = alpha,
                         p.val = p.val) {
  sig.val <-  ifelse(p.val < alpha, 1, 0)
  if (isTRUE(is.data.frame(sig.val)) || isTRUE(is.matrix(sig.val))) {
    pval.power <- apply(sig.val, 2, FUN = sum)
  }
  if (isTRUE(is.vector(sig.val))) {
    pval.power <- sum(sig.val)
  }
  power.parms <- list()

  for (q in 1:length(unlist(pval.power))) {
    power.parms[[q]] <- binom.test(p = 0.05, n = length(unlist(p.val)), x = pval.power[q], 
                                   alternative = "two.sided")
  }
  Power <- vector(length = length(pval.power))
  Lower.95.CI <- vector(length = length(pval.power))
  Upper.95.CI <- vector(length = length(pval.power))
  Alpha <- vector(length = length(pval.power))
  Beta <- vector(length = length(pval.power))
  Names <- vector(length = length(pval.power))
  for (o in 1:length(pval.power)){
    Power[o] = power.parms[[o]]$estimate
    Lower.95.CI[o] = power.parms[[o]]$conf.int[1]
    Upper.95.CI[o] = power.parms[[o]]$conf.int[2]
    Alpha[o] = alpha
    Names[o] = paste0("Arm.", o + 1)
  }
  Beta <- 1 - Power
  power.parms <- data.frame(Power, Lower.95.CI, Upper.95.CI, Alpha, Beta)
  rownames(power.parms) <- Names
  return(power.parms)
}