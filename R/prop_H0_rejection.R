#' Power simulations for cluster-randomized trials: proportion of F-test rejections
#'
#' This function can be called by some of the clusterPower functions to 
#' calculate the proportion of F-test rejections and produce exact confidence intervals.
#' 
#' @author Alexandria C. Sakrejda (\email{acbro0@@umass.edu} and Ken Kleinman (\email{ken.kleinman@@gmail.com})
#' 
#' @param nsim A scalar; the number of simulations chosen by the user.
#' @param alpha A numeric; the user-selected alpha cutoff.
#' @param sig.LRT A logical vector indicating whether or not the model was significant.
#' 
#' @return A dataframe
#' \describe{
#'   \item{Ftest}{Exact confidence intervals produced using \code{binom.test()}}
#' }
#' 
#' @export

prop_H0_rejection <- function (alpha = alpha,
                               nsim = nsim,
                               sig.LRT = sig.LRT) {
  # Proportion of times P(>F)
  LRT.holder.abbrev <- sum(sig.LRT)
  f.test <- binom.test(p = 0.05, n = length(sig.LRT), x = LRT.holder.abbrev)
  Power = f.test$estimate
  Lower.95.CI = f.test$conf.int[1]
  Upper.95.CI = f.test$conf.int[2]
  Beta <- 1 - Power
  Alpha <- alpha
  Ftest <- data.frame(Power, Lower.95.CI, Upper.95.CI, Alpha, Beta)
  num.returned <- data.frame("Converged" = length(sig.LRT), 
                             "Requested" = nsim)
  Ftest <- cbind(Ftest, num.returned)
  return(Ftest)
}

