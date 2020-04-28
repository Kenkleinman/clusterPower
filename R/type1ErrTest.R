#' Power simulations for cluster-randomized trials: trigger a conservative estimate warning.
#'
#' ClusterPower uses maximum likelihood inference, which depends on 
#' asymptotic results. These results hove known limitations in finite 
#' settings, such as when there are small numbers of clusters or small 
#' numbers of observation per cluster. 
#' 
#' @author Alexandria C. Sakrejda (\email{acbro0@@umass.edu} and Ken Kleinman (\email{ken.kleinman@@gmail.com})
#' 
#' @param sigma_sq_ User-selected sigma_sq value.
#' @param sigma_b_sq_ User-selected sigma_b_sq value.
#' @param nsubjects_ User-selected number of observations.
#' 
#' @return If triggered, a warning. Otherwise silent.
#' \describe{
#'   \item{warning}{With the parameters you entered, the power estimate may be conservative.}
#' }
#' 
#' @export


type1ErrTest <- function(sigma_sq_, sigma_b_sq_, nsubjects_) {
  ICC <- createMissingVarianceParam(sigma_sq = sigma_sq_,
                                    sigma_b_sq = sigma_b_sq_)
  nobs <- nsubjects_
  alert <- NA
  clusters <- unlist(lapply(nsubjects_, length))
  for (i in 1:length(nobs)) {
    for (j in 1:length(nobs[[i]])) {
      nobstemp <- nobs[[i]][[j]]
      ICCtemp <- ICC[i]
      clusterstemp <- clusters[i]
      alert[j] <-
        as.vector(ifelse((nobstemp < 5 &
                            ICCtemp < 0.1 & clusterstemp < 20) |
                           (nobstemp < 10 & ICCtemp < 0.05 & clusterstemp < 20) |
                           (nobstemp < 20 & ICCtemp < 0.02 & clusterstemp < 20) |
                           (nobstemp < 50 &
                              ICCtemp < 0.01 & clusterstemp < 10),
                         TRUE,
                         FALSE
        ))
    }
  }
  if (sum(alert) > 0) {
    warning(
      "ClusterPower uses maximum likelihood inference, which depends on asymptotic results.
These results have known limitations in finite settings, such as when there are small numbers of clusters or small numbers of observation per cluster.
With the parameters you entered, the power estimate may be conservative."
    )
  }
}