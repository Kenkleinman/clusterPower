#' Simulation-based power estimation for continuous outcome individually 
#' randomized group treatment trials.
#'
#' This function uses iterative simulations to determine 
#' approximate power for individually randomized group treatment trials with a 
#' normally-distributed outcome of interest. Users can modify a variety of 
#' parameters to suit the simulations to their desired experimental situation. 
#' This function returns the summary power values for each treatment arm.
#' 
#' Runs the power simulation.
#' 
#' Users must specify the desired number of simulations, number of subjects per 
#' cluster, number of clusters per treatment arm, expected absolute difference 
#' between treatments, two of the following: ICC, within-cluster variance, or 
#' between-cluster variance; significance level, analytic method, progress updates, 
#' and simulated data set output may also be specified.
#' 
#' 
#' @param nsim Number of datasets to simulate; accepts integer (required).
#' @param nsubjects Number of subjects per treatment group; accepts either a scalar (equal cluster sizes, both groups), 
#' a vector of length two (equal cluster sizes within groups), or a vector of length \code{sum(nclusters)} 
#' (unequal cluster sizes within groups) (required).
#' @param nclusters Number of clusters per group; accepts single integer or vector of length 2 for unequal number 
#' of clusters per treatment group (required)
#' @param difference Expected absolute treatment effect; accepts numeric (required).
#' At least 2 of the following must be specified:
#' @param ICC Intra-cluster correlation coefficient; accepts a value between 0 - 1
#' @param sigma Within-cluster variance; accepts numeric
#' @param sigma_b Between-cluster variance; accepts numeric
#' If clusters differ between treatment groups, at least 1 of the following 
#' must be specified: ICC2, sigma2.
#' @param ICC2 Intra-cluster correlation coefficient for clusters in TREATMENT group
#' @param sigma2 Within-cluster variance for clusters in TREATMENT group
#' @param sigma_b2 Between-cluster variance for clusters in TREATMENT group
#' @param alpha Significance level; default = 0.05.
#' @param quiet When set to FALSE, displays simulation progress and estimated completion time; default is FALSE.
#' @param all.sim.data Option to output list of all simulated datasets; default = FALSE.
#' 
#' @return A list with the following components:
#' \itemize{
#'   \item Character string indicating total number of simulations and simulation type
#'   \item Number of simulations
#'   \item Data frame with columns "Power" (Estimated statistical power), 
#'                "lower.95.ci" (Lower 95% confidence interval bound), 
#'                "upper.95.ci" (Upper 95% confidence interval bound)
#'   \item Analytic method used for power estimation
#'   \item Significance level
#'   \item Vector containing user-defined cluster sizes
#'   \item Vector containing user-defined number of clusters in each treatment group
#'   \item Data frame reporting ICC for Treatment/Non-Treatment groups
#'   \item Vector containing expected difference between groups based on user inputs
#'   \item Data frame with columns: 
#'                   "Estimate" (Estimate of treatment effect for a given simulation), 
#'                   "Std.err" (Standard error for treatment effect estimate), 
#'                   "Test.statistic" (z-value (for GLMM) or Wald statistic (for GEE)), 
#'                   "p.value", 
#'                   "sig.val" (Is p-value less than alpha?)
#'   \item List of data frames, each containing: 
#'                   "y" (Simulated response value), 
#'                   "trt" (Indicator for treatment group), 
#'                   "clust" (Indicator for cluster)
#'                   }
#' 
#' 
#' @examples 
#' \dontrun{
#' irgtt.normal.sim <- cps.irgtt.normal(nsim = 100, nsubjects = c(100, 10), 
#'                        nclusters = 10, difference = 16,
#'                         sigma = 100, sigma_b2 = 25, alpha = 0.05,
#'                         quiet = FALSE, all.sim.data = FALSE)
#' }
#' @author Alexandria C. Sakrejda (\email{acbro0@@umass.edu}), Alexander R. Bogdan, 
#'   and Ken Kleinman (\email{ken.kleinman@@gmail.com})
#' @export


cps.irgtt.normal <-  function(nsim = NULL, nsubjects = NULL, nclusters = NULL, difference = NULL,
                      ICC = NULL, sigma = NULL, sigma_b = 0,
                      ICC2 = ICC, sigma2 = sigma, sigma_b2 = 0,
                      alpha = 0.05, method = 'glmm', quiet = FALSE,
                      all.sim.data = FALSE){
  if (sigma_b == 0 & sigma_b2 == 0){
    warning("Sigma_b in both arms is 0. This is equivalent to a t-test. Did you mean to 
            enter a sigma_b value for the arm containing clustered observations?")
  }
  if (sigma_b != 0 & sigma_b2 != 0){
    warning("Sigma_b is not zero for either arm. Did you want to use cps.normal()?")
  }
  if (sigma_b != 0 & sigma_b2 == 0){
    stop("Non-clustered group must be entered as the reference group.")
    nclust <- c(nclusters, 1)
  }
  if (sigma_b2 != 0){
    nclust <- c(1, nclusters)
  }
  
  sim = cps.normal(nsim = nsim, nsubjects = nsubjects, nclusters = nclust, 
                   difference = difference, ICC = ICC, ICC2 = ICC2, 
                   sigma = sigma, sigma2 = sigma2, alpha = alpha, 
                   sigma_b = sigma_b, sigma_b2 = sigma_b2,
                   method = "glmm", quiet = quiet, all.sim.data = all.sim.data,
                   irgtt = TRUE)
  return(sim)
}
  
  
  