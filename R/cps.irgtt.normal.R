#' Simulation-based power estimation for continuous outcome individually 
#' randomized group treatment trials.
#'
#' This function uses iterative simulations to determine 
#' approximate power for individually randomized group treatment trials with a 
#' normally-distributed outcome of interest. Users can modify a variety of 
#' parameters to suit the simulations to their desired experimental situation. 
#' This function returns the summary power values for each arm.
#' 
#' Runs the power simulation.
#' 
#' Users must specify the desired number of simulations, number of subjects per 
#' cluster, number of clusters per arm, expected means for the 
#' arm 1 and arm 2 (respectively), two of the following: ICC, 
#' within-cluster variance, or between-cluster variance; significance level, 
#' analytic method, progress updates, and simulated data set output may also be 
#' specified.
#' 
#' 
#' @param nsim Number of datasets to simulate; accepts integer (required).
#' @param nsubjects Number of subjects per cluster in each arm; accepts either a scalar (equal cluster sizes, both groups), 
#' a vector of length two (equal cluster sizes within groups), or a vector of length \code{sum(nclusters)} 
#' (unequal cluster sizes within groups) (required).
#' @param nclusters Number of clusters in the clustered group; accepts a scalar (required)
#' @param mu Expected mean of arm 1; accepts numeric (required).
#' @param mu2 Expected mean of arm 2; accepts numeric (required).
#' @param alpha Significance level; default = 0.05.
#' @param quiet When set to FALSE, displays simulation progress and estimated completion time; default is FALSE.
#' @param allSimData Option to output list of all simulated datasets; default = FALSE.
#' @param seed Option to set seed. Default is NA.
#' @param poorFitOverride Option to override \code{stop()} if more than 25\%
#' of fits fail to converge; default = FALSE.
#' @param lowPowerOverride Option to override \code{stop()} if the power
#' is less than 0.5 after the first 50 simulations and every ten simulations
#' thereafter. On function execution stop, the actual power is printed in the
#' stop message. Default = FALSE. When TRUE, this check is ignored and the
#' calculated power is returned regardless of value.
#' @param timelimitOverride Logical. When FALSE, stops execution if the estimated completion time
#' is more than 2 minutes. Defaults to TRUE.
#' 
#' At least 2 of the following must be specified:
#' @param ICC Intra-cluster correlation coefficient; accepts a value between 0 - 1
#' @param sigma_sq Within-cluster variance; accepts numeric
#' @param sigma_b_sq Between-cluster variance; defaults to 0. Accepts numeric.
#' If clusters differ between arms, at least 1 of the following 
#' must be specified: ICC2, sigma_sq2.
#' @param ICC2 Intra-cluster correlation coefficient for clusters in arm 2
#' @param sigma_sq2 Within-cluster variance for clusters in arm 2
#' @param sigma_b_sq2 Between-cluster variance for clusters in arm 2.
#' @param nofit Option to skip model fitting and analysis and return the simulated data.
#' Defaults to \code{FALSE}.
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
#'   \item Vector containing expected group means based on user inputs
#'   \item Data frame with columns: 
#'                   "Estimate" (Estimate of treatment effect for a given simulation), 
#'                   "Std.err" (Standard error for treatment effect estimate), 
#'                   "Test.statistic" (z-value (for GLMM) or Wald statistic (for GEE)), 
#'                   "p.value", 
#'                   "sig.val" (Is p-value less than alpha?)
#'   \item If \code{allSimData = TRUE}, a list of data frames, each containing: 
#'                   "y" (Simulated response value), 
#'                   "trt" (Indicator for treatment group), 
#'                   "clust" (Indicator for cluster)
#'                   }
#' If \code{nofit = T}, a data frame of the simulated data sets, containing:
#' \itemize{
#'   \item "arm" (Indicator for treatment arm)
#'   \item "cluster" (Indicator for cluster)
#'   \item "y1" ... "yn" (Simulated response value for each of the \code{nsim} data sets).
#'   }
#' 
#' @examples 
#' \dontrun{
#' irgtt.normal.sim <- cps.irgtt.normal(nsim = 100, nsubjects = c(100, 10), 
#'                        nclusters = 8, mu = 1.1, mu2 = 1.5,
#'                         sigma_sq = 0.1, sigma_sq2 = 0.2, 
#'                         sigma_b_sq2 = 0.1, alpha = 0.05,
#'                         quiet = FALSE, allSimData = FALSE)
#' }
#' @author Alexandria C. Sakrejda (\email{acbro0@@umass.edu})
#' @author Alexander R. Bogdan
#' @author Ken Kleinman (\email{ken.kleinman@@gmail.com})
#' 
#' @export



cps.irgtt.normal <-
  function(nsim = NULL,
           nsubjects = NULL,
           nclusters = NULL,
           mu = NULL,
           mu2 = NULL,
           ICC = NULL,
           sigma_sq = NULL,
           sigma_b_sq = 0,
           ICC2 = ICC,
           sigma_sq2 = sigma_sq,
           sigma_b_sq2 = 0,
           alpha = 0.05,
           quiet = FALSE,
           allSimData = FALSE,
           nofit = FALSE,
           seed = NA,
           poorFitOverride = FALSE,
           lowPowerOverride = FALSE,
           timelimitOverride = TRUE) {
    if (sigma_b_sq == 0 & sigma_b_sq2 == 0) {
      warning(
        "sigma_b_sq in both arms is 0. This is equivalent to a t-test. Did you mean to enter a sigma_b_sq value for the arm containing clustered observations?"
      )
    }
    if (sigma_b_sq != 0 & sigma_b_sq2 != 0) {
      warning("sigma_b_sq is not zero for either arm. Did you want to use cps.normal()?")
    }
    if (sigma_b_sq != 0 & sigma_b_sq2 == 0) {
      stop("Non-clustered group must be the reference group.")
    }
    if (sigma_b_sq2 != 0) {
      nclust <- c(1, nclusters)
    }
    
    sim = cps.normal(
      nsim = nsim,
      nsubjects = nsubjects,
      nclusters = nclust,
      mu = mu,
      mu2 = mu2,
      ICC = ICC,
      ICC2 = ICC2,
      sigma_sq = sigma_sq,
      sigma_sq2 = sigma_sq2,
      alpha = alpha,
      sigma_b_sq = sigma_b_sq,
      sigma_b_sq2 = sigma_b_sq2,
      method = "glmm",
      quiet = quiet,
      allSimData = allSimData,
      lowPowerOverride = lowPowerOverride,
      timelimitOverride = timelimitOverride,
      nofit = nofit,
      seed = seed,
      irgtt = TRUE,
      poorFitOverride = poorFitOverride
    )
    return(sim)
  }


