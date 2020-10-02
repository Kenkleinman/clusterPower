#' Power simulations for cluster-randomized trials: Individually randomized
#' group treatment trial designs, count outcome.
#'
#' This function utilizes iterative simulations to determine 
#' approximate power for cluster-randomized controlled trials. Users 
#' can modify a variety of parameters to suit the simulations to their
#' desired experimental situation.
#' 
#' Runs the power simulation for count outcomes.
#' 
#' Users must specify the desired number of simulations, number of subjects per 
#' cluster, number of clusters per arm, between-cluster variance, 
#' two of the following: expected count in arm 1 (no clusters), expected count 
#' in arm 2 (clustered arm), expected difference in counts between arms; significance level, 
#' analytic method, and whether or not progress updates should be displayed 
#' while the function is running.
#' 
#' 
#' @param nsim Number of datasets to simulate; accepts integer (required).
#' @param nsubjects Number of subjects per cluster; accepts integer (required). 
#' @param nclusters Number of clusters in the arm; accepts integer (required). 
#' Arm 1 cluster size defaults to 1.
#' 
#' At least 2 of the following 3 arguments must be specified:
#' @param c1 Expected outcome count in arm 1
#' @param c2 Expected outcome count in arm 2
#' @param sigma_b_sq Between-cluster variance; defaults to 0. Accepts numeric.
#' If between cluster variances differ between arms, the following must 
#' also be specified:
#' @param sigma_b_sq2 Between-cluster variance for clusters in arm 2
#' @param family Distribution from which responses are simulated. Accepts Poisson 
#' ('poisson') or negative binomial ('neg.binom') (required); default = 'poisson'
#' @param analysis Family used for regression; currently only applicable for GLMM. 
#' Accepts 'poisson' or 'neg.binom' (required); default = 'poisson'
#' @param negBinomSize Only used when generating simulated data from the 
#' negative binomial (family = 'neg.binom'), this is the target for number of 
#' successful trials, or the dispersion parameter (the shape parameter of the gamma 
#' mixing distribution). Must be strictly positive but need not be integer. 
#' Defaults to 1.
#' @param alpha Significance level. Default = 0.05.
#' @param quiet When set to FALSE, displays simulation progress and estimated 
#' completion time. Default = FALSE.
#' @param allSimData Option to output list of all simulated datasets. Default = FALSE.
#' @param poorFitOverride Option to override \code{stop()} if more than 25\%
#' of fits fail to converge; default = FALSE.
#' @param lowPowerOverride Option to override \code{stop()} if the power
#' is less than 0.5 after the first 50 simulations and every ten simulations
#' thereafter. On function execution stop, the actual power is printed in the
#' stop message. Default = FALSE. When TRUE, this check is ignored and the
#' calculated power is returned regardless of value.
#' @param timelimitOverride Logical. When FALSE, stops execution if the estimated completion time
#' is more than 2 minutes. Defaults to TRUE.
#' @param nofit Option to skip model fitting and analysis and return the simulated data.
#' Defaults to \code{FALSE}.
#' @param seed Option to set seed. Default is NA.
#' @param opt Option to fit with a different optimizer (using the package \code{optimx}). 
#' Defaults to L-BFGS-B.
#' 
#' @return A list with the following components
#' \itemize{
#'   \item Character string indicating total number of simulations, distribution of 
#'   simulated data, and regression family
#'   \item Number of simulations
#'   \item Data frame with columns "Power" (Estimated statistical power), 
#'                "lower.95.ci" (Lower 95% confidence interval bound), 
#'                "upper.95.ci" (Upper 95% confidence interval bound)
#'   \item Analytic method used for power estimation
#'   \item Data frame containing families for distribution and analysis of simulated data
#'   \item Significance level
#'   \item Vector containing user-defined cluster sizes
#'   \item Vector containing user-defined number of clusters
#'   \item Data frame reporting between-cluster variances for both arms
#'   \item Vector containing expected counts and risk ratios based on user inputs
#'   \item Data frame with columns: 
#'                   "Estimate" (Estimate of treatment effect for a given simulation), 
#'                   "Std.err" (Standard error for treatment effect estimate), 
#'                   "Test.statistic" (z-value (for GLMM) or Wald statistic (for GEE)), 
#'                   "p.value",
#'                   "sig.val" (Is p-value less than alpha?)
#'   \item If \code{allSimData = TRUE}, a list of data frames, each containing: 
#'                   "y" (Simulated response value), 
#'                   "trt" (Indicator for arm), 
#'                   "clust" (Indicator for cluster)
#' }
#' If \code{nofit = T}, a data frame of the simulated data sets, containing:
#' \itemize{
#'   \item "arm" (Indicator for treatment arm)
#'   \item "cluster" (Indicator for cluster)
#'   \item "y1" ... "yn" (Simulated response value for each of the \code{nsim} data sets).
#'   }
#'   
#' @author Alexandria C. Sakrejda (\email{acbro0@@umass.edu}
#' @author Alexander R. Bogdan 
#' @author Ken Kleinman (\email{ken.kleinman@@gmail.com})
#'
#' 
#' @examples 
#' \dontrun{
#' irgtt.count.sim <- cps.irgtt.count(nsim = 100, nsubjects = c(500, 10), nclusters = 500, 
#'                              c1 = 85, c2 = 450, sigma_b_sq2 = 0.25, 
#'                              family = 'poisson', analysis = 'poisson',
#'                              alpha = 0.05, quiet = FALSE, allSimData = FALSE)
#' }
#'
#' @export


# Define function
cps.irgtt.count <-
  function(nsim = NULL,
           nsubjects = NULL,
           nclusters = NULL,
           c1 = NULL,
           c2 = NULL,
           sigma_b_sq = 0,
           sigma_b_sq2 = 0,
           family = 'poisson',
           analysis = 'poisson',
           negBinomSize = 1,
           alpha = 0.05,
           quiet = FALSE,
           allSimData = FALSE,
           poorFitOverride = FALSE,
           lowPowerOverride = FALSE, 
           timelimitOverride = TRUE,
           nofit = FALSE,
           seed = NA,
           opt = "L-BFGS-B") {
    if (sigma_b_sq == 0 & sigma_b_sq2 == 0) {
      warning(
        "sigma_b_sq in both arms is 0. This is equivalent to a t-test. Did you mean to
            enter a sigma_b_sq value for the arm containing clustered observations?"
      )
    }
    if (sigma_b_sq != 0 & sigma_b_sq2 != 0) {
      warning("sigma_b_sq is not zero for either arm. Did you want to use cps.count()?")
    }
    if (sigma_b_sq != 0 & sigma_b_sq2 == 0) {
      stop("Non-clustered group must be entered as the reference group.")
      nclust <- c(nclusters, 1)
    }
    if (sigma_b_sq2 != 0) {
      nclust <- c(1, nclusters)
    }
    
    sim <-
      cps.count(
        nsim = nsim,
        nsubjects = nsubjects,
        nclusters = nclusters,
        c1 = c1,
        c2 = c2,
        sigma_b_sq = sigma_b_sq,
        sigma_b_sq2 = sigma_b_sq2,
        family = family,
        analysis = analysis,
        negBinomSize = negBinomSize,
        method = 'glmm',
        alpha = alpha,
        quiet = quiet,
        allSimData = allSimData,
        poorFitOverride = poorFitOverride,
        lowPowerOverride = lowPowerOverride, 
        timelimitOverride = timelimitOverride,
        nofit = nofit,
        seed = seed,
        irgtt = TRUE,
        opt = opt
      )
    return(sim)
  }
