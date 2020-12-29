#' Simulation-based power estimation for binary outcome individually 
#' randomized group treatment trials.
#' 
#' This function utilizes iterative simulations to determine 
#' approximate power for individually randomized group treatment trials. Users 
#' can modify a variety of parameters to suit the simulations to their
#' desired experimental situation.
#' 
#' Runs the power simulation for binary outcomes.
#' 
#' Users must specify the desired number of simulations, number of subjects per 
#' cluster, number of clusters per arm, two of the following three terms: 
#' expected probability of outcome in arm 1, expected probability of 
#' outcome in arm 2, expected difference in probabilities between groups; 
#' significance level, progress updates, 
#' and simulated data set output may also be specified.
#' 
#' @param nsim Number of datasets to simulate; accepts integer (required).
#' @param nsubjects Number of subjects per cluster in the clustered (arm 2) 
#' group; accepts integer (required). 
#' @param nclusters Number of clusters per arm; accepts integer (required).
#' @param p1 Expected probability of outcome in arm 1 (required)
#' @param p2 Expected probability of outcome in arm 2 (required)
#' @param sigma_b_sq Between-cluster variance; defaults to 0. Accepts numeric.
#' @param sigma_b_sq2 Between-cluster variance for clusters in arm 2.
#' @param alpha Significance level; default = 0.05
#' @param quiet When set to FALSE, displays simulation progress and estimated 
#' completion time, default is TRUE.
#' @param allSimData Option to output list of all simulated datasets; default = FALSE.
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
#'  
#' @return A list with the following components
#' \itemize{
#'   \item Character string indicating total number of simulations, simulation type, 
#'   and number of convergent models
#'   \item Number of simulations
#'   \item Data frame with columns "Power" (Estimated statistical power), 
#'   "lower.95.ci" (Lower 95% confidence interval bound), "upper.95.ci" (Upper 95% 
#'   confidence interval bound)
#'   \item Analytic method used for power estimation
#'   \item Significance level
#'   \item Vector containing user-defined cluster sizes
#'   \item Vector containing user-defined number of clusters
#'   \item Data frame reporting sigma_b_sq for each group
#'   \item Vector containing expected difference in probabilities based on user inputs
#'   \item Data frame with columns: "Estimate" (Estimate of treatment effect for a given 
#'   simulation), "Std.err" (Standard error for treatment effect estimate), 
#'   "Test.statistic" (z-value (for GLMM) or 
#'   Wald statistic (for GEE)), "p.value", "converge" (Did simulated model converge?), 
#'   "sig.val" (Is p-value less than alpha?)
#'   \item If \code{allSimData = TRUE}, a list of data frames, each containing: 
#'   "y" (Simulated response value), 
#'   "trt" (Indicator for arm), "clust" (Indicator for cluster)
#'   \item List of warning messages produced by non-convergent models; 
#'   Includes model number for cross-referencing against \code{model.estimates}
#' }
#' If \code{nofit = T}, a data frame of the simulated data sets, containing:
#' \itemize{
#'   \item "arm" (Indicator for treatment arm)
#'   \item "cluster" (Indicator for cluster)
#'   \item "y1" ... "yn" (Simulated response value for each of the \code{nsim} data sets).
#'   }
#' 
#' @references Snjiders, T. & Bosker, R. Multilevel Analysis: an Introduction to Basic and 
#' Advanced Multilevel Modelling. London, 1999: Sage.
#' @references Eldridge, S., Ukoumunne, O. & Carlin, J. The Intra-Cluster Correlation 
#' Coefficient in Cluster Randomized Trials: 
#' A Review of Definitions. International Statistical Review (2009), 77, 3, 378-394. 
#' doi: 10.1111/j.1751-5823.2009.00092.x
#' 
#' @examples 
#' \dontrun{
#' irgtt.binary.sim <- cps.irgtt.binary(nsim = 100, nsubjects = c(150, 30), 
#'                         nclusters = 10, p1 = 0.44,
#'                         p2 = 0.2, sigma_b_sq2 = 1, alpha = 0.05, 
#'                         allSimData = FALSE)
#' }
#' @author Alexandria C. Sakrejda (\email{acbro0@@umass.edu})
#' @author Alexander R. Bogdan
#' @author Ken Kleinman (\email{ken.kleinman@@gmail.com})
#' 
#' @export

cps.irgtt.binary <-
  function(nsim = NULL,
           nsubjects = NULL,
           nclusters = NULL,
           p1 = NULL,
           p2 = NULL,
           sigma_b_sq = 0,
           sigma_b_sq2 = 0,
           alpha = 0.05,
           quiet = TRUE,
           allSimData = FALSE,
           poorFitOverride = FALSE,
           lowPowerOverride = FALSE, 
           timelimitOverride = TRUE,
           nofit = FALSE,
           seed = NA) {
    
    if (length(nclusters) > 2) {
      stop(
        "NCLUSTERS can only be a vector of length 1 (no clusters in control group)."
      )
    }
    
    if (sigma_b_sq == 0 & sigma_b_sq2 == 0) {
      warning(
        "sigma_b_sq in both arms is 0. This is equivalent to a t-test. Did you mean to
            enter a sigma_b_sq value for the arm containing clustered observations?"
      )
    }
    if (sigma_b_sq != 0 & sigma_b_sq2 != 0) {
      warning("sigma_b_sq is not zero for either arm. Did you want to use cps.binary()?")
    }
    if (sigma_b_sq != 0 & sigma_b_sq2 == 0) {
      stop("Non-clustered group must be entered as the reference group.")
      nclust <- c(nclusters, 1)
    }
    if (sigma_b_sq2 != 0) {
      nclust <- c(1, nclusters)
    }
    
    sim <-
      cps.binary(
        nsim = nsim,
        nsubjects = nsubjects,
        nclusters = nclust,
        p1 = p1,
        p2 = p2,
        sigma_b_sq = sigma_b_sq,
        sigma_b_sq2 = sigma_b_sq2,
        alpha = alpha,
        method = 'glmm',
        quiet = quiet,
        allSimData = allSimData,
        poorFitOverride = poorFitOverride,
        lowPowerOverride = lowPowerOverride, 
        timelimitOverride = timelimitOverride,
        nofit = nofit,
        seed = seed,
        irgtt = TRUE
      )
    return(sim)
  }
