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
#' cluster, number of clusters per treatment arm, two of the following three terms: 
#' expected probability of outcome in non-treatment group, expected probability of 
#' outcome in treatment group, expected difference in probabilities between groups
#' ; significance level, analytic method, progress updates, 
#' and simulated data set output may also be specified.
#' 
#' The following equations are used to estimate intra-cluster correltation coefficients:
#' P_h: \deqn{ICC = \frac{\sigma_{b}}{\sigma_{b} + \pi^{2}/3}}
#' P_c: \deqn{ICC = \frac{P(Y_{ij} = 1, Y_{ih} = 1) - \pi_{j}\pi_{h}}{\sqrt{\pi_{j}(1 - \pi_{j})\pi_{h}(1 - \pi_{h})}}}
#' P_lmer: \deqn{ICC = \frac{\sigma_{b}}{\sigma_{b} + \sigma_{w}}}
#' 
#' @param nsim Number of datasets to simulate; accepts integer (required).
#' @param nsubjects Number of subjects per cluster in the clustered (treatment) group; accepts integer (required). 
#' @param nclusters Number of clusters per treatment group; accepts integer (required).
#' @param p1 Expected probability of outcome in non-treatment group (required)
#' @param p2 Expected probability of outcome in treatment group (required)
#' @param sigma_b Between-cluster variance; defaults to 0. Accepts numeric.
#' @param sigma_b2 Between-cluster variance for clusters in TREATMENT group
#' @param alpha Significance level; default = 0.05
#' @param method Analytical method, either Generalized Linear Mixed Effects Model (GLMM) or Generalized Estimating Equation (GEE). Accepts c('glmm', 'gee') (required); default = 'glmm'.
#' @param quiet When set to FALSE, displays simulation progress and estimated completion time, default is TRUE.
#' @param all.sim.data Option to output list of all simulated datasets; default = FALSE
#' @param seed Option to set seed. Default is NA.
#'  
#' @return A list with the following components
#' \itemize{
#'   \item Character string indicating total number of simulations, simulation type, and number of convergent models
#'   \item Number of simulations
#'   \item Data frame with columns "Power" (Estimated statistical power), 
#'   "lower.95.ci" (Lower 95% confidence interval bound), "upper.95.ci" (Upper 95% confidence interval bound)
#'   \item Analytic method used for power estimation
#'   \item Significance level
#'   \item Vector containing user-defined cluster sizes
#'   \item Vector containing user-defined number of clusters
#'   \item Data frame reporting sigma_b for each group
#'   \item Vector containing expected difference in probabilities based on user inputs
#'   \item Data frame containing three estimates of ICC
#'   \item Data frame with columns: "Estimate" (Estimate of treatment effect for a given simulation), 
#'   "Std.err" (Standard error for treatment effect estimate), "Test.statistic" (z-value (for GLMM) or 
#'   Wald statistic (for GEE)), "p.value", "converge" (Did simulated model converge?), 
#'   "sig.val" (Is p-value less than alpha?)
#'   \item List of data frames, each containing: "y" (Simulated response value), 
#'   "trt" (Indicator for treatment group), "clust" (Indicator for cluster)
#'   \item List of warning messages produced by non-convergent models; 
#'   Includes model number for cross-referencing against \code{model.estimates}
#' }
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
#' irgtt.binary.sim <- cps.irgtt.binary(nsim = 100, nsubjects = 30, nclusters = 10, p1 = 0.5,
#'                         p2 = 0.2, sigma_b = 0, sigma_b2 = 1, alpha = 0.05, 
#'                         all.sim.data = FALSE)
#' }
#' @author Alexandria C. Sakrejda (\email{acbro0@@umass.edu}), Alexander R. Bogdan, 
#'   and Ken Kleinman (\email{ken.kleinman@@gmail.com})
#' @export

#FIXME add irgtt equation in cps.binary
# Define function
cps.irgtt.binary <- function(nsim = NULL, nsubjects = NULL, nclusters = NULL,
                      p1 = NULL, p2 = NULL, sigma_b = 0, sigma_b2 = 0, 
                      alpha = 0.05,
                      quiet = TRUE, all.sim.data = FALSE, seed = NA){
  
  if (sigma_b == 0 & sigma_b2 == 0){
    warning("Sigma_b in both arms is 0. This is equivalent to a t-test. Did you mean to 
            enter a sigma_b value for the arm containing clustered observations?")
  }
  if (sigma_b != 0 & sigma_b2 != 0){
    warning("Sigma_b is not zero for either arm. Did you want to use cps.binary()?")
  }
  if (sigma_b != 0 & sigma_b2 == 0){
    stop("Non-clustered group must be entered as the reference group.")
    nclust <- c(nclusters, 1)
  }
  if (sigma_b2 != 0){
    nclust <- c(1, nclusters)
  }
  
  sim <- cps.binary(nsim = nsim, nsubjects = nsubjects, nclusters = nclust,
                        p1 = p1, p2 = p2, or1 = NULL, or2 = NULL, or.diff = NULL, 
                        sigma_b = sigma_b, sigma_b2 = sigma_b2, alpha = alpha, method = 'glmm', 
                        quiet = quiet, all.sim.data = all.sim.data, seed = seed, irgtt = TRUE)
  return(sim)
  }