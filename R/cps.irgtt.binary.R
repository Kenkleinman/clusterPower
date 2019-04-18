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
#' @param nsubjects Number of subjects per cluster; accepts integer (required). 
#' @param nclusters Number of clusters per treatment group; accepts integer (required).
#' At least 2 of the following 3 arguments must be specified when using expected probabilities:
#' @param p1 Expected probability of outcome in non-treatment group
#' @param p2 Expected probability of outcome in treatment group
#' @param p.diff Expected difference in probability of outcome between groups, defined as p.diff = p1 - p2
#' At least 2 of the following 3 arguments must be specified when using expected odds ratios:
#' @param or1 Expected odds ratio for outcome in non-treatment group
#' @param or2 Expected odds ratio for outcome in treatment group
#' @param or.diff Expected difference in odds ratio for outcome between groups, defined as or.diff = or1 - or2
#' @param sigma_b Between-cluster variance; if sigma_b2 is not specified, 
#' between cluster variances are assumed to be equal for both groups. Accepts numeric.
#' If between cluster variances differ between treatment groups, sigma_b2 must also be specified:
#' @param sigma_b2 Between-cluster variance for clusters in TREATMENT group
#' @param alpha Significance level; default = 0.05
#' @param method Analytical method, either Generalized Linear Mixed Effects Model (GLMM) or Generalized Estimating Equation (GEE). Accepts c('glmm', 'gee') (required); default = 'glmm'.
#' @param quiet When set to FALSE, displays simulation progress and estimated completion time, default is TRUE.
#' @param all.sim.data Option to output list of all simulated datasets; default = FALSE
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
#' @references Elridge, S., Ukoumunne, O. & Carlin, J. The Intra-Cluster Correlation 
#' Coefficient in Cluster Randomized Trials: 
#' A Review of Definitions. International Statistical Review (2009), 77, 3, 378-394. 
#' doi: 10.1111/j.1751-5823.2009.00092.x
#' 
#' @examples 
#' \dontrun{
#' binary.sim = cps.binary(nsim = 100, nsubjects = 50, nclusters = 6, p1 = 0.4,
#'                         p2 = 0.2, sigma_b = 100, alpha = 0.05, method = 'glmm', 
#'                         all.sim.data = FALSE)
#' }
#' @author Alexandria C. Sakrejda (\email{acbro0@@umass.edu}), Alexander R. Bogdan, 
#'   and Ken Kleinman (\email{ken.kleinman@@gmail.com})
#' @export

# Define function
cps.irgtt.binary <-  function(nsim = NULL, nsubjects = NULL, nclusters = NULL, p.diff = NULL,
                      p1 = NULL, p2 = NULL, sigma_b = NULL, sigma_b2 = NULL, 
                      alpha = 0.05, method = 'glmm', 
                      quiet = TRUE, all.sim.data = FALSE, seed = NA){
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


