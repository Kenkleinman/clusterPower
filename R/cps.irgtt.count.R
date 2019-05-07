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
#' cluster, number of clusters per treatment arm, between-cluster variance, 
#' two of the following: expected count in non-treatment group, expected count 
#' in treatment group, difference in counts between groups; significance level, 
#' analytic method, and whether or not progress updates should be displayed 
#' while the function is running.
#' 
#' 
#' @param nsim Number of datasets to simulate; accepts integer (required).
#' @param nsubjects Number of subjects per cluster; accepts integer (required). 
#' @param nclusters Number of clusters per treatment group; accepts integer (required).
#' At least 2 of the following 3 arguments must be specified:
#' @param c1 Expected outcome count in non-treatment group
#' @param c2 Expected outcome count in treatment group
#' @param c.diff Expected difference in outcome count between groups, defined as c.diff = c1 - c2
#' @param sigma_b Between-cluster variance; if sigma_b2 is not specified, 
#' between cluster variances are assumed to be equal between groups. Accepts numeric
#' If between cluster variances differ between treatment groups, the following must also be specified:
#' @param sigma_b2 Between-cluster variance for clusters in TREATMENT group
#' @param family Distribution from which responses are simulated. Accepts Poisson ('poisson') or negative binomial ('neg.binom') (required); default = 'poisson'
#' @param analysis Family used for regression; currently only applicable for GLMM. Accepts 'poisson' or 'neg.binom' (required); default = 'poisson'
#' @param alpha Significance level. Default = 0.05.
#' @param quiet When set to FALSE, displays simulation progress and estimated completion time. Default = FALSE.
#' @param all.sim.data Option to output list of all simulated datasets. Default = FALSE
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
#'   \item Data frame reporting between-cluster variances for Treatment/Non-Treatment groups
#'   \item Vector containing expected counts and risk ratios based on user inputs
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
#' }
#' @author Alexander R. Bogdan 
#' @author Alexandria C. Sakrejda (\email{acbro0@@umass.edu}
#' @author Ken Kleinman (\email{ken.kleinman@@gmail.com})
#'
#' 
#' @examples 
#' \dontrun{
#' count.sim <- cps.irgtt.count(nsim = 100, nsubjects = c(400, 50), nclusters = c(1, 60), c1 = 1000,
#'                       c2 = 2500, sigma_b2 = 1, family = 'poisson', analysis = 'poisson',
#'                       alpha = 0.05, quiet = FALSE, all.sim.data = TRUE)
#' }
#'
#' @export

# Define function
cps.irgtt.count <- function(nsim = NULL, nsubjects = NULL, nclusters = NULL, c1 = NULL, c2 = NULL, 
                     c.diff = NULL, sigma_b = 0, sigma_b2 = 0, family = 'poisson', 
                     analysis = 'poisson', alpha = 0.05, quiet = FALSE, 
                     all.sim.data = FALSE){
  
  if (sigma_b == 0 & sigma_b2 == 0){
    warning("Sigma_b in both arms is 0. This is equivalent to a t-test. Did you mean to 
            enter a sigma_b value for the arm containing clustered observations?")
  }
  if (sigma_b != 0 & sigma_b2 != 0){
    warning("Sigma_b is not zero for either arm. Did you want to use cps.count()?")
  }
  if (sigma_b != 0 & sigma_b2 == 0){
    stop("Non-clustered group must be entered as the reference group.")
    nclust <- c(nclusters, 1)
  }
  if (sigma_b2 != 0){
    nclust <- c(1, nclusters)
  }
  
  sim <- cps.count(nsim = nsim, nsubjects = nsubjects, nclusters = nclusters, c1 = c1, c2 = c2, 
                   c.diff = c.diff, sigma_b = sigma_b, sigma_b2 = sigma_b2, family = family, 
                   analysis = analysis, method = 'glmm', alpha = alpha, quiet = quiet, 
                   all.sim.data = all.sim.data, irgtt = TRUE)
  return(sim)
  }
