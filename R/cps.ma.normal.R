# Validation functions
validateVariance <- function(x){
  warning("FIXME: not actually validating variance yet")
}

is.wholenumber = function(x, tol = .Machine$double.eps^0.5)  abs(x - round(x)) < tol


##FIXME: TO DO
# 1. update the return values
# 2. update the example/ man text
# 3. input validation
# 5. make validate nsubjects fxn in validation file
# 9. write some usage examples
# 10. debug
# 11. testthat tests


#' Power simulations for cluster-randomized trials: Simple Designs, Continuous Outcome.
#'
#' This set of functions utilize iterative simulations to determine 
#' approximate power for cluster-randomized controlled trials. Users 
#' can modify a variety of parameters to suit the simulations to their
#' desired experimental situation.
#' 
#' Runs the power simulation.
#' 
#' Users must specify the desired number of simulations, number of subjects per 
#' cluster, number of clusters per treatment arm, expected absolute difference 
#' between treatments, two of the following: ICC, within-cluster variance, or 
#' between-cluster variance; significance level, analytic method, progress updates, 
#' and simulated data set output may also be specified.
#' 
#' # Return enough data to save the simulated datasets
#' # FIXME: update this
#' 
#' @author Alexandria C. Sakrejda
#' @author Alexander R. Bogdan
#' 
#' @examples # FIXME: update this too
#' \dontrun{
#' Put usage example here
#' }
#'
#' @param nsim Number of datasets to simulate; accepts integer (required).
#' @param nsubjects Number of subjects per treatment group; accepts a list with one entry per arm. 
#' Each entry is a vector containing the number of subjects per cluster (required).
#' @param means Expected absolute treatment effect for each arm; accepts a vector of length \code{narms} (required).
#' @param sigma Within-cluster variance; accepts a vector of length \code{narms} (required).
#' @param sigma_b Between-cluster variance; accepts a vector of length \code{narms} (required).
#' @param alpha Significance level; default = 0.05.
#' @param all.sim.data Option to output list of all simulated datasets; default = FALSE.
#' 
#' @return A list with the following components
#' \describe{
#'   \item{overview}{Character string indicating total number of simulations and simulation type}
#'   \item{nsim}{Number of simulations}
#'   \item{power}{Data frame with columns "Power" (Estimated statistical power), 
#'                "lower.95.ci" (Lower 95% confidence interval bound), 
#'                "upper.95.ci" (Upper 95% confidence interval bound)}
#'   \item{method}{Analytic method used for power estimation}
#'   \item{alpha}{Significance level}
#'   \item{cluster.sizes}{Vector containing user-defined cluster sizes}
#'   \item{n.clusters}{Vector containing user-defined number of clusters in each treatment group}
#'   \item{variance.parms}{Data frame reporting ICC for Treatment/Non-Treatment groups}
#'   \item{inputs}{Vector containing expected difference between groups based on user inputs}
#'   \item{model.estimates}{Data frame with columns: 
#'                   "Estimate" (Estimate of treatment effect for a given simulation), 
#'                   "Std.err" (Standard error for treatment effect estimate), 
#'                   "Test.statistic" (z-value (for GLMM) or Wald statistic (for GEE)), 
#'                   "p.value", 
#'                   "sig.val" (Is p-value less than alpha?)}
#'   \item{sim.data}{List of data frames, each containing: 
#'                   "y" (Simulated response value), 
#'                   "trt" (Indicator for treatment group), 
#'                   "clust" (Indicator for cluster)}
#' }
#' 
#' #' @examples 
#' \dontrun{
#' 
#' nsubjects.example <- list(c(20,20,20,25), c(15, 20, 20, 21), c(17, 20, 21))
#' means.example <- c(30, 21, 53)
#' sigma.example <- c(1, 1, 0.9)
#' sigma_b.example <- c(0.1, 0.15, 0.1)
#' 
#' multi.cps.normal <- cps.ma.normal(nsim = 2, nsubjects = nsubjects.example, 
#'                        means = means.example, sigma = sigma.example, 
#'                        sigma_b = sigma_b.example, alpha = 0.05,
#'                        quiet = FALSE, ICC=NULL, method = 'glmm', 
#'                        all.sim.data = FALSE,
#'                        seed = NULL)
#' }
#' 
#' @export
#' 
 cps.ma.normal <- function(nsim = 1000, nsubjects = NA, 
                           narms = NA, nclusters = NA,
                        means = NA, sigma = NA, 
                        sigma_b = NA, alpha = 0.05,
                        quiet = FALSE, ICC=NULL, method = 'glmm', 
                        all.sim.data = FALSE, seed = 123){

  # supplies sigma or sigma_b if user supplies ICC
  if (exists("sigma", mode = "object")==FALSE){
    sigma <- createMissingVarianceParam(sigma = sigma, 
                                        sigma_b = sigma_b, ICC = ICC)
  }
  if (exists("sigma_b", mode = "object")==FALSE){
    sigma_b <- createMissingVarianceParam(sigma = sigma, 
                                        sigma_b = sigma_b, ICC = ICC)
  }
   
   # creates nsubjects structure if nclusters and nsubjects are scalar
   if (length(nclusters)==1){
     nclusters <- rep(nclusters, narms)
   }
   if (length(nsubjects)==1){
     nsubjects.scalar <- nsubjects
     nsubjects <- list()
     for (i in 1:narms){
       nsubjects[[i]] <- rep(nsubjects.scalar, nclusters[i])
     }
   }
  
   # run the simulations 
   cps.out <- cps.ma.normal.internal(nsim = nsim, nsubjects = nsubjects, 
                      means = means, sigma = sigma, 
                      sigma_b = sigma_b, alpha = alpha, 
                      quiet = quiet, method = method, 
                      all.sim.data = all.sim.data,
                      seed = seed)
   
 }
   
 # Organize the output
 Estimates = matrix(NA, nrow = narms*nsim, ncol = narms)
 std.error = matrix(NA, nrow = narms*nsim, ncol = narms)
 t.val = matrix(NA, nrow = narms*nsim, ncol = narms)
 p.val = matrix(NA, nrow = narms*nsim, ncol = narms)
 
 for (i in 1:nsim){
   Estimates[i,] <- normal.ma.rct$estimates[[i]]$coefficients[,1]
   std.error[i,] <- normal.ma.rct$estimates[[i]]$coefficients[,2]
   t.val[i,] <- normal.ma.rct$estimates[[i]]$coefficients[,3]
   p.val[i,] = 2 * stats::pt(-abs(normal.ma.rct$estimates[[i]]$coefficients[,3]), 
                             df = sum(nclusters) - 2)
 }
 
 # Calculate and store power estimate & confidence intervals
 pval.power = sum(cps.model.est[, 'sig.val']) / nrow(cps.model.est)
 power.parms = data.frame(Power = round(pval.power, 3),
                          Lower.95.CI = round(pval.power - abs(stats::qnorm(alpha / 2)) * 
                                                sqrt((pval.power * (1 - pval.power)) / nsim), 3),
                          Upper.95.CI = round(pval.power + abs(stats::qnorm(alpha / 2)) * 
                                                sqrt((pval.power * (1 - pval.power)) / nsim), 3))
 
   
   
   
   
    
  ## Output objects
    
  # Store simulation output in data frame
  cps.model.est = data.frame(Estimate = as.vector(unlist(est.vector)),
                             Std.err = as.vector(unlist(se.vector)),
                             Test.statistic = as.vector(unlist(stat.vector)),
                             p.value = as.vector(unlist(pval.vector)))
  cps.model.est[, 'sig.val'] = ifelse(cps.model.est[, 'p.value'] < alpha, 1, 0)
  
  # Calculate and store power estimate & confidence intervals
  pval.power = sum(cps.model.est[, 'sig.val']) / nrow(cps.model.est)
  power.parms = data.frame(Power = round(pval.power, 3),
                           Lower.95.CI = round(pval.power - abs(stats::qnorm(alpha / 2)) * sqrt((pval.power * (1 - pval.power)) / nsim), 3),
                           Upper.95.CI = round(pval.power + abs(stats::qnorm(alpha / 2)) * sqrt((pval.power * (1 - pval.power)) / nsim), 3))
  
  # Create list containing all output (class 'crtpwr') and return
  complete.output = list("power" = power.parms,
                         "pval.power" = pval.power,
                         "model.estimates" = cps.model.est, 
                         "sim.data" = simulated.datasets)
  return(complete.output)
}
