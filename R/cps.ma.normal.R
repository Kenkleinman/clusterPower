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
#user must supply either


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
#' @param method Analytical method, either Generalized Linear Mixed Effects Model (GLMM) or 
#' Generalized Estimating Equation (GEE). Accepts c('glmm', 'gee') (required); default = 'glmm'.
#' @param quiet When set to FALSE, displays simulation progress and estimated completion time; default is FALSE.
#' @param seed Option to set.seed. Default is NULL.
#' 
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
 
cps.ma.normal <- function(nsim = 1000, nsubjects = NULL, 
                           narms = NULL, nclusters = NULL,
                        means = NULL, sigma = NULL, 
                        sigma_b = NULL, alpha = 0.05,
                        quiet = FALSE, ICC=NULL, method = 'glmm', 
                        all.sim.data = FALSE, seed = 123){

  # supplies sigma or sigma_b if user supplies ICC
  if (exists("ICC", mode = "object")==TRUE){
  if (exists("sigma", mode = "object")==FALSE){
    sigma <- createMissingVarianceParam(sigma_b = sigma_b, ICC = ICC)
  }
  if (exists("sigma_b", mode = "object")==FALSE){
    sigma_b <- createMissingVarianceParam(sigma = sigma, ICC = ICC)
  }
  }
   
  # create narms and nclusters if not supplied by the user
  if (exists("narms", mode = "object")==FALSE){
    narms <- length(nsubjects)
  }
  if (exists("nclusters", mode = "object")==FALSE){
    nclusters <- length(unlist(nsubjects))
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
   normal.ma.rct <- cps.ma.normal.internal(nsim = nsim, nsubjects = nsubjects, 
                      means = means, sigma = sigma, 
                      sigma_b = sigma_b, alpha = alpha, 
                      quiet = quiet, method = method, 
                      all.sim.data = all.sim.data,
                      seed = seed)
   
   models <- normal.ma.rct[[1]]

 # Organize the output
 Estimates = matrix(NA, nrow = nsim, ncol = narms)
 std.error = matrix(NA, nrow = nsim, ncol = narms)
 t.val = matrix(NA, nrow = nsim, ncol = narms)
 p.val = matrix(NA, nrow = nsim, ncol = narms)
 
 for (i in 1:nsim){
   Estimates[i,] <- models[[i]][[10]][,1]
   std.error[i,] <- models[[i]][[10]][,2]
   t.val[i,] <- models[[i]][[10]][,3]
   p.val[i,] <- 2 * stats::pt(-abs(models[[i]][[10]][,3]), 
                             df = sum(nclusters) - 2)
 }
 
 # Organize the row/col names for the output
 keep.names <- rownames(models[[1]][[10]])
 keep.names[grepl(keep.names, pattern = "(Intercept)")==TRUE] <- "intercept"
 
 names.Est <- rep(NA, length(narms))
 names.st.err <- rep(NA, length(narms))
 names.tval <- rep(NA, length(narms))
 names.pval <- rep(NA, length(narms))
 names.power <- rep(NA, length(narms))
 
 for (i in 1:length(keep.names)){
   names.Est[i] <- paste(keep.names[i], ".Estimate", sep="")
   names.st.err[i] <- paste(keep.names[i], ".Std.Err", sep="")
   names.tval[i] <- paste(keep.names[i], ".tval", sep="")
   names.pval[i] <- paste(keep.names[i], ".pval", sep="")
   names.power[i] <- paste(keep.names[i], ".power", sep="")
 }
 colnames(Estimates) <- names.Est
 colnames(std.error) <- names.st.err
 colnames(t.val) <- names.tval
 colnames(p.val) <- names.pval
 
 # Calculate and store power estimate & confidence intervals
 sig.val <-  ifelse(p.val < alpha, 1, 0)
 pval.power <- apply (sig.val, 2, FUN=function(x) {sum(x, na.rm=TRUE)/nsim})
 power.parms <-  data.frame(Power = round(pval.power, 3),
                          Lower.95.CI = round(pval.power - abs(stats::qnorm(alpha / 2)) * 
                                                sqrt((pval.power * (1 - pval.power)) / nsim), 3),
                          Upper.95.CI = round(pval.power + abs(stats::qnorm(alpha / 2)) * 
                                                sqrt((pval.power * (1 - pval.power)) / nsim), 3))
 rownames(power.parms) <- names.power
 power.parms$pvalue <- pval.power
 
  ## Output objects
    
  # Store simulation output in data frame
 ma.model.est <-  data.frame(Estimates, std.error, t.val, p.val)
 
  # Create list containing all output (class 'crtpwr') and return
 if(all.sim.data == TRUE){
 complete.output <-  list("power" <-  power.parms,
                         "model.estimates" <- ma.model.est, 
                         "sim.data" <-  normal.ma.rct[[2]], 
                         "failed.to.converge" <-  normal.ma.rct[[3]])
 } else {
 complete.output <-  list("power" <-  power.parms,
                            "model.estimates" <- ma.model.est,
                          "proportion.failed.to.converge" <- normal.ma.rct[[2]])
 }
  return(complete.output)
}
