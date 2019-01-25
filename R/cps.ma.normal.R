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
#user must supply either nsubjects or narms and nclusters
# 12. format gee output


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
#' @param singular.fit.override Option to override \code{stop()} if more than 25% of fits fail to converge; default = FALSE 
#'  
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
#'   \item{model.estimates}{Data frame with columns corresponding to each arm with the suffixes as 
#'   follows: 
#'                   ".Estimate" (Estimate of treatment effect for a given simulation), 
#'                   "Std.Err" (Standard error for treatment effect estimate), 
#'                   ".tval" (for GLMM) | ".wald" (for GEE), 
#'                   ".pval"
#'   \item{sim.data}{List of \code{nsim} data frames, each containing: 
#'                   "y" (Simulated response value), 
#'                   "trt" (Indicator for treatment group), 
#'                   "clust" (Indicator for cluster)}
#'   \item{proportion.failed.to.converge} {Character string containing the percent of
#'   \code{nsim} in which the glmm fit was singular, produced only when method == "glmm" & 
#'   all.sim.data==FALSE}
#'   \item{failed.to.converge}{Vector containing of length \code{nsim} denoting whether 
#'   or not a simulation glmm fit was singular, produced only when method == "glmm" & 
#'   all.sim.data==TRUE}
#'          
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
#'                        seed = NULL, 
#'                        singular.fit.override = FALSE)
#' }
#' 
#' @export
#' 
 
cps.ma.normal <- function(nsim = 1000, nsubjects = NULL, 
                           narms = NULL, nclusters = NULL,
                        means = NULL, sigma = NULL, 
                        sigma_b = NULL, alpha = 0.05,
                        quiet = FALSE, ICC=NULL, method = 'glmm', 
                        all.sim.data = FALSE, seed = 123, 
                        singular.fit.override = FALSE){

  # input object validation steps
  if (exists("nsubjects", mode = "any")==FALSE){
    stop("nsubjects must be specified. See ?cps.ma.normal for help.")
  }
  if (length(nsubjects)==1 & exists("nclusters", mode = "numeric")==FALSE){
    stop("When nsubjects is scalar, user must supply nclusters (clusters per arm)")
  }
  if (length(nsubjects)==1 & length(nclusters)==1 & 
      exists("narms", mode = "numeric")==FALSE){
    stop("User must provide narms when nsubjects and nclusters are both scalar.")
  }
  
  # create narms and nclusters if not provided directly by user
  if (exists("nsubjects", mode = "list")==TRUE){
    # create narms and nclusters if not supplied by the user
    if (exists("narms", mode = "numeric")==FALSE){
      narms <- length(nsubjects)
    }
    if (exists("nclusters", mode = "numeric")==FALSE){
      nclusters <- vapply(nsubjects, length, 0)
    }
  }
  if(length(nclusters)==1 & exists("nsubjects", mode = "list")==FALSE){
    nclusters <- rep(nclusters, narms)
  }
  if(length(nclusters)>1 & length(nsubjects)==1){
    narms <- length(nclusters)
  }

  # Create nsubjects structure from narms and nclusters when nsubjects is scalar
  if (length(nsubjects)==1){
    str.nsubjects <- lapply(nclusters, function(x) rep(nsubjects, x))
  } else {
    str.nsubjects <- nsubjects
  }
  
  # supplies sigma or sigma_b if user supplies ICC
  if (exists("ICC", mode = "numeric")==TRUE){
  if (exists("sigma", mode = "numeric")==FALSE){
    sigma <- createMissingVarianceParam(sigma_b = sigma_b, ICC = ICC)
  }
  if (exists("sigma_b", mode = "numeric")==FALSE){
    sigma_b <- createMissingVarianceParam(sigma = sigma, ICC = ICC)
  }
  }

   # run the simulations 
   normal.ma.rct <- cps.ma.normal.internal(nsim = nsim, 
                                           str.nsubjects = str.nsubjects, 
                                           means = means, sigma = sigma, 
                                           sigma_b = sigma_b, alpha = alpha, 
                                           quiet = quiet, method = method, 
                                           all.sim.data = all.sim.data,
                                           seed = seed,
                                           singular.fit.override = singular.fit.override)
   
   models <- normal.ma.rct[[1]]
   
#Organize output for GLMM
 if(method=="glmm"){
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
 
 # Organize the row/col names for the model estimates output
   keep.names <- rownames(models[[1]][[10]])
   keep.names[grepl(keep.names, pattern = "(Intercept)")==TRUE] <- "intercept"
 
   names.Est <- rep(NA, narms)
   names.st.err <- rep(NA, narms)
   names.tval <- rep(NA, narms)
   names.pval <- rep(NA, narms)
   names.power <- rep(NA, narms)
 
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
   
   # Organize the LRT output
   LRT.holder <- matrix(unlist(normal.ma.rct[[2]]), ncol=4, nrow=nsim, 
                     byrow=TRUE, 
                     dimnames = list(seq(1:nsim), 
                            c("Df", "Sum Sq", "Mean Sq", "F value")))
   
 
 # Calculate and store power estimate & confidence intervals
   sig.val <-  ifelse(p.val < alpha, 1, 0)
   pval.power <- apply (sig.val, 2, FUN=function(x) {sum(x, na.rm=TRUE)/nsim})
   power.parms <-  data.frame(Power = round(pval.power, 3),
                            Lower.95.CI = round(pval.power - abs(stats::qnorm(alpha / 2)) * 
                                                  sqrt((pval.power * (1 - pval.power)) / nsim), 3),
                            Upper.95.CI = round(pval.power + abs(stats::qnorm(alpha / 2)) * 
                                                  sqrt((pval.power * (1 - pval.power)) / nsim), 3))
   rownames(power.parms) <- names.power
    
  # Store simulation output in data frame
   ma.model.est <-  data.frame(Estimates, std.error, t.val, p.val)
   
   ## Output objects for GLMM

   # Create list containing all output (class 'crtpwr') and return
   if(all.sim.data == TRUE){
     complete.output <-  list("power" <-  power.parms,
                              "model.estimates" <- ma.model.est, 
                              "overall.sig" <- LRT.holder,
                              "sim.data" <-  normal.ma.rct[[3]], 
                              "failed.to.converge" <-  normal.ma.rct[[4]])
   } else {
     complete.output <-  list("power" <-  power.parms,
                              "model.estimates" <- ma.model.est,
                              "overall.sig" <- LRT.holder,
                              "proportion.failed.to.converge" <- normal.ma.rct[[3]])
   }
   return(complete.output)
 }
   #Organize output for GEE method
   if (method=="gee"){
     # Organize the output
     Estimates = matrix(NA, nrow = nsim, ncol = narms)
     std.error = matrix(NA, nrow = nsim, ncol = narms)
     Wald = matrix(NA, nrow = nsim, ncol = narms)
     Pr = matrix(NA, nrow = nsim, ncol = narms)
     
     for (i in 1:nsim){
       Estimates[i,] <- models[[i]]$coefficients[,1]
       std.error[i,] <- models[[i]]$coefficients[,2]
       Wald[i,] <- models[[i]]$coefficients[,3]
       Pr[i,] <- models[[i]]$coefficients[,4]
     }
     
     # Organize the row/col names for the output
     keep.names <- rownames(models[[1]]$coefficients)
     keep.names[grepl(keep.names, pattern = "(Intercept)")==TRUE] <- "intercept"
     
     names.Est <- rep(NA, length(narms))
     names.st.err <- rep(NA, length(narms))
     names.wald <- rep(NA, length(narms))
     names.pval <- rep(NA, length(narms))
     names.power <- rep(NA, length(narms))
     
     for (i in 1:length(keep.names)){
       names.Est[i] <- paste(keep.names[i], ".Estimate", sep="")
       names.st.err[i] <- paste(keep.names[i], ".Std.Err", sep="")
       names.wald[i] <- paste(keep.names[i], ".wald", sep="")
       names.pval[i] <- paste(keep.names[i], ".pval", sep="")
       names.power[i] <- paste(keep.names[i], ".power", sep="")
     }
     colnames(Estimates) <- names.Est
     colnames(std.error) <- names.st.err
     colnames(Wald) <- names.wald
     colnames(Pr) <- names.pval
     
     # Organize the LRT output
     LRT.holder <- matrix(unlist(normal.ma.rct[[2]]), ncol=3, nrow=nsim, 
                          byrow=TRUE, 
                          dimnames = list(seq(1:nsim), c("Df", "X2", "P(>|Chi|)")))
     
     # Calculate and store power estimate & confidence intervals
     sig.val <-  ifelse(Pr < alpha, 1, 0)
     pval.power <- apply (sig.val, 2, FUN=function(x) {sum(x, na.rm=TRUE)/nsim})
     power.parms <-  data.frame(Power = round(pval.power, 3),
                                Lower.95.CI = round(pval.power - abs(stats::qnorm(alpha / 2)) * 
                                                      sqrt((pval.power * (1 - pval.power)) / nsim), 3),
                                Upper.95.CI = round(pval.power + abs(stats::qnorm(alpha / 2)) * 
                                                      sqrt((pval.power * (1 - pval.power)) / nsim), 3))
     rownames(power.parms) <- names.power
     power.parms$pvalue <- pval.power
     
     # Store GEE simulation output in data frame
     ma.model.est <-  data.frame(Estimates, std.error, Wald, Pr)
   
   ## Output objects for GEE
   
     # Create list containing all output (class 'crtpwr') and return
     if(all.sim.data == TRUE){
       complete.output <-  list("power" <-  power.parms,
                                "model.estimates" <- ma.model.est, 
                                "overall.sig" <- LRT.holder,
                                "sim.data" <-  normal.ma.rct[[3]])
     } else {
       complete.output <-  list("power" <-  power.parms,
                                "model.estimates" <- ma.model.est,
                                "overall.sig" <- LRT.holder)
     }
     return(complete.output)
   }
}

