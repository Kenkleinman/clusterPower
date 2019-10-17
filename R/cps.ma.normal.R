#' Simulation-based power estimation for continuous outcome multi-arm 
#' cluster-randomized trials.
#'
#' This function uses iterative simulations to determine 
#' approximate power for multi-arm cluster-randomized controlled trials with a 
#' normally-distributed outcome of interest. Users can modify a variety of 
#' parameters to suit the simulations to their desired experimental situation. 
#' This function returns the summary power values for each treatment arm.
#' 
#' Users must specify the desired number of simulations, the group/arm means, 
#' and two of the following: ICC, within-cluster variance, or between-cluster 
#' variance. Significance level, analytic method, progress updates, 
#' poor/singular fit override, and simulated data set output may also be 
#' specified. This function validates the user's input and passes the necessary 
#' arguments to an internal function, which performs the simulations. The 
#' internal function can be called directly by the user to return the fitted 
#' models rather than the power summaries (see \code{?cps.ma.normal.internal}
#' for details).
#' 
#' Users must also supply the number of arms, the subjects per 
#' cluster, and the number of clusters per treatment arm. For a balanced design, 
#' users can provide these values with the arguments \code{narms}, 
#' \code{nsubjects}, and \code{nclusters}, respectively. For unbalanced 
#' designs, the user may provide a list of vectors with one vector per arm,
#' with each vector containing the number of subjects per cluster. See the 
#' examples provided below for a demonstration of the various input options.
#' 
#' @param narms Integer value representing the number of arms. 
#' @param nclusters An integer or vector of integers representing the number 
#' of clusters in each arm.
#' @param ICC The intra-cluster correlation coefficient
#' @param nsim Number of datasets to simulate; accepts integer (required).
#' @param nsubjects Number of subjects per cluster (required); accepts an 
#' integer if all are equal and \code{narms} and \code{nclusters} are provided. 
#' Alternately, the user can supply a list with one entry per arm if the 
#' cluster sizes are the same within the arm, or, if they are not the same 
#' within the arms, the user can supply a list of vectors where each vector 
#' represents an arm and each entry in the vector is the number of subjects 
#' per cluster.
#' @param means Expected absolute treatment effect for each arm; accepts a 
#' vector of length \code{narms} (required).
#' @param sigma_sq Within-cluster variance; accepts a vector of length 
#' \code{narms} (required).
#' @param sigma_b_sq Between-cluster variance; accepts a vector of length 
#' \code{narms} (required).
#' @param alpha Significance level; default = 0.05.
#' @param all.sim.data Option to output list of all simulated datasets; 
#' default = FALSE.
#' @param method Analytical method, either Generalized Linear Mixed Effects 
#' Model (GLMM) or Generalized Estimating Equation (GEE). Accepts c('glmm', 
#' 'gee') (required); default = 'glmm'.
#' @param multi.p.method A string indicating the method to use for adjusting 
#' p-values for multiple comparisons. Choose one of "holm", "hochberg", 
#' "hommel", "bonferroni", "BH", "BY", "fdr", or "none" to leave p-values 
#' unadjusted. The default isc"bonferroni". See \code{?p.adjust} for additional 
#' details.
#' @param quiet When set to FALSE, displays simulation progress and estimated 
#' completion time; default is FALSE.
#' @param seed Option to set.seed. Default is NULL.
#' @param cores a string ("all") or numeric value indicating the number of cores to be 
#' used for parallel computing. 
#' @param poor.fit.override Option to override \code{stop()} if more than 25\% 
#' of fits fail to converge; default = FALSE. 
#' @param low.power.override Option to override \code{stop()} if the power 
#' is less than 0.5 after the first 50 simulations and every ten simulations
#' thereafter. On function execution stop, the actual power is printed in the 
#' stop message. Default = FALSE. When TRUE, this check is ignored and the 
#' calculated power is returned regardless of value. 
#' @param tdist Logical; use t-distribution instead of normal distribution 
#' for simulation values, default = FALSE.
#' @param return.all.models Logical; Returns all of the fitted models, the simulated data,
#' the overall model comparisons, and the convergence report vector. This is equivalent
#' to the output of cps.ma.normal.internal(). See ?cps.ma.normal.internal() for details.
#' @param opt Option to fit with a different optimizer (using the package \textit{optimx}). Default is 'optim'.
#' @return A list with the following components:
#' \describe{
#'   \item{power}{
#'   Data frame with columns "Power" (Estimated statistical power), 
#'                "lower.95.ci" (Lower 95\% confidence interval bound), 
#'                "upper.95.ci" (Upper 95\% confidence interval bound).
#'                }
#'   \item{model.estimates}{
#'   Produced only when all.sim.data=TRUE, data frame with columns 
#'   corresponding to each arm with the suffixes as follows: 
#'                   ".Estimate" (Estimate of treatment effect for a given 
#'                   simulation), 
#'                   "Std.Err" (Standard error for treatment effect estimate), 
#'                   ".tval" (for GLMM) | ".wald" (for GEE), 
#'                   ".pval"
#'                   }
#'   \item{overall.power.table}{
#'   Produced only when all.sim.data=TRUE, table of F-test (when 
#'   method="glmm") or chi-squared (when method="gee") significance test 
#'   results.
#'   }
#'   \item{overall.power.summary}{
#'   Overall power of model compared to H0.
#'   }
#'   \item{simulated.data}{
#'   List of \code{nsim} data frames, each containing: 
#'                   "y" (Simulated response value), 
#'                   "trt" (Indicator for treatment group), 
#'                   "clust" (Indicator for cluster).
#'                   }
#'   \item{model.fit.warning.percent}{
#'   Character string containing the percent of \code{nsim} in which the 
#'   glmm fit was singular or failed to converge, produced only when 
#'   method == "glmm" & all.sim.data==FALSE.
#'   }
#'   \item{model.fit.warning.incidence}{
#'   Vector of length \code{nsim} denoting whether 
#'   or not a simulation glmm fit triggered a "singular fit" 
#'   or "non-convergence" error, produced only when 
#'   method == "glmm" & all.sim.data==TRUE.
#'   }
#'   }
#'          
#' @examples 
#' \dontrun{
#' nsubjects.example <- list(c(20,20,20,25), c(15, 20, 20, 21), c(17, 20, 21))
#' means.example <- c(22, 21, 21.5)
#' sigma_sq.example <- c(1, 1, 0.9)
#' sigma_b_sq.example <- c(0.1, 0.15, 0.1)
#' 
#' multi.cps.normal.unbal <- cps.ma.normal(nsim = 100, nsubjects = nsubjects.example, 
#'                        means = means.example, sigma_sq = sigma_sq.example, 
#'                        sigma_b_sq = sigma_b_sq.example, alpha = 0.05,
#'                        quiet = FALSE, ICC=NULL, method = 'glmm', 
#'                        all.sim.data = FALSE,
#'                        seed = 123, cores = "all",
#'                        poor.fit.override = FALSE)
#'                        
#'  multi.cps.normal <- cps.ma.normal(nsim = 100, narms = 3, 
#'                                    nclusters = c (10,11,10), nsubjects = 100,
#'                                    means = c(21, 21, 21.4),
#'                                    sigma_sq = c(1,1,.9), 
#'                                    sigma_b_sq = c(.1,.15,.1), alpha = 0.05,
#'                                    quiet = FALSE, ICC=NULL, method = 'glmm',
#'                                    all.sim.data = FALSE, seed = 123,
#'                                    poor.fit.override = TRUE, cores="all")
#' }
#' multi.cps.normal.simple <- cps.ma.normal(nsim = 100, narms = 3,
#'                                   nclusters = 10, nsubjects = 25, 
#'                                   means = c(22.1, 21, 22.5),
#'                                   sigma_sq = 1, 
#'                                   sigma_b_sq = 1, alpha = 0.05,
#'                                   quiet = FALSE, ICC=NULL, method = 'glmm',
#'                                   all.sim.data = FALSE, seed = 123,
#'                                   poor.fit.override = TRUE, cores="all")
#' 
#' @author Alexandria C. Sakrejda (\email{acbro0@@umass.edu}), Alexander R. Bogdan, 
#'   and Ken Kleinman (\email{ken.kleinman@@gmail.com})
#' 
#' @export
#' 
 
cps.ma.normal <- function(nsim = 1000, nsubjects = NULL, 
                           narms = NULL, nclusters = NULL,
                        means = NULL, sigma_sq = NULL, 
                        sigma_b_sq = NULL, alpha = 0.05,
                        quiet = FALSE, ICC = NULL, method = 'glmm', 
                        multi.p.method = "bonferroni",
                        all.sim.data = FALSE, seed = NA, 
                        cores=NULL,
                        poor.fit.override = FALSE, 
                        low.power.override = FALSE, 
                        tdist=FALSE,
                        return.all.models = FALSE,
                        opt = "optim"){

  # Create wholenumber function
  is.wholenumber = function(x, tol = .Machine$double.eps^0.5)  abs(x - round(x)) < tol
  
  # create proportion of F-test rejections fxn
  prop_H0_rejection <- function (alpha=alpha, nsim=nsim, LRT.holder.abbrev=LRT.holder.abbrev){
    print(paste("Proportion of F-test rejections = ", 
        round(LRT.holder.abbrev, 3), ", CI:",
        round(LRT.holder.abbrev - abs(stats::qnorm(alpha / 2)) * 
                sqrt((LRT.holder.abbrev * (1 - LRT.holder.abbrev)) / nsim), 3), ", ", 
        round(LRT.holder.abbrev + abs(stats::qnorm(alpha / 2)) * 
                sqrt((LRT.holder.abbrev * (1 - LRT.holder.abbrev)) / nsim), 3), ".", sep=""))
  }

  # create narms and nclusters if not provided directly by user
  if (isTRUE(is.list(nsubjects))) {
    # create narms and nclusters if not supplied by the user
    if (is.null(narms)) {
      narms <- length(nsubjects)
    }
    if (is.null(nclusters)) {
      nclusters <- vapply(nsubjects, length, 0)
    }
  }
  if (length(nclusters) == 1 & !isTRUE(is.list(nsubjects))) {
    nclusters <- rep(nclusters, narms)
  }
  if (length(nclusters) > 1 & length(nsubjects) == 1) {
    narms <- length(nclusters)
  }
    
  # input validation steps
  if(!is.wholenumber(nsim) || nsim < 1 || length(nsim)>1){
    stop("nsim must be a positive integer of length 1.")
  }
  if (is.null(nsubjects)){
    stop("nsubjects must be specified. See ?cps.ma.normal for help.")
  }
  if (length(nsubjects) == 1 & !isTRUE(is.numeric(nclusters))) {
    stop("When nsubjects is scalar, user must supply nclusters (clusters per arm)")
  }
  if (length(nsubjects) == 1 & length(nclusters)==1 & 
      !isTRUE(is.list(narms))) {
    stop("User must provide narms when nsubjects and nclusters are both scalar.")
  }

  validateVariance(dist="norm", 
                   difference=means, alpha=alpha, ICC=ICC, sigma_sq=sigma_sq, 
                   sigma_b_sq=sigma_b_sq, ICC2=NA, sigma_sq2=NA, 
                   sigma_b_sq2=NA, method=method, quiet=quiet, 
                   all.sim.data=all.sim.data, 
                   poor.fit.override=poor.fit.override)

  # nclusters must be positive whole numbers
  if (sum(is.wholenumber(nclusters)==FALSE)!=0 || sum(unlist(nclusters) < 1)!=0){
    stop("nclusters must be postive integer values.")
  }
  #FIXME change this to the following conditions:
  # obs/cluster      ICC       #clusters
  #   <5            <.1          <20
  #   <10           <.05         <20
  #   <20           <.02         <20
  #   <50           <.01         <10
  
  if (nclusters < 20 & sum((createMissingVarianceParam(sigma_b_sq = sigma_b_sq, sigma_sq = sigma_sq) < 0.05) != 0)){
    warning("WARNING: Type 1 error rate increases when ICC and nclusters are small. True power may be lower than estimates.")
  }
  # nsubjects must be positive whole numbers
  if (sum(is.wholenumber(unlist(nsubjects))==FALSE)!=0 || sum(unlist(nsubjects)< 1)!=0){
    stop("nsubjects must be positive integer values.")
  }
  # Create nsubjects structure from narms and nclusters when nsubjects is scalar
  if (length(nsubjects)==1){
    str.nsubjects <- lapply(nclusters, function(x) rep(nsubjects, x))
  } else {
    str.nsubjects <- nsubjects
  }
  
  # allows for means, sigma_sq, sigma_b_sq, and ICC to be entered as scalar
  if (length(sigma_sq)==1){
    sigma_sq <- rep(sigma_sq, narms)
  }
  if (length(sigma_b_sq)==1){
    sigma_b_sq <- rep(sigma_b_sq, narms)
  }
  if (length(ICC)==1){
    ICC <- rep(ICC, narms)
  }
  if (length(means)==1){
    means <- rep(means, narms)
  }

  # supplies sigma_sq or sigma_b_sq if user supplies ICC
  if (length(ICC)!=0){
  if (length(sigma_sq)==0){
    sigma_sq <- createMissingVarianceParam(sigma_b_sq = sigma_b_sq, ICC = ICC)
  }
  if (length(sigma_b_sq)==0){
    sigma_b_sq <- createMissingVarianceParam(sigma_sq = sigma_sq, ICC = ICC)
  }
  }
  
  if (length(sigma_sq)!=narms){
    stop("Length of variance parameters (sigma_sq, sigma_b_sq, ICC) 
         must equal narms, or be provided as a scalar if sigma_sq for all arms are equal.")
  }
  
   # run the simulations 
   normal.ma.rct <- cps.ma.normal.internal(nsim = nsim, 
                                           str.nsubjects = str.nsubjects, 
                                           means = means, sigma_sq = sigma_sq, 
                                           sigma_b_sq = sigma_b_sq, alpha = alpha, 
                                           quiet = quiet, method = method, 
                                           all.sim.data = all.sim.data,
                                           seed = seed,
                                           cores = cores,
                                           poor.fit.override = poor.fit.override,
                                           low.power.override = low.power.override,
                                           tdist = tdist,
                                           opt = opt)
   
   models <- normal.ma.rct[[1]]
   
#Organize output for GLMM
 if(method=="glmm"){
   Estimates = matrix(NA, nrow = nsim, ncol = narms)
   std.error = matrix(NA, nrow = nsim, ncol = narms)
   t.val = matrix(NA, nrow = nsim, ncol = narms)
   p.val = matrix(NA, nrow = nsim, ncol = narms)
 
  if(max(sigma_sq) != min(sigma_sq)) {
    for (i in 1:nsim) {
      Estimates[i,] <- models[[i]][20][[1]][,1]
      std.error[i,] <- models[[i]][20][[1]][,2]
      t.val[i,] <- models[[i]][20][[1]][,4]
      p.val[i,] <- models[[i]][20][[1]][,5]
    }
    keep.names <- rownames(models[[1]][20][[1]])
  } else {
   for (i in 1:nsim){
     Estimates[i,] <- models[[i]][[10]][,1]
     std.error[i,] <- models[[i]][[10]][,2]
     t.val[i,] <- models[[i]][[10]][,4]
     p.val[i,] <- models[[i]][[10]][,5]
   }
    keep.names <- rownames(models[[1]][[10]])
   }
 
 # Organize the row/col names for the model estimates output
 
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
   LRT.holder <- matrix(unlist(normal.ma.rct[[2]]), ncol=6, nrow=nsim, 
                        byrow=TRUE, 
                        dimnames = list(seq(1:nsim), 
                                        c("Sum Sq", "Mean Sq", "NumDF", "DenDF", "F value", "P(>F)")))
   
   # Proportion of times P(>F)
   sig.LRT <-  ifelse(LRT.holder[,6] < alpha, 1, 0)
   LRT.holder.abbrev <- sum(sig.LRT)/nsim
   
 
 # Calculate and store power estimate & confidence intervals
   sig.val <-  ifelse(p.val < alpha, 1, 0)
   pval.power <- apply (sig.val, 2, FUN=function(x) {sum(x, na.rm=TRUE)/nsim})
  # FIXME use binom.test
   power.parms <-  data.frame(Power = round(pval.power, 3),
                          Lower.95.CI = round(pval.power - abs(stats::qnorm(alpha / 2)) * 
                                                  sqrt((pval.power * (1 - pval.power)) / nsim), 3),
                            Upper.95.CI = round(pval.power + abs(stats::qnorm(alpha / 2)) * 
                                                  sqrt((pval.power * (1 - pval.power)) / nsim), 3))
   rownames(power.parms) <- names.power
    
  # Store simulation output in data frame
   ma.model.est <-  data.frame(Estimates, std.error, t.val, p.val)
   ma.model.est <- ma.model.est[, -grep('.*ntercept.*', names(ma.model.est))] 
   
   ## Output objects for GLMM

   # Create list containing all output (class 'crtpwr') and return
   if(all.sim.data == TRUE){
     complete.output <-  list("power" <-  power.parms[-1,],
                              "model.estimates" <-  ma.model.est, 
                              "overall.power" <- LRT.holder,
                              "overall.power2" <- prop_H0_rejection(alpha=alpha, nsim=nsim, LRT.holder.abbrev=LRT.holder.abbrev),
                              "sim.data" <-  normal.ma.rct[[3]], 
                              "failed.to.converge" <-  normal.ma.rct[[4]])
   } 
   if (return.all.models == TRUE) {
     complete.output <-  list("power" <-  power.parms[-1,],
                              "model.estimates" <-  ma.model.est, 
                              "overall.power" <- LRT.holder,
                              "overall.power2" <- prop_H0_rejection(alpha=alpha, nsim=nsim, LRT.holder.abbrev=LRT.holder.abbrev),
                              "all.models" <-  normal.ma.rct)
     } else {
     complete.output <-  list("power" <-  power.parms[-1,],
                              "overall.power" <- prop_H0_rejection(alpha=alpha, nsim=nsim, LRT.holder.abbrev=LRT.holder.abbrev),
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
       Pr[i,] <- p.adjust(models[[i]]$coefficients[,4], method = multi.p.method)
     }
     
     # Organize the row/col names for the output
     keep.names <- rownames(models[[1]]$coefficients)
     
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
                          dimnames = list(seq(1:nsim), 
                                          c("Df", "X2", "P(>|Chi|)")))
     
     # Proportion of times P(>F)
     sig.LRT <-  ifelse(LRT.holder[,3] < alpha, 1, 0)
     LRT.holder.abbrev <- sum(sig.LRT)/nsim
     
     # Calculate and store power estimate & confidence intervals
     sig.val <-  ifelse(Pr < alpha, 1, 0)
     pval.power <- apply (sig.val, 2, FUN=function(x) {sum(x, na.rm=TRUE)/nsim})
     power.parms <-  data.frame(Power = round(pval.power, 3),
                                Lower.95.CI = round(pval.power - abs(stats::qnorm(alpha / 2)) * 
                                                      sqrt((pval.power * (1 - pval.power)) / nsim), 3),
                                Upper.95.CI = round(pval.power + abs(stats::qnorm(alpha / 2)) * 
                                                      sqrt((pval.power * (1 - pval.power)) / nsim), 3))
     rownames(power.parms) <- names.power
     
     # Store GEE simulation output in data frame
     ma.model.est <-  data.frame(Estimates, std.error, Wald, Pr)
     ma.model.est <- ma.model.est[, -grep('.*ntercept.*', names(ma.model.est))] 
   
   ## Output objects for GEE
   
     # Create list containing all output (class 'crtpwr') and return
     if(all.sim.data == TRUE){
       complete.output <-  list("power" <-  power.parms[-1,],
                                "model.estimates" <-  ma.model.est, 
                                "overall.power" <- LRT.holder,
                                "overall.power2" <- prop_H0_rejection(
                                  alpha=alpha, nsim=nsim, 
                                  LRT.holder.abbrev=LRT.holder.abbrev),
                                "sim.data" <-  normal.ma.rct[[3]])
     } else {
       complete.output <-  list("power" <-  power.parms[-1,],
                                "model.estimates" <-  ma.model.est, 
                                "overall.power" <- LRT.holder,
                                "overall.power2" <- prop_H0_rejection(
                                  alpha=alpha, nsim=nsim, 
                                  LRT.holder.abbrev=LRT.holder.abbrev))
     }
     return(complete.output)
   }
}

