
##FIXME: TO DO
# 11. testthat tests


#' Power simulations for cluster-randomized trials: Multi-Arm Designs, Binary Outcome.
#'
#' This set of functions utilize iterative simulations to determine 
#' approximate power for multi-arm cluster-randomized controlled trials. Users 
#' can modify a variety of parameters to suit the simulations to their
#' desired experimental situation.
#' 
#' 
#' Users must specify the desired number of simulations, number of subjects per 
#' cluster, number of clusters per treatment arm, group probs, two of the following: ICC, within-cluster variance, or 
#' between-cluster variance. Significance level, analytic method, progress updates, poor/singular fit override,
#' and simulated data set output may also be specified. This function validates the user's input 
#' and passes the necessary arguments to \code{cps.ma.normal.internal}, which performs the simulations.
#' 
#' @author Alexandria C. Sakrejda (\email{acbro0@@umass.edu}, Alexander R. Bogdan, and Ken Kleinman (\email{ken.kleinman@@gmail.com})
#'
#' @param nsim Number of datasets to simulate; accepts integer (required).
#' @param nsubjects Number of subjects per treatment group; accepts a list with one entry per arm. 
#' Each entry is a vector containing the number of subjects per cluster (required).
#' @param probs Expected absolute treatment effect probabilities for each arm; accepts a scalar or a vector of length \code{narms} (required).
#' @param sigma_b_sqrd Between-cluster variance; accepts a vector of length \code{narms} (required).
#' @param alpha Significance level; default = 0.05.
#' @param all.sim.data Option to output list of all simulated datasets; default = FALSE.
#' @param method Analytical method, either Generalized Linear Mixed Effects Model (GLMM) or 
#' Generalized Estimating Equation (GEE). Accepts c('glmm', 'gee') (required); default = 'glmm'.
#' @param quiet When set to FALSE, displays simulation progress and estimated completion time; default is FALSE.
#' @param seed Option to set.seed. Default is NULL.
#' @param poor.fit.override Option to override \code{stop()} if more than 25% of fits fail to converge or 
#' power<0.5 after 50 iterations; default = FALSE.
#' @param overall.power Logical value indicating whether the user would like to return the overall p-value. 
#' The default is FALSE. This option uses \code{pbkrtest::PBmodcomp}, which can take a long time and 
#' provides an approximation based on parametric bootstrapping. There is no reliable alternative method for 
#' binomial outcomes, so proceed with caution if you choose to obtain estimates for overall power.
#' @param cores a string ("all") or scalar numeric value indicating the number of cores to be used for 
#' parallel computing. When this option is set to NULL, no parallel computing is used. Default = NULL.
#'  
#' @return A list with the following components
#' \describe{
#'   \item{power}{Data frame with columns "Power" (Estimated statistical power), 
#'                "lower.95.ci" (Lower 95% confidence interval bound), 
#'                "upper.95.ci" (Upper 95% confidence interval bound)}
#'   \item{model.estimates}{Produced only when all.sim.data=TRUE, data frame with columns corresponding 
#'   to each arm with the suffixes as follows: 
#'                   ".Estimate" (Estimate of treatment effect for a given simulation), 
#'                   "Std.Err" (Standard error for treatment effect estimate), 
#'                   ".tval" (for GLMM) | ".wald" (for GEE), 
#'                   ".pval"
#'   \item{overall.power}{Produced only when all.sim.data=TRUE, table of F-test (when method="glmm") or 
#'   chi^{2} (when method="gee") significance test results.
#'   \item{overall.power2}{Overall power of model compared to H0.}
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
#' probs.example <- c(0.30, 0.21, 0.53)
#' sigma_b_sqrd.example <- c(25, 25, 120)
#' 
#' bin.ma.rct <- cps.ma.binary(nsim = 10, nsubjects = nsubjects.example, 
#'                                      probs = probs.example,
#'                                      sigma_b_sqrd = sigma_b_sqrd.example, alpha = 0.05,
#'                                      quiet = FALSE, method = 'gee', 
#'                                      all.sim.data = FALSE, seed = 123)
#'}
#' 
#' @export
#' 

cps.ma.binary <- function(nsim = 1000, nsubjects = NULL, 
                          narms = NULL, nclusters = NULL,
                          probs = NULL, sigma_b_sqrd = NULL, 
                          alpha = 0.05,
                          quiet = FALSE, ICC=NULL, method = 'glmm', 
                          all.sim.data = FALSE, seed = 123, 
                          cores=NULL,
                          overall.power=FALSE,
                          poor.fit.override = FALSE){

  # Create wholenumber function
  is.wholenumber = function(x, tol = .Machine$double.eps^0.5)  abs(x - round(x)) < tol
 
  # create proportion of F-test rejections fxn
  prop_H0_rejection <- function (alpha=alpha, nsim=nsim, LRT.holder.abbrev=LRT.holder.abbrev, test="F"){
    print(paste("Proportion of ", test, "significance-test rejections = ", 
                round(LRT.holder.abbrev, 3), ", CI:",
                round(LRT.holder.abbrev - abs(stats::qnorm(alpha / 2)) * 
                        sqrt((LRT.holder.abbrev * (1 - LRT.holder.abbrev)) / nsim), 3), ", ", 
                round(LRT.holder.abbrev + abs(stats::qnorm(alpha / 2)) * 
                        sqrt((LRT.holder.abbrev * (1 - LRT.holder.abbrev)) / nsim), 3), ".", sep=""))
  }
   
  # input validation steps
  if(!is.wholenumber(nsim) || nsim < 1 || length(nsim)>1){
    stop("nsim must be a positive integer of length 1.")
  }
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
 
  if(length(nclusters)==1 & (exists("nsubjects", mode = "list")==FALSE)){
    nclusters <- rep(nclusters, narms)
  }
  if(length(nclusters)>1 & length(nsubjects)==1){
    narms <- length(nclusters)
  }
 
  # nclusters must be whole numbers
  if (sum(is.wholenumber(nclusters)==FALSE)!=0 || nclusters < 1){
    stop("nclusters must be postive integer values.")
  }
 
  # nsubjects must be whole numbers
  if (sum(is.wholenumber(unlist(nsubjects))==FALSE)!=0 || unlist(nsubjects) < 1){
    stop("nsubjects must be positive integer values.")
  }

  # Create nsubjects structure from narms and nclusters when nsubjects is scalar
  if (length(nsubjects)==1){
    str.nsubjects <- lapply(nclusters, function(x) rep(nsubjects, x))
  } else {
    str.nsubjects <- nsubjects
  }
 
  # allows for probs, sigma_b_sqrd to be entered as scalar
  if (length(sigma_b_sqrd)==1){
    sigma_b_sqrd <- rep(sigma_b_sqrd, narms)
  }
  if (length(probs)==1){
    probs <- rep(probs, narms)
  }
   
  if (length(sigma_b_sqrd)!=narms){
    stop("Length of variance parameters sigma_b_sqrd
         must equal narms, or be provided as a scalar 
         if sigma_b_sqrd for all arms are equal.")
  }
  
  # run the simulations 
  binary.ma.rct <- cps.ma.binary.internal(nsim = nsim, 
                                          str.nsubjects = str.nsubjects, 
                                          probs = probs,
                                          sigma_b_sqrd = sigma_b_sqrd, alpha = alpha, 
                                          quiet = quiet, method = method, 
                                          all.sim.data = all.sim.data,
                                          seed = seed,
                                          cores = cores,
                                          overall.power=FALSE,
                                          poor.fit.override = poor.fit.override)
  
  models <- binary.ma.rct[[1]]
  print("pass 1")
  #Organize output for GLMM
  if(method=="glmm"){
    Estimates = matrix(NA, nrow = nsim, ncol = narms)
    std.error = matrix(NA, nrow = nsim, ncol = narms)
    z.val = matrix(NA, nrow = nsim, ncol = narms)
    p.val = matrix(NA, nrow = nsim, ncol = narms)
    
    for (i in 1:nsim){
      Estimates[i,] <- models[[i]][[10]][,1]
      std.error[i,] <- models[[i]][[10]][,2]
      z.val[i,] <- models[[i]][[10]][,3]
      p.val[i,] <- models[[i]][[10]][,4]
    }
    
    # Organize the row/col names for the model estimates output
    keep.names <- rownames(models[[1]][[10]])
    
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
    LRT.holder <- matrix(binary.ma.rct[[2]][[1]][[1]][[3]], ncol=2, nrow=nsim, 
                         byrow=TRUE, 
                         dimnames = list(seq(1:nsim), 
                                         rownames(bin.ma.rct[[2]][[1]][[1]])))
    
    # Proportion of times P(>F)
    sig.LRT <-  ifelse(LRT.holder[,1] < alpha, 1, 0)
    LRT.holder.abbrev <- sum(sig.LRT)/nsim
    sig.PBtest <-  ifelse(LRT.holder[,2] < alpha, 1, 0)
    PBtest.holder.abbrev <- sum(sig.PBtest)/nsim
    
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
    ma.model.est <- ma.model.est[, -grep('.*ntercept.*', names(ma.model.est))] 
    
    ## Output objects for GLMM
    
    # Create list containing all output (class 'crtpwr') and return
    if(all.sim.data == TRUE){
      complete.output <-  list("power" <-  power.parms[-1,],
                               "model.estimates" <-  ma.model.est, 
                               "overall.power" <- LRT.holder,
                               "overall.power2" <- paste(prop_H0_rejection(alpha=alpha, nsim=nsim, 
                                                                     LRT.holder.abbrev=LRT.holder.abbrev, test="LRT"),
                                                         "; ", prop_H0_rejection(alpha=alpha, nsim=nsim, 
                                                                                 LRT.holder.abbrev=PBtest.holder.abbrev, test="Bootstrap"),
                                                         sep=""),
                               "sim.data" <-  binary.ma.rct[[3]], 
                               "failed.to.converge" <-  binary.ma.rct[[4]])
    } else {
      complete.output <-  list("power" <-  power.parms[-1,],
                               "overall.power" <- paste(prop_H0_rejection(alpha=alpha, nsim=nsim, 
                                                                    LRT.holder.abbrev=LRT.holder.abbrev, test="LRT"), 
                                                        "; ", prop_H0_rejection(alpha=alpha, nsim=nsim, 
                                                                                LRT.holder.abbrev=PBtest.holder.abbrev, test="Bootstrap"),
                                                        sep=""),
                               "proportion.failed.to.converge" <- binary.ma.rct[[3]])
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
    LRT.holder <- matrix(unlist(binary.ma.rct[[2]]), ncol=3, nrow=nsim, 
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
                               "overall.power2" <- prop_H0_rejection(alpha=alpha, nsim=nsim, LRT.holder.abbrev=LRT.holder.abbrev),
                               "sim.data" <-  binary.ma.rct[[3]])
    } else {
      complete.output <-  list("power" <-  power.parms[-1,],
                               "overall.power" <- prop_H0_rejection(alpha=alpha, nsim=nsim, LRT.holder.abbrev=LRT.holder.abbrev))
    }
    return(complete.output)
  }
  }

