# Validation functions
validateVariance <- function(x){
  warning("FIXME: not actually validating variance yet")
}

is.wholenumber = function(x, tol = .Machine$double.eps^0.5)  abs(x - round(x)) < tol


##FIXME: TO DO
# 1. update the return values
# 2. update the example/ man text
# 3. input validation
# 4. make validateVariance fxn in "validation" file
# 5. make validate nsubjects fxn in validation file
# 6. make a wrapper function that takes ICC or sigma/sigma_b, nclusters?, narms?, formats output
# 7. make a seperate function for taking ICC, sigma, sigma_b
# 9. write some usage examples
# 10. debug
# 11. testthat tests
# 12. add output element with the pairwise arm comparison p-values.  
# 13. Make sure man page notes that the responsibility for correcting for multiple testing lies with the user.
# 14. Must be able to set the seed on the simulation methods.
# 15. set.seed() option in the wrapper


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
#' @param method Analytical method, either Generalized Linear Mixed Effects Model (GLMM) or 
#' Generalized Estimating Equation (GEE). Accepts c('glmm', 'gee') (required); default = 'glmm'.
#' @param quiet When set to FALSE, displays simulation progress and estimated completion time; default is FALSE.
#' @param all.sim.data Option to output list of all simulated datasets; default = FALSE.
#' 
#' @return A list with the following components
#' \describe{
#'   \item{power}{Data frame with columns "Power" (Estimated statistical power), 
#'                "lower.95.ci" (Lower 95% confidence interval bound), 
#'                "upper.95.ci" (Upper 95% confidence interval bound)}
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
#' @examples 
#' \dontrun{
#' 
#' nsubjects.example <- list(c(20,20,20,25), c(15, 20, 20, 21), c(17, 20, 21))
#' means.example <- c(30, 21, 53)
#' sigma.example <- c(1, 1, 0.9)
#' sigma_b.example <- c(0.1, 0.15, 0.1)
#' 
#' normal.ma.rct <- cps.ma.normal.internal(nsim = 100, nsubjects = nsubjects.example, 
#'                                        means = means.example, sigma = sigma.example, 
#'                                        sigma_b = sigma_b.example, alpha = 0.05, 
#'                                        quiet = FALSE, method = 'glmm', 
#'                                        all.sim.data = FALSE)
#' }
#' 
#' @export

cps.ma.normal.internal = function(nsim = NULL, nsubjects = NULL,
                      means = NULL, sigma = NULL, sigma_b = NULL,
                      alpha = 0.05,
                      all.sim.data = FALSE){

  # Create vectors to collect iteration-specific values
  est.vector = vector(mode = "numeric", length = nsim)
  se.vector = vector(mode = "numeric", length = nsim)
  stat.vector = vector(mode = "numeric", length = nsim)
  pval.vector = vector(mode = "numeric", length = nsim)
  simulated.datasets = list()
  
  # Create NCLUSTERS, NARMS, from NSUBJECTS
  narms = length(nsubjects)
  nclusters = sapply(nsubjects, length)
  
  # FIXME: Validate NSIM, NSUBJECTS
  min1.warning = " must be an integer greater than or equal to 1"
  if(!is.wholenumber(nsim) || nsim < 1){
    stop(paste0("NSIM", min1.warning))
  }

  # FIXME: Validate MEANS, ALPHA
  if(!is.numeric(alpha) || alpha < 0 || alpha > 1){
    stop("ALPHA must be a numeric value between 0 - 1")
  }
  
  # Validate SIGMA, SIGMA_B
  validateVariance(sigma)
  validateVariance(sigma_b)
  
  # Validate ALL.SIM.DATA
  if(!is.logical(all.sim.data)){
    stop("ALL.SIM.DATA must be either TRUE (Output all simulated data sets) or FALSE (No simulated data output")
  }
  
  # Create indicators for treatment group & cluster
  trt = list()
  for (arm in 1:length(nsubjects)){
    trt[[arm]] = list()
    for (cluster in 1:length(nsubjects[[arm]])){
      trt[[arm]][[cluster]] = rep(arm, nsubjects[[arm]][[cluster]])
    }
  }
  clust = list()
  for (arm in 1:length(nsubjects)){
  clust[[arm]] = list()
    for (cluster in 1:length(nsubjects[[arm]])){
    clust[[arm]][[cluster]] = rep(cluster, nsubjects[[arm]][[cluster]])
    }
  }
  
  # Create a container for the simulated.dataset output
  sim.dat = data.frame(y = NA, trt = unlist(trt), clust = unlist(clust))
  
  # Create simulation loop
  for(i in 1:nsim){
    # Generate between-cluster effects for non-treatment and treatment
    randint = mapply(function(nc, s, mu) stats::rnorm(nc, mean = mu, sd = sqrt(s)), 
                                                      nc = nclusters, s = sigma_b, 
                                                      mu = means)
    # Create y-value
    y.bclust = vector(mode = "list", length = narms)
    y.wclust = vector(mode = "list", length = narms)
    y = vector(mode = "list", length = narms)
    for (j in 1:narms){
      y.bclust[[j]] = sapply(1:nclusters[j], function(x) rep(randint[[j]][x], length.out = nsubjects[[j]][x]))
      y.wclust[[j]] = lapply(nsubjects[1:nclusters[j]], function(x) stats:: rnorm(x, mean = randint[j], sd = sqrt(sigma[j])))
      y[[j]] = y.bclust[[j]] + y.wclust[[j]]
    }
    # Create data frame for simulated dataset
    sim.dat[["y"]] = unlist(y)
    
    # Fit GLMM (lmer)
    if(method == 'glmm'){
      my.mod = lme4::lmer(y ~ trt + (1|clust), data = sim.dat)
      glmm.values = summary(my.mod)$coefficient
      p.val = 2 * stats::pt(-abs(glmm.values['trt', 't value']), df = sum(nclusters) - 2)
      est.vector[i] = glmm.values['trt', 'Estimate']
      se.vector[i] = glmm.values['trt', 'Std. Error']
      stat.vector[i] = glmm.values['trt', 't value']
      pval.vector[i] = p.val
      simulated.datasets[[i]] = sim.dat
  }
  
    # FIXME: Not set up for multi-arm RCT yet
  # Fit GEE (geeglm)
  #if(method == 'gee'){
  #  sim.dat = dplyr::arrange(sim.dat, clust)
  #  my.mod = geepack::geeglm(y ~ trt, data = sim.dat,
  #                           id = clust, corstr = "exchangeable")
  #  gee.values = summary(my.mod)$coefficients
  #  est.vector = append(est.vector, gee.values['trt', 'Estimate'])
  #  se.vector = append(se.vector, gee.values['trt', 'Std.err'])
  #  stat.vector = append(stat.vector, gee.values['trt', 'Wald'])
  #  pval.vector = append(pval.vector, gee.values['trt', 'Pr(>|W|)'])
  #}
  
  # Update simulation progress information
  if(quiet == FALSE){
    if(i == 1){
      avg.iter.time = as.numeric(difftime(Sys.time(), start.time, units = 'secs'))
      time.est = avg.iter.time * (nsim - 1) / 60
      hr.est = time.est %/% 60
      min.est = round(time.est %% 60, 0)
      message(paste0('Begin simulations :: Start Time: ', Sys.time(), 
                     ' :: Estimated completion time: ', hr.est, 'Hr:', min.est, 'Min'))
    }
    
    # Iterate progress bar
    prog.bar$update(i / nsim)
    Sys.sleep(1/100)
    
    if(i == nsim){
      total.est = as.numeric(difftime(Sys.time(), start.time, units = 'secs'))
      hr.est = total.est %/% 3600
      min.est = total.est %/% 60
      sec.est = round(total.est %% 60, 0)
      message(paste0("Simulations Complete! Time Completed: ", Sys.time(), 
                     "\nTotal Runtime: ", hr.est, 'Hr:', min.est, 'Min:', sec.est, 'Sec'))
    }
  }
    
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
