#' Power simulations for cluster-randomized trials: calculating and validating ICC, 
#' sigma, and sigma_b inputs
#'
#' Usually called from within a function, createMissingVarianceParam takes 2 of 3 
#' arguments (ICC, sigma, or sigma_b) and returns the missing argument. It gives 
#' an error if one of the inputs is missing or not specified properly. Note that 
#' it is required that the user specify at least 2 of the arguments.
#' 
#' @author Alexandria C. Sakrejda

#' @param sigma Within-cluster variance; accepts a vector of length \code{narms}.
#' @param sigma_b Between-cluster variance; accepts a vector of length \code{narms}.
#' @param ICC Intra-cluster correlation coefficient; accepts a vector of length \code{narms}
# with values between 0 - 1.

#' @return A vector of length \code{narms} 
#' \describe{
#'   \item{missing.arg}{The argument that the user did not specify is returned as a character 
#'   vector with the first item being the name of the}
#' }
#' 
#' @examples 
#' \dontrun{
#' 
#' ICC <- createMissingVarianceParam(sigma = c(1, 1, 0.9), sigma_b = c(0.1, 0.15, 0.1))
#' 
#' @export

















validateVariance <- function(x){
  if (exists(nsim)){ 
# validate nsim
    min1.warning = " must be an integer greater than or equal to 1"
    if(!is.wholenumber(nsim) || nsim < 1){
      stop(paste0("nsim", min1.warning))
    }
  }
  # validate nsubjects
  if (exists(nsubjects)){
    if(!is.numeric(unlist(nsubjects)) || length(unlist(nsubjects)[unlist(nsubjects)<0])!=0 || 
       !is.wholenumber(sum(unlist(nsubjects)))){
      stop("nsubjects must be postive integers")
    }
  }
    # FIXME: Validate MEANS, ALPHA
    if (exists(alpha)){
      if(!is.numeric(alpha) || alpha < 0 || alpha > 1){
      stop("ALPHA must be a numeric value between 0 - 1")
        }
    }
  if (exists(means)){
    if(!is.numeric(means) || ##FIXME, not sure how to proceed until wrapper is written
  }
  
  
  # Validate METHOD, QUIET, ALL.SIM.DATA
  if(!is.element(method, c('glmm', 'gee'))){
    stop("METHOD must be either 'glmm' (Generalized Linear Mixed Model)
         or 'gee'(Generalized Estimating Equation)")
  }
  if(!is.logical(quiet)){
    stop("QUIET must be either TRUE (No progress information shown) or FALSE (Progress information shown)")
  }
  if(!is.logical(all.sim.data)){
    stop("ALL.SIM.DATA must be either TRUE (Output all simulated data sets) or FALSE (No simulated data output")
  }
  
  
  
  
  
  
  
    # Validate SIGMA, SIGMA_B
    validateVariance(sigma)
    validateVariance(sigma_b)
    
    # Validate ALL.SIM.DATA
    if(!is.logical(all.sim.data)){
      stop("ALL.SIM.DATA must be either TRUE (Output all simulated data sets) or FALSE (No simulated data output")
    }
  }
  warning("FIXME: not actually validating variance yet")
}









cps.ma.normal.internal = function(nsim = NULL, nsubjects = NULL,
                                  means = NULL, sigma = NULL, sigma_b = NULL,
                                  alpha = 0.05,
                                  quiet = FALSE, method = 'glmm', 
                                  all.sim.data = FALSE){
  

  
  # Create NCLUSTERS, NARMS, from NSUBJECTS
  narms = length(nsubjects)
  nclusters = sapply(nsubjects, length)
 
  

  
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

  
  is.wholenumber = function(x, tol = .Machine$double.eps^0.5)  abs(x - round(x)) < tol
}