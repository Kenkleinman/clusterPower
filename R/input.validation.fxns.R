


## Create missing variance parameters
createMissingVarianceParam <- function(ICC = NULL, 
                                        sigma = NULL,
                                        sigma_b = NULL){
  parm1.arg.list = list(ICC, sigma, sigma_b)
  parm1.args = unlist(lapply(parm1.arg.list, is.null))
  if(sum(parm1.args) > 1){
    stop("At least two of the following terms must be specified: ICC, sigma, sigma_b")
  }
  if(!is.null(c(ICC, sigma)) && is.null(sigma_b)){
    sigma_b = ICC * sigma / (1 - ICC)
    return(sigma_b)
    }
  if(!is.null(c(ICC, sigma_b)) && is.null(sigma)){
    sigma = sigma_b / ICC - sigma_b
    return(sigma)
    }
  if(!is.null(c(sigma, sigma_b)) && is.null(ICC)){
    ICC = sigma_b / (sigma_b + sigma)
    return(ICC)
  }
  if(sum(parm1.args) == 0 && ICC != sigma_b / (sigma_b + sigma)){
    stop("At least one of the following terms has been misspecified: ICC, sigma, sigma_b")
  }
}
















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