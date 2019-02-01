
validateVariance <- function(difference=means, alpha=alpha, ICC=ICC, sigma=sigma, 
                             sigma_b=sigma_b, ICC2=ICC, sigma2=sigma, 
                             sigma_b2=sigma_b, method=method, quiet=quiet, 
                             all.sim.data=all.sim.data){
  # Validate DIFFERENCE, ALPHA
  min0.warning = " must be a numeric value greater than 0"
  if(!is.numeric(difference) || difference < 0){
    stop("DIFFERENCE", min0.warning)
  }
  if(!is.numeric(alpha) || alpha < 0 || alpha > 1){
    stop("ALPHA must be a numeric value between 0 - 1")
  }

  # Validate ICC, SIGMA, SIGMA_B, ICC2, SIGMA2, SIGMA_B2
  parm1.arg.list = list(ICC, sigma, sigma_b)
  parm1.args = unlist(lapply(parm1.arg.list, is.null))
  if(sum(parm1.args) > 1){
    stop("At least two of the following terms must be specified: ICC, sigma, sigma_b")
  }
  if(sum(parm1.args) == 0 && ICC != sigma_b / (sigma_b + sigma)){
    stop("At least one of the following terms has been misspecified: ICC, sigma, sigma_b")
  }
  parm2.arg.list = list(ICC2, sigma2, sigma_b2)
  parm2.args = unlist(lapply(parm2.arg.list, is.null))
  if(sum(parm2.args) > 1 && sum(parm2.args) != 3){
    stop("At least two of the following terms must be provided to simulate treatment-specific
         variances: ICC2, sigma2, sigma_b2")
  }
  if(sum(parm2.args) == 0 && ICC2 != sigma_b2 / (sigma_b2 + sigma2)){
    stop("At least one of the following terms has been misspecified: ICC2, sigma2, sigma_b2")
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
}