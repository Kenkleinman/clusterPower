#' Power simulations for cluster-randomized trials: validating ICC, 
#' sigma, and sigma_b inputs
#'
#' Usually called from within a function, validateVariance takes variance parameters
#' and validates those inputs. This function gives an error if one of the inputs 
#' is missing or not specified properly. 
#' 
#' @author Alexandria C. Sakrejda (\email{acbro0@@umass.edu}, Alexander R. Bogdan, and Ken Kleinman (\email{ken.kleinman@@gmail.com})

#' @param ICC Intra-cluster correlation coefficient; accepts a vector of length \code{narms}
#' with values between 0 - 1.
#' @param difference Expected absolute treatment effect; accepts numeric.
#' @param sigma Within-cluster variance; accepts numeric
#' @param sigma_b Between-cluster variance; accepts numeric
#' @param probs Outcome probabilities per arm. Accepts a numeric vector or a scalar if all arms are equal.
#' @param ICC2 Intra-cluster correlation coefficient for clusters in TREATMENT group
#' @param sigma2 Within-cluster variance for clusters in TREATMENT group
#' @param sigma_b2 Between-cluster variance for clusters in TREATMENT group
#' @param alpha Significance level; default = 0.05.
#' @param method Analytical method, either Generalized Linear Mixed Effects Model (GLMM) or 
#' Generalized Estimating Equation (GEE). Accepts c('glmm', 'gee') (required); default = 'glmm'.
#' @param quiet Logical. When set to FALSE, displays simulation progress and estimated completion time; default is FALSE.
#' @param all.sim.data Logical. Option to output list of all simulated datasets; default = FALSE.
#' @param cores Integer number of cores the user would like to use for parallel computing or the string "all". 
#' 
#' @return A vector of length \code{narms} 
#' \describe{
#'   \item{errors}{Stops function execution if validation fails for any component.}
#' }
#' 
#' @examples 
#' \dontrun{
#' 
#' var <- validateVariance(difference=30, alpha=0.05, ICC=0.2, sigma=100,
#'   sigma_b=25, ICC2=0.5, sigma2=120, sigma_b2=120, method=c("glmm", "gee"), 
#'   quiet=FALSE, all.sim.data=FALSE, poor.fit.override=TRUE)
#' 
#' @export


validateVariance <- function(difference=means, alpha=alpha, ICC=ICC, sigma=sigma_sq, 
                             sigma_b=sigma_b_sq, ICC2=NA, sigma2=NA, 
                             sigma_b2=NA, method=method, quiet=quiet, 
                             all.sim.data=all.sim.data, 
                             poor.fit.override=poor.fit.override, 
                             cores=NA,
                             probs=NA){
  # Validate DIFFERENCE, ALPHA
  min0.warning = " must be a numeric value greater than 0"
  if(!is.numeric(difference) || difference < 0){
    stop("difference", min0.warning)
  }
  if(!is.numeric(alpha) || alpha < 0 || alpha > 1){
    stop("alpha must be a numeric value between 0 - 1")
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
  
  if (!is.na(probs)){
    if(is.numeric(probs)==FALSE){
      stop("'probs' must be a numeric scalar or vector.")
    }
  }
  
  if (!is.na(ICC2) | !is.na(sigma2) | !is.na(sigma_b2)){
    parm2.arg.list = list(ICC2, sigma2, sigma_b2)
    parm2.args = unlist(lapply(parm2.arg.list, is.null))
    if(sum(parm2.args) > 1 && sum(parm2.args) != 3){
      stop("At least two of the following terms must be provided to simulate treatment-specific
           variances: ICC2, sigma2, sigma_b2")
    }
    if(sum(parm2.args) == 0 && ICC2 != sigma_b2 / (sigma_b2 + sigma2)){
      stop("At least one of the following terms has been misspecified: ICC2, sigma2, sigma_b2")
    }
  }

  # Validate METHOD, QUIET, ALL.SIM.DATA
  if(!is.element(method, c('glmm', 'gee'))){
    stop("method must be either 'glmm' (Generalized Linear Mixed Model)
         or 'gee'(Generalized Estimating Equation)")
  }
  if(!is.logical(quiet)){
    stop("quiet must be either TRUE (No progress information shown) or FALSE (Progress information shown)")
  }
  if(!is.logical(all.sim.data)){
    stop("all.sim.data must be either TRUE (Output all simulated data sets) or FALSE (No simulated data output")
  }
  if(!is.logical(poor.fit.override)){
    stop("poor.fit.override must be either FALSE (Stop simulations if estimated power is <0.5 or 
         more than 25% of fits do not converge) or TRUE (Continue simulations).")
  }
  if (!is.na(cores)){
    if(!is.integer(cores) & cores!="all"){
      stop("cores must be an integer, or 'all'.")
    }
    if (is.integer(cores)==TRUE){
      if (cores>parallel::detectCores()){
        stop("'cores' exceeds the number of available cores.")
      }
    }
  }
}

