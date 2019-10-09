#' Power simulations for cluster-randomized trials: validating ICC, 
#' sigma_sq, and sigma_b_sq inputs.
#'
#' Usually called from within a simulation-producing function, this function 
#' takes variance parameters provided by the user and performs a series of 
#' validation checks. This function gives an error and stops function execution 
#' if any required input is missing, is the wrong type, or not specified 
#' properly. 
#' 
#' @author Alexandria C. Sakrejda (\email{acbro0@@umass.edu}, and Ken Kleinman (\email{ken.kleinman@@gmail.com})

#' @param ICC Intra-cluster correlation coefficient; accepts a vector of 
#' length \code{narms} with values between 0 - 1.
#' @param dist String indicating outcome distribution. Takes one of "norm", 
#' "bin", or "count".
#' @param difference Expected absolute treatment effect; accepts numeric.
#' @param sigma_sq Within-cluster variance; accepts numeric
#' @param sigma_b_sq Between-cluster variance; accepts numeric
#' @param probs Outcome probabilities per arm. Accepts a numeric vector or a 
#' scalar if all arms are equal.
#' @param ICC2 Intra-cluster correlation coefficient for clusters in the 
#' treatment group (used in 2-arm simulation methods only).
#' @param sigma_sq2 Within-cluster variance for clusters in treatment group
#' (used in 2-arm simulation methods only).
#' @param sigma_b_sq2 Between-cluster variance for clusters in treatment group 
#' (used in 2-arm simulation methods only).
#' @param alpha Significance level; default = 0.05.
#' @param method Analytical method, either Generalized Linear Mixed Effects 
#' Model (GLMM) or Generalized Estimating Equation (GEE). Accepts c('glmm', 
#' 'gee') (required); default = 'glmm'.
#' @param quiet Logical. When set to FALSE, displays simulation progress and 
#' estimated completion time; default is FALSE.
#' @param all.sim.data Logical. Option to output list of all simulated 
#' datasets; default = FALSE.
#' @param cores Integer number of cores the user would like to use for 
#' parallel computing or the string "all". 
#' 
#' @return A vector of length \code{narms} 
#' \describe{
#'   \item{errors}{Stops function execution if validation fails for any 
#'   component, but no return if input has no errors.}
#' }
#' 
#' @examples 
#' \dontrun{
#' 
#' validateVariance(dist="norm", difference=30, alpha=0.05, ICC=0.2, 
#'                   sigma_sq=100, sigma_b_sq=25, ICC2=0.5, sigma_sq2=120, 
#'                   sigma_b_sq2=120, method=c("glmm", "gee"), quiet=FALSE, 
#'                   all.sim.data=FALSE, poor.fit.override=TRUE)
#' 
#' }
#' @export


validateVariance <- function(dist=NULL, difference=NULL, alpha=alpha, 
                             ICC=ICC, sigma_sq=sigma_sq, 
                             sigma_b_sq=sigma_b_sq, ICC2=NA, sigma_sq2=NA, 
                             sigma_b_sq2=NA, method=method, quiet=quiet, 
                             all.sim.data=all.sim.data, 
                             poor.fit.override=poor.fit.override, 
                             cores=NA,
                             probs=NA){
  if(!is.element(dist, c('norm', 'bin', 'count'))){
    stop("dist must be specified as 'norm', 'bin', or 'count'.")
  }
  if (dist=="norm"){
  # Validate DIFFERENCE, ALPHA
  min0.warning = " must be a numeric value"
  if(!is.numeric(difference)){
    stop("difference or means", min0.warning)
  }

  # Validate ICC, sigma_sq, sigma_b_sq, ICC2, sigma_sq2, sigma_b_sq2
  if(length(difference)==1){
    parm1.arg.list = list(ICC, sigma_sq, sigma_b_sq)
    parm1.args = unlist(lapply(parm1.arg.list, is.null))
    if(sum(parm1.args) > 1){
      stop("At least two of the following terms must be specified: ICC, sigma_sq, sigma_b_sq")
    }
    if(sum(parm1.args) == 0 && ICC != sigma_b_sq / (sigma_b_sq + sigma_sq)){
      stop("At least one of the following terms has been misspecified: ICC, sigma_sq, sigma_b_sq")
    }
  }
  if (!is.na(ICC2) | !is.na(sigma_sq2) | !is.na(sigma_b_sq2)){
    parm2.arg.list = list(ICC2, sigma_sq2, sigma_b_sq2)
    parm2.args = unlist(lapply(parm2.arg.list, is.null))
    if(sum(parm2.args) > 1 && sum(parm2.args) != 3){
      stop("At least two of the following terms must be provided to simulate treatment-specific
           variances: ICC2, sigma_sq2, sigma_b_sq2")
    }
    if(sum(parm2.args) == 0 && ICC2 != sigma_b_sq2 / (sigma_b_sq2 + sigma_sq2)){
      stop("At least one of the following terms has been misspecified: ICC2, sigma_sq2, sigma_b_sq2")
    }
  }
  }
  
  if (dist=="bin"){
      if(any(is.numeric(probs)==FALSE) | any(is.na(probs))){
        stop("probs must be a scalar or numeric vector.")
      }
      if (any(probs >=1)){
          stop("probabilities must be less than 1.")
        }
      if (any(probs <=0)){
          stop("probabilities must be greater than zero.")
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
  #validate alpha
  if(!is.numeric(alpha) || alpha < 0 || alpha > 1){
    stop("alpha must be a numeric value between 0 - 1")
  }
  #validate cores
    if(!is.numeric(cores) & cores!="all" & !is.na(cores)){
      stop("cores must be NA, an integer, or 'all'.")
    }
    if (is.numeric(cores)==TRUE){
      if (cores>parallel::detectCores()){
        stop("'cores' exceeds the number of cores detected.")
      }
    }
}

