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
#'   vector with the first item being the name of the returned values.}
#' }
#' 
#' @examples 
#' \dontrun{
#' 
#' ICC <- createMissingVarianceParam(sigma = c(1, 1, 0.9), sigma_b = c(0.1, 0.15, 0.1))
#' 
#' @export


## Create missing variance parameters
createMissingVarianceParam <- function(ICC = NULL, 
                                       sigma = NULL,
                                       sigma_b = NULL){
  parm1.arg.list <-  list(ICC, sigma, sigma_b)
  parm1.args <-  unlist(lapply(parm1.arg.list, is.null))
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

