#' Power simulations for cluster-randomized trials: calculating ICC, 
#' sigma_sqrd, and sigma_b_sqrd inputs
#'
#' Usually called from within a function, createMissingVarianceParam takes 2 of 3 
#' arguments (ICC, sigma_sqrd, or sigma_b_sqrd) and returns the missing argument. It gives 
#' an error if one of the inputs is missing or not specified properly. Note that 
#' it is required that the user specify at least 2 of the arguments.
#' 
#' @author Alexandria C. Sakrejda (\email{acbro0@@umass.edu}, Alexander R. Bogdan, and Ken Kleinman (\email{ken.kleinman@@gmail.com})
#' 
#' @param sigma_sqrd Within-cluster variance; accepts a vector of length \code{narms}.
#' @param sigma_b_sqrd Between-cluster variance; accepts a vector of length \code{narms}.
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
#' ICC <- createMissingVarianceParam(sigma_sqrd = c(1, 1, 0.9), sigma_b_sqrd = c(0.1, 0.15, 0.1))
#' 
#' @export


## Create missing variance parameters
createMissingVarianceParam <- function(ICC = NULL, 
                                       sigma_sqrd = NULL,
                                       sigma_b_sqrd = NULL){
  parm1.arg.list <-  list(ICC, sigma_sqrd, sigma_b_sqrd)
  parm1.args <-  unlist(lapply(parm1.arg.list, is.null))
  if(sum(parm1.args) > 1){
    stop("At least two of the following terms must be specified: ICC, sigma_sqrd, sigma_b_sqrd")
  }
  if(!is.null(c(ICC, sigma_sqrd)) && is.null(sigma_b_sqrd)){
    sigma_b_sqrd = ICC * sigma_sqrd / (1 - ICC)
    return(sigma_b_sqrd)
  }
  if(!is.null(c(ICC, sigma_b_sqrd)) && is.null(sigma_sqrd)){
    sigma_sqrd = sigma_b_sqrd / ICC - sigma_b_sqrd
    return(sigma_sqrd)
  }
  if(!is.null(c(sigma_sqrd, sigma_b_sqrd)) && is.null(ICC)){
    ICC = sigma_b_sqrd / (sigma_b_sqrd + sigma_sqrd)
    return(ICC)
  }
  if(sum(parm1.args) == 0 && ICC != sigma_b_sqrd / (sigma_b_sqrd + sigma_sqrd)){
    stop("At least one of the following terms has been misspecified: ICC, sigma_sqrd, sigma_b_sqrd")
  }
}

