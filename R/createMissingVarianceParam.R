#' Calculates ICC, sigma_sq, or sigma_b_sq from the other 2 variables.
#'
#' Usually called from within another function, this function takes 2 of 3 
#' arguments (ICC, sigma_sq, or sigma_b_sq) and returns the missing argument. 
#' Note that it is required that the user specify at least 2 of the arguments.
#' 
#' @param sigma_sq Within-cluster variance; accepts a vector of length \code{narms}.
#' @param sigma_b_sq Between-cluster variance; accepts a vector of length \code{narms}.
#' @param ICC Intra-cluster correlation coefficient; accepts a vector of length \code{narms} with values between 0 - 1.
#' 
#' @return A vector of length \code{narms} containing the argument that the user did not specify.
#' @author Alexandria C. Sakrejda (\email{acbro0@@umass.edu} and Ken Kleinman (\email{ken.kleinman@@gmail.com})
#' @examples 
#' \dontrun{
#' ICC <- createMissingVarianceParam(sigma_sq = c(1, 1, 0.9), sigma_b_sq = c(0.1, 0.15, 0.1))
#' }
#' 
#' @export


## Create missing variance parameters
createMissingVarianceParam <- function(ICC = NULL, 
                                       sigma_sq = NULL,
                                       sigma_b_sq = NULL){
  
  parm1.arg.list <-  list(ICC, sigma_sq, sigma_b_sq)
  parm1.args <-  unlist(lapply(parm1.arg.list, is.null))
  if(sum(parm1.args) > 1){
    stop("At least two of the following terms must be specified: ICC, sigma_sq, sigma_b_sq")
  }
  if(!is.null(c(ICC, sigma_sq)) && is.null(sigma_b_sq)){
    sigma_b_sq = ICC * sigma_sq / (1 - ICC)
    return(sigma_b_sq)
  }
  if(!is.null(c(ICC, sigma_b_sq)) && is.null(sigma_sq)){
    sigma_sq = sigma_b_sq / ICC - sigma_b_sq
    return(sigma_sq)
  }
  if(!is.null(c(sigma_sq, sigma_b_sq)) && is.null(ICC)){
    ICC = sigma_b_sq / (sigma_b_sq + sigma_sq)
    return(ICC)
  }
  if(sum(parm1.args) == 0 && ICC != sigma_b_sq / (sigma_b_sq + sigma_sq)){
    stop("At least one of the following terms has been misspecified: ICC, sigma_sq, sigma_b_sq")
  }
}

