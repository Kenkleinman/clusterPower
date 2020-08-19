#' Power simulations for cluster-randomized trials: check if whole number.
#' 
#' 
#' @author Alexandria C. Sakrejda (\email{acbro0@@umass.edu} and Ken Kleinman (\email{ken.kleinman@@gmail.com})
#' 
#' @param x The numeric value to be tested.
#' @param tol The machine tolerance. Defaults to \code{.Machine$double.eps ^ 0.5}.
#' 
#' @return Logical.
#' @export

is.wholenumber <- function(x, tol = .Machine$double.eps ^ 0.5){
  abs(x - round(x)) < tol
}