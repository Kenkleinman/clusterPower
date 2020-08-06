#' Updates the Shiny app options if the function arguments change
#'
#'
#' @author Alexandria C. Sakrejda (\email{acbro0@@umass.edu}
#'
#' @param fxnName String. The name of the function from which arguments should
#' be drawn.
#'
#' @return List of numericInput expressions.
#' @export



argMatch <- function(fxnName) {
  require("clusterPower")
  
# compare the function arguments to the objects above  
  temp <- dplyr::intersect(fxnArgs(), names(formals(fxnName)))
  
# store the expressions in a list 
  holder <- list()
  for (i in 1:length(temp)) {
    if (!is.null(temp[i])) {
      holder[[i]] <- eval(parse(text = temp[i]))
    } else {
      next()
    }
  }
  return(holder)
}