#' Updates the Shiny app options if the function arguments change
#'
#'
#' @author Alexandria C. Sakrejda (\email{acbro0@@umass.edu}
#'
#' @param fxnName String. The name of the function from which arguments ahould
#' be drawn.
#'
#' @return List of numericInput expressions.
#' @export



argMatch <- function(fxnName) {
  library(shiny)
  library("clusterPower")
  
  sigma_sq <- numericInput("sigma_sq", "sigma_sq", value = 0.01)
  sigma_b_sq <-
    numericInput("sigma_b_sq", "sigma_b_sq", value = 0.01)
  CV <- numericInput("CV", "CV", value = 0.01)
  d <- numericInput("d", "Means Difference", value = 0.01)
  ICC <-
    numericInput("ICC", "Intracluster Correlation Coefficient (ICC)", value = 0.01)
  vart <-
    numericInput("vart", "Total Variation of the Outcome", value = 0.01)
  seed <- numericInput("seed", , value = 0.01)
  
  temp <- dplyr::intersect(objects(), names(formals(fxnName)))
  
  holder <- list()
  for (i in 1:length(temp)) {
    holder[[i]] <- eval(parse(text = temp[i]))
  }
  return(holder)
}