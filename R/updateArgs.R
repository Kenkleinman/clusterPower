#' Concatenates the Shiny app options if the function arguments change
#'
#'
#' @author Alexandria C. Sakrejda (\email{acbro0@@umass.edu}
#'
#' @param fxnName String. The name of the function from which arguments should
#' be drawn.
#'
#' @return List of numericInput expressions.
#' @export


updateArgs <- function(fxnName) {
  argNames <-
    c("nsubjects",
      "nclusters",
      "alpha",
      dplyr::intersect(c(fxnArgs(), "power"), names(formals(fxnName))))
  holder <- as.vector(NA, mode = "character")
  arghelper <- function(argname) {
    x <- paste0(argname, "=input$", argname)
    return(x)
  }
  for (i in 1:length(argNames)) {
    holder[i] <- arghelper(argNames[i])
  }
  return(cat(holder, sep = ", "))
}
