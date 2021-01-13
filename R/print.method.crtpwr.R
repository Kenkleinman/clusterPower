#' A print method for crtpwr objects
#' 
#' @param x A crtpwr object
#' 
#' @param ... Other arguments passed to the function 
#'

print.crtpwr <- function(x, ...) {
  cat(paste0("\n", x[['overview']], "\n"))
  cat(paste0("\nPower Estimate (alpha = ", x[['alpha']], "):\n"))
  print(x[['power']], row.names = FALSE)
  cat("\n")
}
