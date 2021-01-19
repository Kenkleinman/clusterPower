#' A summary method for 2-arm simulation crtpwr objects (crtpwr)
#' 
#' @param object A crtpwr object
#' @param ... other arguments passed to the summary method
#'  
#' 
#' @noRd
setMethod(
  "summary",
  signature(object = "crtpwr"),
  definition = function(object, ...)
  {
    cat(paste0("\n", object[['overview']], "\n"))
    cat(paste0("\nPower Estimate (alpha = ", object[['alpha']], "):\n"))
    print(object[['power']], row.names = FALSE)
    cat(paste("\nMethod:", object[['method']], "\n"))
    cat("\nVariance Parameters:\n")
    print(object[['variance.parms']])
    cat("\nClusters:\n")
    print(object[['n.clusters']])
    cat("\nObservations:\n")
    print(object[['cluster.sizes']])
    cat("\nConvergence:")
    print(table(object[['convergence']]))
  }
)
