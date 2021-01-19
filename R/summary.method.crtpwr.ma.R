#' A summary method for multi-arm crtpwr objects (crtpwr.ma)
#' 
#' @param object A crtpwr.ma object
#' @param ... other arguments passed to the summary method
#'  
#' @noRd
setMethod(
  "summary",
  signature(object = "crtpwr.ma"),
  definition = function(object, ...)
  {
    cat(paste0("\n", object[['overview']], "\n"))
    cat(paste0("\nPower Estimate (alpha = ", object[['alpha']], "):\n"))
    print(object[['power']], row.names = FALSE)
    cat(paste("\nMethod:", object[['method']], "\n"))
    cat("\nOverall Study Power Estimate:\n")
    print(object[['overall.power2']])
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

