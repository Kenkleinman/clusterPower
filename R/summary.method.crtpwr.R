#' @export

setMethod(
  "summary",
  signature(object = "crtpwr"),
  definition = function(object, ...)
  {
    cat(paste0("\n", object[['overview']], "\n"))
    cat(paste0("\nPower Estimate (alpha = ", object[['alpha']], "):\n"))
    print(object[['power']], row.names = FALSE)
    #cat("\nType II Error Probability:\n")
    #print(object[['beta']])
    cat(paste("\nMethod:", object[['method']], "\n"))
    if (isTRUE(exists("LRT.holder"))) {
      cat("\nOverall Study Power Estimate:\n")
      print(object[['overall.power2']])
    }
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
#setClass('crtpwr')
