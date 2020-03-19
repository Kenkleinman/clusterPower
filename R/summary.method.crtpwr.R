#' @export

setMethod("summary",
    signature(object = "crtpwr"),
    definition = function (object, ...) 
    {
      cat(paste0("\n", object[['overview']], "\n"))
      cat(paste0("\nPower Estimate (alpha = ", object[['alpha']], "):\n"))
      print(object[['power']], row.names = FALSE)
      cat(paste("\nMethod:", object[['method']], "\n"))
      cat("\nVariance Parameters:\n")
      print(object[['variance.parms']])
      cat("\nCluster Sizes:\n")
  #    if(length(object[['cluster.sizes']][[1]]) < 10 & length(object[['cluster.sizes']][[2]]) < 10){
  #      cat(paste0("  Non-Treatment (n = ", object[['n.clusters']][1], "): ", format(object[['cluster.sizes']][[1]]), "\n"))
  #      cat(paste0("  Treatment (n = ", object[['n.clusters']][2], "): ", format(object[['cluster.sizes']][[2]]), "\n"))
  #    }else{
  #      cat("   Non-Treatment (n =", object[['n.clusters']][1], "):", format(object[['cluster.sizes']][[1]][1:10]), "...\n", sep = " ") 
  #      cat("   Treatment (n =", object[['n.clusters']][2], "):", format(object[['cluster.sizes']][[2]][1:10]), "...\n", sep = " ") 
  #    }
      cat("\nArms:")
      print(table(object[['narms']]))
      cat("\nClusters:")
      print(table(object[['n.clusters']]))
      print(table(object[['cluster.sizes']]))
      cat("\nArms:")
      print(table(object[['narms']]))
      cat("\nConvergence:")
      print(table(object[['convergence']]))
    }
)
#setClass('crtpwr')
