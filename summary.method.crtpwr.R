setMethod("summary",
    signature(object = "crtpwr"),
    function (object, ...) 
    {
      cat(paste0("\n", object[['overview']], "\n"))
      cat(paste0("\nPower Estimate (alpha = ", object[['alpha']], "):\n"))
      print(object[['power']], row.names = FALSE)
      cat(paste0("\nExpected Difference: ", object[['inputs']], "\n"))
      cat("\nVariance Parameters:\n")
      print(object[['variance.parms']])
      cat("\nCluster Sizes:\n")
      cat(paste0("  Non-Treatment (n = ", object[['n.clusters']][1], "): ", format(object[['cluster.sizes']][1]), "\n"))
      cat(paste0("  Treatment (n = ", object[['n.clusters']][2], "): ", format(object[['cluster.sizes']][2]), "\n"))
      cat(paste("\nMethod:", switch(object[['method']], glmm = "Generalized Linear Mixed Model", gee = "Generalized Estimating Equation"), "\n"))
      cat("\n")
    }
)
setClass('crtpwr')
