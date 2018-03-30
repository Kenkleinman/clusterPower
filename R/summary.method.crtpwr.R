setMethod("summary",
    signature(object = "crtpwr"),
    function (object, ...) 
    {
      cat(paste0("\n", object[['overview']], "\n"))
      cat(paste0("\nPower Estimate (alpha = ", object[['alpha']], "):\n"))
      print(object[['power']], row.names = FALSE)
      cat("\nInput Values:\n")
      print(object[["inputs"]])
      cat("\nVariance Parameters:\n")
      print(object[['variance.parms']])
      cat("\nCluster Sizes:\n")
      if(length(x[['cluster.sizes']][[1]]) < 10 & length(x[['cluster.sizes']][[1]]) < 10){
        cat(paste0("  Non-Treatment (n = ", x[['n.clusters']][1], "): ", format(x[['cluster.sizes']][[1]]), "\n"))
        cat(paste0("  Treatment (n = ", x[['n.clusters']][2], "): ", format(x[['cluster.sizes']][[2]]), "\n"))
      }else{
        cat("   Non-Treatment (n =", x[['n.clusters']][1], "):", format(x[['cluster.sizes']][[1]][1:10]), "...\n", sep = " ") 
        cat("   Treatment (n =", x[['n.clusters']][2], "):", format(x[['cluster.sizes']][[2]][1:10]), "...\n", sep = " ") 
        #cat(paste0("  Non-Treatment (n = ", x[['n.clusters']][1], "): ", format(x[['cluster.sizes']][[1]][1:10]), " ...\n"))
        #cat(paste0("  Treatment (n = ", x[['n.clusters']][2], "): ", format(x[['cluster.sizes']][[2]][1:10]), "...\n")) 
      }
      cat(paste("\nMethod:", switch(object[['method']], glmm = "Generalized Linear Mixed Model", gee = "Generalized Estimating Equation"), "\n"))
      cat("\n")
    }
)
setClass('crtpwr')
