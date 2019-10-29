
optimizerSearch <- function(model){
  require("lme4")
  diff_optims <- allFit(model, maxfun = 1e5, parallel = 'multicore', ncpus = parallel::detectCores())
  is.OK <- sapply(diff_optims, is, "merMod")
  diff_optims.OK <- diff_optims[is.OK]
  if (isTRUE(class(model) == "lmerMod")){
    convergence_results <- lapply(diff_optims.OK,function(x) x@optinfo$conv$lme4$messages)
  }
  if (isTRUE(class(model) == "glmerMod")){
    convergence_results <- lapply(diff_optims.OK,function(x) x@optinfo$conv$lme4$messages)
  }
  working_indices <- sapply(convergence_results, is.null)
  if(sum(working_indices)==0){
    print("No algorithms from allFit converged. Default to ")
    print("You may still be able to use the results, but proceed with extreme caution.")
    if (isTRUE(class(model) == "lmerMod")){
      goodopt <- "nloptwrap"
    } 
    if (isTRUE(class(model) == "glmerMod")){
      goodopt <- "bobyqa"
    }
    } else {
    goodopt <- names(working_indices[working_indices==TRUE][1])
    }
  return(goodopt)
}

