#' Power simulations for cluster-randomized trials: find an optimizer.
#'
#' This function can be called by some of the clusterPower functions to search
#' for an optimizer that does not produce convergence errors.
#'
#' @author Alexandria C. Sakrejda (\email{acbro0@@umass.edu}, and Joshua Nugent.

#' @param model a lme4 or glmer model
#'
#' @return A string
#' \describe{
#'   \item{goodopt}{The name of the first optimizer tested that showed convergence.}
#' }
#'
#' @export


optimizerSearch <- function(model) {
  if (isTRUE(missing(model))) {
    message("Error: model does not exist.")
  }
  diff_optims <- try(lme4::allFit(
    model,
    maxfun = 1e5,
    parallel = 'multicore',
    ncpus = parallel::detectCores()
  ))
  if (class(diff_optims) == "try-error") {
    stop("optimizerSearch is not compatible with this model structure.")
  }
  is.OK <- sapply(diff_optims, is, "merMod")
  diff_optims.OK <- diff_optims[is.OK]
  convergence_results <- lapply(diff_optims.OK,
                                function(x)
                                  x@optinfo$conv$lme4$messages)
  working_indices <- sapply(convergence_results, is.null)
  if (isTRUE(sum(unlist(working_indices), na.rm = TRUE) == 0)) {
    print("No algorithms from allFit converged. Default to nloptwrap")
    print("You may still be able to use the results, but proceed with extreme caution.")
    goodopt <- "nloptwrap"
  } else {
    goodopt <- names(working_indices[working_indices == TRUE][1])
  }
  return(goodopt)
}
