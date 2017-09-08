#' Shiny app for running 
#'
#' This function runs a Shiny app that employs the \code{crtpwr.2mean},
#' \code{crtpwr.2prop}, and \code{crtpwr.2rate}. The app can do calculations
#' for several sets of parameters and produces graphs.
#' 
#' @section Authors:
#' Jonathan Moyer (\email{jon.moyer@@gmail.com})

#' @export
runExample <- function(){
  appDir <- system.file("shiny-examples", "app", package = "clusterPower")
  if (appDir == "") {
    stop("Could not find example directory. Try re-installing `mypackage`.", call. = FALSE)
  }
  shiny::runApp(appDir, display.mode = "normal")
}