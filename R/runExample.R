#' Run a Shiny app for power analysis functions
#'
#' This function runs the clusterPower Shiny apps. 
#' 
#' @param appname which app should be launched. Choices are either "main" or "analytic".
#' Default is "main".
#' 
#' @author Jonathan Moyer (\email{jon.moyer@@gmail.com})
#' @author Alexandria Sakrejda (\email{acbro0@@umass.edu})
#' 
#' @noRd
runExample <- function(appname = "main") {
  # find the examples
  validappnames <- list.files(system.file("shiny-examples", package = "clusterPower"))
  
  validappnamesMsg <-
    paste0(
      "Valid examples are: '",
      paste(validappnames, collapse = "', '"),
      "'")
  
  # find and launch the app
  appDir <- system.file("shiny-examples", appname, package = "clusterPower")
  shiny::runApp(appDir, display.mode = "normal")
}