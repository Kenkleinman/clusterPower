#' Run a Shiny app for power analysis functions
#'
#' This function runs the clusterPower Shiny apps. 
#' 
#' @author Jonathan Moyer (\email{jon.moyer@@gmail.com})
#' @author Alexandria Sakrejda (\email{acbro0@@umass.edu})
#' 
#' @export
runExample <- function(appname = "main") {
  # find the examples
  validappnames <- list.files(system.file("shiny-examples", package = "clusterPower"))
  
  validappnamesMsg <-
    paste0(
      "Valid examples are: '",
      paste(validappnames, collapse = "', '"),
      "'")
  
  # possible errors
 # if (missing(appname) || !nzchar(appname) ||
#      !appname %in% validappnames) {
#    stop(
#      'Please run `runExample()` with a valid example app as an argument.\n',
#      validappnamesMsg,
#      call. = FALSE)
#  }
  
  # find and launch the app
  appDir <- system.file("shiny-examples", appname, package = "clusterPower")
  shiny::runApp(appDir, display.mode = "normal")
}