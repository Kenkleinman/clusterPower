#' Power calculations for simple cluster randomized trials, count outcome
#'
#' Compute the power of a simple cluster randomized trial with a count outcome,
#' or determine parameters to obtain a target power.
#'
#' @section Authors:
#' Jonathan Moyer (\email{jon.moyer@@gmail.com})
#'
#' @param alpha The level of significance of the test, the probability of a
#'   Type I error.
#' @param power The power of the test, 1 minus the probability of a Type II
#'   error.
#' @param m The number of clusters per condition. It must be greater than 1.
#' @param py The number of person-years of observation per cluster.
#' @param r1 The expected mean event rate per unit time in the treatment group.
#' @param r2 The mean event rate per unit time in the control group.
#' @param cvb The between-cluster coefficient of variation.
#' @param tol Numerical tolerance used in root finding. The default provides
#'   at least four significant digits.
#' @return The computed argument.
#' #' #' @examples 
#' # Find the number of clusters per condition needed for a trial with alpha = .05, 
#' # power = 0.8, 10 person-years per cluster, rate
#' in condition 1 of .1 and condition 2 of .2, and cvb = 0.1.
#' crtpwr.2rate(py=10, r1=.1, r2=.2, cvb=.1)
#' # 
#' # The result, showimg m of greater than 24, suggests 25 clusters per condition should be used.
#' @export

crtpwr.2rate<- function(alpha = 0.05, power = 0.80,
                         m = NA, py = NA,
                         r1 = NA, r2 = NA,
                         cvb = NA,
                         tol = .Machine$double.eps^0.25){
  
  if(!is.na(m) && m <= 1) {
    stop("'m' must be greater than 1.")
  }
  
  needlist <- list(alpha, power, m, py, r1, r2, cvb)
  neednames <- c("alpha", "power", "m", "py", "r1", "r2", "cvb")
  needind <- which(unlist(lapply(needlist, is.na))) # find NA index
  
  if (length(needind) != 1) {
    stop("Exactly one of 'alpha', 'power', 'm', 'py', 'r1', 'r2', or 'cvb' must be NA.")
  }
  
  target <- neednames[needind]
  
  pwr <- quote({
    DEFF <- 1 + cvb^2*(r1^2 + r2^2)*py/(r1 + r2)
    zcrit <- qnorm(alpha/2, lower.tail = FALSE)
    vard <- (r1 + r2)*DEFF/py
    pnorm(sqrt((m - 1)*(r1 - r2)^2/vard) - zcrit, lower.tail = TRUE)
  })
  
  # calculate alpha
  if (is.na(alpha)) {
    alpha <- stats::uniroot(function(alpha) eval(pwr) - power,
                     interval = c(1e-10, 1 - 1e-10),
                     tol = tol, extendInt = "yes")$root
  }
  
  # calculate power
  if (is.na(power)) {
    power <- eval(pwr)
  }
  
  # calculate m
  if (is.na(m)) {
    m <- stats::uniroot(function(m) eval(pwr) - power,
                 interval = c(2 + 1e-10, 1e+07),
                 tol = tol)$root
  }
  
  # calculate py
  if (is.na(py)) {
    py <- stats::uniroot(function(py) eval(pwr) - power,
                  interval = c(1e-10, 1e+07),
                  tol = tol, extendInt = "upX")$root
  }
  
  # calculate r1
  if (is.na(r1)) {
      r1 <- stats::uniroot(function(r1) eval(pwr) - power,
                    interval = c(1e-7, 1e7),
                    tol = tol, extendInt = "yes")$root
  }
  
  # calculate r2
  if (is.na(r2)) {
    r2 <- stats::uniroot(function(r2) eval(pwr) - power,
                  interval = c(1e-7, 1e7),
                  tol = tol, extendInt = "yes")$root
  }
  
  # calculate cvb
  if (is.na(cvb)) {
    cvb <- stats::uniroot(function(cvb) eval(pwr) - power,
                  interval = c(1e-7, 1e+07),
                  tol = tol, extendInt = "downX")$root
  }

  structure(get(target), names = target)

}
