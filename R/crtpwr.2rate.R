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
#' @param l1 The expected mean event rate per unit time in the treatment group.
#' @param l2 The mean event rate per unit time in the control group.
#' @param cvb The between-cluster coefficient of variation.
#' @param tol Numerical tolerance used in root finding. The default provides
#'   at least four significant digits.
#' @return The computed argument.
#' @export

crtpwr.2prop <- function(alpha = 0.05, power = 0.80,
                         m = NA, py = NA,
                         l1 = NA, l2 = NA,
                         cvb = NA,
                         tol = .Machine$double.eps^0.25){
  
  if(!is.na(m) && m <= 1) {
    stop("'m' must be greater than 1.")
  }
  
  needlist <- list(alpha, power, m, n, py, l1, l2, cvb)
  neednames <- c("alpha", "power", "m", "py", "l1", "l2", "cvb")
  needind <- which(unlist(lapply(needlist, is.na))) # find NA index
  
  if (length(needind) != 1) {
    stop("Exactly one of 'alpha', 'power', 'm', 'py', 'l1', 'l2', or 'cvb' must be NA.")
  }
  
  target <- neednames[needind]
  
  p.body <- quote({
    DEFF <- 1 + cvb^2*(l1^2 + l2^2)*py/(l1 + l2)
    zcrit <- qnorm(alpha/2, lower.tail = FALSE)
    vard <- (l1 + l2)*DEFF/py
    pnorm(sqrt((m - 1)*(l1 - l2)^2/vard) - zcrit, lower.tail = FALSE)
  })
  
  # calculate alpha
  if (is.na(alpha)) {
    alpha <- uniroot(function(alpha) eval(p.body) - power,
                     interval = c(1e-10, 1 - 1e-10),
                     tol = tol, extendInt = "yes")$root
  }
  
  # calculate power
  if (is.na(power)) {
    power <- eval(p.body)
  }
  
  # calculate m
  if (is.na(m)) {
    m <- uniroot(function(m) eval(p.body) - power,
                 interval = c(2 + 1e-10, 1e+07),
                 tol = tol)$root
  }
  
  # calculate py
  if (is.na(n)) {
    py <- uniroot(function(py) eval(p.body) - power,
                  interval = c(1e-10, 1e+07),
                  tol = tol, extendInt = "upX")$root
  }
  
  # calculate l1
  if (is.na(l1)) {
      l1 <- uniroot(function(l1) eval(p.body) - power,
                    interval = c(1e-7, 1e7),
                    tol = tol, extendInt = "yes")$root
  }
  
  # calculate l2
  if (is.na(l2)) {
    l1 <- uniroot(function(l2) eval(p.body) - power,
                  interval = c(1e-7, 1e7),
                  tol = tol, extendInt = "yes")$root
  }
  
  # calculate cvb
  if (is.na(cvb)) {
    cv <- uniroot(function(cvb) eval(p.body) - power,
                  interval = c(1e-7, 1e+07),
                  tol = tol, extendInt = "downX")$root
  }

  structure(get(target), names = target)
  
  # method <- paste("Clustered two-sample proportion power calculation: ", target, sep = "")
  # note <- "'m' is the number of clusters in each group and 'n' is the number of individuals in each cluster."
  # structure(list(m = m, n = n, cv = cv,
  #                p1 = p1, p1dec = p1dec, p1inc = p1inc,
  #                p2 = p2, p2dec = p2dec, p2inc = p2inc,
  #                icc = icc, alpha = alpha, power = power,
  #                note = note, method = method),
  #           class = "power.htest")
  
}
