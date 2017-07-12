#' Power calculations for simple cluster randomized trials, continuous outcome
#'
#' Compute the power of a simple cluster randomized trial with a continuous outcome,
#' or determine parameters to obtain a target power.
#'
#' Exactly one of \code{alpha}, \code{power}, \code{m}, \code{n},
#'   \code{nsd}, \code{d}, \code{icc}, and \code{varw}  must be passed as \code{NULL}.
#'   Note that \code{alpha}, \code{power}, and \code{nsd} have non-\code{NULL}
#'   defaults, so if those are the parameters of interest they must be
#'   explicitly passed as \code{NULL}.
#'
#' @section Note:
#'   'uniroot' is used to solve power equation for unknowns, so you may see
#'   errors from it, notably about inability to bracket the root when
#'   invalid arguments are given.
#'
#' @author Jonathan Moyer (\email{jon.moyer@@gmail.com})
#'
#' @param alpha The level of significance of the test, the probability of a
#'   Type I error.
#' @param power The power of the test, 1 minus the probability of a Type II
#'   error.
#' @param m The number of clusters per condition. It must be greater than 1.
#' @param n The mean of the cluster sizes.
#' @param nsd The standard deviation of the cluster sizes.
#' @param d The difference in condition means.
#' @param icc The intraclass correlation.
#' @param varw The within-cluster variation.
#' @param tol Numerical tolerance used in root finding. The default provides
#'   at least four significant digits.
#' @return The computed argument.
#' @references Eldridge SM, Ashby D, Kerry S. Sample size for cluster randomized
#'   trials: effect of coefficient of variation of cluster size and analysis method.
#'  Int J Epidemiol. 2006 Oct;35(5):1292-300.
#' @export

crtpwr.2mean <- function(alpha = 0.05, power = 0.80, m = NULL,
                         n = NULL, nsd = 0,
                         d = NULL, icc = NULL,
                         varw = NULL,
                         tol = .Machine$double.eps^0.25){
  
  if(!is.null(m) && m <= 1) {
    stop("'m' must be greater than 1.")
  }

  # list of needed inputs
  needlist <- list(alpha, power, m, n, nsd, d, icc, varw)
  neednames <- c("alpha", "power", "m", "n", "nsd", "d", "icc", "varw")
  needind <- which(unlist(lapply(needlist, is.null)))
  # check to see that exactly one needed param is null
  
  if (length(needind) != 1) {
    neederror = "Exactly one of 'alpha', 'power', 'm', 'n', 'nsd', 'd', 'icc' and 'varw' must be NULL."
    stop(neederror)
  } 
  
  target <- neednames[needind]
  
  # evaluate power
  p.body <- quote({
    
    # design effect
    DEFF <- 1 + ((((nsd/n)^2 + 1)*n) - 1)*icc
    
    tcrit <- qt(alpha/2, 2*(m - 1), lower.tail = FALSE)
    
    ncp <- sqrt(m*n/(2*DEFF)) * d/sqrt(varw)
    
    pt(tcrit, 2*(m - 1), ncp, lower.tail = FALSE) +
      pt(-tcrit, 2*(m - 1), ncp, lower.tail = TRUE)
  })
  
  # calculate alpha
  if (is.null(alpha)) {
    alpha <- uniroot(function(alpha) eval(p.body) - power,
                     interval = c(1e-10, 1 - 1e-10),
                     tol = tol, extendInt = "yes")$root
  }
  
  # calculate power
  if (is.null(power)) {
    power <- eval(p.body)
  }
  
  # calculate m
  if (is.null(m)) {
    m <- uniroot(function(m) eval(p.body) - power,
                 interval = c(2 + 1e-10, 1e+07),
                 tol = tol, extendInt = "upX")$root
  }
  
  # calculate n
  if (is.null(n)) {
    n <- uniroot(function(n) eval(p.body) - power,
                 interval = c(2 + 1e-10, 1e+07),
                 tol = tol, extendInt = "upX")$root
  }
  
  # calculate nsd
  if (is.null(nsd)) {
    nsd <- uniroot(function(nsd) eval(p.body) - power,
                   interval = c(1e-10, 1e+07),
                   tol = tol, extendInt = "downX")$root
  }

  # calculate d
  if (is.null(d)) {
    d <- uniroot(function(d) eval(p.body) - power,
                 interval = c(1e-07, 1e+07),
                 tol = tol, extendInt = "upX")$root
  }

  # calculate icc
  if (is.null(icc)){
    icc <- uniroot(function(icc) eval(p.body) - power,
                   interval = c(1e-07, 1e+07),
                   tol = tol, extendInt = "downX")$root
  }
  
  # calculate varw
  if (is.null(varw)) {
    varw <- uniroot(function(varw) eval(p.body) - power,
                    interval = c(1e-07, 1e+07),
                    tol = tol, extendInt = "downX")$root
  }

  method <- paste("Clustered two-sample t-test power calculation: ", target, sep = "")
  note <- "'m' is the number of clusters in each group and 'n' is the number of individuals in each cluster."
  structure(list(alpha = alpha, power = power, m = m, n = n, nsd = nsd, d = d,
                 icc = icc, varw = varw, note = note, method = method),
            class = "power.htest")
  
}
