#' Power calculations for difference-in-difference cluster randomized trials, continuous outcome
#'
#' Compute the power of a difference-in-difference cluster randomized trial design with a continuous outcome,
#' or determine parameters to obtain a target power.
#'
#' Exactly one of \code{alpha}, \code{power}, \code{m}, \code{n},
#'   \code{d}, \code{icc}, \code{rho_c}, \code{rho_s}, and \code{varw}
#'   must be passed as \code{NA}. Note that \code{alpha} and\code{power}
#'   have non-\code{NA} defaults, so if those are the parameters of 
#'   interest they must be explicitly passed as \code{NA}.
#'   
#' If \code{n} is a vector the values, \code{m} will be recalculated
#'    using the values in \code{n}. 
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
#' @param n The mean of the cluster sizes, or a vector of cluster sizes for one arm.
#' @param d The difference in mean change between conditions (i.e. "difference-in-difference").
#' @param icc The intraclass correlation.
#' @param rho_c The correlation between baseline and post-test outcomes at the
#'   cluster level. This value can be used in both cross-sectional and cohort
#'   designs. If this quantity is unknown, a value of 0 is a conservative estimate.
#' @param rho_s The correlation between baseline and post-test outcomes at the
#'   subject level. This should be used for a cohort design or a mixture of cohort
#'   and cross-sectional designs. In a purely cross-sectional design (baseline subjects
#'   are completely different from post-test subjects), this value should be 0.
#' @param varw The within-cluster variation.
#' @param tol Numerical tolerance used in root finding. The default provides
#'   at least four significant digits.
#' @return The computed argument.
#' @examples 
#' # Find the number of clusters per condition needed for a trial with alpha = 0.05, 
#' # power = 0.80, n = 100, d = 0.50 units, icc = 0.05, rho_c = 0.50, rho_s = 0.70,
#' # and varw = 1 unit.
#' crtpwr.2meanD(n = 100 , d = 0.5, icc = 0.05, rho_c = 0.50, rho_s = 0.70, varw = 1)
#' # 
#' # The result, m = 4.683358, suggests 5 clusters per condition should be used.
#' 
#' @references Rutterford C, Copas A, Eldridge S. (2015) Methods for sample size
#'   determination in cluster randomized trials. Int J Epidemiol. 44(3):1051-1067.
#' @references Teerenstra S, Eldridge S, Graff M, de Hoop E, Borm, GF. (2012) A simple
#'   sample size formula for analysis of covariance in cluster randomized trials.
#'   Statist Med. 31:2169-2178 
#'   
#' @export

crtpwr.2meanD <- function(alpha = 0.05, power = 0.80, m = NA,
                          n = NA, d = NA, icc = NA,
                          rho_c = NA, rho_s = NA,
                          varw = NA,
                          tol = .Machine$double.eps^0.25){
  
  # if n is a vector, 
  if(length(n) > 1){
    m <- length(n)
  }
  
  if(!is.na(m) && m <= 1) {
    stop("'m' must be greater than 1.")
  }
  
  # list of needed inputs
  needlist <- list(alpha, power, m, n, d, icc, rho_c, rho_s, varw)
  neednames <- c("alpha", "power", "m", "n", "d", "icc",
                 "rho_c", "rho_s", "varw")
  needind <- which(unlist(lapply(needlist, is.na)))
  # check to see that exactly one needed param is NA
  
  if (length(needind) != 1) {
    neederror = "Exactly one of 'alpha', 'power', 'm', 'n', 'd', 'icc', 'rho_c', 'rho_s', and 'varw' must be NA."
    stop(neederror)
  } 
  
  target <- neednames[needind]
  
  # evaluate power
  pwr <- quote({
    
    # variance inflation
    DEFF <- 1 + (n - 1)*icc
    r <- (n*icc/DEFF)*rho_c + ((1 - icc)/DEFF)*rho_s
    # VIF <- DEFF*(1 - r^2)
    VIF <- DEFF*2*(1 - r)

    tcrit <- qt(alpha/2, 2*(m - 1), lower.tail = FALSE)
    
    ncp <- sqrt(m*n/(2*VIF)) * abs(d)/sqrt(varw)
    
    pt(tcrit, 2*(m - 1), ncp, lower.tail = FALSE)#+
    #pt(-tcrit, 2*(m - 1), ncp, lower.tail = TRUE)
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
                        tol = tol, extendInt = "upX")$root
  }
  
  # calculate n
  if (is.na(n)) {
    n <- stats::uniroot(function(n) eval(pwr) - power,
                        interval = c(2 + 1e-10, 1e+07),
                        tol = tol, extendInt = "upX")$root
  }
  
  # calculate d
  if (is.na(d)) {
    d <- stats::uniroot(function(d) eval(pwr) - power,
                        interval = c(1e-07, 1e+07),
                        tol = tol, extendInt = "upX")$root
  }
  
  # calculate icc
  if (is.na(icc)){
    icc <- stats::uniroot(function(icc) eval(pwr) - power,
                          interval = c(1e-07, 1e+07),
                          tol = tol, extendInt = "downX")$root
  }
  
  # calculate rho_c
  if (is.na(rho_c)){
    rho_c <- stats::uniroot(function(rho_c) eval(pwr) - power,
                          interval = c(1e-07, 1e+07),
                          tol = tol, extendInt = "upX")$root
  }
  
  # calculate rho_s
  if (is.na(rho_s)){
    rho_s <- stats::uniroot(function(rho_s) eval(pwr) - power,
                          interval = c(1e-07, 1e+07),
                          tol = tol, extendInt = "upX")$root
  }
  
  # calculate varw
  if (is.na(varw)) {
    varw <- stats::uniroot(function(varw) eval(pwr) - power,
                           interval = c(1e-07, 1e+07),
                           tol = tol, extendInt = "downX")$root
  }
  
  structure(get(target), names = target)
}
