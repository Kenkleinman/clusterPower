#' Power calculations for  difference-in-difference cluster randomized trials, dichotomous outcome
#'
#' Compute the power of a difference-in-difference cluster randomized trial design with a 
#' binary outcome, or determine parameters to obtain a target power.
#' 
#' Exactly one of \code{alpha}, \code{power}, \code{nclusters}, \code{nsubjects},
#'   \code{p}, \code{d}, \code{ICC}, \code{rho_c}, and \code{rho_s} must be passed as \code{NA}.
#'   Note that \code{alpha} and \code{power} have non-\code{NA}
#'   defaults, so if those are the parameters of interest they must be
#'   explicitly passed as \code{NA}.
#'
#' @section Authors:
#' Jonathan Moyer (\email{jon.moyer@@gmail.com}), Ken Kleinman (\email{ken.kleinman@@gmail.com})
#' 
#' @section Note:
#'   This function was inspired by work from Stephane Champely (pwr.t.test) and
#'   Peter Dalgaard (power.t.test). As with those functions, 'uniroot' is used to
#'   solve power equation for unknowns, so you may see
#'   errors from it, notably about inability to bracket the root when
#'   invalid arguments are given.
#'   
#' @param alpha The level of significance of the test, the probability of a
#'   Type I error.
#' @param power The power of the test, 1 minus the probability of a Type II
#'   error.
#' @param nclusters The number of clusters per condition. It must be greater than 1.
#' @param nsubjects The mean of the cluster sizes.
#' @param p The expected mean proportion at the post-test, averaged across both arms.
#' @param d The expected absolute difference.
#' @param ICC The intraclass correlation.
#' @param rho_c The correlation between baseline and post-test outcomes at the
#'   cluster level. This value can be used in both cross-sectional and cohort
#'   designs. If this quantity is unknown, a value of 0 is a conservative estimate.
#' @param rho_s The correlation between baseline and post-test outcomes at the
#'   subject level. This should be used for a cohort design or a mixture of cohort
#'   and cross-sectional designs. In a purely cross-sectional design (baseline subjects
#'   are completely different from post-test subjects), this value should be 0.
#' @param tol Numerical tolerance used in root finding. The default provides
#'   at least four significant digits.
#' @return The computed argument.
#' @examples 
#' # Find the number of clusters per condition needed for a trial with alpha = .05, 
#' # power = 0.8, 50 observations per cluster, expected mean post-test proportion of .50,
#' # expected difference of .1, ICC = 0.05, cluster level correlation of 0.3, and subject level 
#' # correlation of 0.7.
#' cpa.did.binary(nsubjects=50 ,p=.5, d=.1, ICC=.05, rho_c=.3, rho_s=.7)
#' # 
#' # The result, showimg nclusters of greater than 32, suggests 33 clusters per 
#' # condition should be used.
#' 
#' @references Murray D. Design and Analysis of Group-Randomized Trials. New York, NY: Oxford
#' University Press; 1998.
#' 
#' @export


cpa.did.binary <- function(alpha = 0.05,
                           power = 0.80,
                           nclusters = NA,
                           nsubjects = NA,
                           p = NA,
                           d = NA,
                           ICC = NA,
                           rho_c = NA,
                           rho_s = NA,
                           tol = .Machine$double.eps ^ 0.25) {
  if (!is.na(nclusters) && nclusters <= 1) {
    stop("'nclusters' must be greater than 1.")
  }
  
  # check to see that difference is positive
  if (!is.na(d)) {
    if (d <= 0) {
      stop("'d' must be greater than 0.")
    }
  }
  
  # check to see if d is greater than p or 1 - p
  if (!is.na(p) && !is.na(d)) {
    if (d > p | d > 1 - p) {
      stop("'d' cannot be greater than 'p' or 1 - 'p'.")
    }
  }
  
  needlist <-
    list(alpha, power, nclusters, nsubjects, p, d, ICC, rho_c, rho_s)
  neednames <-
    c("alpha",
      "power",
      "nclusters",
      "nsubjects",
      "p",
      "d",
      "ICC",
      "rho_c",
      "rho_s")
  needind <- which(unlist(lapply(needlist, is.na))) # find NA index
  
  if (length(needind) != 1) {
    stop(
      "Exactly one of 'alpha', 'power', 'nclusters', 'nsubjects', 'p', 'd', 'ICC', 'rho_c', or 'rho_s' must be NA."
    )
  }
  
  target <- neednames[needind]
  
  pwr <- quote({
    tcrit <- qt(alpha / 2, 2 * (nclusters - 1), lower.tail = FALSE)
    
    # variance is given by:
    # between cluster: varb = p*(1-p)*ICC
    # within cluster: varw = p*(1-p)*(1 - ICC)
    # 2*2*(p*(1-p)*(1 - ICC)*(1 - rho_s) + nsubjects*p*(1-p)*ICC*(1 - rho_c))/(nsubjects*nclusters)
    varb <- p * (1 - p) * ICC
    varw <- p * (1 - p) * (1 - ICC)
    ncp <-
      d / sqrt(2 * 2 * (varw * (1 - rho_s) + nsubjects * varb * (1 - rho_c)) /
                 (nsubjects * nclusters))
    
    pt(tcrit, 2 * (nclusters - 1), ncp, lower.tail = FALSE)
  })
  
  # calculate alpha
  if (is.na(alpha)) {
    alpha <- stats::uniroot(function(alpha)
      eval(pwr) - power,
      interval = c(1e-10, 1 - 1e-10),
      tol = tol)$root
  }
  
  # calculate power
  if (is.na(power)) {
    power <- eval(pwr)
  }
  
  # calculate nclusters
  if (is.na(nclusters)) {
    nclusters <- stats::uniroot(
      function(nclusters)
        eval(pwr) - power,
      interval = c(2 + 1e-10, 1e+07),
      tol = tol,
      extendInt = "upX"
    )$root
  }
  
  # calculate nsubjects
  if (is.na(nsubjects)) {
    nsubjects <- stats::uniroot(
      function(nsubjects)
        eval(pwr) - power,
      interval = c(2 + 1e-10, 1e+07),
      tol = tol,
      extendInt = "upX"
    )$root
  }
  
  # calculate p
  if (is.na(p)) {
    p <- stats::uniroot(function(p)
      eval(pwr) - power,
      interval = c(1e-7, 1 - 1e-7),
      tol = tol)$root
  }
  
  # calculate d
  if (is.na(d)) {
    d <- stats::uniroot(function(d)
      eval(pwr) - power,
      interval = c(1e-7, 1 - 1e-7),
      tol = tol)$root
  }
  
  # calculate ICC
  if (is.na(ICC)) {
    ICC <- stats::uniroot(function(ICC)
      eval(pwr) - power,
      interval = c(1e-07, 1 - 1e-7),
      tol = tol)$root
  }
  
  # calculate rho_c
  if (is.na(rho_c)) {
    rho_c <- stats::uniroot(function(rho_c)
      eval(pwr) - power,
      interval = c(1e-07, 1 - 1e-7),
      tol = tol)$root
  }
  
  # calculate rho_s
  if (is.na(rho_s)) {
    rho_s <- stats::uniroot(function(rho_s)
      eval(pwr) - power,
      interval = c(1e-07, 1 - 1e-7),
      tol = tol)$root
  }
  
  # # calculate covdf
  # if (is.na(covdf)){
  #   covdf <- stats::uniroot(function(covdf) eval(pwr) - power,
  #                           interval = c(1e-07, 1e+07),
  #                           tol = tol, extendInt = "upX")$root
  # }
  #
  # # calculate pvar_c
  # if (is.na(pvar_c)){
  #   pvar_c <- stats::uniroot(function(pvar_c) eval(pwr) - power,
  #                           interval = c(1e-07, 1 - 1e-7),
  #                           tol = tol)$root
  # }
  #
  # # calculate pvar_s
  # if (is.na(pvar_s)){
  #   pvar_s <- stats::uniroot(function(pvar_s) eval(pwr) - power,
  #                           interval = c(1e-07, 1 - 1e-7),
  #                           tol = tol)$root
  # }
  
  structure(get(target), names = target)
  
}
