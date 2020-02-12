#' Power calculations for multi-arm cluster randomized trials, continuous outcome
#'
#' Compute the power of the overall F-test for a multi-arm cluster randomized trial with a continuous
#' outcome, or determine parameters to obtain a target power.
#'
#' Exactly one of \code{alpha}, \code{power}, \code{narms}, \code{nclusters},
#'   \code{nsubjects}, \code{vara}, \code{varc}, and \code{vare}  must be passed as \code{NA}.
#'   Note that \code{alpha} and \code{power} have non-\code{NA}
#'   defaults, so if those are the parameters of interest they must be
#'   explicitly passed as \code{NA}.
#'   
#' Assuming a balanced design, the between-arm variance \eqn{\sigma_a^2} (corresponding to
#'   the function argument \code{vara}) can be estimated using the formula: 
#'   
#'   \deqn{\sigma_a^2 = \sum\limits_{i=1}^{n_a}(\mu_i - \mu)^2/(n_a-1)}
#'   
#'   where \eqn{n_a} is the number of arms, \eqn{\mu_i} is the estimate of the \eqn{i}-th arm
#'   mean, and \eqn{\mu} is the estimate of the overall mean of the outcome. This 
#'   variance can be computed in R using the \code{var} function and a vector of arm means.
#'   For example, suppose the estimated means for a three-arm trial were 74, 80, and 86 Then the
#'   estimate of the between-arm variance could be computed with \code{var(c(74,80,86))}, 
#'   yielding a value of 36.
#'
#' @section Note:
#'   This function was inspired by work from Stephane Champely (pwr.t.test),
#'   Peter Dalgaard (power.t.test), and Claus Ekstrom (power.anova.test). As 
#'   with those functions, 'uniroot' is used to solve power equation for
#'   unknowns, so you may see errors from it, notably about inability to
#'   bracket the root when invalid arguments are given.
#'
#' @section Authors:
#' Jonathan Moyer (\email{jon.moyer@@gmail.com}), Ken Kleinman (\email{ken.kleinman@@gmail.com})
#'
#' @param alpha The level of significance of the test, the probability of a
#'   Type I error.
#' @param power The power of the test, 1 minus the probability of a Type II
#'   error.
#' @param narms The number of independent arms (conditions). It must be greater than 2.
#' @param nclusters The number of clusters per arm. It must be greater than 1.
#' @param nsubjects The cluster size.
#' @param vara The between-arm variance.
#' @param varc The between-cluster variance.
#' @param vare The within-cluster variance.
#' @param tol Numerical tolerance used in root finding. The default provides
#'   at least four significant digits.
#' @return The computed argument.
#' @examples 
#' # Suppose we are planning a multi-arm trial composed of a control arm and 
#' # two treatment arms. It is known that each arm will contain 5 clusters. We
#' # wish to know the minimum number of subjects per cluster necessary to
#' # attain 80% power at a 5% level of significance. A pilot study was used to
#' # determine estimates of the between-arm variance, the between-cluster 
#' # variance, and the within-cluster variance. The observed means of each arm
#' # in the pilot study were 74, 80, and 86, so the  between-arm variance is 36.
#' # As discussed in the "Details" section above, this can be calculated using
#' # the command var(c(74,80,86)). The within-cluster and between-cluster
#' # standard deviations were observed to be 8 and 3, respectively. This means
#' # the within-cluster and between-cluster variances are 64 and 9, respectively.
#' # These values are entered into the function as follows:
#' 
#' crtpwr.nmean(narms=3,nclusters=5,vara=36,varc=9,vare=64)
#' # 
#' # The result, showing nsubjects of greater than 20, suggests 21 subjects per 
#' # cluster should be used.
#' @references Murray DM. Design and Analysis of Group-Randomized Trials. New York,
#'   NY: Oxford University Press; 1998.
#' @export

crtpwr.nmean <- function(alpha = 0.05, power = 0.80,
                         narms = NA, nclusters = NA, nsubjects = NA,
                         vara = NA, varc = NA, vare = NA, 
                         tol = .Machine$double.eps^0.25){
  
  .Deprecated("cpa.ma.normal")
  
  # if(!is.na(narms) && narms <= 2) {
  #   stop("'narms' must be greater than 2.")
  # }
  
  if(!is.na(nclusters) && nclusters <= 1) {
    stop("'nclusters' must be greater than 1.")
  }
  
  if(!is.na(nsubjects) && nsubjects <= 1) {
    stop("'nsubjects' must be greater than 1.")
  }
  
  # list of needed inputs
  needlist <- list(alpha, power, narms, nclusters, nsubjects, vara, varc, vare)
  neednames <- c("alpha", "power", "narms", "nclusters", "nsubjects", "vara", "varc", "vare")
  needind <- which(unlist(lapply(needlist, is.na)))
  # check to see that exactly one needed param is NA
  
  if (length(needind) != 1) {
    neederror = "Exactly one of 'alpha', 'power', 'narms', 'nclusters', 'nsubjects', 'vara', 'varb', and 'varc' must be NA."
    stop(neederror)
  } 
  
  target <- neednames[needind]

  # evaluate power
  pwr <- quote({
    ndf <- narms - 1
    ddf <- narms*(nclusters-1)
    fcrit <- qf(alpha/2, ndf, ddf, lower.tail = FALSE)
    ncp <- (vare+nsubjects*varc+nsubjects*nclusters*vara)/(vare+nsubjects*varc)
    pf(fcrit, ndf, ddf, ncp, lower.tail = FALSE) 
  })
  
  # calculate alpha
  if (is.na(alpha)) {
    alpha <- stats::uniroot(function(alpha) eval(pwr) - power,
                            interval = c(1e-10, 1 - 1e-10),
                            tol = tol)$root
  }
  
  # calculate power
  if (is.na(power)) {
    power <- eval(pwr)
  }
  
  # calculate narms
  if (is.na(narms)) {
    narms <- stats::uniroot(function(narms) eval(pwr) - power,
                                interval = c(2 + 1e-10, 1e+07),
                                tol = tol, extendInt = "upX")$root
  }
  
  # calculate nclusters
  if (is.na(nclusters)) {
    nclusters <- stats::uniroot(function(nclusters) eval(pwr) - power,
                                interval = c(2 + 1e-10, 1e+07),
                                tol = tol, extendInt = "upX")$root
  }
  
  # calculate nsubjects
  if (is.na(nsubjects)) {
    nsubjects <- stats::uniroot(function(nsubjects) eval(pwr) - power,
                                interval = c(2 + 1e-10, 1e+07),
                                tol = tol, extendInt = "upX")$root
  }
  
  # calculate vara
  if (is.na(vara)){
    vara <- stats::uniroot(function(vara) eval(pwr) - power,
                           interval = c(1e-07, 1e+07),
                           tol = tol, extendInt = "downX")$root
  }
  
  # calculate varc
  if (is.na(varc)){
    varc <- stats::uniroot(function(varc) eval(pwr) - power,
                          interval = c(1e-07, 1e+07),
                          tol = tol, extendInt = "downX")$root
  }
  
  # calculate vare
  if (is.na(vare)) {
    vare <- stats::uniroot(function(vare) eval(pwr) - power,
                           interval = c(1e-07, 1e+07),
                           tol = tol, extendInt = "downX")$root
  }
  
  structure(get(target), names = target)

}
