#' Power calculations for individually randomized group treatment trials, binary outcome
#'
#' Compute the power of an individually randomized group treatment trial (IRGTT) design with a binary outcome,
#' or determine parameters to obtain a target power.
#'
#' Exactly one of \code{alpha}, \code{power}, \code{nclusters}, \code{nsubjects},
#'   \code{ncontrols}, \code{ICC}, \code{p2}, and \code{p1}
#'   must be passed as \code{NA}. Note that \code{alpha} and \code{power}
#'   have non-\code{NA} defaults, so if those are the parameters of 
#'   interest they must be explicitly passed as \code{NA}.
#'
#' @section Note:
#'   This function was inspired by work from Stephane Champely (pwr.t.test) and
#'   Peter Dalgaard (power.t.test). As with those functions, 'uniroot' is used to
#'   solve power equation for unknowns, so you may see
#'   errors from it, notably about inability to bracket the root when
#'   invalid arguments are given.
#'   
#' @section Authors:
#' Jonathan Moyer (\email{jon.moyer@@gmail.com}), Ken Kleinman (\email{ken.kleinman@@gmail.com})
#'
#' @param alpha The level of significance of the test, the probability of a
#'   Type I error.
#' @param power The power of the test, 1 minus the probability of a Type II
#'   error.
#' @param nclusters The number of clusters in the intervention arm.
#' @param nsubjects The number of subjects in each cluster in the intervention arm.
#' @param ncontrols The number of subjects in the control arm.
#' @param ICC The intracluster correlation coefficient, the correlation in outcome measurements between
#'   two individuals from the same cluster in the intervention arm.
#' @param p2 The expected probability of the outcome in the intervention arm.
#' @param p1 The expected probability of the outcome in the control arm.
#' @param decrease Whether or not the proportion in the intervention arm is expected to be
#'   less than the proportion in the control arm. If TRUE it is assumed p2 < p1, while FALSE implies
#'   p2 > p1.
#' @param tol Numerical tolerance used in root finding. The default provides
#'   at least four significant digits.
#' @return The computed argument.
#' @examples 
#' # Find the required number of subjects per intervention cluster an IRGTT with alpha = 0.05,
#' # power = 0.80, nclusters = 23, ncontrols = 146, ICC = 0.05, p2 = 0.397, and p1 = 0.243.
#' 
#' cpa.irgtt.binary(nclusters=23, ncontrols = 146, 
#'   ICC = 0.05, p2 = 0.397, p1 = 0.243, decrease = FALSE)
#' 
#' # 
#' # The result, nsubjects = 7.96624, suggests 8 subjects per cluster 
#' # in the intervention arm should be recruited.
#' # This means that the total number of subjects in the 
#' # study is nclusters * nsubjects + ncontrols = 23 * 8 + 146 = 330.
#' 
#' @references Moerbeek, M. and Wong, W. K. (2008) Sample size formulae for trials comparing
#' group and individual treatments in a multilevel model. Statist. Med., 27:2850-2864. 
#' doi: 10.1002/sim.3115.
#' 
#' @export

cpa.irgtt.binary <- function(alpha = 0.05, power = 0.80, nclusters = NA,
                             nsubjects = NA, ncontrols = NA, 
                             ICC = NA, p2 = NA, p1 = NA,
                             decrease = TRUE,
                             tol = .Machine$double.eps^0.25){
  
  # list of needed inputs
  needlist <- list(alpha, power, nclusters, nsubjects, ncontrols, ICC, p2, p1)
  neednames <- c("alpha", "power", "nclusters", "nsubjects", "ncontrols", "ICC", "p2", "p1")
  needind <- which(unlist(lapply(needlist, is.na)))
  # check to see that exactly one needed param is NA
  
  if (length(needind) != 1) {
    neederror = "Exactly one of 'alpha', 'power', 'nclusters', 'nsubjects', 'ncontrols', 'ICC', 'p2', and 'p1' must be NA."
    stop(neederror)
  } 
  
  target <- neednames[needind]
  
  # evaluate power
  pwr <- quote({
    
    # design effect in intervention arm
    DE <- (nsubjects - 1)*ICC + 1
    
    # variance of treatment effect d
    vard <- p2*(1-p2)*DE/(nclusters*nsubjects) + p1*(1-p1)/ncontrols
    
    zcrit <- qnorm(alpha/2, lower.tail = FALSE)
    zstat <- abs(p2-p1)/sqrt(vard)
    pnorm(zcrit, zstat, lower.tail = FALSE)
    
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
  
  # calculate ncontrols
  if (is.na(ncontrols)) {
    ncontrols <- stats::uniroot(function(ncontrols) eval(pwr) - power,
                                interval = c(2 + 1e-10, 1e+07),
                                tol = tol, extendInt = "upX")$root
  }
  
  # calculate ICC
  if (is.na(ICC)) {
    ICC <- stats::uniroot(function(ICC) eval(pwr) - power,
                          interval = c(1e-10, 1 - 1e-10),
                          tol = tol)$root
  }
  
  # calculate p2
  if (is.na(p2)) {
    if(decrease){
      p2 <- stats::uniroot(function(p2) eval(pwr) - power,
                            interval = c(1e-10, p1),
                            tol = tol)$root
      
    } else {
      p2 <- stats::uniroot(function(p2) eval(pwr) - power,
                            interval = c(p1, 1 - 1e-10),
                            tol = tol, extendInt = "yes")$root
    }
  }
  
  # calculate p1
  if (is.na(p1)) {
    if(decrease){
      p1 <- stats::uniroot(function(p1) eval(pwr) - power,
                            interval = c(p2, 1 - 1e-10),
                            tol = tol)$root
    } else {
      p1 <- stats::uniroot(function(p1) eval(pwr) - power,
                            interval = c(1e-10, p2),
                            tol = tol)$root
    }
    
  }
  
  structure(get(target), names = target)
}

