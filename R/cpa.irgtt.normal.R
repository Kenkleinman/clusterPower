#' Power calculations for individually randomized group treatment trials, continuous outcome
#'
#' Compute the power of an individually randomized group treatment trial (IRGTT) design with a continuous outcome,
#' or determine parameters to obtain a target power.
#'
#' Exactly one of \code{alpha}, \code{power}, \code{nclusters}, \code{nsubjects},
#'   \code{ncontrols}, \code{d}, \code{varu}, \code{varei}, and \code{varr}
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
#' @param d The expected treatment effect.
#' @param varu The variance of the cluster level random effect for clusters in the intervention arm.
#' @param varei The variance of the subject level random error for individuals in the intervention arm.
#' @param varr The variance of the subject level random error for individuals in the control arm.
#' @param tol Numerical tolerance used in root finding. The default provides
#'   at least four significant digits.
#' @return The computed argument.
#' @examples 
#' # Find the required number of control subjects for an IRGTT with alpha = 0.05, power = 0.80,
#' # nclusters = 10, nsubjects = 10, d = 0.5 units, varu = 0.1, varei = 0.9, varr = 1.
#' cpa.irgtt.normal(nclusters=10, nsubjects = 10, 
#'   d = 0.5, varu = 0.1, varei = 0.9, varr = 1)
#' # 
#' # The result, ncontrols = 77.81084, suggests 78 subjects in the control arm should be recruited.
#' # This means that the total number of subjects in the study is nclusters*nsubjects + ncontrols = 10*10 + 78 = 178.
#' 
#' @references Moerbeek, M. and Wong, W. K. (2008) Sample size formulae for trials comparing
#' group and individual treatments in a multilevel model. Statist. Med., 27:2850-2864. 
#' doi: 10.1002/sim.3115.
#' 
#' @export


cpa.irgtt.normal <-
  function(alpha = 0.05,
           power = 0.80,
           nclusters = NA,
           nsubjects = NA,
           ncontrols = NA,
           d = NA,
           varu = NA,
           varei = NA,
           varr = NA,
           tol = .Machine$double.eps ^ 0.25) {
    # list of needed inputs
    needlist <-
      list(alpha,
           power,
           nclusters,
           nsubjects,
           ncontrols,
           d,
           varu,
           varei,
           varr)
    neednames <-
      c("alpha",
        "power",
        "nclusters",
        "nsubjects",
        "ncontrols",
        "d",
        "varu",
        "varei",
        "varr")
    needind <- which(unlist(lapply(needlist, is.na)))
    # check to see that exactly one needed param is NA
    
    if (length(needind) != 1) {
      neederror = "Exactly one of 'alpha', 'power', 'nclusters', 'nsubjects', 'ncontrols', 'd', 'varu', 'varei', and 'varr' must be NA."
      stop(neederror)
    }
    
    target <- neednames[needind]
    
    # evaluate power
    pwr <- quote({
      # variance of treatment effect d
      vard <-
        (nsubjects * varu + varei) / (nclusters * nsubjects) + varr / ncontrols
      
      zcrit <- qnorm(alpha / 2, lower.tail = FALSE)
      zstat <- abs(d) / sqrt(vard)
      pnorm(zcrit, zstat, lower.tail = FALSE)
      
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
    
    # calculate ncontrols
    if (is.na(ncontrols)) {
      ncontrols <- stats::uniroot(
        function(ncontrols)
          eval(pwr) - power,
        interval = c(2 + 1e-10, 1e+07),
        tol = tol,
        extendInt = "upX"
      )$root
    }
    
    # calculate d
    if (is.na(d)) {
      d <- stats::uniroot(
        function(d)
          eval(pwr) - power,
        interval = c(1e-07, 1e+07),
        tol = tol,
        extendInt = "upX"
      )$root
    }
    
    # calculate varu
    if (is.na(varu)) {
      varu <- stats::uniroot(
        function(varu)
          eval(pwr) - power,
        interval = c(1e-07, 1e+07),
        tol = tol,
        extendInt = "downX"
      )$root
    }
    
    # calculate varei
    if (is.na(varei)) {
      varei <- stats::uniroot(
        function(varei)
          eval(pwr) - power,
        interval = c(1e-07, 1e+07),
        tol = tol,
        extendInt = "downX"
      )$root
    }
    
    # calculate varr
    if (is.na(varr)) {
      varr <- stats::uniroot(
        function(varr)
          eval(pwr) - power,
        interval = c(1e-07, 1e+07),
        tol = tol,
        extendInt = "downX"
      )$root
    }
    
    structure(get(target), names = target)
  }
