#' Power calculations for simple cluster randomized trials, continuous outcome
#'
#' Compute the power, number of clusters needed, number of subjects per cluster 
#' needed, or other key parameters, for a simple parallel cluster randomized 
#' trial with a continuous outcome.
#'
#' Exactly one of \code{alpha}, \code{power}, \code{nclusters}, \code{nsubjects},
#'   \code{CV}, and \code{d}  must be passed as \code{NA}.
#'   Note that \code{alpha}, \code{power}, and \code{CV} have non-\code{NA}
#'   defaults, so if those are the parameters of interest they must be
#'   explicitly passed as \code{NA}. The user must supply sufficient variance 
#'   parameters to produce values for both ICC and vart by providing 2 of the 
#'   following: \code{ICC}, \code{vart}, \code{sigma_b_sq}, or \code{sigma_sq}.
#'
#' If \code{nsubjects} is a vector the values, \code{nclusters} and \code{CV} will be recalculated
#'    using the values in \code{nsubjects}. If \code{nsubjects} is a vector and \code{method} is
#'    "taylor", the exact relative efficiency will be calculated as described in
#'    van Breukelen et al (2007).
#'    
#' @section Note:
#'   This function was inspired by work from Stephane Champely (pwr.t.test) and
#'   Peter Dalgaard (power.t.test). As with those functions, 'uniroot' is used to
#'   solve power equation for unknowns, so you may see
#'   errors from it, notably about inability to bracket the root when
#'   invalid arguments are given. This generally means that no solution exists for which the 
#'   omitted parameter and the supplied parameters fulfill the equation.  In particular, the desired
#'   power may not be acheiveable with any number of subjects or clusters.
#'
#' @section Testing details:
#' This function has been verified against reference values from the NIH's GRT
#' Sample Size Calculator, PASS11, \code{CRTsize::n4means}, and
#' \code{clusterPower::cps.normal}.
#'
#'
#' @section Authors:
#' Jonathan Moyer (\email{jon.moyer@@gmail.com}), 
#' Alexandria C. Sakrejda (\email{acbro0@@umass.edu}),
#' and Ken Kleinman (\email{ken.kleinman@@gmail.com})
#'
#'
#' @param alpha The level of significance of the test, the probability of a
#'   Type I error.
#' @param power The power of the test, 1 minus the probability of a Type II
#'   error.
#' @param nclusters The number of clusters per condition. It must be greater than 1.
#' @param nsubjects The mean of the cluster sizes, or a vector of cluster sizes for one arm.
#' @param sigma_sq Within-cluster variance.
#' @param sigma_b_sq Between-cluster variance.
#' @param CV The coefficient of variation of the cluster sizes. When \code{CV} = 0,
#'   the clusters all have the same size.
#' @param d The difference in condition means \eqn{|\beta_1|}
#' @param ICC The intraclass correlation \eqn{\sigma_b^2 / (\sigma_b^2 + \sigma^2)}.
#' Accepts a numeric between 0-1.
#' @param vart The total variation of the outcome (the sum of within- and 
#' between-cluster variation) \eqn{\sigma_b^2 + \sigma^2}.
#' @param method The method for calculating variance inflation due to unequal cluster
#'   sizes. Either a method based on Taylor approximation of relative efficiency
#'   ("taylor"), or weighting by cluster size ("weighted"). Default is \code{"taylor"}.
#' @param tol Numerical tolerance used in root finding. The default provides
#'   at least four significant digits.
#' 
#' @return The computed value of the NA parameter (from among \code{alpha}, \code{power}, \code{nclusters}, \code{nsubjects},
#'   \code{CV}, and \code{d})needed to satisfy the power and 
#' sample size equation.
#'
#' @examples 
#' # Find the number of clusters per condition needed for a trial with alpha = .05, 
#' # power = 0.8, 10 observations per cluster, no variation in cluster size, a difference 
#' # of 1 unit,  ICC = 0.1 and a variance of five units.
#' 
#' cpa.normal(nsubjects=10, d=1, ICC=.1, vart=5)
#'  
#' # The result, showing nclusters of greater than 15, suggests 16 clusters per 
#' # condition should be used.
#' 
#' # Find the power achieved with 16 clusters, 10 subjects per cluster,
#' # difference between condition of 1 unit, ICC = .1, and total variance
#' # of 5 units
#' 
#' cpa.normal(power = NA, nclusters= 16, nsubjects=10 ,d=1, sigma_b_sq=.5, vart=5)
#' 
#' # The result shows the power is 0.801766.
#' 
#' @references Eldridge SM, Ukoumunne OC, Carlin JB. (2009) The Intra-Cluster Correlation
#'   Coefficient in Cluster Randomized Trials: A Review of Definitions. Int Stat Rev.
#'   77: 378-394.
#' @references Eldridge SM, Ashby D, Kerry S. (2006) Sample size for cluster randomized
#'   trials: effect of coefficient of variation of cluster size and analysis method.
#'   Int J Epidemiol. 35(5):1292-300.
#' @references van Breukelen GJP, Candel MJJM, Berger MPF. (2007) Relative efficiency of
#'   unequal versus equal cluster sizes in cluster randomized and multicentre trials.
#'   Statist Med. 26:2589-2603.
#'   
#' @export


cpa.normal <- function(alpha = 0.05,
                       power = 0.80,
                       nclusters = NA,
                       nsubjects = NA,
                       sigma_sq = NA,
                       sigma_b_sq = NA,
                       CV = 0,
                       d = NA,
                       ICC = NA,
                       vart = NA,
                       method = c("taylor", "weighted"),
                       tol = .Machine$double.eps ^ 0.25) {
  method <- match.arg(method)
  
  # if nsubjects is a vector,
  if (length(nsubjects) > 1) {
    nvec <- nsubjects
    nsubjects <- mean(nvec)
    nsd <- stats::sd(nvec)
    CV <- nsd / nsubjects
    nclusters <- length(nvec)
  }
  
  if (!is.na(nclusters) && nclusters <= 1) {
    stop("'nclusters' must be greater than 1.")
  }
  
  # check for sufficient variance parameters
  varname <-
    c("vart",
      "ICC",
      "sigma_sq",
      "sigma_b_sq"
    )
  varind <- which(is.na(varname))
  if (length(varind) != 2) {
    varerror = "At least 2 of ICC, vart, sigma_b,and sigma_b_sq must be supplied."
    stop(varerror)
  }
  
  if (!is.na(sigma_b_sq) && is.na(vart)) {
    vart <- (sigma_b_sq / ICC)
  }
  if (!is.na(sigma_sq) && is.na(vart)) {
    vart <- (sigma_sq / (1 - ICC))
  }
  if (!is.na(sigma_b_sq) && is.na(ICC)) {
    ICC <- (sigma_b_sq / vart)
  }
  if (!is.na(sigma_sq) && is.na(ICC)) {
    ICC <- (1 - (sigma_sq / vart))
  }
  if (is.na(vart) && is.na(ICC)) {
    ICC <- ( (sigma_b_sq / (sigma_b_sq + sigma_sq) ) )
    vart <- (sigma_b_sq + sigma_sq)
  }
  
  # list of needed inputs
  needlist <-
    list(alpha, power, nclusters, nsubjects, CV, d)
  neednames <-
    c("alpha",
      "power",
      "nclusters",
      "nsubjects",
      "CV",
      "d"
      )
  needind <- which(unlist(lapply(needlist, is.na)))
  # check to see that exactly one needed param is NA
  
  if (length(needind) != 1) {
    neederror = "Exactly one of 'alpha', 'power', 'nclusters', 'nsubjects', 'CV', and 'd' must be NA."
    stop(neederror)
  }
  
  target <- neednames[needind]
  
  # evaluate power
  pwr <- quote({
    # variance inflation
    # if nvec exists, calcuate exact relative efficiency
    if (exists("nvec")) {
      if (method == "taylor") {
        a <- (1 - ICC) / ICC
        DEFF <- 1 + (nsubjects - 1) * ICC
        RE <-
          ((nsubjects + a) / nsubjects) * (sum((nvec / (nvec + a))) / nclusters) # exact relative efficiency
        VIF <- DEFF * RE
      } else{
        VIF <- 1 + ((CV ^ 2 + 1) * nsubjects - 1) * ICC
      }
    } else if (!is.na(nsubjects)) {
      if (method == "taylor") {
        DEFF <- 1 + (nsubjects - 1) * ICC
        L <- nsubjects * ICC / DEFF
        REt <- 1 / (1 - CV ^ 2 * L * (1 - L)) # taylor approximation
        VIF <- DEFF * REt
      } else {
        VIF <- 1 + ((CV ^ 2 + 1) * nsubjects - 1) * ICC
      }
    }
    
    tcrit <- qt(alpha / 2, 2 * (nclusters - 1), lower.tail = FALSE)
    
    ncp <-
      sqrt(nclusters * nsubjects / (2 * VIF)) * abs(d) / sqrt(vart)
    
    pt(tcrit, 2 * (nclusters - 1), ncp, lower.tail = FALSE) #+ pt(-tcrit, 2*(nclusters - 1), ncp, lower.tail = TRUE)
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
  
  # calculate CV
  if (is.na(CV)) {
    CV <- stats::uniroot(
      function(CV)
        eval(pwr) - power,
      interval = c(1e-10, 1e+07),
      tol = tol,
      extendInt = "downX"
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
  
  # calculate ICC
  if (is.na(ICC)) {
    ICC <- stats::uniroot(function(ICC)
      eval(pwr) - power,
      interval = c(1e-07, 1 - 1e-07),
      tol = tol)$root
  }
  
  # calculate vart
  if (is.na(vart)) {
    vart <- stats::uniroot(
      function(vart)
        eval(pwr) - power,
      interval = c(1e-07, 1e+07),
      tol = tol,
      extendInt = "downX"
    )$root
  }
  
  structure(get(target), names = target)
}
