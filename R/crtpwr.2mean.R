#' Power calculations for simple cluster randomized trials, continuous outcome
#'
#' Compute the power of a simple cluster randomized trial with a continuous outcome,
#' or determine parameters to obtain a target power.
#'
#' Exactly one of \code{alpha}, \code{power}, \code{m}, \code{n},
#'   \code{cv}, \code{d}, \code{icc}, and \code{varw}  must be passed as \code{NA}.
#'   Note that \code{alpha}, \code{power}, and \code{cv} have non-\code{NA}
#'   defaults, so if those are the parameters of interest they must be
#'   explicitly passed as \code{NA}.
#'   
#' If \code{n} is a vector the values, \code{m} and \code{cv} will be recalculated
#'    using the values in \code{n}. If \code{n} is a vector and \code{method} is
#'    "taylor", the exact relative efficiency will be calculated as described in
#'    van Breukelen et al (2007).
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
#' @param m The number of clusters per condition. It must be greater than 1.
#' @param n The mean of the cluster sizes, or a vector of cluster sizes for one arm.
#' @param cv The coefficient of variation of the cluster sizes. When \code{cv} = 0,
#'   the clusters all have the same size.
#' @param d The difference in condition means.
#' @param icc The intraclass correlation.
#' @param vart The total variation of the outcome (the sum of within- and between-cluster variation).
#' @param method The method for calculating variance inflation due to unequal cluster
#'   sizes. Either a method based on Taylor approximation of relative efficiency 
#'   ("taylor"), or weighting by cluster size ("weighted")
#' @param tol Numerical tolerance used in root finding. The default provides
#'   at least four significant digits.
#' @return The computed argument.
#' @examples 
#' # Find the number of clusters per condition needed for a trial with alpha = .05, 
#' # power = 0.8, 10 observations per cluster, no variation in cluster size, a difference 
#' # of 1 unit,  icc = 0.1 and   a variance of five units.
#' crtpwr.2mean(n=10 ,d=1, icc=.1, varw=5)
#' # 
#' # The result, showimg m of greater than 15, suggests 16 clusters per condition should be used.
#' @references Eldridge SM, Ukoumunne OC, Carlin JB. (2009) The Intra-Cluster Correlation
#'   Coefficient in Cluster Randomized Trials: A Review of Definitions. Int Stat Rev. 
#'   77: 378-394.
#' @references Eldridge SM, Ashby D, Kerry S. (2006) Sample size for cluster randomized
#'   trials: effect of coefficient of variation of cluster size and analysis method.
#'   Int J Epidemiol. 35(5):1292-300.
#' @references van Breukelen GJP, Candel MJJM, Berger MPF. (2007) Relative efficiency of
#'   unequal versus equal cluster sizes in cluster randomized and multicentre trials.
#'   Statist Med. 26:2589-2603.  
#' @export

crtpwr.2mean <- function(alpha = 0.05, power = 0.80, m = NA,
                         n = NA, cv = 0,
                         d = NA, icc = NA,
                         vart = NA,
                         method = c("taylor", "weighted"),
                         tol = .Machine$double.eps^0.25){
  
  method <- match.arg(method)
  
  # if n is a vector, 
  if(length(n) > 1){
    nvec <- n
    n <- mean(nvec)
    nsd <- stats::sd(nvec)
    cv <- nsd/n
    m <- length(nvec)
  }
  
  if(!is.na(m) && m <= 1) {
    stop("'m' must be greater than 1.")
  }

  
  # list of needed inputs
  needlist <- list(alpha, power, m, n, cv, d, icc, vart)
  neednames <- c("alpha", "power", "m", "n", "cv", "d", "icc", "vart")
  needind <- which(unlist(lapply(needlist, is.na)))
  # check to see that exactly one needed param is NA
  
  if (length(needind) != 1) {
    neederror = "Exactly one of 'alpha', 'power', 'm', 'n', 'cv', 'd', 'icc' and 'vart' must be NA."
    stop(neederror)
  } 
  
  target <- neednames[needind]
  
  # evaluate power
  pwr <- quote({
  
    # variance inflation
    # if nvec exists, calcuate exact relative efficiency
    if (exists("nvec")) {
      if(method == "taylor"){
        a <- (1 - icc)/icc
        DEFF <- 1 + (n - 1)*icc
        RE <- ((n + a)/n)*(sum((nvec/(nvec+a)))/m) # exact relative efficiency
        VIF <- DEFF*RE
      } else{
        VIF <- 1 + ((cv^2 + 1)*n - 1)*icc
      }
    } else if(!is.na(n)){
      if(method == "taylor"){
        DEFF <- 1 + (n - 1)*icc
        L <- n*icc/DEFF
        REt <- 1/(1 - cv^2*L*(1 - L)) # taylor approximation
        VIF <- DEFF*REt
      } else {
        VIF <- 1 + ((cv^2 + 1)*n - 1)*icc
      }
    }
    
    tcrit <- qt(alpha/2, 2*(m - 1), lower.tail = FALSE)
    
    ncp <- sqrt(m*n/(2*VIF)) * abs(d)/sqrt(vart)
    
    pt(tcrit, 2*(m - 1), ncp, lower.tail = FALSE) #+ pt(-tcrit, 2*(m - 1), ncp, lower.tail = TRUE)
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
  
  # calculate cv
  if (is.na(cv)) {
    cv <- stats::uniroot(function(cv) eval(pwr) - power,
                  interval = c(1e-10, 1e+07),
                  tol = tol, extendInt = "downX")$root
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
                   interval = c(1e-07, 1 - 1e-07),
                   tol = tol)$root
  }
  
  # calculate vart
  if (is.na(vart)) {
    varw <- stats::uniroot(function(vart) eval(pwr) - power,
                    interval = c(1e-07, 1e+07),
                    tol = tol, extendInt = "downX")$root
  }
  
  structure(get(target), names = target)
  
  # method <- paste("Clustered two-sample t-test power calculation: ", target, sep = "")
  # note <- "'m' is the number of clusters in each group and 'n' is the number of individuals in each cluster."
  # structure(list(alpha = alpha, power = power, m = m, n = n, cv = cv, d = d,
  #                icc = icc, varw = varw, note = note, method = method),
  #           class = "power.htest")
  
}
