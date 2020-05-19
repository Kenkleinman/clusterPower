#' Analytic power calculations for parallel arm cluster-randomized trials with count outcomes
#'
#' Compute the power, number of clusters needed, number of subjects per cluster 
#' needed, or other key parameters, for a simple parallel cluster randomized 
#' trial with a count outcome.
#'
#' Exactly one of \code{alpha}, \code{power}, \code{nclusters}, \code{nsubjects},
#'   \code{r1}, \code{r2}, and \code{CVB} must be passed as \code{NA}.
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
#'   invalid arguments are given.  This generally means that no solution exists for which the 
#'   omitted parameter and the supplied parameters fulfill the equation.  In particular, the desired
#'   power may not be acheiveable with any number of subjects or clusters.
#'
#' @param alpha The level of significance of the test, the probability of a
#'   Type I error.
#' @param power The power of the test, 1 minus the probability of a Type II
#'   error.
#' @param nclusters The number of clusters per condition. It must be greater than 1
#' @param nsubjects The number of units of person-time of observation per cluster
#' @param r1 The expected mean event rate per unit time in one of the conditions
#' @param r2 The expected mean event rate per unit time in the other condition
#' @param CVB The between-cluster coefficient of variation
#' @param r1inc Logical indicating if r1 is expected to be greater than r2. This is
#'   only important to specify if one of r1 or r2 is NA.
#' @param tol Numerical tolerance used in root finding. The default provides
#'   at least four significant digits.
#'   
#' @return #' @return The computed value of the NA parameter (from among \code{alpha}, \code{power}, \code{nclusters}, \code{nsubjects},
#'   \code{CV}, and \code{d})needed to satisfy the power and 
#' sample size equation.
#'
#' @examples 
#' # Find the number of clusters per condition needed for a trial with alpha = 0.05, 
#' # power = 0.80, 10 person-years per cluster, rate in condition 1 of 0.10 
#' # and condition 2 of 0.20, and CVB = 0.10.
#' 
#' cpa.count(nsubjects=10, r1=0.10, r2=0.20, CVB=0.10)
#' 
#' # 
#' # The result, showimg nclusters of greater than 24, suggests 25 clusters per
#' # condition should be used.
#' 
#' Find the largest CVB compatible with 80% power when there are 25 clusters, 10
#' subject-units of time per cluster, and a rate of 0.1 and 0.2 in each condition.  
#' 
#' cpa.count(nsubjects=10, nclusters= 25,r1=0.10, r2=0.20, CVB=NA)
#' 
#' Results show that CVB as high as 0.107 can still yield power this high.
#' 
#' @references Donner A, Klar N. Design and Analysis of Cluster Randomization Trials in Health Research. Chichester, UK; 2009.
#' 
#' @references Hayes JR, Moulton LH. Cluster Randomized Trials. Boca Raton, FL: CRC Press; 2009.
#' @export

cpa.count<- function(alpha = 0.05, power = 0.80,
                        nclusters = NA, nsubjects = NA,
                        r1 = NA, r2 = NA,
                        CVB = NA, r1inc = TRUE,
                        tol = .Machine$double.eps^0.25){
  
  if(!is.na(nclusters) && nclusters <= 1) {
    stop("'nclusters' must be greater than 1.")
  }
  
  needlist <- list(alpha, power, nclusters, nsubjects, r1, r2, CVB)
  neednames <- c("alpha", "power", "nclusters", "nsubjects", "r1", "r2", "CVB")
  needind <- which(unlist(lapply(needlist, is.na))) # find NA index
  
  if (length(needind) != 1) {
    stop("Exactly one of 'alpha', 'power', 'nclusters', 'nsubjects', 'r1', 'r2', or 'CVB' must be NA.")
  }
  
  target <- neednames[needind]
  
  pwr <- quote({
    IF <- 1 + CVB^2*(r1^2 + r2^2)*nsubjects/(r1 + r2)
    zcrit <- qnorm(alpha/2, lower.tail = FALSE)
    vard <- (r1 + r2)*IF/nsubjects
    pnorm(sqrt((nclusters - 1)*(r1 - r2)^2/vard) - zcrit, lower.tail = TRUE)
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
                                tol = tol)$root
  }
  
  # calculate nsubjects
  if (is.na(nsubjects)) {
    nsubjects <- stats::uniroot(function(nsubjects) eval(pwr) - power,
                         interval = c(1e-10, 1e+07),
                         tol = tol, extendInt = "upX")$root
  }
  
  # calculate r1
  if (is.na(r1)) {
    if(r1inc){
      r1 <- stats::uniroot(function(r1) eval(pwr) - power,
                           interval = c(r2 + 1e-7, 1 - 1e-7),
                           tol = tol, extendInt = "yes")$root
    } else {
      r1 <- stats::uniroot(function(r1) eval(pwr) - power,
                           interval = c(1e-7, r2 - 1e-7),
                           tol = tol, extendInt = "yes")$root
    }
  }
  
  # calculate r2
  if (is.na(r2)) {
    if(r1inc){
      r2 <- stats::uniroot(function(r2) eval(pwr) - power,
                           interval = c(1e-7, r1 - 1e-7),
                           tol = tol, extendInt = "yes")$root
      
    } else {
      r2 <- stats::uniroot(function(r2) eval(pwr) - power,
                           interval = c(r1 + 1e-7, 1 - 1e-7),
                           tol = tol, extendInt = "yes")$root
    }
  }
  
  # calculate CVB
  if (is.na(CVB)) {
    CVB <- stats::uniroot(function(CVB) eval(pwr) - power,
                          interval = c(1e-7, 1e+07),
                          tol = tol, extendInt = "downX")$root
  }
  
  structure(get(target), names = target)
  
}