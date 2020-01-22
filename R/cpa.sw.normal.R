#' Power calculations for stepped wedge cluster randomized trials, continuous outcome
#'
#' Compute the power of a stepped wedge cluster randomized trial design with a continuous outcome,
#' or determine parameters to obtain a target power.
#'
#' Exactly one of \code{alpha}, \code{power}, \code{nclusters}, \code{nsubjects},
#'   \code{ntimes}, \code{d}, \code{icc}, \code{rho_c}, \code{rho_s}, and \code{vart}
#'   must be passed as \code{NA}. Note that \code{alpha} and\code{power}
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
#' @param nclusters The number of clusters switching to the intervention condition at each time point.
#' @param nsubjects The number of subjects in each cluster.
#' @param ntimes The number of time points in the trial (not including baseline).
#' @param d The expected treatment effect.
#' @param icc The intra-cluster correlation, the correlation in outcome measurements between
#'   two individuals from the same cluster.
#' @param rho_c The cluster autocorrelation, the correlation between two population means
#'   from the same cluster at different times. This value can be used in both cross-sectional and cohort
#'   designs.
#' @param rho_s The individual autocorrelation, the correlation between two outcome measurements
#'   in the same individual at different times. In a purely cross-sectional design (new subjects are
#'   obtained at each time point), this value should be 0.
#' @param vart The total variation of the outcome (the sum of variances from cluster and individual
#'   level random effects).
#' @param tol Numerical tolerance used in root finding. The default provides
#'   at least four significant digits.
#' @return The computed argument.
#' @examples 
#' # Find the required number of clusters switching to intervention at each time point for a trial 
#' # with alpha = 0.05, power = 0.80, nsubjects = 50, ntimes = 5, d = 1.5 units, icc = 0.2,
#' # rho_c = 0.80, rho_s = 0.50, and vart = 16 square-units.
#' cpa.sw.normal(nsubjects = 50, ntimes = 5, d = 1.5, icc = 0.2, rho_c = 0.80, rho_s = 0.50, vart = 16)
#' # 
#' # The result, nclusters = 2.134918, suggests 3 clusters switching per time point should be used. This
#' # means that the total number of clusters in the study is nclusters*ntimes = 3*5 = 15.
#' 
#' @references Hooper, R., Teerenstra, S., Hoop, E., and Eldridge, S. (2016)
#'   Sample size calculation for stepped wedge and other longitudinal cluster
#'   randomised trials. Statist. Med., 35: 4718-4728. doi: 10.1002/sim.7028.
#'   
#' @references Hussey, M. and Hughes, J. (2007) Design and analysis of stepped wedge cluster
#'   randomized trials. Contemp. Clin. Trials, 28: 182-191. doi: 10.1016/j.cct.2006.05.007.
#'   
#' @export

cpa.sw.normal <- function(alpha = 0.05, power = 0.80, nclusters = NA,
                          nsubjects = NA, ntimes = NA, d = NA, icc = NA,
                          rho_c = NA, rho_s = NA,
                          vart = NA,
                          tol = .Machine$double.eps^0.25){
  
  # list of needed inputs
  needlist <- list(alpha, power, nclusters, nsubjects, ntimes, d, icc, rho_c, rho_s, vart)
  neednames <- c("alpha", "power", "nclusters", "nsubjects", "ntimes", "d", "icc",
                 "rho_c", "rho_s", "vart")
  needind <- which(unlist(lapply(needlist, is.na)))
  # check to see that exactly one needed param is NA
  
  if (length(needind) != 1) {
    neederror = "Exactly one of 'alpha', 'power', 'nclusters', 'nsubjects', 'ntimes', 'd', 'icc', 'rho_c', 'rho_s', and 'vart' must be NA."
    stop(neederror)
  } 
  
  target <- neednames[needind]
  
  # evaluate power
  pwr <- quote({
    
    # variance inflation
    # DEFFc is design effect due to clustering
    DEFFc <- 1 + (nsubjects - 1)*icc 
    
    r <- (nsubjects*icc*rho_c + (1 - icc)*rho_s)/DEFFc
    
    # DEFFr is design effect due to repeated assessment
    DEFFr <- 3*ntimes*(1 - r)*(1 + ntimes*r)/( (ntimes^2 - 1)*(2 + ntimes*r) )
    
    VIF <-DEFFr*DEFFc
    
    # zcrit <- qnorm(alpha/2, lower.tail = FALSE)
    # zstat <- abs(d)/sqrt( 4*vart/(nclusters*nsubjects*ntimes)*VIF )
    # pnorm(zcrit, zstat, lower.tail = FALSE)
    
    tcrit <- qt(alpha/2, ntimes*nclusters - (ntimes + 2), lower.tail = FALSE)
    ncp <- abs(d)/sqrt( 4*vart/(nclusters*nsubjects*ntimes)*VIF )
    pt(tcrit, ntimes*nclusters - (ntimes + 2), ncp, lower.tail = FALSE)#+
    #pt(-tcrit, 2*(nclusters - 1), ncp, lower.tail = TRUE)
    
    
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
  
  # calculate ntimes
  if (is.na(ntimes)) {
    ntimes <- stats::uniroot(function(ntimes) eval(pwr) - power,
                             interval = c(1 + 1e-10, 1e+07),
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
                          interval = c(1e-07, 1 - 1e-07),
                          tol = tol)$root
  }
  
  # calculate rho_c
  if (is.na(rho_c)){
    rho_c <- stats::uniroot(function(rho_c) eval(pwr) - power,
                            interval = c(1e-07, 1 - 1e-07),
                            tol = tol)$root
  }
  
  # calculate rho_s
  if (is.na(rho_s)){
    rho_s <- stats::uniroot(function(rho_s) eval(pwr) - power,
                            interval = c(1e-07, 1 - 1e-07),
                            tol = tol)$root
  }
  
  # calculate vart
  if (is.na(vart)) {
    vart <- stats::uniroot(function(vart) eval(pwr) - power,
                           interval = c(1e-07, 1e+07),
                           tol = tol, extendInt = "downX")$root
  }
  
  structure(get(target), names = target)
}

