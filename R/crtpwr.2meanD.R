#' Power calculations for difference-in-difference cluster randomized trials, continuous outcome
#'
#' Compute the power of a difference-in-difference cluster randomized trial design with a continuous outcome,
#' or determine parameters to obtain a target power.
#'
#' Exactly one of \code{alpha}, \code{power}, \code{nclusters}, \code{nsubjects},
#'   \code{d}, \code{icc}, \code{rho_c}, \code{rho_s}, and \code{vart}
#'   must be passed as \code{NA}. Note that \code{alpha} and\code{power}
#'   have non-\code{NA} defaults, so if those are the parameters of 
#'   interest they must be explicitly passed as \code{NA}.
#'   
#' If \code{nsubjects} is a vector the values, \code{nclusters} will be recalculated
#'    using the values in \code{nsubjects}. 
#'
#' @section Note:
#'   This function was inspired by work from Stephane Champely (pwr.t.test) and
#'   Peter Dalgaard (power.t.test). As with those functions, 'uniroot' is used to
#'   solve power equation for unknowns, so you may see
#'   errors from it, notably about inability to bracket the root when
#'   invalid arguments are given.
#'
#' @section Authors:
#' @section Authors:
#' Jonathan Moyer (\email{jon.moyer@@gmail.com}), Ken Kleinman (\email{ken.kleinman@@gmail.com})
#'
#' @param alpha The level of significance of the test, the probability of a
#'   Type I error.
#' @param power The power of the test, 1 minus the probability of a Type II
#'   error.
#' @param nclusters The number of clusters per condition. It must be greater than 1.
#' @param nsubjects The mean of the cluster sizes, or a vector of cluster sizes for one arm.
#' @param d The difference in mean change between conditions (i.e. "difference-in-difference").
#' @param icc The intraclass correlation.
#' @param rho_c The correlation between baseline and post-test outcomes at the
#'   cluster level. This value can be used in both cross-sectional and cohort
#'   designs. If this quantity is unknown, a value of 0 is a conservative estimate.
#' @param rho_s The correlation between baseline and post-test outcomes at the
#'   subject level. This should be used for a cohort design or a mixture of cohort
#'   and cross-sectional designs. In a purely cross-sectional design (baseline subjects
#'   are completely different from post-test subjects), this value should be 0.
#' @param vart The total variation of the outcome (the sum of within- and between-cluster variation).
#' @param tol Numerical tolerance used in root finding. The default provides
#'   at least four significant digits.
#' @return The computed argument.
#' @examples 
#' # Find the number of clusters per condition needed for a trial with alpha = 0.05, 
#' # power = 0.80, nsubjects = 100, d = 0.50 units, icc = 0.05, rho_c = 0.50, rho_s = 0.70,
#' # and vart = 1 unit.
#' crtpwr.2meanD(nsubjects = 100 , d = 0.5, icc = 0.05, rho_c = 0.50, rho_s = 0.70, vart = 1)
#' # 
#' # The result, nclusters = 4.683358, suggests 5 clusters per condition should be used.
#' 
#' @references Rutterford C, Copas A, Eldridge S. (2015) Methods for sample size
#'   determination in cluster randomized trials. Int J Epidemiol. 44(3):1051-1067.
#' @references Teerenstra S, Eldridge S, Graff M, de Hoop E, Borm, GF. (2012) A simple
#'   sample size formula for analysis of covariance in cluster randomized trials.
#'   Statist Med. 31:2169-2178 
#'   
#' @export

crtpwr.2meanD <- function(alpha = 0.05, power = 0.80, nclusters = NA,
                          nsubjects = NA, d = NA, icc = NA,
                          rho_c = NA, rho_s = NA,
                          vart = NA,
                          tol = .Machine$double.eps^0.25){
  
  # if nsubjects is a vector, 
  if(length(nsubjects) > 1){
    nclusters <- length(nsubjects)
  }
  
  if(!is.na(nclusters) && nclusters <= 1) {
    stop("'nclusters' must be greater than 1.")
  }
  
  # list of needed inputs
  needlist <- list(alpha, power, nclusters, nsubjects, d, icc, rho_c, rho_s, vart)
  neednames <- c("alpha", "power", "nclusters", "nsubjects", "d", "icc",
                 "rho_c", "rho_s", "vart")
  needind <- which(unlist(lapply(needlist, is.na)))
  # check to see that exactly one needed param is NA
  
  if (length(needind) != 1) {
    neederror = "Exactly one of 'alpha', 'power', 'nclusters', 'nsubjects', 'd', 'icc', 'rho_c', 'rho_s', and 'vart' must be NA."
    stop(neederror)
  } 
  
  target <- neednames[needind]
  
  # evaluate power
  pwr <- quote({
    
    # variance inflation
    DEFF <- 1 + (nsubjects - 1)*icc
    r <- (nsubjects*icc/DEFF)*rho_c + ((1 - icc)/DEFF)*rho_s
    VIF <- DEFF*(1 - r^2)

    tcrit <- qt(alpha/2, 2*(nclusters - 1), lower.tail = FALSE)
    
    ncp <- sqrt(nclusters*nsubjects/(2*VIF)) * abs(d)/sqrt(vart)
    
    pt(tcrit, 2*(nclusters - 1), ncp, lower.tail = FALSE)#+
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
