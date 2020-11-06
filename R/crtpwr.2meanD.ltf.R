#' Power calculations for difference-in-difference cluster randomized trials accounting for 
#' loss to and gain at follow-up, continuous outcome
#'
#' Compute the power of a difference-in-difference cluster 
#' randomized trial design with a continuous outcome,
#' accounting for potential loss to or gain at follow-up, or 
#' determine parameters to obtain a target power. 
#'
#' Exactly one of \code{alpha}, \code{power}, \code{nclusters}, \code{nsubjects},
#'   \code{d}, \code{icc}, \code{rho_c}, \code{rho_s}, and \code{vart}
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
#' @section Authors:
#' Jonathan Moyer (\email{jon.moyer@@gmail.com}), Ken Kleinman (\email{ken.kleinman@@gmail.com})
#'
#' @param alpha The level of significance of the test, the probability of a
#'   Type I error.
#' @param power The power of the test, 1 minus the probability of a Type II
#'   error.
#' @param nclusters The number of clusters per condition at baseline. It must be greater than 1.
#' @param nsubjects The number of subjects per cluster at baseline.
#' @param d The difference in mean change between conditions (i.e. "difference-in-difference").
#' @param icc The intraclass correlation.
#' @param rho_c The correlation between baseline and post-test outcomes at the
#'   cluster level. This value can be used in both cross-sectional and cohort
#'   designs. If this quantity is unknown, a value of 0 is a conservative estimate.
#' @param rho_s The correlation between baseline and post-test outcomes at the
#'   subject level assuming the a cohort design with no loss or gain to follow-up. Entering a value
#'   of 0 for this parameter is equivalent to a cross-sectional design and automatically sets \code{ltf_0}, \code{ltf_1},
#'   \code{gtf_0}, and \code{gtf_1} all to 1. 
#' @param vart The total variation of the outcome (the sum of within- and between-cluster variation).
#' @param ltf A length 2 vector indicating the proportion lost to follow-up per cluster in terms of \code{nsubjects}
#'   in the control and treatment arms. The first element corresponds to the proportion lost in the control
#'   group, while the second element correponds to the proportion lost in the treatment group. Entering a scalar
#'   for \code{ltf} will set both elements of the vector to that value. If \code{rho_s} is 0 then both elements 
#'   will be set to 1.
#' @param gtf A length 2 vector indicating the proportion gained at follow-up per cluster in terms of \code{nsubjects}
#'   in the control and treatment arms. The first element corresponds to the proportion gained in the control
#'   group, while the second element correponds to the proportion gained in the treatment group. By default, \code{gtf}
#'   is set to \code{ltf}. Entering a scalar for \code{gtf} will set both elements of the vector to that value. 
#'   If \code{rho_s} is 0 then both elements will be set to 1.
#' @param tol Numerical tolerance used in root finding. The default provides
#'   at least four significant digits.
#' @return The computed argument.
#' @examples 
#' # Find the number of clusters per condition needed for a trial with alpha = 0.05, 
#' # power = 0.80, nsubjects = 20, d = 0.50 units, icc = 0.05, rho_c = 0.50, rho_s = 0.70,
#' # and vart = 1 square unit if 50 percent of subjects in each cluster are lost to follow-up 
#' # and replaced.
#' crtpwr.2meanD.ltf(nsubjects = 100 , d = 0.5, icc = 0.05, 
#'   rho_c = 0.50, rho_s = 0.70, vart = 1,ltf=0.5)
#' # 
#' # The result, nclusters = 8.099772, suggests 9 clusters per condition should be used.
#' 
#' @references Rutterford C, Copas A, Eldridge S. (2015) Methods for sample size
#'   determination in cluster randomized trials. Int J Epidemiol. 44(3):1051-1067.
#' @references Teerenstra S, Eldridge S, Graff M, de Hoop E, Borm, GF. (2012) A simple
#'   sample size formula for analysis of covariance in cluster randomized trials.
#'   Statist Med. 31:2169-2178 
#'   
#' @export
#' 

crtpwr.2meanD.ltf <- function(alpha = 0.05, power = 0.80, nclusters = NA,
                              nsubjects = NA, d = NA, icc = NA,
                              rho_c = NA, rho_s = NA, vart = NA, 
                              ltf = c(0,0), gtf = ltf,
                              tol = .Machine$double.eps^0.25){
  
  if(!is.na(nclusters) && nclusters <= 1) {
    stop("'nclusters' must be greater than 1.")
  }
  
  if(!is.na(nsubjects) && nsubjects <= 1) {
    stop("'nsubjects' must be greater than 1.")
  }
  
  # rho_s = 0 corresponds to a cross-sectional design in which ltfs and gtfs equal 1
  if(!is.na(rho_s) && rho_s == 0) {
    ltf <- gtf <- c(1,1)
  }
  
  if(length(ltf) > 2 | length(gtf) > 2){
    stop("'ltf' and 'gtf' cannot contain have length greater than 2.")
  }
  
  if(length(ltf) == 1){
    ltf <- c(ltf,ltf)
  }
  
  if(length(gtf) == 1){
    gtf <- c(gtf,gtf)
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
    
    varc <- icc*rho_c*vart
    varct <- icc*(1-rho_c)*vart
    vars <- (1-icc)*rho_s*vart
    varst <- (1-icc)*(1-rho_s)*vart
    
    # ltf[1], gtf[1] are loss, gain to follow-up in control arm
    # ltf[2], gtf[2] are loss, gain to follow-up in treatment arm
    eta <- (1-ltf[1])/(1-ltf[1]+gtf[1]) + (1-ltf[2])/(1-ltf[2]+gtf[2]) - 1
    rho_s_star <- rho_s - 0.25*( 1/(1-ltf[1]+gtf[1]) + 1/(1-ltf[2]+gtf[2]) - 2*(eta*vars+varst)/(vars+varst) )
    
    vardid <- 4*( varct/nclusters + (1-rho_s_star)*(vars+varst)/(nclusters*nsubjects) )
    
    tcrit <- qt(alpha/2, 2*(nclusters - 1), lower.tail = FALSE)
    
    ncp <- abs(d)/sqrt(vardid)
    
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