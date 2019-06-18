#' Power calculations for difference-in-difference cluster randomized trials, continuous outcome
#'
#' Wrapper function that uses crtpwr.2meanD() to compute the power of a 
#' difference-in-difference cluster randomized trial design with a 
#' continuous outcome, or determine parameters to obtain a target power. 
#' See ?crtpwr.2meanD for additional details.
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
#' Alexandria Sakrejda, Jonathan Moyer (\email{jon.moyer@@gmail.com}), Ken Kleinman (\email{ken.kleinman@@gmail.com})
#'
#' @param alpha The level of significance of the test, the probability of a
#'   Type I error.
#' @param power The power of the test, 1 minus the probability of a Type II
#'   error.
#' @param n The number of clusters per condition. It must be greater than 1.
#' @param m The mean of the cluster sizes.
#' @param did The difference in mean change between conditions (i.e. "difference-in-difference").
#' @param ICC The intraclass correlation.
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
#' # power = 0.80, m = 100, did = 0.50 units, ICC = 0.05, rho_c = 0.50, rho_s = 0.70,
#' # and vart = 1 unit.
#' cpa.did.normal(m = 100 , did = 0.5, ICC = 0.05, rho_c = 0.50, rho_s = 0.70, vart = 1)
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

cpa.did.normal <- function(alpha = 0.05, power = 0.80, n = NA,
                          m = NA, did = NA, ICC = NA,
                          rho_c = NA, rho_s = NA,
                          vart = NA,
                          tol = .Machine$double.eps^0.25){
  
  did.analytic.normal <- crtpwr.2meanD(alpha = alpha, 
                            power = power, nclusters = n,
                            nsubjects = m, d = did, icc = ICC,
                            rho_c = rho_c, rho_s = rho_s,
                            vart = vart,
                            tol = tol)
  
  return(did.analytic.normal)
}