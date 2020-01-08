#' Power calculations for  difference-in-difference cluster randomized trials, dichotomous outcome
#'
#' Wrapper function that uses crtpwr.2propD() to compute the power 
#' of a difference-in-difference cluster randomized trial design 
#' with a binary outcome, or determine parameters to obtain a target 
#' power. For additional details, see ?crtpwr.2propD.
#' 
#' Exactly one of \code{alpha}, \code{power}, \code{n}, \code{m},
#'   \code{p}, \code{did}, \code{ICC}, \code{rho_c}, and \code{rho_s} must be passed as \code{NA}.
#'   Note that \code{alpha} and \code{power} have non-\code{NA}
#'   defaults, so if those are the parameters of interest they must be
#'   explicitly passed as \code{NA}.
#'
#' @section Authors:
#' Alexandria Sakrejda, Jonathan Moyer (\email{jon.moyer@@gmail.com}), Ken Kleinman (\email{ken.kleinman@@gmail.com})
#' 
#' @section Note:
#'   This function was inspired by work from Stephane Champely (pwr.t.test) and
#'   Peter Dalgaard (power.t.test). As with those functions, 'uniroot' is used to
#'   solve power equation for unknowns, so you may see
#'   errors from it, notably about inability to bracket the root when
#'   invalid arguments are given.
#'   
#' @param alpha The level of significance of the test, the probability of a
#'   Type I error.
#' @param power The power of the test, 1 minus the probability of a Type II
#'   error.
#' @param n The number of clusters per condition. It must be greater than 1.
#' @param m The mean of the cluster sizes.
#' @param p The expected mean proportion at the post-test, averaged across treatment and control arms.
#' @param did The expected absolute difference.
#' @param ICC The intraclass correlation.
#' @param rho_c The correlation between baseline and post-test outcomes at the
#'   cluster level. This value can be used in both cross-sectional and cohort
#'   designs. If this quantity is unknown, a value of 0 is a conservative estimate.
#' @param rho_s The correlation between baseline and post-test outcomes at the
#'   subject level. This should be used for a cohort design or a mixture of cohort
#'   and cross-sectional designs. In a purely cross-sectional design (baseline subjects
#'   are completely different from post-test subjects), this value should be 0.
#' @param tol Numerical tolerance used in root finding. The default provides
#'   at least four significant digits.
#' @return The computed argument.
#' @examples 
#' # Find the number of clusters per condition needed for a trial with alpha = .05, 
#' # power = 0.8, 50 observations per cluster, expected mean post-test proportion of .50,
#' # expected difference of .1, ICC = 0.05, cluster level correlation of 0.3, and subject level 
#' # correlation of 0.7.
#' cpa.did.binary(m=50 ,p=.5, did=.1, ICC=.05, rho_c=.3, rho_s=.7)
#' # 
#' # The result, showimg nclusters of greater than 32, suggests 33 clusters per 
#' # condition should be used.
#' 
#' @references Murray D. Design and Analysis of Group-Randomized Trials. New York, NY: Oxford
#' University Press; 1998.
#' 
#' @export

cpa.did.binary <- function(alpha = 0.05, power = 0.80,
                          n = NA, m = NA,
                          p = NA, did = NA, ICC = NA, 
                          rho_c = NA, rho_s = NA,
                          tol = .Machine$double.eps^0.25){
  
  did.analytic.binary <- crtpwr.2propD(alpha = alpha, power = power,
                            nclusters = n, nsubjects = m,
                            p = p, d = did, icc = ICC, 
                            rho_c = rho_c, rho_s = rho_s,
                            tol = tol)
  
  return(did.analytic.binary)
}