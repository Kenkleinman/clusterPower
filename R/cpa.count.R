#' Analytic power calculations for parallel arm cluster-randomized trials with count outcomes
#'
#' Wrapper function that uses crtpwr.2rate to compute the power of a 
#' simple cluster randomized trial with a count outcome,
#' or determine parameters to obtain a target power. For additional
#' details, see ?crtpwr.2rate help page.
#'
#' @section Authors:
#' Jonathan Moyer (\email{jon.moyer@@gmail.com}), Ken Kleinman (\email{ken.kleinman@@gmail.com})
#' 
#' @param alpha The level of significance of the test, the probability of a
#'   Type I error.
#' @param power The power of the test, 1 minus the probability of a Type II
#'   error.
#' @param nclusters The number of clusters per condition. It must be greater than 1.
#' @param py The number of person-years of observation per cluster.
#' @param r1 The expected mean event rate per unit time in the treatment group.
#' @param r2 The mean event rate per unit time in the control group.
#' @param cvb The between-cluster coefficient of variation.
#' @param r1inc Logical indicating if r1 is expected to be greater than r2. This is
#'   only important to specify if one of r1 or r2 is NA.
#' @param tol Numerical tolerance used in root finding. The default provides
#'   at least four significant digits.
#' @return The computed argument.
#' @examples 
#' # Find the number of clusters per condition needed for a trial with alpha = 0.05, 
#' # power = 0.80, 10 person-years per cluster, rate in condition 1 of 0.10 
#' # and condition 2 of 0.20, and CVB = 0.10.
#' cpa.count(m=10, r1=0.10, r2=0.20, CV=0.10)
#' # 
#' # The result, showimg nclusters of greater than 24, suggests 25 clusters per
#' # condition should be used.
#' 
#' @references Donner A, Klar N. Design and Analysis of Cluster Randomization Trials in Health Research. Chichester, UK; 2009.
#' 
#' @references Hayes JR, Moulton LH. Cluster Randomized Trials. Boca Raton, FL: CRC Press; 2009.
#' @export

cpa.count <- function(alpha = 0.05, power = 0.80,
                        n = NA, m = NA,
                        r1 = NA, r2 = NA,
                        CV = NA, r1inc = TRUE,
                        tol = .Machine$double.eps^0.25){
  
  simple.analytic.count <- crtpwr.2rate(alpha = alpha, power = power,
                                        nclusters = n, py = m,
                                        r1 = r1, r2 = r2,
                                        cvb = CV, r1inc = r1inc,
                                        tol = tol)
  return(simple.analytic.count)
}