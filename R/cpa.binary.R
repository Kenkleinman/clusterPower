#' Analytic power calculations for parallel arm cluster-randomized trials with dichotomous outcomes
#'
#' Wrapper function that uses crtpwr.2prop to compute the power of 
#' a simple cluster randomized trial with a binary outcome,
#' or determine parameters to obtain a target power. See ?crtpwr.2prop
#' for additional details.
#' 
#' @section Authors:
#' Alexandria Sakrejda, Jonathan Moyer (\email{jon.moyer@@gmail.com}), 
#' Ken Kleinman (\email{ken.kleinman@@gmail.com})
#' @param alpha The level of significance of the test to be applied, the probability of a
#'   Type I error.
#' @param power The power of the test, 1 minus the probability of a Type II
#'   error.
#' @param n The number of clusters per condition. It must be greater than 1.
#' @param m The mean of the number of observations per cluster (mean cluster size).
#' @param CV The coefficient of variation of the cluster sizes. When \code{cv} = 0,
#'   the clusters all have the same size, or the same number of observations.
#' @param p1 The expected proportion with the outcome in the treatment group.
#' @param p2 The expected proportion with the outcome in the control group.
#' @param ICC The intraclass correlation.
#' @param pooled Logical indicating if pooled standard error should be used.
#' @param p1inc Logical indicating if p1 is expected to be greater than p2.
#' @param tol Numerical tolerance used in root finding. The default provides
#'   at least four significant digits.
#' @return The computed argument.
#' @examples 
#' # Find the number of clusters per condition needed for a trial with alpha = .05, 
#' # power = 0.8, 10 observations per cluster, no variation in cluster size, probability
#' # in condition 1 of .1 and condition 2 of .2, and icc = 0.1.
#' cpa.binary(m=10 ,p1=.1, p2=.2, ICC=.1)
#' # 
#' # The result, showimg nclusters of greater than 37, suggests 38 clusters per 
#' # condition should be used.
#' @references Donner A, Klar N. Design and Analysis of Cluster Randomization
#' Trials in Health Research. London; Arnold; 2000.
#' @references Wu S, Crespi CM, Wong WK. Comparison of Methods for Estimating Intraclass
#' Correlation Coefficient for Binary Responses in Cancer Prevention Cluster Randomized
#' Trials. Contemp Clin Trials. 2012; 33(5): 869-880. doi:10.1016/j.cct.2012.05.004 London: Arnold; 2000.
#' @export

cpa.binary <- function(alpha = 0.05, power = 0.80,
                         n = NA, m = NA, CV = 0,
                         p1 = NA, p2 = NA,
                         ICC = NA, pooled = FALSE,
                         p1inc = TRUE,
                         tol = .Machine$double.eps^0.25){
   
   simple.analytic.binary <- crtpwr.2prop(alpha = alpha, power = power,
                                          nclusters = n, nsubjects = m, 
                                          cv = CV,
                                          p1 = p1, p2 = p2,
                                          icc = ICC, pooled = pooled,
                                          p1inc = p1inc,
                                          tol = tol)
   return(simple.analytic.binary)
}