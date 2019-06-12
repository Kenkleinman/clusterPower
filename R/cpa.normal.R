#' Power calculations for simple cluster randomized trials, continuous outcome
#'
#' Wrapper function that uses crtpwr.2mean to compute the power of 
#' a simple cluster randomized trial with a continuous outcome,
#' or determine parameters to obtain a target power. See ?crtpwr.2mean 
#' for more details.
#'#' @section Authors:
#' Alexandria Sakrejda, Jonathan Moyer (\email{jon.moyer@@gmail.com}), Ken Kleinman (\email{ken.kleinman@@gmail.com})
#'
#' @param alpha The level of significance of the test, the probability of a
#'   Type I error.
#' @param power The power of the test, 1 minus the probability of a Type II
#'   error.
#' @param n The number of clusters per condition. It must be greater than 1.
#' @param m The mean of the cluster sizes, or a vector of cluster sizes for one arm.
#' @param cv The coefficient of variation of the cluster sizes. When \code{cv} = 0,
#'   the clusters all have the same size.
#' @param difference The difference in condition means.
#' @param ICC The intraclass correlation.
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
#' cpa.normal(m=10, difference=1, ICC=.1, vart=5)
#' # 
#' # The result, showimg nclusters of greater than 15, suggests 16 clusters per 
#' # condition should be used.
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

cpa.normal <- function(alpha = 0.05, power = 0.80, n = NA,
                         m = NA, cv = 0,
                         difference = NA, ICC = NA,
                         vart = NA,
                         method = c("taylor", "weighted"),
                         tol = .Machine$double.eps^0.25){
  
  simple.analytic.normal <- crtpwr.2mean(alpha = alpha, power = power, 
                                         nclusters = n,
                                         nsubjects = m, cv = cv,
                                         d = difference, icc = ICC,
                                         vart = vart,
                                         method = method,
                                         tol = tol)
  return(simple.analytic.normal)
}