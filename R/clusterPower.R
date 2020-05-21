#' clusterPower: power analysis and sample size calculations 
#' for cluster randomized trials and related designs.
#' 
#' @docType package
#' @name clusterPower
#' 
#' @description 
#' \loadmathjax
#'
#' The clusterPower package is design for experiments with correlated, A.K.A. clustered observations.
#' It contains many functions for calculating the power, sample size, 
#' and parameters necessary for acheiving desired power and sample size from analytic equations.  It also
#' estimation of power via Monte Carlo methods via very flexible functions for a wide range of 
#' methods for which analytic methods or approximations may not exist.
#' 
#' 
#' @section Introduction:
#' Most of the functions in this package are based on the generalized
#' linear mixed model approach to the analysis of cluster randomized 
#' trials.
#' 
#' For example, for approximately normal outcomes in a simple parallel design, 
#' the data generating model is:
#' 
#' \mjsdeqn{y_{ij}|b_i = \beta_0 + \beta_1 x_{ij} + b_i + e_{ij}}
#' 
#' where \eqn{i} indexes cluster and \eqn{j} indexes subject within
#' cluster.  In this case we assume random effect \mjseqn{b_i}, the cluster
#' means, are distributed Normal \mjseqn{(0, \sigma^2_b)} and the residual error
#' \mjseqn{e_{ij} ~ N(0,\sigma^2)}.  In this special case, we define the
#' "total variance" as \mjseqn{\sigma^2_b + \sigma^2} and the Intracluster
#' Correlation coefficient as \mjseqn{\sigma^2_b)/(\sigma^2_b + \sigma^2)}.  
#' The ICC is useful in some special cases as a simplifying statistic,
#' though it has no natural analogue in generalized linear mixed models
#' when the distribution lacks a variance that is independent of
#' the mean.
#' 
#' For example, for dichotomous outcomes, one covenient data generating model is:
#' 
#' \mjsdeqn{logit(Pr(y_{ij}|b_i = \beta_0 + \beta_1 x_{ij} + b_i}
#' 
#' 
#' 
#' 
NULL