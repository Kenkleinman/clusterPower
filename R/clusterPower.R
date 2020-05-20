#' clusterPower: power analysis and sample size calculations 
#' for cluster randomized trials and related designs.
#' 
#' @docType package
#' @name clusterPower
#' @title A package for calculating power for experiments with clustered data
#' 
#' @section Introduction:
#' Most of the functions in this package are based on the generalized
#' linear mixed model approach to the analysis of cluster randomized 
#' trials.
#' 
#' For example, for approximately normal outcomes in a simple parallel design, 
#' the outcome can be modeled as:
#' 
#' \deqn{y_{ij} = \beta_0 + \beta_1 x_{ij} + b_i + e_{ij}}
#' 
#' where \eqn{i} indexes cluster and \eqn{j} indexes subject within
#' cluster.  In this case we assume random effect \eqn{b_i}, the cluster
#' means, are distributed Normal \eqn{(0, \sigma^2_b)} and the residual error
#' \eqn{e_{ij} ~ N(0,\sigma^2)}.  In this special case, we define the
#' "total variance" as \eqn{\sigma^2_b + \sigma^2} and the Intracluster
#' Correlation coefficient as \eqn{\sigma^2_b)/(\sigma^2_b + \sigma^2)}.  
#' The ICC is useful in some special cases as a simplifying statistic,
#' though it has no natural analogue in generalized linear mixed models
#' when the distribution lacks a variance that is independent of
#' the mean.
#' 
#' 
NULL