% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cps.sw.normal.R
\name{cps.sw.normal}
\alias{cps.sw.normal}
\title{Power simulations for cluster-randomized trials: Stepped Wedge Design, Continuous Outcome.}
\usage{
cps.sw.normal(
  nsim = NULL,
  nsubjects = NULL,
  nclusters = NULL,
  mu0 = 0,
  mu1 = NULL,
  steps = NULL,
  sigma_sq = NULL,
  sigma_b_sq = NULL,
  alpha = 0.05,
  method = "glmm",
  quiet = FALSE,
  allSimData = FALSE,
  poorFitOverride = FALSE,
  lowPowerOverride = FALSE,
  timelimitOverride = TRUE,
  seed = NULL
)
}
\arguments{
\item{nsim}{Number of datasets to simulate; accepts integer (required).}

\item{nsubjects}{Number of subjects per cluster; accepts either a scalar (equal cluster sizes) 
or a vector of length \code{nclusters} (user-defined size for each cluster) (required).}

\item{nclusters}{Number of clusters; accepts non-negative integer scalar (required).}

\item{mu0}{Expected baseline mean; accepts numeric, default 0.  Required..}

\item{mu1}{Expected post-treatment mean; accepts numeric.  Required.}

\item{steps}{Number of crossover steps; a baseline step (all clusters in non-treatment group) is assumed. 
Accepts positive scalar (indicating the total number of steps; clusters per step is obtained by 
\code{nclusters / steps}) or a vector of non-negative integers corresponding either to the number 
of clusters to be crossed over at each time point (e.g c(2,4,4,2); nclusters = 10) or the cumulative 
number of clusters crossed over by a given time point (e.g. c(2,4,8,10); nclusters = 10) (required).}

\item{sigma_sq}{Within-cluster variance; accepts non-negative numeric scalar (indicating equal within-cluster variances for both 
treatment groups) or a vector of length 2 specifying within-cluster variances for the non-treatment and treatment groups, 
respectively (required).}

\item{sigma_b_sq}{Between-cluster variance; accepts non-negative numeric scalar (indicating equal 
between-cluster variances for both treatment groups) or a vector of length 2 specifying treatment-specific 
between-cluster variances (required).}

\item{alpha}{Significance level. Default = 0.05.}

\item{method}{Analytical method, either Generalized Linear Mixed Effects Model (GLMM) or 
Generalized Estimating Equation (GEE). Accepts c('glmm', 'gee') (required); default = 'glmm'.}

\item{quiet}{When set to FALSE, displays simulation progress and estimated completion time; default is FALSE.}

\item{allSimData}{Option to output list of all simulated datasets; default = FALSE.}

\item{poorFitOverride}{Option to override \code{stop()} if more than 25\%
of fits fail to converge; default = FALSE.}

\item{lowPowerOverride}{Option to override \code{stop()} if the power
is less than 0.5 after the first 50 simulations and every ten simulations
thereafter. On function execution stop, the actual power is printed in the
stop message. Default = FALSE. When TRUE, this check is ignored and the
calculated power is returned regardless of value.}

\item{timelimitOverride}{Logical. When FALSE, stops execution if the estimated completion time
is more than 2 minutes. Defaults to TRUE.}

\item{seed}{Option to set.seed. Default is NULL.}
}
\value{
A list with the following components
\itemize{
  \item Character string indicating total number of simulations and simulation type
  \item Number of simulations
  \item Data frame with columns "Power" (Estimated statistical power), 
               "lower.95.ci" (Lower 95% confidence interval bound), 
               "upper.95.ci" (Upper 95% confidence interval bound)
  \item Analytic method used for power estimation
  \item Significance level
  \item Vector containing user-defined cluster sizes
  \item Vector containing user-defined number of clusters
  \item Data frame reporting ICC, within & between cluster variances for 
  Treatment/Non-Treatment groups at each time point
  \item Vector containing expected means for each arm based on user inputs
  \item Data frame containing mean response values for each treatment group at each time point
  \item Matrix showing cluster crossover at each time point
  \item Data frame with columns: 
                  "Estimate" (Estimate of treatment effect for a given simulation), 
                  "Std.err" (Standard error for treatment effect estimate), 
                  "Test.statistic" (z-value (for GLMM) or Wald statistic (for GEE)), 
                  "p.value", 
                  "sig.val" (Is p-value less than alpha?)
  \item If \code{allSimData = TRUE}, a list of data frames, each containing: 
                  "y" (Simulated response value), 
                  "trt" (Indicator for treatment group),
                  "time.point" (Indicator for step; "t1" = time point 0) 
                  "clust" (Indicator for cluster), 
                  "period" (Indicator for at which step a cluster crosses over)
}

If \code{nofit = T}, a data frame of the simulated data sets, containing:

\itemize{
  \item "arm" (Indicator for treatment arm)
  \item "cluster" (Indicator for cluster)
  \item "y1" ... "yn" (Simulated response value for each of the \code{nsim} data sets).
  }
}
\description{
This set of functions utilize iterative simulations to determine 
approximate power for stepped wedge cluster-randomized controlled trials. Users 
can modify a variety of parameters to suit the simulations to their
desired experimental situation.
}
\details{
Runs power simulations for stepped wedge cluster-randomized controlled trials 
with continuous outcome. The stepped wedge trial design is a type of cross-over
design in which clusters change treatments in waves. Initially all the 
clusters recieve the same standard treatment, and at the end of the trial all
of the clusters will be recieving the treatment of interest. More than one 
cluster can change treatments in a wave, but the order in which clusters 
change treatments is randomly determined. The outcome of interest is assessed 
in each cluster during each wave.

Users must specify the desired number of simulations, number of subjects per 
cluster, number of clusters per arm, expected means for each arm, 
within-cluster variance, between-cluster variance, significance level, 
analytic method, progress updates, and simulated data set output may also be 
specified.
}
\examples{

# Estimate power for a trial with 3 steps and 9 clusters in arm 1 (often the 
# standard-of-care or 'control' arm) at the initiation of the study. Those 
# clusters have 14 subjects each, with sigma_b = 1 and sigma_b_sq = 1. We 
# have estimated arm outcome means of 1 and 2.1 in the first and second arms, 
# respectively, and 100 simulated data sets analyzed by the GLMM method. Using seed = 123, 
# the resulting power should be 0.82.

\dontrun{
normal.sw.rct = cps.sw.normal(nsim = 100, nsubjects = 14, nclusters = 9, 
                              mu0 = 1, mu1 = 2.1, steps = 3, sigma_sq = 1, 
                              sigma_b_sq = 1, alpha = 0.05, method = 'glmm', 
                              seed = 123)
}

}
\author{
Alexander R. Bogdan

Alexandria C. Sakrejda (\email{acbro0@umass.edu})

Ken Kleinman (\email{ken.kleinman@gmail.com})
}
