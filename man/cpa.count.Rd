% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cpa.count.R
\name{cpa.count}
\alias{cpa.count}
\title{Analytic power calculations for parallel arm cluster-randomized trials with count outcomes}
\usage{
cpa.count(
  alpha = 0.05,
  power = 0.8,
  nclusters = NA,
  nsubjects = NA,
  r1 = NA,
  r2 = NA,
  CVB = NA,
  r1inc = TRUE,
  tol = .Machine$double.eps^0.25
)
}
\arguments{
\item{alpha}{The level of significance of the test, the probability of a
Type I error.}

\item{power}{The power of the test, 1 minus the probability of a Type II
error.}

\item{nclusters}{The number of clusters per condition. It must be greater than 1}

\item{nsubjects}{The number of units of person-time of observation per cluster}

\item{r1}{The mean event rate per unit time in one of the conditions}

\item{r2}{The mean event rate per unit time in the other condition}

\item{CVB}{The between-cluster coefficient of variation}

\item{r1inc}{Logical indicating if r1 is expected to be greater than r2. This is
only important to specify if one of r1 or r2 is NA.}

\item{tol}{Numerical tolerance used in root finding. The default provides
at least four significant digits.}
}
\value{
The computed value of the NA parameter (among \code{alpha}, \code{power}, \code{nclusters},
  \code{nsubjects},
  \code{r1}, \code{r2} and \code{CVB}) needed to satisfy the power and 
sample size equation.
}
\description{
\loadmathjax

Compute the power, number of clusters needed, number of subjects per cluster 
needed, or other key parameters for a simple parallel cluster randomized 
trial with a count outcome.

Exactly one of \code{alpha}, \code{power}, \code{nclusters}, \code{nsubjects},
  \code{r1}, \code{r2}, and \code{CVB} must be passed as \code{NA}.
  Note that \code{alpha} and \code{power} have non-\code{NA}
  defaults, so if those are the parameters of interest they must be
  explicitly passed as \code{NA}.
}
\section{Authors}{

Jonathan Moyer (\email{jon.moyer@gmail.com}), Ken Kleinman (\email{ken.kleinman@gmail.com})
}

\section{Note}{


 This function implements the approach of Hayes and Bennet (1999).  An estimate for the
 intracluster correlation coefficient (ICC) is used to calculate a design effect that 
 accounts for variance inflation due to clustering. 

 The coefficient of variation \code{CVB} is the variance of the cluster rates divided by the
 mean of the cluster rates.  
 
       
                   
 The CVB refers neither to 
 any natural parameter of a data generating model nor to any function of its parameters. 
 For this reason we do not offer the user a option to input
 the variance between the cluster means.  If you prefer to use that input, we suggest using the
 cps.count function.
     
        
  This function was inspired by work from Stephane Champely (pwr.t.test) and
  Peter Dalgaard (power.t.test). As with those functions, 'uniroot' is used to
  solve power equation for unknowns, so you may see
  errors from it, notably about inability to bracket the root when
  invalid arguments are given.  This generally means that no solution exists for which the 
  omitted parameter and the supplied parameters fulfill the equation.  In particular, the desired
  power may not be achievable with any number of subjects or clusters.
}

\section{Testing details}{

This function has been verified against reference values from 
\code{CRTsize::n4incidence}, and \code{clusterPower::cps.count}.
}

\examples{
# Find the number of clusters per condition needed for a trial with alpha = 0.05, 
# power = 0.80, 10 person-years per cluster, rate in condition 1 of 0.10 
# and condition 2 of 0.20, and CVB = 0.10.

cpa.count(nsubjects=10, r1=0.10, r2=0.20, CVB=0.10)

# The result, showimg nclusters of greater than 24, suggests 25 clusters per
# condition should be used.

# Find the largest CVB compatible with 80\% power when there are 25 clusters, 10
# subject-units of time per cluster, and a rate of 0.1 and 0.2 in each condition.  

cpa.count(nsubjects=10, nclusters= 25,r1=0.10, r2=0.20, CVB=NA)

# Results show that CVB as high as 0.107 can still yield power this high.

}
\references{
Donner A, Klar N. Design and Analysis of Cluster Randomization Trials in Health Research. Chichester, UK; 2009.

Hayes JR, Bennett S. Simple sample size calculation for cluster-randomized trials. International Journal of Epidemiology 1999; 28:319-326

Hayes JR, Moulton LH. Cluster Randomized Trials. Boca Raton, FL: CRC Press; 2009.
}
