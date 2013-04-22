\name{mixed.eff.params}
\alias{mixed.eff.params}

\title{
Calculation of variance in Poisson mixed model setting.
}
\description{
This function is designed to calculate the overall variance for cluster-level outcomes in a mixed-effect Poisson model. Conditional expectation calculations are implemented. 
}
\usage{
mixed.eff.params(pi0, btw.clust.var, Tk)
}
\arguments{
  \item{pi0}{ the baseline cluster-level mean on the scale of the link function }
  \item{btw.clust.var}{ the between-cluster-variance }
  \item{Tk}{ the at-risk time for each cluster  }
}
\details{
mixed.eff.params() is used by the hayes.power.poisson() function to compute the effective coefficient of variation, or k, for a particular study design.
}
\value{ A numeric vector with the following three named elements, in order: ["expectation"] the overall mean of cluster-level outcomes, ["variance"] the overall variance of cluster-level outcomes, ["hayes.k"] the estimated coefficient of variation.
}

\references{
	Reich NG et al.  PLoS ONE.  Empirical Power and Sample Size Calculations for Cluster-Randomized and Cluster-Randomized Crossover Studies. 2012.  \url{http://ow.ly/fEn39}
	Hayes RJ and Bennett S. Int J Epi. Simple sample size calculation for cluster-randomized trials. 1999. \url{http://www.ncbi.nlm.nih.gov/pubmed/10342698}
}

\author{
Nicholas G. Reich
}

\seealso{
\code{\link{hayes.power.poisson}}
}

\examples{
    
mixed.eff.params(pi0=log(1), btw.clust.var=.5, Tk=100)

}

\keyword{ random effects models }
