\name{hayes.power.poisson}
\alias{hayes.power.poisson}

\title{
An implementation of power calculations for cluster-randomized study based on the coefficient of variation.
}

\description{
This function calculates the power for a specified cluster-randomized study based on the methods described by Hayes et al (1999).
}

\usage{
hayes.power.poisson(n.clusters, period.effect, btw.clust.var, at.risk.params, 
                    cluster.size, effect.size, alpha=.05)
}

\arguments{
  \item{n.clusters}{ number of clusters}
  \item{period.effect}{ period effect, on the link scale. See details. }
  \item{btw.clust.var}{the between-cluster variance}
  \item{at.risk.params}{ the expected at-risk time per individual in the study }
  \item{cluster.size}{ the number of individuals in each cluster}
  \item{effect.size}{ effect size, specified on the GLM link scale }
  \item{alpha}{ desired type I error rate}
}

\details{
Calculates, for a specified study design, the power of that study to detect the specified effect size. The model is specified as a Poisson log-linear random effects model (period.effect and btw.clust.var are parameters from the model specified in Reich et al (2012)). Based on this model specification, the coefficient of varation between cluster-level outcomes is calculated using conditional expectation (see mixed.eff.params()) and then the formula from Hayes and Bennett (1999) is implemented.
}

\value{ A numeric vector of length 1, containing the estimated power for the given study specifications.
}

\references{
	Reich NG et al.  PLoS ONE.  Empirical Power and Sample Size Calculations for Cluster-Randomized and Cluster-Randomized Crossover Studies. 2012.  \url{http://ow.ly/fEn39}
	Hayes RJ and Bennett S. Int J Epi. Simple sample size calculation for cluster-randomized trials. 1999. \url{http://www.ncbi.nlm.nih.gov/pubmed/10342698}
}

\author{
Nicholas G. Reich
}

\seealso{
\code{\link{mixed.eff.params}}
}

\examples{
    
hayes.power.poisson(n.clusters=36, period.effect=log(.015), btw.clust.var=0, 
                    at.risk.params=20, cluster.size=20, effect.size=log(.7))

}

\keyword{ random effects models }
\keyword{ power }
