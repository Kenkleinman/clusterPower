\name{fixed.effect}
\alias{fixed.effect}
\alias{fixed.effect.cluster.level}
\alias{random.effect}

\title{
Canned estimation functions for the power simulations.
}
\description{
These functions are designed to be used by the power.sim.XXX functions as the functions which estimate the treatment effect. They fit simple fixed and random effects models and return the estimated treatment effect. These functions are not designed to be called directly by the user.
}
\usage{
fixed.effect(dat, incl.period.effect, outcome.type, alpha)
random.effect(dat, incl.period.effect, outcome.type, alpha)
fixed.effect.cluster.level(dat, incl.period.effect, outcome.type, alpha)
}
\arguments{
  \item{dat}{ observed data as a data.frame with columns named, "y", "trt" and "clust". "per" column is optional if period.var==0. }
  \item{incl.period.effect}{ indicator of whether to include a period effect }
  \item{outcome.type}{  one of "gaussian", "binomial", "poisson" }
  \item{alpha}{ the type 1 error rate }
}
\details{
random.effect() relies on a call to glmer() from the lme4 package.  fixed.effect() relies on a call to glm(). fixed.effect.cluster.level() will save lots of time if you just want to run a cluster-level analysis and you have lots of observations.
}
\value{ A numeric vector with the following three elements, in order: [1] a point estimate for the treatment effect, [2] lower bound of (1-alpha) confidence interval, [3] lower bound of (1-alpha) confidence interval.
}
\author{
Nicholas G. Reich
}

\seealso{
\code{\link{power.sim.normal}}, \code{\link{power.sim.binomial}},\code{\link{power.sim.poisson}}
}


\keyword{ regression }
