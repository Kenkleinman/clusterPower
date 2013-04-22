\name{expit}
\alias{expit}
\alias{logit}

\title{
The expit and logit functions.
}
\description{
The expit and logit functions are useful shortcuts when using logistic regression models.
}
\usage{
expit(x)
logit(p)
}
\arguments{
  \item{x}{a real number}
  \item{p}{a number between 0 and 1, i.e. a probability}
}
\details{
The logit function is defined as logit(p) = log(p)/log(1-p) and can also be described as the log odds of a given probability.  The expit is the inverse of the logit function and is defined as expit(x) = exp(x)/(1+exp(x)).}
\value{
A numeric value, possibly vectorized, depending on the input.
}
\author{
Nicholas G. Reich
}

\examples{
expit(-2)
curve(expit(x), from=-5, to=5)

logit(.5)
curve(logit(x), from=0, to=1)

}
\keyword{ logistic }
