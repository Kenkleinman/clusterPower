#' Power calculations for stepped-wedge trials with a count outcome.
#' 
#' This function uses the \code{SWSamp} package by Gianluca Baio for 
#' estimating power based on analytic formula of Hussey and
#' Hughes (2007) where sample size calculations are based on 
#' an assumption of a normally-distributed outcome.
#' 
#' @param lambda1 Baseline rate for outcome of interest
#' @param RR Estimated relative risk of the intervention
#' @param nclusters Number of clusters
#' @param steps Number of time steps
#' @param nsubjects Average size of each cluster
#' @param ICC Intra-class correlation coefficient (default = 0.01)
#' @param alpha Significance level (default=0.05)
#' @param which.var String character specifying which variance to report.
#' Options are the default value \code{'within'} or \code{'total'}.
#' @param X A design matrix indicating the time
#' at which each of the clusters should switch to the intervention arm. 
#' Default is NULL and this matrix is automatically computed, but can it can 
#' be passed as a user-defined matrix with (nclusters) rows and (steps + 1) columns.
#' @param all.returned.objects Logical. Default = FALSE, indicating that only the
#' estimated power should be returned. When TRUE, all objects (listed below) are 
#' returned.
#' 
#' @return \item{power}{ The resulting power } 
#' 
#' When all.returned.objects = TRUE, returned items also include:
#' \item{sigma.y}{The estimated total (marginal) sd for the outcome} 
#' \item{sigma.e}{The estimated residual sd} 
#' \item{sigma.a}{The resulting cluster-level sd} 
#' \item{setting}{A list including the following values: 
#' - n.clusters = The number of clusters (nclusters)
#' - n.time.points = The number of steps in the SW design (steps)
#' - avg.cluster.size = The average cluster size (nsubjects)
#' - design.matrix = The design matrix for the SWT under consideration }
#' 
#' @author Alexandria C. Sakrejda (\email{acbro0@@umass.edu})
#' @author Ken Kleinman (\email{ken.kleinman@@gmail.com})
#' 
#' @references Baio, G; Copas, A; Ambler, G; Hargreaves, J; Beard, E; and Omar,
#' RZ Sample size calculation for a stepped wedge trial. Trials, 16:354. Aug
#' 2015.
#' 
#' Hussey M and Hughes J. Design and analysis of stepped wedge cluster
#' randomized trials. Contemporary Clinical Trials. 28(2):182-91. Epub 2006 Jul
#' 7. Feb 2007
#' @examples
#' 
#' cpa.sw.count(lambda1 = 1.75, RR = 0.9, nclusters = 21, steps = 6, nsubjects = 30, ICC = 0.01)
#' 
#' @export cpa.sw.count
cpa.sw.count <-
  function(lambda1,
           RR,
           nclusters,
           steps,
           nsubjects,
           ICC = 0.01,
           alpha = 0.05,
           which.var = "within",
           X = NULL,
           all.returned.objects = FALSE) {
    ## Validate user entries
    if (!is.integer(nclusters) ||
        nclusters < 1 ||
        length(nclusters) > 1 ||
        is.na(nclusters)) {
      errorCondition(message = "nclusters must be a positive scalar.")
    }
    if (!is.integer(steps) ||
        steps < 1 ||
        length(steps) > 1 ||
        is.na(steps)) {
      errorCondition(message = "steps must be a positive scalar.")
    }
    if (!is.integer(nsubjects) ||
        nsubjects < 1 ||
        length(nsubjects) > 1 ||
        is.na(nsubjects)) {
      errorCondition(message = "nsubjects must be a positive scalar.")
    }
    if (which.var != "total" ||
        which.var != "within" ||
        is.na(which.var)) {
      errorCondition(message = "which.var must be either 'total' or 'within'.")
    }
    if (!is.logical(all.returned.objects) ||
        is.na(all.returned.objects)) {
      errorCondition(message = "all.returned.objects must be logical.")
    }
    # define Baio's sw.design.mat fxn
    
    sw.design.mat <- function (I, J, H = NULL) 
    {
      if (sum(sapply(list(I, J, H), is.null)) != 1) {
        warning("exactly one of 'I', 'J' and 'H' must be NULL")
      }
      if (is.null(I)) {
        I <- H * J
      }
      if (is.null(J)) {
        J <- I/H
      }
      if (is.null(H)) {
        H <- I/J
      }
      X <- matrix(0, I, (J + 1))
      for (i in 2:(J + 1)) {
        X[1:((i - 1) * H), i] <- 1
      }
      row.names(X) <- sample(1:I, I)
      colnames(X) <- c("Baseline", paste0("Time ", 
                                          1:J))
      return(X)
    }
    
    
    #define Baio's HH.count
    HH.count <-
      function (lambda1,
                RR,
                I,
                J,
                K,
                rho = 0,
                alpha = 0.05,
                which.var = "within",
                X = NULL)
      {
        if (is.null(X)) {
          X <- sw.design.mat(I = I, J = J, H = NULL)
        }
        else {
          row.names(X) <- sample(1:I, I)
          colnames(X) <- c("Baseline", paste0("Time ",
                                              1:J))
        }
        U <- sum(X)
        W <- sum(apply(X, 2, sum) ^ 2)
        V <- sum(apply(X, 1, sum) ^ 2)
        lambda2 <- RR * lambda1
        theta <- abs(lambda1 - lambda2)
        if (which.var == "within") {
          sigma.e = (sqrt(lambda1) + sqrt(lambda2)) / 2
          sigma.a <- sqrt(rho * sigma.e ^ 2 / (1 - rho))
          sigma.y <- sqrt(sigma.e ^ 2 + sigma.a ^ 2)
        }
        if (which.var == "total") {
          sigma.y <- (sqrt(lambda1) + sqrt(lambda2)) / 2
          sigma.a <- sqrt(sigma.y ^ 2 * rho)
          sigma.e <- sqrt(sigma.y ^ 2 - sigma.a ^ 2)
        }
        sigma <- sqrt(sigma.e ^ 2 / K)
        v <- (I * sigma ^ 2 * (sigma ^ 2 + ((J + 1) * sigma.a ^ 2))) / ((I *
                                                                           U - W) * sigma ^
                                                                          2 + (U ^ 2 + I * (J + 1) * U - (J + 1) *
                                                                                 W - I * V) * sigma.a ^
                                                                          2)
        power <- pnorm(theta / sqrt(v) - qnorm(1 - alpha / 2))
        setting <-
          list(
            n.clusters = I,
            n.time.points = J,
            avg.cluster.size = K,
            design.matrix = X
          )
        list(
          power = power,
          lambda1 = lambda1,
          lambda2 = lambda2,
          sigma.y = sigma.y,
          sigma.e = sigma.e,
          sigma.a = sigma.a,
          setting = setting
        )
      }
    
    # execute the function
    o <-
      HH.count(
        lambda1 = lambda1,
        RR = RR,
        I = nclusters,
        J = steps,
        K = nsubjects,
        rho = ICC,
        alpha = alpha,
        which.var = which.var,
        X = X
      )
    if (all.returned.objects == FALSE) {
      o <- o$power
    }
    return(o)
  }