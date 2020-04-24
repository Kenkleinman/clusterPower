###### Define some FORTRAN-calling functions  ########

syminverse <- function(invVar = invVar, 
                       Var = Var, steps = steps){
  o = .Fortran("syminverse", a = as.numeric(invVar), 
               c = as.numeric(Var), n = as.integer(steps + 2))
  return(o)
}

vectorsquare <- function(derlikelihood = as.vector(derlikelihood), steps = steps){
  derlen = steps + 2
  derlikelihood2 = matrix(0, nrow = (steps + 2), ncol = (steps + 2))
  o = .Fortran("vectorsquare", a = as.numeric(derlikelihood), 
               n = as.integer(derlen), c = as.numeric(derlikelihood2))
  return(o)
}

der_likelihood_time <- function(mu = mu, 
                                beta = beta, 
                                gammaobj = gammaobj, 
                                tau2 = tau2, 
                                z0 = z0,
                                z1 = z1, 
                                XX = XX, 
                                steps = steps, 
                                nsubjects = nsubjects, 
                                a = a, 
                                b = b, 
                                mincomp = mincomp, 
                                maxcomp = maxcomp, 
                                GQ = GQ, 
                                t = t, 
                                wts = wts){
  stopifnot(length(gammaobj) == steps)
  stopifnot(length(z0) == steps)
  stopifnot(length(z1) == steps)
  stopifnot(length(XX) == steps)
  stopifnot(length(mincomp) == steps + 2)
  stopifnot(length(maxcomp) == steps + 2)
  stopifnot(length(t) == GQ)
  stopifnot(length(wts) == GQ)
  derlikelihood = rep(0.0, times = (steps + 2))
  prob = 0.0
  o = .Fortran("der_likelihood_time", mu = as.numeric(mu), beta = as.numeric(beta), 
               gammaobj = as.numeric(gammaobj), tau2 = as.numeric(tau2), z0 = as.integer(z0),
               z1 = as.integer(z1), XX = as.integer(XX), JJ = as.integer(steps), 
               KK = as.integer(nsubjects), a = as.numeric(a), b = as.numeric(b), 
               mincomp = as.integer(mincomp), maxcomp = as.integer(maxcomp), 
               GQ = as.integer(GQ), GQX = as.numeric(t), GQW = as.numeric(wts), 
               derlikelihood = as.numeric(derlikelihood), prob = as.numeric(prob))
  return(o)
}  




#' Power simulations for cluster-randomized trials: Stepped Wedge Design, Binary Outcome
#'
#' This function uses a modified Cox method to determine power for stepped 
#' wedge cluster-randomized controlled trials. Users can modify a variety of 
#' parameters to suit their desired experimental situation.
#' 
#' The stepped wedge trial design is a type of cross-over
#' design in which clusters change treatments in waves. Initially all the 
#' clusters recieve the same standard treatment, and at the end of the trial all
#' of the clusters will be recieving the treatment of interest. More than one 
#' cluster can change treatments in a wave, but the order in which clusters 
#' change treatments is randomly determined. The outcome of interest is assessed 
#' in each cluster during each wave.
#' 
#' Users must specify the number of subjects per cluster, number of clusters, 
#' the number of time steps, the baseline effect, the expected treatment effect, 
#' expected absolute difference between treatment 
#' arms, ICC, and time effect.
#' 
#' Fortran for this package was kindly provided by Dr. Zhou and others. 
#' 
#' @param nsubjects Number of subjects per cluster; accepts a scalar. Equal cluster sizes 
#' are assumed (required).
#' @param nclusters Number of clusters; accepts non-negative integer scalar (required).
#' @param steps Number of crossover steps; Accepts positive scalar indicating the total 
#' number of steps (required).
#' @param d Total change over the study period (assume that time effects are linear 
#' across time steps); accepts numeric (required).
#' @param ICC ICC for participants at first time step; accepts numeric (required).
#' @param beta Estimated treatment effect; accepts numeric (required).
#' @param mu Estimated baseline effect; accepts numeric (required).
#' @param tol Machine tolerance. Accepts numeric. Default is 1e-5. 
#' @param GQ Number of quadriture points used in Gaussian Legendre integration; accepts 
#' a scalar. Default is 100.
#' 
#' @return The estimated power.
#' 
#' @examples 
#' \dontrun{
#' holder <- cpa.sw.binary(nclusters = 12, 
#'   steps = 3, 
#'   nsubjects = 90, 
#'   d = -0.05, 
#'   ICC = 0.01, 
#'   beta = -0.5, 
#'   mu = 0.18, 
#'   tol = 1e-5, 
#'   GQ = 100)
#' }
#'
#' @author Alexandria C. Sakrejda (\email{acbro0@@umass.edu})
#' @author Ken Kleinman (\email{ken.kleinman@@gmail.com})
#'
#' @references Zhou X, Liao X, Kunz LM, Normand ST, Wang M, Spiegelman D. A maximum 
#' likelihood approach to power calculations for stepped wedge designs of binary 
#' outcomes. Biostatistics. 2020 Jan 1;21(1):102-121. doi: 10.1093/biostatistics/kxy031
#' @export
cpa.sw.binary <- function(nclusters = NA, 
                            steps = NA, 
                            nsubjects = NA, 
                            d = NA, 
                            ICC = NA, 
                            beta = NA, 
                            mu = NA, 
                            tol = 1e-5, 
                            GQ = 100){
  
  ######## main function code ###################
  
  
  p0 <- vector(mode = "numeric", length = steps)
  gammaobj <- vector(mode = "numeric", length = steps)
  p0[1] <- mu
  p11 <-  mu + beta
  p0stepchange <- d / (steps - 1)
  tau2 <-  ICC / (1 - ICC) * mu * (1 - mu)
  for (i in 2:steps) {
    p0[i] <-  p0[i - 1] + p0stepchange
  }
  gammaobj <- p0 - mu
  
  # mincomp and maxcomp are steps+2 vectors of 0 and 1's, 
  # representing the weights of gammaobj(1),...,gammaobj(steps), mu, beta.
  comp <- rep(0, times = (steps + 2))
  maxcomp <- comp
  mincomp <- comp
  a <-  100 
  b <- -100
  for (i in 1:steps) {
      temp = mu + gammaobj[i]
      if (temp < a) {
        a = temp
        mincomp <- comp
        mincomp[steps + 1] = 1
        mincomp[i] = 1
        }
      if (temp > b) {
        b = temp
        maxcomp <- comp
        maxcomp[steps + 1] = 1
        maxcomp[i] = 1
        }
      temp = mu + beta + gammaobj[i]
      if (temp < a) {
        a = temp
        mincomp <- comp
        mincomp[steps + 1] = 1
        mincomp[steps + 2] = 1
        mincomp[i] = 1
        }
      if (temp > b) {
        b = temp
        maxcomp <- comp
        maxcomp[steps + 1] = 1
        maxcomp[steps + 2] = 1
        maxcomp[i] = 1
        }
  }
  rm(comp)
    
    a <- -a
    b <-  1 - b
    GQholder <- statmod::gauss.quad(GQ, kind = "legendre", alpha = 0, beta = 0)
    t <- GQholder[[1]]
    wts <- GQholder[[2]]
    
    # scale the quadrature formula to a nonstandard interval
    shft <- (a + b) / 2.0
    slp <- (b - a) / 2.0
    st <- slp * t
    t <- shft + st
    wts <- wts * slp
    rm(slp)
    rm(st)

    DD <- nclusters / (steps - 1)   # nclusters is a multiple of (steps-1)
  # assign intervention
    interventionX <- matrix(data = 0, nrow = (steps - 1), ncol = (steps))
    
    for (i in 1:(steps - 1)) {
      for (j in (i + 1):steps) {
        interventionX[i,j] <- 1
      }
    }

    invVar <- matrix(0, nrow = (steps + 2), ncol = (steps + 2))
    
     for (i in 1:(steps - 1)) {
       z0 <- rep(0, times = steps)
       finish <- 0 
      while (finish < 1) {
      XX <- interventionX[i,]
        z1 <- nsubjects - z0

        Dholder <- der_likelihood_time(mu = mu, 
                                             beta = beta, 
                                             gammaobj = gammaobj, 
                                             tau2 = tau2, 
                                             z0 = z0,
                                             z1 = z1, 
                                             XX = XX, 
                                             steps = steps, 
                                             nsubjects = nsubjects, 
                                             a = a, 
                                             b = b, 
                                             mincomp = mincomp, 
                                             maxcomp = maxcomp, 
                                             GQ = GQ, 
                                             t = t, 
                                             wts = wts)
       
        
        prob <- Dholder$prob
        derlikelihood <- Dholder$derlikelihood
        
    VecHolder <- vectorsquare(derlikelihood = derlikelihood, steps = steps)
    derlikelihood2 <- VecHolder$c

    invVar <- invVar + derlikelihood2 * prob

        #update the counter
        finish <- 0
        z0[1] <- z0[1] + 1
        for (p in 1:(steps - 1)) {
          if (z0[p] > nsubjects) {
            z0[p] <- 0
            z0[p + 1] <- z0[p + 1] + 1
            } else {
              break
            }
        }
        
        if (z0[steps] > nsubjects) {finish = 1}
      }
     }

    mat <- matrix(0, nrow = (steps + 2), ncol = (steps + 2))
    mat <- syminverse(invVar = invVar, Var = mat, steps = steps)
    Var <- matrix(mat$c, nrow = (steps + 2), ncol = (steps + 2))

    sebeta <- sqrt(Var[2,2] / DD)
    power <- pnorm(beta/sebeta - 1.959964, FALSE) + pnorm(-beta/sebeta - 1.959964, FALSE)
    return(power)
  }
    

