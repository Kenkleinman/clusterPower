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
#' Users must specify the number of subjects per cluster, number of clusters 
#' per treatment arm, the number of time steps, the baseline effect, the 
#' expected treatment effect, expected absolute difference between treatment 
#' arms, ICC, and significance level.
#' 
#' 
#' @param nsubjects Number of subjects per cluster; accepts either a scalar (equal cluster sizes) 
#' or a vector of length \code{nclusters} (user-defined size for each cluster) (required).
#' @param nclusters Number of clusters; accepts non-negative integer scalar (required).
#' @param steps Number of crossover steps; a baseline step (all clusters in non-treatment group) is assumed. 
#' Accepts positive scalar (indicating the total number of steps; clusters per step is obtained by 
#' \code{nclusters / steps}) or a vector of non-negative integers corresponding either to the number 
#' of clusters to be crossed over at each time point (e.g c(2,4,4,2); nclusters = 10) or the cumulative 
#' number of clusters crossed over by a given time point (e.g. c(2,4,8,10); nclusters = 10) (required).
#' @param alpha Significance level. Default = 0.05.
#' @param quiet When set to FALSE, displays simulation progress and estimated completion time; default is FALSE.
#' 
#' @return The estimated power.
#' 
#' @examples 
#' \dontrun{
#' }
#'
#' @author Alexandria C. Sakrejda (\email{acbro0@@umass.edu})
#' @author Ken Kleinman (\email{ken.kleinman@@gmail.com})
#'
#' @references Zhou X, Liao X, Kunz LM, Normand ST, Wang M, Spiegelman D. A maximum 
#' likelihood approach to power calculations for stepped wedge designs of binary 
#' outcomes. Biostatistics. 2020 Jan 1;21(1):102-121. doi: 10.1093/biostatistics/kxy031
#' 
#' @export
# cpa.sw.binary(nclusters (I), ntimes(J), nsubjects(K), d(delta, p0totalchange), ICC(rho0), beta, mu)



cpa.sw.binary <- function(nclusters = 12, 
                            ntimes = 2, 
                            nsubjects = 40, 
                            d = -10, 
                            ICC = 0.01, 
                            beta = -2, 
                            mu = 0.18, 
                            tol = 1e-5, 
                            GQ = 100){
    
  ###### Define some FORTRAN-calling functions  ########
  
  der_likelihood_time <- function(mu = as.numeric(mu), 
                                  beta = as.numeric(beta), 
                                  gamma = as.numeric(gamma), 
                                  tau2 = as.numeric(tau2), 
                                  z0 = as.integer(z0),
                                  z1 = as.integer(z1), 
                                  XX = as.integer(XX), 
                                  ntimes = as.integer(ntimes), 
                                  nsubjects = as.integer(nsubjects), 
                                  a = as.numeric(a), 
                                  b = as.numeric(b), 
                                  mincomp = as.integer(mincomp), 
                                  maxcomp = as.integer(maxcomp), 
                                  GQ = as.numeric(GQ), 
                                  t = as.numeric(t), 
                                  wts = as.numeric(wts), 
                                  derlikelihood = as.numeric(derlikelihood), 
                                  prob = as.numeric(prob)){
    .Fortran("der_likelihood_time", mu = mu, beta = beta, gamma = gamma, tau2 = tau2, z0 = z0,
        z1 = z1, XX = XX, JJ = ntimes, KK = nsubjects, a = a, b = b, 
        mincomp = mincomp, maxcomp = maxcomp, GQ = GQ, GQX = t, GQW = wts, 
        derlikelihood = derlikelihood, prob = prob)
  return(derlikelihood)
    }
  
  p0 <- vector(mode = "numeric", length = ntimes)
  gamma <- vector(mode = "numeric", length = ntimes)
  p0[1] <- mu
  p11 <-  mu + beta
  p0stepchange <- d / (ntimes - 1)
  tau2 = ICC / (1 - ICC) * mu * (1 - mu)
  for (i in 2:ntimes) {
    p0[i] = p0[i - 1] + p0stepchange
  }
  gamma <- p0 - mu
  
  # mincomp and maxcomp are ntimes+2 vectors of 0 and 1's, 
  # representing the weights of gamma(1),...,gamma(ntimes), mu, beta.
  comp <- rep(0, times = (ntimes + 2))
  maxcomp <- comp
  mincomp <- comp
  a = 100 
  b = -100
  for (i in 1:ntimes) {
      temp = mu + gamma[i]
      if (temp < a) {
        a = temp
        mincomp <- comp
        mincomp[ntimes + 1] = 1
        mincomp[i] = 1
        }
      if (temp > b) {
        b = temp
        maxcomp <- comp
        maxcomp[ntimes + 1] = 1
        maxcomp[i] = 1
        }
      temp = mu + beta + gamma[i]
      if (temp < a) {
        a = temp
        mincomp <- comp
        mincomp[ntimes + 1] = 1
        mincomp[ntimes + 2] = 1
        mincomp[i] = 1
        }
      if (temp > b) {
        b = temp
        maxcomp <- comp
        maxcomp[ntimes + 1] = 1
        maxcomp[ntimes + 2] = 1
        maxcomp[i] = 1
        }
  }
  rm(comp)
    
    a = -a
    b = 1 - b
    GQholder <- statmod::gauss.quad(GQ, kind = "legendre", alpha = 0, beta = 0)
    t <- GQholder[[1]]
    wts <- GQholder[[2]]
    
    # scale the quadrature formula to a nonstandard interval
    shft = (a + b) / 2.0
    slp = (b - a) / 2.0
    st <- slp * t
    t = shft + st
    wts = wts * slp
    rm(slp)
    rm(st)

  #power = LinearPower_time(mu, beta, gamma, tau2, nclusters, ntimes, nsubjects, a, b, mincomp, maxcomp, GQ, GQX, GQW)
    DD <-  nclusters / (ntimes - 1)   # nclusters is a multiple of (ntimes-1)
  # assign intervention
    interventionX <- matrix(data = 0, nrow = (ntimes), ncol = (ntimes - 1))
    
    for (i in 1:(ntimes - 1)) {
      for (j in (i + 1):ntimes) {
        interventionX[j,i] <- 1
      }
    }
    
    invVar = matrix(0, nrow = (ntimes + 2), ncol = (ntimes + 2))
    derlikelihood = as.vector(0.0)
    
     for (n in 1:(ntimes - 1)) {
       z0 = rep(0, times = ntimes)
       finish = 0  #need this?
      while (isTRUE(finish < 1)) { #THIS IS SLOW
      XX <- interventionX[,n]
        z1 = nsubjects - z0
   # call der_likelihood_time(mu,beta,gamma,tau2, z0, z1, X(i,:), ntimes, nsubjects, a, b, &
   #                            mincomp, maxcomp, GQ, GQX, GQW, derlikelihood, prob)
    prob <- 0.0
    derlikelihood <- 0.0

        derlikelihood <- der_likelihood_time(mu = as.numeric(mu), 
                                             beta = as.numeric(beta), 
                                             gamma = as.numeric(gamma), 
                                             tau2 = as.numeric(tau2), 
                                             z0 = as.integer(z0),
                                             z1 = as.integer(z1), 
                                             XX = as.integer(XX), 
                                             ntimes = as.integer(ntimes), 
                                             nsubjects = as.integer(nsubjects), 
                                             a = as.numeric(a), 
                                             b = as.numeric(b), 
                                             mincomp = as.integer(mincomp), 
                                             maxcomp = as.integer(maxcomp), 
                                             GQ = as.numeric(GQ), 
                                             t = as.numeric(t), 
                                             wts = as.numeric(wts), 
                                             derlikelihood = as.numeric(derlikelihood), 
                                             prob = as.numeric(prob))
        
        
    #call vectorsquare(derlikelihood, ntimes+2, derlikelihood2)
          derlen <- length(derlikelihood)
          mat2 <- matrix(rep(derlikelihood, times = derlen) , nrow = derlen, 
            ncol = derlen, byrow = TRUE)
          derlikelihood2 <- mat2 * derlikelihood
          rm(mat2)
          rm(derlen)

    #call linearpower_time
          
          holder <- list()
          for (h in 1:(40^2)){
          holder[[h]] <- append(holder, invVar)
          }
          
        invVar = invVar + derlikelihood2 * prob
        
   # finish = updatez(z0, ntimes, nsubjects)
        finish = 0
        z0[1] = z0[1] + 1
        for (p in 1:(ntimes - 1)) {
          if (z0[p] > nsubjects) {
            z0[p] = 0
            z0[p + 1] = z0[p + 1] + 1
            } else {
              break
            }
        }
        print(c(z0,z1))
        
        if (z0[ntimes] > nsubjects) {finish = 1}
      }
     }
    browser()
return(list(derlikelihood, derlikelihood2, invVar, prob))
  }
    
    ######################################
    ##### DEBUGGED TO HERE  ##############
    ######################################

