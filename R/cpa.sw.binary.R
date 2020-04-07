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
#@export
# cpa.sw.binary(nclusters (I), ntimes(J), nsubjects(K), d(delta, p0totalchange), ICC(rho0), beta, mu)


cpa.sw.binary <- function(nclusters, ntimes, nsubjects, d, ICC, beta, mu, 
                          tol = 1e-5, GQ = 100){
  ###delete this later
  nclusters = 12
  ntimes = 3
  nsubjects = 90
  mu = 0.18
  beta = -0.05
  ICC = 0.01
  d = -0.01
  GQ = 100
  tol = 1e-5
  ######
  
  
  p0 <- vector(mode = "numeric", length = ntimes)
  gamma <- vector(mode = "numeric", length = ntimes)
  p0[1] <- mu
  p11 <-  mu + beta
  p0stepchange <- d / (ntimes - 1)
  tau2 = ICC / (1 - ICC) * mu * (1 - mu)
  gamma[1] <- 0.0
  for (i in 2:ntimes) {
    p0[i] = p0[i - 1] + p0stepchange
    gamma[i] <- p0[i] - mu
  }
  
  # mincomp and maxcomp are ntimes+2 vectors of 0 and 1's, 
  # representing the weights of gamma(1),...,gamma(ntimes), mu, beta.
  comp <- rep(0, times = (ntimes + 2))
  mincomp <- comp
  maxcomp <- comp
  
  if (p0totalchange > tol || p0totalchange < -tol) {
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
#####################################  
 #   Debugged to here
#####################################
    
    
  #power = LinearPower_time(mu, beta, gamma, tau2, nclusters, ntimes, nsubjects, a, b, mincomp, maxcomp, GQ, GQX, GQW)
    DD <-  nclusters / (ntimes - 1)   # nclusters is a multiple of (ntimes-1)
  # assign intervention
    interventionX <- matrix(data = 1, nrow = (ntimes - 1), ncol = (ntimes - 1))
    invVar = 0.0
    z0 = as.vector(0)
    finish = 0  #need this?
    for (i in 1:(ntimes - 1)) {
      while (isTRUE(finish < 1)) {
        z1 = nsubjects - z0
   # call der_likelihood_time(mu,beta,gamma,tau2, z0, z1, X(i,:), ntimes, nsubjects, a, b, &
   #                            mincomp, maxcomp, GQ, GQX, GQW, derlikelihood, prob)
    
        likelihoodf_denom = 0.0
        likelihoodf_numer = 0.0
        likelihoodf_denomb2 = 0.0
        derlikelihood_mu = 0.0
        derlikelihood_beta = 0.0
        derlikelihood_tau2 = 0.0
        derlikelihood_gamma = as.vector(0.0)
        prob = 0.0
        for (i in 1:GQ) {
          x = GQX[i]
          exx = dexp(-0.5 * x * x / tau2)
          ff = 1.0
          ffprob = 1.0
          ff_mu = 0.0
          ff_beta = 0.0
            for (j in 1:ntimes) {
              ff1 = mu + beta * XX[j] + gamma[j] + x  #WTF is XX???
              ff0 = 1 - ff1
              ff = ff * (ff0^z0[j]) * (ff1^z1[j])
    # since GQX are not at limits, we ignore the cases where ff0=0 or ff1=0
              temp = z1[j] / ff1 - z0[j] / ff0
              ff_mu = ff_mu + temp
              ff_beta = ff_beta + temp * XX[j]  ##WTF is XX????
              k = j - 1
              if (k > 0) {
                ff_gamma[k] = temp # where is ff_gamma specified?
                }
              ff01 = ff0 * ff1
    # compute binomial
              if (z0[j] < z1[j]) {
                ffprob = ffprob * ff1^(z1[j] - z0[j])
                for (k in 0:(z0[j] - 1)) {
                  ffprob = ffprob * (nsubjects - k) / (z0[j] - k) * ff01
                  }
                } else {      
                  ffprob = ffprob * ff0^(z0[j] - z1[j])
                  for (k in 0:(z1[j] - 1)) {
                    ffprob = ffprob * (nsubjects - k) / (z1[j] - k) * ff01
                    }
                }
              }
      # compute the possible combinations of (z0,z1)
          prob = prob + GQW[i] * ffprob * exx
          likelihoodf_denom = likelihoodf_denom + GQW[i] * exx
          likelihoodf_numer = likelihoodf_numer + GQW[i] * ff * exx
          likelihoodf_denomb2 = likelihoodf_denomb2 + GQW[i] * x * x * exx
          derlikelihood_mu = derlikelihood_mu + GQW[i] * ff * ff_mu * exx
          derlikelihood_beta = derlikelihood_beta + GQW[i] * ff * ff_beta * exx
          derlikelihood_gamma = derlikelihood_gamma + GQW[i] * ff * ff_gamma * exx
          derlikelihood_tau2 = derlikelihood_tau2 + GQW[i] * ff * x * x * exx
          }
    # calculate f(a)exp(-0.5*a*a)
        eaa = dexp(-0.5 * a * a / tau2)
        ff = 1.0
        for (i in 1:ntimes) {
          ff1 = mu + beta * XX[j] + gamma[j] + a  #still haven't specified XX
          ff0 = 1 - ff1
          ff = ff * (ff0^z0[j]) * (ff1^z1[j])
          }
        faeaa = ff * eaa
    # calculate f(b)exp(-0.5*b*b)
        ebb = dexp(-0.5 * b * b / tau2)
        ff = 1.0
        for (i in 1:ntimes) {
          ff1 = mu + beta * XX[j] + gamma[j] + b
          ff0 = 1 - ff1
          ff = ff * (ff0^z0[j]) * (ff1^z1[j])
          }
        fbebb = ff * ebb
    
    # calculate derlikelihood_mu
        derlikelihood_mu = derlikelihood_mu + faeaa * mincomp[ntimes + 1] - fbebb * maxcomp[ntimes + 1]
        derlikelihood_mu = derlikelihood_mu / likelihoodf_numer - (eaa * mincomp[ntimes + 1] - ebb * maxcomp[ntimes + 1]) / likelihoodf_denom
    # calculate derlikelihood_beta
        derlikelihood_beta = derlikelihood_beta + faeaa * mincomp[ntimes + 2] - fbebb * maxcomp[ntimes + 2]
        derlikelihood_beta = derlikelihood_beta / likelihoodf_numer - (eaa * mincomp[ntimes + 2] - ebb * maxcomp[ntimes + 2]) / likelihoodf_denom
    # calculate derlikelihood_gamma
        for (j in 2:ntimes) {
          k = j - 1
          derlikelihood_gamma[k] = derlikelihood_gamma[k] + faeaa * mincomp[j] - fbebb * maxcomp[j]  #is this a vector now?
          derlikelihood_gamma[k] = derlikelihood_gamma[k] / likelihoodf_numer - (eaa * mincomp[j] - ebb * maxcomp[j]) / likelihoodf_denom #is this a vector now?
          }
    
    # calculate derlikelihood_tau2
        derlikelihood_tau2 = 0.5 * (derlikelihood_tau2 / likelihoodf_numer - likelihoodf_denomb2 / likelihoodf_denom) / tau2 / tau2
        prob = prob / likelihoodf_denom
        derlikelihood[1] = derlikelihood_mu
        derlikelihood[2] = derlikelihood_beta
        derlikelihood[3:(ntimes + 1)] = derlikelihood_gamma
        derlikelihood[ntimes + 2] = derlikelihood_tau2
  
    #call vectorsquare(derlikelihood, ntimes+2, derlikelihood2)
        mat <- as.matrix(nrow = length(derlikelihood), ncol = length(derlikelihood))
        for (i in 1:ntimes - 1) {
          mat[i,i] = derlikelihood[i] * derlikelihood[i]
          for (j in (i + 1):(ntimes + 2)) {
            mat[i,j] = derlikelihood[i] * derlikelihood[j]
            mat[j,i] = mat[i,j]
            }
          }
        mat[(ntimes + 2),(ntimes + 2)] = derlikelihood[(ntimes + 2)] * derlikelihood[(ntimes + 2)]
    
    #call linearpower_time
        derlikelihood2 <- mat
        invVar = invVar + derlikelihood2 * prob
   # finish = updatez(z0, ntimes, nsubjects)
   
        z0[1] = z0[1] + 1
        for (j in 1:(ntimes - 1)) {
          if (z0[j] > nsubjects) {
            z0[j] = 0
            z0[j + 1] = z0[j + 1] + 1
            } else {
              break
            }
          }
        if (z0[ntimes] > nsubjects) {finish = 1}
        }
      }
  
  #call syminverse(invVar,Var,ntimes+2)
    k = 0
    aa <- as.vector(0)
    for (i in 1:(ntimes + 2)) {
      for (j in 1:i) {
        k = k + 1
        aa[k] = invVar[i,j]  #wait, this is a matrix?
        }
      }
  
  #call syminv(aa, n, cc, nullty, ifault)
  #syminv ( a, n, c, nullty, ifault )
    n = ntimes
   #
  #  Compute the Cholesky factorization of A.
  #  The result is stored in C.
   #
    nn = (n * (n + 1)) / 2
  #call cholesky ( a, n, nn, c, nullty, ifault )
   # cholesky ( a, n, nn, u, nullty, ifault )
    cholmat <- chol(a) #this doesn't make sense, originally c
  #back to synverse
    k = 0
    for (i in 1:n) {
      for (j in 1:i - 1) {
        k = k + 1
        c[i,j] = cc[k]
        c[j,i] = cc[k]
        }
      k = k + 1
      c[i,i] = cc[k]
      }
  #end synverse
    sebeta = sqrt(Var[2,2] / DD)
    power = alnorm(beta / sebeta - 1.959964, upper) + #is alnorm==dnorm?
      alnorm(-beta / sebeta - 1.959964, upper)
    } else {
      if (beta > 0) {
        a = -mu
        b = 1 - mu - beta
      } else {
        a = -mu - beta
        b = 1 - mu
      }
    
    #call legendre_handle (GQ, a, b, GQX, GQW)
    # Gaussian Legendre will not take two limits, a and b
      GQholder <- statmod::gauss.quad(GQ, kind = "legendre", alpha = alpha, beta = beta)
      GQX <- GQholder[[1]]
      GQW <- GQholder[[2]]
    
    #power = LinearPower_notime(mu, beta, tau2, nclusters, ntimes, nsubjects, a, b, GQ, GQX, GQW)

      NI = ntimes * nsubjects
      DD = nclusters / (ntimes - 1)
    
      z0 = 0
      for (j in 1:(ntimes - 1)) {
        z0[j] = j * nsubjects
        }
      z1 = NI - z0
      h11 = 0.0
      h12 = 0.0
      h13 = 0.0
      h22 = 0.0
      h23 = 0.0
      h33 = 0.0
      for (i in 1:(ntimes - 1)) {
        for (z00 in 0:z0[i]) {
          z01 = z0[i] - z00
          for (z10 in (0:z1[i])) {
            z11 = z1[i] - z10
    #call der_likelihood_notime(mu, beta, tau2, z00, z01, z10, z11, GQ, GQX, GQW, &
    #                             derlikelihood_mu, derlikelihood_beta, derlikelihood_tau2, prob)
            z0 = z00 + z01
            z1 = z10 + z11
            derlikelihood_mu = 0.0
            derlikelihood_beta = 0.0
            derlikelihood_tau2 = 0.0
            likelihoodf_denom = 0.0
            likelihoodf_denomb2 = 0.0
            likelihoodf_numer = 0.0
            prob = 0.0
            for (i in 1:GQ) {
              x = GQX[i]
              ff = 1.0
              ffprob = 1.0
              ff01 = mu + x
              ff00 = 1 - ff01
              ff11 = mu + beta + x
              ff10 = 1 - ff11
              exx = dexp(-0.5 * x * x / tau2)
              ff = (ff00^z00) * (ff01^z01) * (ff10^z10) * (ff11^z11)
              likelihoodf_numer = likelihoodf_numer + GQW[i] * ff * exx
              likelihoodf_denom = likelihoodf_denom + GQW[i] * exx
              likelihoodf_denomb2 = likelihoodf_denomb2 + GQW[i] * x * x * exx
          
        #  ! since GQX are not at limits, we ignore the cases where ff00=0 or ff01=0 or ff10=0 or ff11=0
              ff1 = ff * (z01 / ff01 - z00 / ff00 + z11 / ff11 - z10/ff10)
              derlikelihood_mu = derlikelihood_mu + GQW[i] * ff1 * exx
              ff1 = ff * (z11 / ff11 - z10 / ff10)
              derlikelihood_beta = derlikelihood_beta + GQW[i] * ff1 * exx
              ff1 = ff * x * x
              derlikelihood_tau2 = derlikelihood_tau2 + GQW[i] * ff1 * exx
          
          #! prob
              ff0prob = ff00 * ff01
              ff1prob = ff10 * ff11
              if (z00 < z01) {
                ffprob = ffprob * ff01^(z01 - z00)
                for (k in 0:(z00 - 1)) {
                  ffprob = ffprob * (z0 - k) / (z00 - k) * ff0prob
                  }
                } else {
                  ffprob = ffprob * ff00^(z00 - z01)
                  for (k in 0:(z01 - 1)) {
                    ffprob = ffprob * (z0 - k) / (z01 - k) * ff0prob
                    }
                  }
              if (z10 < z11) {
                ffprob = ffprob * ff11^(z11 - z10)
                for (k in 0:(z10 - 1)) {
                  ffprob = ffprob * (z1 - k) / (z10 - k) * ff1prob
                  }
                } else {
                  ffprob = ffprob * ff10^(z10 - z11)
                  for (k in 0:(z11 - 1)) {
                    ffprob = ffprob * (z1 - k) / (z11 - k) * ff1prob
                    }
                  }
              prob = prob + GQW[i] * ffprob * exx
              }
            if (beta >= 0) {
          #! we don't consider the case of beta = 0
              exx1 = dexp(-0.5 * mu * mu / tau2)
              exx2 = dexp(-0.5 * (1 - mu - beta) * (1 - mu - beta) / tau2)
              if (z01 == 0) {
                derlikelihood_mu = derlikelihood_mu + ((1 - beta)^z10) * (beta^z11) * exx1
                }
              if (z10 == 0) {
                temp = ((1 - beta)^z01) * (beta^z00) * exx2
                derlikelihood_mu = derlikelihood_mu - temp
                derlikelihood_beta = derlikelihood_beta - temp
                }
              derlikelihood_mu = derlikelihood_mu / likelihoodf_numer - (exx1 - exx2) / likelihoodf_denom
              derlikelihood_beta = derlikelihood_beta / likelihoodf_numer + exx2 / likelihoodf_denom
              } else {
     #   ! beta < 0
                exx1 = dexp(-0.5 * (mu + beta) * (mu + beta) / tau2)
                exx2 = dexp(-0.5 * (1 - mu) * (1 - mu) / tau2)
                if (z00 == 0) {
                  derlikelihood_mu = derlikelihood_mu - ((-beta)^z10) * ((1 + beta)^z11) * exx2
                  }
                if (z11 == 0) {
                  temp = ((-beta)^z01) * ((1 + beta)^z00) * exx1
                  derlikelihood_mu = derlikelihood_mu + temp
                  derlikelihood_beta = derlikelihood_beta + temp
                  }
                derlikelihood_mu = derlikelihood_mu / likelihoodf_numer - (exx1 - exx2)/likelihoodf_denom
                derlikelihood_beta = derlikelihood_beta / likelihoodf_numer + exx1 / likelihoodf_denom
                }
            derlikelihood_tau2 = 0.5 * (derlikelihood_tau2 / likelihoodf_denom - likelihoodf_denomb2 / likelihoodf_denom) / tau2 / tau2
            prob = prob / likelihoodf_denom

            h11 = h11 + derlikelihood_mu * derlikelihood_mu * prob
            h22 = h22 + derlikelihood_beta * derlikelihood_beta * prob
            h33 = h33 + derlikelihood_tau2 * derlikelihood_tau2 * prob
            h12 = h12 + derlikelihood_mu * derlikelihood_beta * prob
            h13 = h13 + derlikelihood_mu * derlikelihood_tau2 * prob
            h23 = h23 + derlikelihood_beta * derlikelihood_tau2 * prob
            }
          }
        }
      sebeta = sqrt(abs((h33 * h11 - h13 * h13) / (h11 * h22 * h33 + 2.0 * h12 * h23 * h13 - h13 * h13 * h22 - h12 * h12 * h33 - h23 * h23 * h11)) / DD)
      power = alnorm(beta / sebeta - 1.959964, upper) + alnorm(-beta / sebeta - 1.959964, upper)
      }
  return(power)
}

