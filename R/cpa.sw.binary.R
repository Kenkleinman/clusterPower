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
#' @return A list with the following components
#' \itemize{
#'   \item Character string indicating total number of simulations and simulation type
#'   \item 
#' }
#' 
#' @examples 
#' \dontrun{
#' }
#'
#' @author Alexandria C. Sakrejda (\email{acbro0@@umass.edu})
#' @author Ken Kleinman (\email{ken.kleinman@@gmail.com})
#'
#' 
#@export
# cpa.sw.binary(nclusters (I), ntimes(J), nsubjects(K), d(delta, p0totalchange), ICC(rho0), beta, mu)


cpa.sw.binary <- function(nclusters, ntimes, nsubjects, d, ICC, beta, mu, 
                          tol = 1e-5, GQ = 100){
  p0 <- vector(mode = "numeric", length = ntimes)
  gamma <- vector(mode = "numeric", length = ntimes)
  p0[1] <-  mu + beta
  p0stepchange <- d / (ntimes - 1)
  tau2 = ICC / (1 - ICC) * mu * (1 - mu)
  gamma[1] <- 0.0
  for (i in 2:ntimes){
    p0[i] = p0[i - 1] + p0stepchange
    gamma[i] <- p0[i] - mu
  }
  
  if (p0totalchange > tol || p0totalchange < -tol){
  a = 100 
  b = -100
  for(i in 1:ntimes){
    temp = mu + gamma[i]
     if (temp < a){
       a = temp
       mincomp = 0
       mincomp[ntimes+1] = 1
       mincomp[i] = 1
     }
    if (temp > b){
      b = temp
      maxcomp = 0
      maxcomp[ntimes+1] = 1
      maxcomp[i] = 1
    }
    temp = mu + beta + gamma[i]
    if (temp < a){
      a = temp
      mincomp = 0
      mincomp[ntimes + 1] = 1
      mincomp[ntimes + 2] = 1
      mincomp[i] = 1
      }
    if (temp > b){
      b = temp
      maxcomp = 0
      maxcomp[ntimes + 1] = 1
      maxcomp[ntimes + 2] = 1
      maxcomp[i] = 1
    }
  }
  a = -a
  b = 1-b
  
  GQholder <- statmod::gauss.quad(GQ, kind = "legendre", alpha = 0, beta = 0)
  GQX <- GQholder[[1]]
  GQW <- GQholder[[2]]
  #power = LinearPower_time(mu, beta, gamma, tau2, II, JJ, KK, a, b, mincomp, maxcomp, GQ, GQX, GQW)
  DD <-  nclusters/(ntimes-1)   # II is a multiple of (JJ-1)
  # assign intervention
  interventionX <- matrix(data = 1, nrow = (ntimes-1), ncol = (ntimes-1))

  
  # compute the possible combinations of (z0,z1)
  # II : number of clusters
  # JJ : number of steps
  # KK : number of subjects per cluster per step
  invVar = 0.0
  for (i in 1:(ntimes - 1)){
    z0 = 0
    finish = 0
    do while (finish<1)
      z1 = nsubjects - z0
    
   # call der_likelihood_time(mu,beta,gamma,tau2, z0, z1, X(i,:), JJ, KK, a, b, &
   #                            mincomp, maxcomp, GQ, GQX, GQW, derlikelihood, prob)
    
    likelihoodf_denom = 0.0
    likelihoodf_numer = 0.0
    likelihoodf_denomb2 = 0.0
    derlikelihood_mu = 0.0
    derlikelihood_beta = 0.0
    derlikelihood_tau2 = 0.0
    derlikelihood_gamma = 0.0
    
    prob = 0.0
    
    for (i in 1:GQ){
      x = GQX[i]
      exx = dexp(-0.5 * x * x / tau2)
      ff = 1.0
      ffprob = 1.0
      ff_mu = 0.0
      ff_beta = 0.0
    
      for (j in 1:ntimes){
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
        if (z0[j] < z1[j]){
          ffprob = ffprob * ff1^(z1[j]-z0[j])
          for (k in 0:(z0[j]-1)){
            ffprob = ffprob * (nsubjects - k) / (z0[j] - k) * ff01
          }
        } else {      
          ffprob = ffprob * ff0^(z0[j] - z1[j])
          for (k in 0:(z1[j]-1)){
            ffprob = ffprob * (nsubjects - k) / (z1[j] - k) * ff01
          }
        }
      }
      
      
      
      
      
      
    prob = prob + GQW(i)*ffprob*exx
    likelihoodf_denom = likelihoodf_denom + GQW(i)*exx
    likelihoodf_numer = likelihoodf_numer + GQW(i)*ff*exx
    likelihoodf_denomb2 = likelihoodf_denomb2 + GQW(i)*x*x*exx
    
    derlikelihood_mu = derlikelihood_mu + GQW(i)*ff*ff_mu*exx
    derlikelihood_beta = derlikelihood_beta + GQW(i)*ff*ff_beta*exx
    derlikelihood_gamma = derlikelihood_gamma + GQW(i)*ff*ff_gamma*exx
    derlikelihood_tau2= derlikelihood_tau2 + GQW(i)*ff*x*x*exx
    }
    enddo
    
    ! calculate f(a)exp(-0.5*a*a)
    eaa = dexp(-0.5d0*a*a/tau2)
    ff = 1.0d0
    do j=1,JJ
    ff1 = mu+beta*XX(j)+gamma(j)+a
    ff0 = 1-ff1
    ff = ff*(ff0**z0(j))*(ff1**z1(j))
    end do
    faeaa = ff*eaa
    ! calculate f(b)exp(-0.5*b*b)
    ebb = dexp(-0.5d0*b*b/tau2)
    ff = 1.0d0
    do j=1,JJ
    ff1 = mu+beta*XX(j)+gamma(j)+b
    ff0 = 1-ff1
    ff = ff*(ff0**z0(j))*(ff1**z1(j))
    end do
    fbebb = ff*ebb
    
    ! calculate derlikelihood_mu
    derlikelihood_mu = derlikelihood_mu + faeaa*dble(mincomp(JJ+1)) &
      - fbebb*dble(maxcomp(JJ+1))
    derlikelihood_mu = derlikelihood_mu / likelihoodf_numer - &
      (eaa*dble(mincomp(JJ+1)) - ebb*dble(maxcomp(JJ+1))) / likelihoodf_denom
    ! calculate derlikelihood_beta
    derlikelihood_beta = derlikelihood_beta + faeaa*dble(mincomp(JJ+2)) &
      - fbebb*dble(maxcomp(JJ+2))
    derlikelihood_beta = derlikelihood_beta / likelihoodf_numer - &
      (eaa*dble(mincomp(JJ+2)) - ebb*dble(maxcomp(JJ+2))) / likelihoodf_denom
    ! calculate derlikelihood_gamma
    do j=2,JJ
    k = j-1
    derlikelihood_gamma(k) = derlikelihood_gamma(k) + &
      faeaa*dble(mincomp(j)) - fbebb*dble(maxcomp(j))
    derlikelihood_gamma(k) = derlikelihood_gamma(k)/likelihoodf_numer &
      -(eaa*dble(mincomp(j))-ebb*dble(maxcomp(j)))/likelihoodf_denom
    end do
    
    ! calculate derlikelihood_tau2
    derlikelihood_tau2 = 0.5d0*(derlikelihood_tau2/likelihoodf_numer- &
                                  likelihoodf_denomb2/likelihoodf_denom)/tau2/tau2
    prob = prob/likelihoodf_denom
    
    derlikelihood(1) = derlikelihood_mu
    derlikelihood(2) = derlikelihood_beta
    derlikelihood(3:(JJ+1)) = derlikelihood_gamma
    derlikelihood(JJ+2) = derlikelihood_tau2
    
    
    call vectorsquare(derlikelihood, JJ+2, derlikelihood2)
    invVar = invVar + derlikelihood2 * prob
    finish = updatez(z0, JJ, KK)
  }
  
  end do
  end do
  
  
  call syminverse(invVar,Var,JJ+2)
  sebeta = sqrt(Var(2,2)/DD)
  power = alnorm(beta/sebeta-1.959964d0,upper) + alnorm(-beta/sebeta-1.959964d0,upper)
  
  
  
  
  } else

  
  
}

