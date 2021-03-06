#' Power simulations for cluster-randomized trials: Stepped Wedge Design, Binary Outcome
#'
#' This function uses a modified Cox method to determine power for stepped
#' wedge cluster-randomized controlled trials. Users can modify a variety of
#' parameters to suit their desired experimental situation.
#'
#' The stepped wedge trial design is a type of cross-over
#' design in which clusters change treatments in waves. Initially all the
#' clusters receive the same standard treatment, and at the end of the trial all
#' of the clusters will be receiving the treatment of interest. More than one
#' cluster can change treatments in a wave, but the order in which clusters
#' change treatments is randomly determined. The outcome of interest is assessed
#' in each cluster during each wave.
#'
#' Users must specify the number of subjects per cluster, number of clusters,
#' the number of time steps, the baseline effect, the expected treatment effect,
#' expected absolute difference between treatment
#' arms, ICC, and time effect.
#'
#' @param nsubjects Number of subjects per cluster; accepts a scalar. Equal cluster sizes
#' are assumed (required).
#' 
#' @param nclusters Number of clusters; accepts non-negative integer scalar (required).
#' 
#' @param alpha Significance level (default=0.05).
#' 
#' @param steps Number of crossover steps; Accepts positive scalar indicating the total
#' number of steps, NOT including the baseline (required).
#' 
#' @param timeEffect Expected time effect over the entire study period (assumed to be linear
#' across time steps); accepts numeric (required). Default = 0 (no time effects).
#' 
#' @param ICC Intracluster correlation coefficient as defined by Hussey and Hughes (2007) 
#' for participants at first time step; accepts numeric (required). 
#' 
#' @param p1 Estimated treatment effect; accepts numeric (required).
#' 
#' @param p0 Estimated baseline effect; accepts numeric (required).
#' 
#' @param tol Machine tolerance. Accepts numeric. Default is 1e-5.
#' 
#' @param GQ Number of quadriture points used in Gaussian Legendre integration; accepts
#' a scalar. Default is 100.
#' 
#' @param quiet Suppresses the progress bar; logical. Default is FALSE.
#'
#' @return The estimated power.
#'
#' @examples
#' 
#' # Estimate power for a trial with 3 steps and 9 clusters at the 
#' # initiation of the study. Those 
#' # clusters have 14 subjects each with no time effects. 
#' # We estimated arm outcome proportions of 
#' # 0.2 (pre-treatment) and 0.31 (post-treatment) and intracluster 
#' # correlation coefficient (ICC) of 0.05. 
#' # The resulting power should be 0.7992842.
#' 
#' \dontrun{
#' sw.bin <- cpa.sw.binary(nclusters = 9,
#'   steps = 3,
#'   nsubjects = 14,
#'   timeEffect = 0,
#'   ICC = 0.05,
#'   p1 = 0.31,
#'   p0 = 0.2)
#' }
#'
#' @author Alexandria C. Sakrejda (\email{acbro0@@umass.edu})
#' @author Ken Kleinman (\email{ken.kleinman@@gmail.com})
#'
#' @references Zhou X, Liao X, Kunz LM, Normand ST, Wang M, Spiegelman D. A maximum
#' likelihood approach to power calculations for stepped wedge designs of binary
#' outcomes. Biostatistics. 2020 Jan 1;21(1):102-121. doi: 10.1093/biostatistics/kxy031
#' @references Hussey, MA AND Hughes, JP. (2007). Design and analysis of stepped wedge 
#' cluster randomized trials. Contemporary Clinical Trials 28, 182–191.
#' @note Much of the FORTRAN code for this package was kindly provided by Dr. Zhou.
#' 
#' @export

cpa.sw.binary <- function(nclusters = NA,
                          steps = NA,
                          nsubjects = NA,
                          alpha = 0.05,
                          timeEffect = 0,
                          ICC = NA,
                          p0 = NA,
                          p1 = NA,
                          tol = 1e-5,
                          GQ = 100,
                          quiet = FALSE) {
  ###### Define some FORTRAN-calling functions  ########
  
  steps = steps + 1 # to include the baseline
  
  if (p1 < 0 | p0 < 0) {
    stop("Proportions (p1 & p2) cannot be negative.")
  }
  
  beta <- p1 - p0
  
  syminverse <- function(invVar = invVar,
                         Var = Var,
                         steps = steps) {
    o = .Fortran(
      "syminverse",
      a = as.numeric(invVar),
      c = as.numeric(Var),
      n = as.integer(steps + 2)
    )
    return(o)
  }
  
  vectorsquare <-
    function(derlikelihood = as.vector(derlikelihood),
             steps = steps) {
      derlen = steps + 2
      derlikelihood2 = matrix(0, nrow = (steps + 2), ncol = (steps + 2))
      o = .Fortran(
        "vectorsquare",
        a = as.numeric(derlikelihood),
        n = as.integer(derlen),
        c = as.numeric(derlikelihood2)
      )
      return(o)
    }
  
  der_likelihood_time <- function(mu = p0,
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
                                  wts = wts) {
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
    o = .Fortran(
      "der_likelihood_time",
      mu = as.numeric(mu),
      beta = as.numeric(beta),
      gammaobj = as.numeric(gammaobj),
      tau2 = as.numeric(tau2),
      z0 = as.integer(z0),
      z1 = as.integer(z1),
      XX = as.integer(XX),
      JJ = as.integer(steps),
      KK = as.integer(nsubjects),
      a = as.numeric(a),
      b = as.numeric(b),
      mincomp = as.integer(mincomp),
      maxcomp = as.integer(maxcomp),
      GQ = as.integer(GQ),
      GQX = as.numeric(t),
      GQW = as.numeric(wts),
      derlikelihood = as.numeric(derlikelihood),
      prob = as.numeric(prob)
    )
    return(o)
  }
  
  legendre_handle <- function(order = GQ,
                              a = a,
                              b = b) {
    x <- rep(0, times = order)
    w <- rep(0, times = order)
    o = .Fortran(
      "legendre_handle",
      order = as.integer(order),
      a = as.numeric(a),
      b = as.numeric(b),
      x = as.numeric(x),
      w = as.numeric(w)
    )
    return(o)
  }
  
  computeparameter <- function(JJ = steps,
                               mu = p0,
                               beta = beta,
                               p0 = p0,
                               p11 = p11,
                               rho0 = ICC) {
    tau2 <- 0.0
    gamma <- rep(0, times = JJ)
    o = .Fortran(
      "computeparameter",
      JJ = as.integer(JJ),
      mu = as.numeric(mu),
      beta = as.numeric(beta),
      gamma = as.numeric(gamma),
      tau2 = as.numeric(tau2),
      p0 = as.numeric(p0),
      p11 = as.numeric(p11),
      rho0 = as.numeric(rho0)
    )
    return(o)
  }
  
  LinearPower_notime_subroutine <-
    function(mu = p0,
             beta = beta,
             tau2 = tau2,
             II = II,
             JJ = JJ,
             KK = KK,
             a = a,
             b = b,
             GQ = GQ,
             GQX = GQX,
             GQW = GQW) {
      power = 0.0
      o = .Fortran(
        "LinearPower_notime_subroutine",
        mu = as.numeric(mu),
        beta = as.numeric(beta),
        tau2 = as.numeric(tau2),
        II = as.integer(II),
        JJ = as.integer(JJ),
        KK = as.integer(KK),
        a = as.numeric(a),
        b = as.numeric(b),
        GQ = as.integer(GQ),
        GQX = as.numeric(GQX),
        GQW = as.numeric(GQW),
        power = as.numeric(power)
      )
      return(o)
    }
  
  ###############################################
  ######## main function code ###################
  ###############################################
  
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
  if (!is.integer(GQ) ||
      GQ < 1 ||
      length(GQ) > 1 ||
      is.na(GQ)) {
    errorCondition(message = "GQ must be a positive scalar.")
  }
  if (is.na(timeEffect) ||
      is.na(ICC) || is.na(beta) || is.na(p0) || is.na(tol)) {
    errorCondition("User must provide a value for timeEffect, ICC, beta, p0, and tol. See documentation for details.")
  }
  if (!is.logical(quiet)) {
    errorCondition("Provide a logical for quiet.")
  }
  
  #### cpa.sw.binary R code ##################
  
  #  Update progress information
  if (quiet == FALSE) {
    message(paste0('Begin calculations :: Start Time: ', Sys.time()))
  }
  
  p0 <- rep(p0, times = steps)
  p11 <-  p0[1] + beta
  p0stepchange <- timeEffect / (steps - 1)
  for (i in 2:steps) {
    p0[i] <-  p0[i - 1] + p0stepchange
  }
  parholder <- computeparameter(
    JJ = steps,
    mu = p0,
    beta = beta,
    p0 = p0,
    p11 = p11,
    rho0 = ICC
  )
  tau2 <- parholder$tau2
  gammaobj <- parholder$gamma
  # mincomp and maxcomp are steps+2 vectors of 0 and 1's,
  # representing the weights of gammaobj(1),...,gammaobj(steps), p0, beta.
  comp <- rep(0, times = (steps + 2))
  maxcomp <- comp
  mincomp <- comp
  
  ## Set start.time for progress iterator & initialize progress bar
  if (quiet == FALSE) {
    start.time = Sys.time()
    prog.bar =  progress::progress_bar$new(
      format = "(:spin) [:bar] :percent eta :eta",
      total = (nsubjects + 1),
      clear = FALSE,
      width = 100
    )
    prog.bar$tick(0)
  }
  
  if (timeEffect > tol || timeEffect < -tol) {
    a <-  100
    b <- -100
    
    for (i in 1:steps) {
      temp = p0[i] + gammaobj[i]
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
      temp = p0[i] + beta + gammaobj[i]
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
    quadholder <- legendre_handle(order = GQ, a = a, b = b)
    t <- quadholder$x
    wts <- quadholder$w
    
    DD <-
      nclusters / (steps - 1)   # nclusters is a multiple of (steps-1)
    # assign intervention
    interventionX <-
      matrix(data = 0,
             nrow = (steps - 1),
             ncol = (steps))
    
    for (i in 1:(steps - 1)) {
      for (j in (i + 1):steps) {
        interventionX[i, j] <- 1
      }
    }
    invVar <- matrix(0, nrow = (steps + 2), ncol = (steps + 2))
    
    h <- 0
    
    for (i in 1:(steps - 1)) {
      z0 <- rep(0, times = steps)
      finish <- 0
      h <- h
      while (finish < 1) {
        h <- h + 1
        XX <- interventionX[i,]
        z1 <- nsubjects - z0
        
        Dholder <- der_likelihood_time(
          mu = p0,
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
          wts = wts
        )
        prob <- Dholder$prob
        derlikelihood <- Dholder$derlikelihood
        VecHolder <-
          vectorsquare(derlikelihood = derlikelihood, steps = steps)
        derlikelihood2 <- VecHolder$c
        invVar <- invVar + derlikelihood2 * prob
        
        # Iterate progress bar
        totaliter <- ((nsubjects + 1) ^ steps * (steps - 1))
        if (quiet == FALSE) {
          hspecial <- (h / 100000)
          if (hspecial == as.integer(hspecial) ||
              hspecial == totaliter ||
              h == 10) {
            prog.bar$update(h / totaliter)
            Sys.sleep(1 / 100)
          }
        }
        
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
        
        if (z0[steps] > nsubjects) {
          finish = 1
        }
      }
    }
    
    mat <- matrix(0, nrow = (steps + 2), ncol = (steps + 2))
    mat <- syminverse(invVar = invVar,
                      Var = mat,
                      steps = steps)
    Var <- matrix(
      mat$c,
      nrow = (steps + 2),
      ncol = (steps + 2),
      byrow = FALSE
    )
    sebeta <- sqrt(Var[2, 2] / DD)
    z_half_alpha <- qnorm(1 - alpha / 2) 
    power <- pnorm(beta / sebeta - z_half_alpha, lower.tail = TRUE) +
      pnorm(-beta / sebeta - z_half_alpha, lower.tail = TRUE)
    
  } else {
    if (beta > 0) {
      a = -p0
      b = 1 - p0 - beta
    } else {
      a = -p0 - beta
      b = 1 - p0
    }
    quadholder <- legendre_handle(order = GQ,
                                  a = a,
                                  b = b)
    t <- quadholder$x
    wts <- quadholder$w
    
    Linpower <-
      LinearPower_notime_subroutine(
        mu = p0,
        beta = beta,
        tau2 = tau2,
        II = nclusters,
        JJ = steps,
        KK = nsubjects,
        a = a,
        b = b,
        GQ = GQ,
        GQX = t,
        GQW = wts
      )
    power <- Linpower$power
  }
  if (quiet == FALSE) {
    ## show elapsed time
    total.est = as.numeric(difftime(Sys.time(), start.time, units = 'secs'))
    hr.est = total.est %/% 3600
    min.est = total.est %/% 60
    sec.est = round(total.est %% 60, 0)
    message(
      paste0(
        "Time Completed: ",
        Sys.time(),
        "\nTotal Runtime: ",
        hr.est,
        'Hr:',
        min.est,
        'Min:',
        sec.est,
        'Sec'
      )
    )
  }
  names(power) <- "power"
  return(power)
}
