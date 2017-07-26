#' Power calculations for simple cluster randomized trials, binary outcome
#'
#' Compute the power of a simple cluster randomized trial with a binary outcome,
#' or determine parameters to obtain a target power.
#'
#' @section Authors:
#' Jonathan Moyer (\email{jon.moyer@@gmail.com})
#'
#' @param alpha The level of significance of the test, the probability of a
#'   Type I error.
#' @param power The power of the test, 1 minus the probability of a Type II
#'   error.
#' @param m The number of clusters per condition. It must be greater than 1.
#' @param n The mean of the cluster sizes.
#' @param cv The coefficient of variation of the cluster sizes. When \code{cv} = 0,
#'   the clusters all have the same size.
#' @param p1 The expected proportion in the treatment group.
#' @param p2 The proportion in the control group.
#' @param icc The intraclass correlation.
#' @param pooled Logical indicating if pooled standard error should be used.
#' @param p1inc Logical indicating if p1 is expected to be greater than p2.
#' @param tol Numerical tolerance used in root finding. The default provides
#'   at least four significant digits.
#' @return The computed argument.
#' @export

crtpwr.2prop <- function(alpha = 0.05, power = 0.80,
                         m = NA, n = NA, cv = 0,
                         p1 = NA, p2 = NA,
                         icc = NA, pooled = FALSE,
                         p1inc = TRUE,
                         tol = .Machine$double.eps^0.25){
  
  if(!is.na(m) && m <= 1) {
    stop("'m' must be greater than 1.")
  }
  
  needlist <- list(alpha, power, m, n, cv, p1, p2, icc)
  neednames <- c("alpha", "power", "m", "n", "cv", "p1", "p2", "icc")
  needind <- which(unlist(lapply(needlist, is.na))) # find NA index
  
  if (length(needind) != 1) {
    stop("Exactly one of 'alpha', 'power', 'm', 'n', 'cv', 'p1', 'p2', or 'icc' must be NA.")
  }
  
  target <- neednames[needind]
  
  p.body <- quote({
    DEFF <- 1 + ((cv^2 + 1)*n - 1)*icc
    if (pooled) {
      p <- (p1+p2)/2
      sdd <- sqrt(p*(1 - p)*2*DEFF/(m*n))
    } else {
      sdd <- sqrt((p1*(1-p1) + p2*(1-p2))*DEFF/(m*n))
    }
    zcrit <- qnorm(alpha/2, lower.tail = FALSE)
    pnorm(abs(p1 - p2)/sdd - zcrit, lower.tail = TRUE)# +
    #pnorm(-zcrit - (p1 - p2)/sdd, lower.tail = TRUE)
  })
  
  # calculate alpha
  if (is.na(alpha)) {
    alpha <- uniroot(function(alpha) eval(p.body) - power,
                     interval = c(1e-10, 1 - 1e-10),
                     tol = tol, extendInt = "yes")$root
  }
  
  # calculate power
  if (is.na(power)) {
    power <- eval(p.body)
  }
  
  # calculate m
  if (is.na(m)) {
    m <- uniroot(function(m) eval(p.body) - power,
                 interval = c(2 + 1e-10, 1e+07),
                 tol = tol)$root
  }
  
  # calculate p1
  if (is.na(p1)) {
    if(p1inc){
      p1 <- uniroot(function(p1) eval(p.body) - power,
                    interval = c(p2 + 1e-7, 1 - 1e-7),
                    tol = tol, extendInt = "yes")$root
    } else {
      p1 <- uniroot(function(p1) eval(p.body) - power,
                    interval = c(1e-7, p2 - 1e-7),
                    tol = tol, extendInt = "yes")$root
    }
  }
  
  # calculate p2
  if (is.na(p2)) {
    if(p1inc){
      p2 <- uniroot(function(p2) eval(p.body) - power,
                    interval = c(1e-7, p1 - 1e-7),
                    tol = tol, extendInt = "yes")$root
      
    } else {
      p2 <- uniroot(function(p2) eval(p.body) - power,
                    interval = c(p1 + 1e-7, 1 - 1e-7),
                    tol = tol, extendInt = "yes")$root
    }
  }
  
  # calculate n
  if (is.na(n)) {
    n <- uniroot(function(n) eval(p.body) - power,
                 interval = c(2 + 1e-10, 1e+07),
                 tol = tol, extendInt = "upX")$root
  }
  
  # calculate cv
  if (is.na(cv)) {
    
    cv <- uniroot(function(cv) eval(p.body) - power,
                  interval = c(1e-7, 1e+07),
                  tol = tol, extendInt = "downX")$root
  }
  
  # calculate icc
  if (is.na(icc)){
    icc <- uniroot(function(icc) eval(p.body) - power,
                   interval = c(1e-07, 1 - 1e-7),
                   tol = tol, extendInt = "downX")$root
  }
  
  structure(get(target), names = target)
  
  # method <- paste("Clustered two-sample proportion power calculation: ", target, sep = "")
  # note <- "'m' is the number of clusters in each group and 'n' is the number of individuals in each cluster."
  # structure(list(m = m, n = n, cv = cv,
  #                p1 = p1, p1dec = p1dec, p1inc = p1inc,
  #                p2 = p2, p2dec = p2dec, p2inc = p2inc,
  #                icc = icc, alpha = alpha, power = power,
  #                note = note, method = method),
  #           class = "power.htest")
  
}
