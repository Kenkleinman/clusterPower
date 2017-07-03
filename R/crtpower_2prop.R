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
#' @param n The mean of the cluster sizes, or a vector of cluster sizes with
#'   length equal to twice \code{m}.
#' @param cv The coefficient of variation of the cluster sizes. When \code{cv} = 0,
#'   the clusters all have the same size.
#' @param p1 The expected proportion in the treatment group.
#' @param p2 The proportion in the control group.
#' @param icc The intraclass correlation.
#' @param pooled Logical indicating if pooled standard error should be used.
#' @param tol Numerical tolerance used in root finding. The default provides
#'   at least four significant digits.
#' @return The computed argument.

crtpower_2prop <- function(alpha = 0.05, power = 0.80,
                           m = NULL, n = NULL, cv = NULL,
                           p1 = NULL, p2 = NULL,
                           icc = NULL, pooled = FALSE,
                           tol = .Machine$double.eps^0.25){

  if(!is.null(m) && m <= 1) {
    stop("'m' must be greater than 1.")
  }

  # if n is a vector of cluster sizes, calculate mean and cv of cluster sizes
  if (length(n) > 1) {
    if (!is.null(m) && length(n) != 2*m) { # use && to evaluate !is.null(m) first
      stop("length(n) is not equal to 2*m. Enter a vector of the correct length,
           or enter one number for mean cluster size.")
    }
    nsd <- sd(n) # find sd of cluster sizes
    n <- mean(n) # find mean cluster size
    cv <- nsd/n  # find coeffient of variation
  }

  needlist <- list(alpha, power, m, n, cv, p1, p2, icc)
  neednames <- c("alpha", "power", "m", "n", "cv", "p1", "p2", "icc")
  needind <- which(unlist(lapply(needlist, is.null))) # find null index

  if (length(needind) != 1) {
    stop("Exactly one of 'alpha', 'power', 'm', 'n', 'cv', 'p1', 'p2', or 'icc' must be NULL.")
  }

  # p1dec corresponds to the case when p1 < p2
  # p1inc corresponds to the case when p2 > p1
  # these values will be calculated if p1 is null
  p1dec <- NULL
  p1inc <- NULL

  # p2dec corresponds to the case when p1 < p2
  # p2inc corresponds to the case when p2 > p1
  # these values will be calculated if p2 is null
  p2dec <- NULL
  p2inc <- NULL

  p.body <- quote({
    DEFF <- 1 + ((cv^2 + 1)*n - 1)*icc
    if (pooled) {
      p <- (p1+p2)/2
      sdd <- sqrt(p*(1 - p)*2*DEFF/(m*n))
    } else {
      sdd <- sqrt((p1*(1-p1) + p2*(1-p2))*DEFF/(m*n))
    }
    zcrit <- qnorm(alpha/2, lower.tail = FALSE)
    pnorm(zcrit - abs(p1 - p2)/sdd, lower.tail = FALSE)# +
      #pnorm(-zcrit - (p1 - p2)/sdd, lower.tail = TRUE)
  })


  # calculate alpha
  if (is.null(alpha)) {
    alpha <- uniroot(function(alpha) eval(p.body) - power,
                     interval = c(1e-10, 1 - 1e-10),
                     tol = tol, extendInt = "yes")$root
    #return(alpha)
  }

  # calculate power
  if (is.null(power)) {
    power <- eval(p.body)
    #return(power)
  }

  # calculate m
  if (is.null(m)) {
    m <- uniroot(function(m) eval(p.body) - power,
                 interval = c(2 + 1e-10, 1e+07),
                 tol = tol, extendInt = "upX")$root
    #return(m)
  }

  # calculate p1
  if (is.null(p1)) {
    p1dec <- uniroot(function(p1) eval(p.body) - power,
                    interval = c(0.0001, p2 - 0.0001),
                    tol = tol, extendInt = "yes")$root

    p1inc <- uniroot(function(p1) eval(p.body) - power,
                     interval = c(p2 + 0.0001, 0.9999),
                     tol = tol, extendInt = "yes")$root
  }

  # calculate p2
  if (is.null(p2)) {
    p2dec <- uniroot(function(p2) eval(p.body) - power,
                 interval = c(p1 + 0.0001, 0.999),
                 tol = tol, extendInt = "yes")$root

    p2inc <- uniroot(function(p2) eval(p.body) - power,
                  interval = c(0.0001, p1 - 0.0001),
                  tol = tol, extendInt = "yes")$root

  }

  # calculate n
  if (is.null(n)) {
    n <- uniroot(function(n) eval(p.body) - power,
                 interval = c(2 + 1e-10, 1e+07),
                 tol = tol, extendInt = "upX")$root
    #return(n)
  }

  # calculate cv
  if (is.null(cv)) {

    cv <- uniroot(function(cv) eval(p.body) - power,
                  interval = c(1e-7, 1e+07),
                  tol = tol, extendInt = "downX")$root
    #return(cv)
  }

  # calculate icc
  # if icc is null but varw, varb not null
  if (is.null(icc)){
    icc <- uniroot(function(icc) eval(p.body) - power,
                   interval = c(1e-07, 0.9999999),
                   tol = tol, extendInt = "downX")$root
  }

  method <- "Clustered two-sample proportion power calculation"
  note <- "'m' is the number of clusters in each group and 'n' is the number of individuals in each cluster."
  structure(list(m = m, n = n, cv = cv,
                 p1 = p1, p1dec = p1dec, p1inc = p1inc,
                 p2 = p2, p2dec = p2dec, p2inc = p2inc,
                 icc = icc, alpha = alpha, power = power,
                 note = note, method = method),
            class = "power.htest")

}
