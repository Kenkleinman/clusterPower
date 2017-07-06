#' Power calculations for simple cluster randomized trials, continuous outcome
#'
#' Compute the power of a simple cluster randomized trial with a continuous outcome,
#' or determine parameters to obtain a target power.
#'
#' Exactly one of \code{alpha}, \code{power}, \code{d}, \code{ICC}, \code{m},
#'   \code{n}, and \code{cv} must be passed as \code{NULL}. Note that
#'   \code{alpha}, \code{power}, and \code{cv} have non-\code{NULL}
#'   defaults, so if those are the parameters of interest they must be
#'   explicitly passed as \code{NULL}.
#'
#' @section Note:
#'   'uniroot' is used to solve power equation for unknowns, so you may see
#'   errors from it, notably about inability to bracket the root when
#'   invalid arguments are given.
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
#' @param nsd The standard deviation of the cluster sizes.
#' @param cv The coefficient of variation of the cluster sizes. When \code{cv} = 0,
#'   the clusters all have the same size.
#' @param d The difference in condition means.
#' @param icc The intraclass correlation.
#' @param varw The within-cluster variation.
#' @param varb The between-cluster variation.
#' @param tol Numerical tolerance used in root finding. The default provides
#'   at least four significant digits.
#' @return The computed argument.

crtpower_2mean <- function(alpha = 0.05, power = 0.80, m = NULL,
                            n = NULL, nsd = 0, cv = 0,
                            d = NULL, icc = NULL,
                            varw = NULL, varb = NULL,
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

  # if cv = 0, set nsd = 0
  # if nsd = 0, set cv = 0
  if(!is.null(cv) && cv == 0) {
    nsd <- 0
  } else if (!is.null(nsd) && nsd == 0) {
    cv <- 0
  }

  # checking set of n, nsd, and cv
  nlist <- list(n, nsd, cv)
  nnames <- c("n","nsd","cv")
  nind <- which(unlist(lapply(nlist, is.null))) # find null index
  # check to make sure that both n and cv are not NULL
  # if only one of n and cv is specified it will be assumed that the user wants other in the pair
  if (is.null(n) & is.null(cv)){
    stop("At least one of 'n' and 'cv' must not be NULL.")
  }

  # if none of n, nsd, and cv are null, double check to make sure their equality holds
  if (length(nind) == 0) {
    if(cv != nsd/n) {
      message("Recalcuating cv so that cv = nsd/n")
      cv <- nsd/n
    }
  }

  # if one of n, nsd, and cv is null, calculate it
  # skip this if nsd and cv are both 0
  if (length(nind) == 1 && !is.null(nsd) && nsd != 0) {
    assign(nnames[nind], calc_n(nind, n, nsd, cv))
    # # if num_null also is 1, then return the value just calculated
    # if (all_null == 1) {
    #   return(get(nnames[nind]))
    # }
  }

  # checking set of icc, varw, varb
  icclist <- list(icc, varw, varb)
  iccnames <- c("icc","varw","varb")
  iccind <- which(unlist(lapply(icclist, is.null))) # find null index
  # # check to make sure that at least two of icc, varw, and varb is not null
  # # if only one of n and cv is specified it will be assumed that the user wants other in the pair
  # #   the user can calculate nsd from these two if so desired
  # if (length(iccind) > 1){
  #   stop("At least two of 'icc', 'varw', and 'varb' must not be NULL.")
  # }

  # if none of icc, varw, and varb are null, double check to make sure their equality holds
  # if (length(iccind) == 0) {
  #   if(icc != varb/(varw + varb)) {
  #     message("Recalcuating icc so that icc = varb/(varw + varb)")
  #     icc <- varb/(varw + varb)
  #   }
  # }

  # if one of icc, varw, and varb is null, calculate it
  if (length(iccind) == 1) {
    assign(iccnames[iccind], calc_icc(iccind, icc, varw, varb))
    # # if num_null also is 1, then return the value just calculated
    # if (all_null == 1) {
    #   return(get(iccnames[iccind]))
    # }
  }

  # list of needed inputs
  needlist <- list(alpha, power, m, d)
  needind <- which(unlist(lapply(needlist, is.null)))
  # check to see that exactly one needed param is null

  if (length(needind) > 1) {
    neederror = "More than one of 'alpha', 'power', 'm', and 'd' are NULL. Only up to one of the four may be NULL."
  } else {
    neederror = NULL
  }

  if (length(nind) == 3) {
    nerror <- "All three of 'n', 'nsd', and 'cv' are NULL. Only up to two of the three may be NULL."
  } else {
    nerror <- NULL
  }

  if (length(iccind) == 3) {
    iccerror <- "All three of 'icc', 'varw', and 'varb' are NULL. Only up to two of the three may be NULL."
  } else {
    iccerror <- NULL
  }

  if(length(nind) == 2 & length(iccind) == 2){
    botherror <- "Two of the set (n, nsd, cv) and two of the set (icc, varw, varb) are NULL. At least one set must have at most one NULL."
  } else {
    botherror <- NULL
  }

  errorlist <- list(neederror, nerror, iccerror, botherror)
  errorind <- which(!unlist(lapply(errorlist, is.null)))
  if (length(errorind) > 0){
    errormsg <- paste(unlist(errorlist),collapse = "\n  ")
    stop(errormsg)
  }

  # create call to evaluate power

  # if two of n, nsd, and cv are NULL:
  #   if nsd is not null, it will be assumed the user wants n
  #   if n is not null, it will be assumed the user wants cv
  #   if cv is not null, it will be assumed the user wants n
  # if two of icc, varw, and varb are NULL:
  #   if varb is not null, it will be assumed the user wants varw
  #   if varw is not null, it will be assumed the user wants icc
  #   if icc is not null, it will be assumed the user wants varw

  # DEFF gets updated so define inside p.body
  if (length(nind) == 2 | length(iccind) == 2) {

    p.body <- quote({

      DEFF <- getDEFF(n, nsd, cv, varw, varb, icc)

      qu <- qt(alpha/2, 2*(m - 1), lower.tail = FALSE)

      ncp <- sqrt(m*n/(2*DEFF)) * d/sqrt(varw)

      pt(qu, 2*(m - 1), ncp, lower.tail = FALSE) +
        pt(-qu, 2*(m - 1), ncp, lower.tail = TRUE)
    })

    # else if one of n, cv, icc, and varw are missing
    # DEFF gets updated so define inside p.body
  } else if (is.null(n) | is.null(cv) | is.null(icc) | is.null(varw)) {

    p.body <- quote({

      DEFF <- 1 + (((cv^2 + 1)*n) - 1)*icc

      qu <- qt(alpha/2, 2*(m - 1), lower.tail = FALSE)

      ncp <- sqrt(m*n/(2*DEFF)) * d/sqrt(varw)

      pt(qu, 2*(m - 1), ncp, lower.tail = FALSE) +
        pt(-qu, 2*(m - 1), ncp, lower.tail = TRUE)
    })

    # else if one of alpha, power, m, and d is NULL
    # DEFF doesn't get updated so define outside of p.body
  } else {

    DEFF <- 1 + (((cv^2 + 1)*n) - 1)*icc

    p.body <- quote({

      qu <- qt(alpha/2, 2*(m - 1), lower.tail = FALSE)

      ncp <- sqrt(m*n/(2*DEFF)) * d/sqrt(varb*(1-icc)/icc)

      pt(qu, 2*(m - 1), ncp, lower.tail = FALSE) +
        pt(-qu, 2*(m - 1), ncp, lower.tail = TRUE)
    })
  }

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

  # calculate d
  if (is.null(d)) {
    d <- uniroot(function(d) eval(p.body) - power,
                 interval = c(1e-07, 1e+07),
                 tol = tol, extendInt = "upX")$root
    #return(d)
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

    if(!is.null(n) & !is.null(nsd)) {
      cv <- nsd/n
    } else {
      cv <- uniroot(function(cv) eval(p.body) - power,
                    interval = c(1e-7, 1e+07),
                    tol = tol, extendInt = "downX")$root
    }

    #return(cv)
  }

  # calculate nsd
  # putting it after n, cv because those two would have been determined by now
  if (is.null(nsd)) {
    nsd <- calc_n(2, n, nsd, cv)
  }

  # calculate varw
  if (is.null(varw)) {
    varw <- uniroot(function(varw) eval(p.body) - power,
                    interval = c(1e-07, 1e+07),
                    tol = tol, extendInt = "downX")$root
    #return(varw)
  }

  # calculate icc
  # if icc is null but varw, varb not null
  if (is.null(icc)){

    if(!is.null(varw) & !is.null(varb)) {
      icc <- varb/(varw + varb)
    } else {
      icc <- uniroot(function(icc) eval(p.body) - power,
                     interval = c(1e-07, 1e+07),
                     tol = tol, extendInt = "downX")$root
    }
  }

  # calculate varb
  # icc and varw should have been determined by now
  if (is.null(varb)) {
    varb <- calc_icc(3, icc, varw, varb)
  }

  method <- "Clustered two-sample t-test power calculation"
  note <- "'m' is the number of clusters in each group and 'n' is the number of individuals in each cluster."
  structure(list(m = m, n = n, nsd = nsd, cv = cv, alpha = alpha, power = power, d = d,
                 icc = icc, varw = varw, varb = varb, note = note, method = method),
            class = "power.htest")

}
