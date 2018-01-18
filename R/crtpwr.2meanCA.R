#' Power calculations for simple cluster randomized trials, continuous outcome with covariate adjustment
#'
#' Compute the power of a simple cluster randomized trial with a continuous outcome after adjusting for covariates,
#' or determine parameters to obtain a target power.
#'
#' Exactly one of \code{alpha}, \code{power}, \code{m}, \code{n},
#'   \code{d}, \code{icc}, \code{varw}, \code{covdf}, and \code{rho_b} must be passed as \code{NA}.
#'   Note that \code{alpha} and \code{power} have non-\code{NA}
#'   defaults, so if those are the parameters of interest they must be
#'   explicitly passed as \code{NA}.
#'   
#' If \code{n} is a vector the values, \code{m} will be recalculated
#'    using the values in \code{n}.
#'
#' @section Note:
#'   This function was inspired by work from Stephane Champely (pwr.t.test) and
#'   Peter Dalgaard (power.t.test). As with those functions, 'uniroot' is used to
#'   solve power equation for unknowns, so you may see
#'   errors from it, notably about inability to bracket the root when
#'   invalid arguments are given.
#'
#' @section Authors:
#' Jonathan Moyer (\email{jon.moyer@@gmail.com})
#' Ken Kleinman (\email(ken.kleinman@@gmail.com))
#'
#' @param alpha The level of significance of the test, the probability of a
#'   Type I error.
#' @param power The power of the test, 1 minus the probability of a Type II
#'   error.
#' @param m The number of clusters per condition. It must be greater than 1.
#' @param n The mean of the cluster sizes, or a vector of cluster sizes for one arm.
#' @param d The difference in condition means.
#' @param icc The intraclass correlation.
#' @param varw The within-cluster variation.
#' @param covdf The degrees of freedom used by the cluster-level covariates. If \code{covdf} is 0,
#'   the value of \code{rho_b} must be 0 as well.
#' @param rho_b The between-cluster residual correlation between the outcome and the covariates.
#' @param tol Numerical tolerance used in root finding. The default provides
#'   at least four significant digits.
#' @return The computed argument.
#' @examples 
#' # Find the number of clusters per condition needed for a trial with alpha = 0.05, 
#' # power = 0.8, 10 observations per cluster, a difference of 1 unit,  icc = 0.1,
#' # a variance of 5 units, adjusting for 1 covariate, and a between-cluster correlation of 0.7:
#' crtpwr.2meanCA(n=10 ,d=1, icc=.1, varw=5, covdf=1, rho_b=.7)
#' # 
#' # The result, showimg m of greater than 12, suggests 13 clusters per condition should be used.
#' @references Moerbeek M, Teerenstra S. (2016) Power Analysis of Trials with Multilevel Data. 
#'   Boca Raton: CRC Press; 2016: 84-87
#' @export

crtpwr.2meanCA <- function(alpha = 0.05, power = 0.80, m = NA,
                           n = NA, d = NA, icc = NA,
                           varw = NA, covdf = NA, rho_b = NA,
                           tol = .Machine$double.eps^0.25){
  
  # if n is a vector, 
  if(length(n) > 1){
    nvec <- n
    n <- mean(nvec)
    nsd <- stats::sd(nvec)
    cv <- nsd/n
    m <- length(nvec)
  }
  
  if(!is.na(m) && m <= 1) {
    stop("'m' must be greater than 1.")
  }
  
  if(covdf == 0 && rho_b != 0){
    stop("Because 'covdf' is 0, 'rho_b' must be 0.")
  }
  
  # list of needed inputs
  needlist <- list(alpha, power, m, n, d, icc, varw, covdf, rho_b)
  neednames <- c("alpha", "power", "m", "n", "d", "icc", "varw", "covdf", "rho_b")
  needind <- which(unlist(lapply(needlist, is.na)))
  # check to see that exactly one needed param is NA
  
  if (length(needind) != 1) {
    neederror = "Exactly one of 'alpha', 'power', 'm', 'n', 'd', 'icc', 'varw', 'covdf', and 'rho_b' must be NA."
    stop(neederror)
  } 
  
  target <- neednames[needind]
  
  # evaluate power
  pwr <- quote({
    
    # design effect
    DEFF <- 1 + (n - 1)*icc*(1 - rho_b^2)
    
    # # variance inflation
    # # if nvec exists, calcuate exact relative efficiency
    # if (exists("nvec")) {
    #   if(method == "taylor"){
    #     a <- (1 - icc)/icc
    #     RE <- ((n + a)/n)*(sum((nvec/(nvec+a)))/m) # exact relative efficiency
    #     VIF <- DEFF*RE
    #   } else{
    #     VIF <- 1 + ((cv^2 + 1)*n - 1)*icc
    #   }
    # } else if(!is.na(n)){
    #   if(method == "taylor"){
    #     L <- n*icc/DEFF
    #     REt <- 1/(1 - cv^2*L*(1 - L)) # taylor approximation
    #     VIF <- DEFF*REt
    #   } else {
    #     VIF <- 1 + ((cv^2 + 1)*n - 1)*icc
    #   }
    # }
    
    tcrit <- qt(alpha/2, 2*(m - 1) - covdf, lower.tail = FALSE)
    
    ncp <- sqrt(m*n/(2*DEFF)) * abs(d)/sqrt(varw)
    
    pt(tcrit, 2*(m - 1) - covdf, ncp, lower.tail = FALSE) #+ pt(-tcrit, 2*(m - 1), ncp, lower.tail = TRUE)
  })
  
  # calculate alpha
  if (is.na(alpha)) {
    alpha <- stats::uniroot(function(alpha) eval(pwr) - power,
                            interval = c(1e-10, 1 - 1e-10),
                            tol = tol)$root
  }
  
  # calculate power
  if (is.na(power)) {
    power <- eval(pwr)
  }
  
  # calculate m
  if (is.na(m)) {
    m <- stats::uniroot(function(m) eval(pwr) - power,
                        interval = c(2 + covdf + 1e-10, 1e+7),
                        tol = tol, extendInt = "upX")$root
  }
  
  # calculate n
  if (is.na(n)) {
    n <- stats::uniroot(function(n) eval(pwr) - power,
                        interval = c(2 + 1e-10, 1e+07),
                        tol = tol, extendInt = "upX")$root
  }
  
  # calculate d
  if (is.na(d)) {
    d <- stats::uniroot(function(d) eval(pwr) - power,
                        interval = c(1e-10, 1e+07),
                        tol = tol, extendInt = "upX")$root
  }
  
  # calculate icc
  if (is.na(icc)){
    icc <- stats::uniroot(function(icc) eval(pwr) - power,
                          interval = c(1e-10, 1e+07),
                          tol = tol, extendInt = "downX")$root
  }
  
  # calculate varw
  if (is.na(varw)) {
    varw <- stats::uniroot(function(varw) eval(pwr) - power,
                           interval = c(1e-10, 1e+07),
                           tol = tol, extendInt = "downX")$root
  }
  
  # calculate covdf
  if (is.na(covdf)) {
    covdf <- stats::uniroot(function(covdf) eval(pwr) - power,
                            interval = c(1e-10, 1e+07),
                            tol = tol, extendInt = "downX")$root
  }
  
  # calculate rho_b
  if (is.na(rho_b)) {
    rho_b <- stats::uniroot(function(rho_b) eval(pwr) - power,
                            interval = c(1e-10, 1 - 1e-7),
                            tol = tol, extendInt = "upX")$root
  }
  
  structure(get(target), names = target)
  
  # method <- paste("Clustered two-sample t-test power calculation: ", target, sep = "")
  # note <- "'m' is the number of clusters in each group and 'n' is the number of individuals in each cluster."
  # structure(list(alpha = alpha, power = power, m = m, n = n, cv = cv, d = d,
  #                icc = icc, varw = varw, note = note, method = method),
  #           class = "power.htest")
  
}
