#' Power calculations for simple cluster randomized trials with matching, binary outcome
#'
#' Compute the power of a simple cluster randomized trial with a binary outcome with cluster-level matching,
#' or determine parameters to obtain a target power.
#'
#' @section Authors:
#' Jonathan Moyer (\email{jon.moyer@@gmail.com}), Ken Kleinman (\email{ken.kleinman@@gmail.com})
#'
#' @section Note:
#'   This function was inspired by work from Stephane Champely (pwr.t.test) and
#'   Peter Dalgaard (power.t.test). As with those functions, 'uniroot' is used to
#'   solve power equation for unknowns, so you may see
#'   errors from it, notably about inability to bracket the root when
#'   invalid arguments are given.
#'
#' @param alpha The level of significance of the test, the probability of a
#'   Type I error.
#' @param power The power of the test, 1 minus the probability of a Type II
#'   error.
#' @param m The number of clusters per condition. It must be greater than 1.
#' @param n The mean of the cluster sizes.
#' @param p1 The expected proportion in the treatment group.
#' @param p2 The proportion in the control group.
#' @param cvm The coefficient of variation in the outcome within matched clusters.
#' @param p1inc Logical indicating if p1 is expected to be greater than p2.
#' @param tol Numerical tolerance used in root finding. The default provides
#'   at least four significant digits.
#' @return The computed argument.
#' @examples 
#' # Find the number of clusters per condition needed for a trial with alpha = 0.05, 
#' # power = 0.80, 80 observations per probability
#' # in condition 1 of 0.80 and condition 2 of 0.60, within-pair coefficient of variation equal to 0.2.
#' crtpwr.2propM(n=80 ,p1=0.80, p2=0.60, cvm = 0.2)
#' # 
#' # The result, showimg m of greater than 10, suggests 11 clusters per condition should be used.
#' 
#' @references Hayes JR, Moulton LH. Cluster Randomized Trials. Boca Raton, FL: CRC Press; 2009.
#' 
#' @export

crtpwr.2propM <- function(alpha = 0.05, power = 0.80,
                         m = NA, n = NA, p1 = NA, p2 = NA,
                         cvm = 0, p1inc = TRUE,
                         tol = .Machine$double.eps^0.25){
  
  if(!is.na(m) && m <= 1) {
    stop("'m' must be greater than 1.")
  }
  
  needlist <- list(alpha, power, m, n, p1, p2, cvm)
  neednames <- c("alpha", "power", "m", "n", "p1", "p2", "cvm")
  needind <- which(unlist(lapply(needlist, is.na))) # find NA index
  
  if (length(needind) != 1) {
    stop("Exactly one of 'alpha', 'power', 'm', 'n', 'p1', 'p2', or 'cvm' must be NA.")
  }
  
  target <- neednames[needind]
  
  pwr <- quote({
    zcrit <- qnorm(alpha/2, lower.tail = FALSE)
    sdd <- sqrt((p1*(1-p1) + p2*(1-p2))/n + cvm^2*(p1^2 + p2^2))
    pnorm(sqrt((m - 2)*(p1 -p2)^2)/sdd - zcrit, lower.tail = TRUE)# +
    #pnorm(-zcrit - (p1 - p2)/sdd, lower.tail = TRUE)
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
                        interval = c(2 + 1e-10, 1e+07),
                        tol = tol)$root
  }
  
  # calculate p1
  if (is.na(p1)) {
    if(p1inc){
      p1 <- stats::uniroot(function(p1) eval(pwr) - power,
                           interval = c(p2 + 1e-7, 1 - 1e-7),
                           tol = tol)$root
    } else {
      p1 <- stats::uniroot(function(p1) eval(pwr) - power,
                           interval = c(1e-7, p2 - 1e-7),
                           tol = tol)$root
    }
  }
  
  # calculate p2
  if (is.na(p2)) {
    if(p1inc){
      p2 <- stats::uniroot(function(p2) eval(pwr) - power,
                           interval = c(1e-7, p1 - 1e-7),
                           tol = tol)$root
      
    } else {
      p2 <- stats::uniroot(function(p2) eval(pwr) - power,
                           interval = c(p1 + 1e-7, 1 - 1e-7),
                           tol = tol)$root
    }
  }
  
  # calculate n
  if (is.na(n)) {
    n <- stats::uniroot(function(n) eval(pwr) - power,
                        interval = c(2 + 1e-10, 1e+07),
                        tol = tol, extendInt = "upX")$root
  }
  
  # calculate cvm
  if (is.na(cvm)) {
    
    cvm <- stats::uniroot(function(cvm) eval(pwr) - power,
                         interval = c(1e-7, 1e+07),
                         tol = tol, extendInt = "downX")$root
  }
  
  structure(get(target), names = target)
  
}
