# Helper functions for the app

# "safe" versions of functions to catch errors
crtpwr.2mean.safe <- function(alpha,power,m,n,cv,d,icc,varw,method){
  # make safe version
  fun <- safely(crtpwr.2mean, otherwise = NA)
  # store result
  res <- fun(alpha,power,m,n,cv,d,icc,varw,method)
  # if res$error NULL, set to NA, otherwise set to message
  if(is.null(res$error)){
    res$error = 'None'
  } else {
    res$error <- res$error$message
  }
  res
}

crtpwr.2meanD.safe <- function(alpha,power,m,n,d,icc,rho_c,rho_s,varw){
  # make safe version
  fun <- safely(crtpwr.2meanD, otherwise = NA)
  # store result
  res <- fun(alpha,power,m,n,d,icc,rho_c,rho_s,varw)
  # if res$error NULL, set to NA, otherwise set to message
  if(is.null(res$error)){
    res$error = 'None'
  } else {
    res$error <- res$error$message
  }
  res
}

crtpwr.2meanM.safe <- function(alpha,power,m,n,d,icc,varw,rho_m){
  # make safe version
  fun <- safely(crtpwr.2meanM, otherwise = NA)
  # store result
  res <- fun(alpha,power,m,n,d,icc,varw,rho_m)
  # if res$error NULL, set to NA, otherwise set to message
  if(is.null(res$error)){
    res$error = 'None'
  } else {
    res$error <- res$error$message
  }
  res
}


crtpwr.2prop.safe <- function(alpha,power,m,n,cv,p1,p2,icc,pooled,p1inc){
  # make safe version
  fun <- safely(crtpwr.2prop, otherwise = NA)
  # store result
  res <- fun(alpha,power,m,n,cv,p1,p2,icc,pooled,p1inc)
  # if res$error NULL, set to NA, otherwise set to message
  if(is.null(res$error)){
    res$error = 'None'
  } else {
    res$error <- res$error$message
  }
  res
}

crtpwr.2rate.safe <- function(alpha,power,m,py,r1,r2,cvb,r1inc){
  # make safe version
  fun <- safely(crtpwr.2rate, otherwise = NA)
  # store result
  res <- fun(alpha,power,m,py,r1,r2,cvb,r1inc)
  # if res$error NULL, set to NA, otherwise set to message
  if(is.null(res$error)){
    res$error = 'None'
  } else {
    res$error <- res$error$message
  }
  res
}

# function to shorten error messages
shorten_error <- function(err, target){
  ifelse(err == "did not succeed extending the interval endpoints for f(lower) * f(upper) <= 0",
         paste0("No '", target, "' found for given inputs."),
         err)
}

# function to convert app input into numeric vectors
make_sequence <- function(x){
  # trim any spaces before and after the supplied string
  x <- str_trim(x)
  if(x == ""){
    # if x is an empty string, just return NA
    as.numeric(x)
  } else if(str_detect(x, "^[0-9.]")){
    # if the initial character of x is a digit or period, use the supplied
    # numbers to create a vector of possibilities
    temp <- as.numeric(str_split(x,"[^0-9.]",simplify=TRUE))
    temp[!is.na(temp)]
  } else {
    # if the initial character of x is not a digit or period, use the three
    # supplied numbers to create a sequence 'from x to y by z'
    temp <- as.numeric(str_split(x,"[^0-9.]",simplify=TRUE))
    temp <- temp[!is.na(temp)]
    seq(temp[1], temp[2], by = temp[3])
  }
}