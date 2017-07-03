#------------------------------------------------------------------------------

# Miscellaneous functions

#------------------------------------------------------------------------------
#------------------------------------------------------------------------------

# function to get design effect

# INPUT:
# n is the mean cluster size
# ICC is the intraclass correlation
# cv is the coefficient of variation

# OUTPUT:
# the design effect of the study
# getDEFF <- function(n, ICC, cv) {
#   1 + (((cv^2 + 1)*n) - 1)*ICC
# }

# getDEFF <- function(n, cv, icc) {
#     1 + (((cv^2 + 1)*n) - 1)*icc
# }

getDEFF <- function(n, nsd, cv, varw, varb, icc) {

  # if n or cv is not null, use following formula for first part
  if(!is.null(nsd)) {
    deffstring1 <- "1 + ((((nsd/n)^2 + 1)*n) - 1)"
  } else {
    deffstring1 <- "1 + (((cv^2 + 1)*n) - 1)"
  }

  if (!is.null(varb)){
    deffstring2 <- "varb/(varw+varb)"
  } else {
    deffstring2 <- "icc"
  }
  deffstring <- paste(deffstring1,deffstring2,sep="*")
  eval(parse(text=deffstring))
}
#------------------------------------------------------------------------------

# function to calculate one of n, n_sd, or cv if other two are not null

calc_n <- function(ind, n, nsd, cv){
  # input index of list element that's NULL
  # 1 --> n = nsd/cv
  # 2 --> nsd = n*cv
  # 3 --> cv = nsd/n
  switch(ind,
         "1" = nsd/cv,
         "2" = n*cv,
         "3" = nsd/n)
}

#------------------------------------------------------------------------------

# function to calculate one of ICC, varw, varb if other two are not null
calc_icc <- function(ind, icc, varw, varb){
  # input index of list element that's NULL
  # 1 --> ICC
  # 2 --> varw
  # 3 --> varb
  switch(ind,
         "1" = varb/(varb+varw),
         "2" = varb*(1-icc)/icc,
         "3" = varw*icc/(1-icc))
}

#------------------------------------------------------------------------------

calc_p <- function(ind, p1, p2, d){
  # input index of list element that's NULL
  # 1 --> p1
  # 2 --> p2
  # 3 --> p3
  switch(ind,
         "1" = d + p2,
         "2" = p1 - d,
         "3" = p1 - p2)
}
