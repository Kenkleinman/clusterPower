print.crtpwr = function (x, ...) {
  cat(paste0("\n", x[['overview']], "\n"))
  cat(paste0("\nPower Estimate (alpha = ", x[['alpha']], "):\n"))
  print(x[['power']], row.names = FALSE)
  cat(paste("\nMethod:", switch(x[['method']], glmm = "Generalized Linear Mixed Model", gee = "Generalized Estimating Equation"), "\n"))
  cat("\n")
}
setClass('crtpwr')
