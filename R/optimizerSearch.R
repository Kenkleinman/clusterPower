


sim <- function(nsub = 2, nclust = 3, sigma_residual2 = 1,  sig_b2 = c(.3, 10), betas = c(0,1)){
  narm <- 2
  n <- nsub * nclust * narm
  y <- rep_len(NA, n)
  arm <- rep(0:(narm-1), each = nsub*nclust)
  arm_factor <- as.factor(arm)
  clustid <- rep(1:(nclust*narm), each = nsub)
  clustRElist <- rnorm(narm*nclust, mean = 0, sd = rep(sqrt(sig_b2), each = nclust))
  clustRE <- rep(clustRElist, each = nsub)
  sig_b2 <- rep(sig_b2, each = nclust*nsub)
  error <- rnorm(n, mean = 0, sd = sqrt(sigma_residual2))
  beta <- rep(betas, each = nclust*nsub)
  y <- beta + clustRE + error
  output <- cbind.data.frame(arm, arm_factor, clustid, sig_b2, clustRE, error, beta, y)
  return(output)
}
set.seed(2)
dat <- sim(nsub = 100, nclust = 100, sigma_residual2 = .2,  sig_b2 = c(.3, 1), betas = c(0,1))

model <- lmer(y ~ arm + (1 + arm | clustid), data = dat)





diff_optims <- allFit(model, maxfun = 1e5, parallel = 'multicore', ncpus = detectCores())
is.OK <- sapply(diff_optims, is, "merMod")
diff_optims.OK <- diff_optims[is.OK]
lapply(diff_optims.OK,function(x) x@optinfo$conv$lme4$messages)

convergence_results <- lapply(diff_optims.OK,function(x) x@optinfo$conv$lme4$messages)
working_indices <- sapply(convergence_results, is.null)

if(sum(working_indices)==0){
  print("No algorithms from allFit converged.")
  print("You may still be able to use the results, but proceed with extreme caution.")
  first_fit <- NULL
} else {
  goodopt <- names(working_indices[working_indices==TRUE][1])
}

model2 <- lmer(y ~ arm + (1 + arm | clustid), data = dat, 
               control = lmerControl(optimizer = goodopt))
