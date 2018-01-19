# code to create crt did data with dropouts, replace dropouts, test model fit
library(lme4)


# X1 is number of groups (0 = control, 1 = treatment)
# X2 is number of time points (0 = pre, 1 = post)
# j, k, m are number of clusters per group, number of subjects per cluster,
#   and number of observations per subject.
# B0, B1, B2, B3 are fixed effects for intercept, group, time, and group*time
#   interaction, respectively.
# sdb1, sdb2, sdb3, sdb4, sde are variances for cluster intercepts,
#   cluster*time slopes, subject intercepts, subject*time slopes, and
#   measurement error.
# ltf is proportion of cluster size lost to follow up (considered the same across clusters)

generate_did_data <- function(j, k, m,
                              B0, B1, B2, B3,
                              sdb1, sdb2, sdb3, sdb4, sde,
                              ltf = 0, 
                              X1 = 0:1, X2 = 0:1){
  
  # determine number lost to follow up, use to determine number of all subjects
  kloss <- ceiling(k*ltf)
  kmax <- k + kloss
  
  data <- expand.grid(X1 = X1, X2 = X2,
                      B0 = B0, B1 = B1, B2 = B2, B3 = B3,
                      j = 1:2, k = 1:kmax, m = 1:m, stringsAsFactors = FALSE)
  
  data$clust <- paste(data$X1, data$j, sep=".")
  data$subj <- paste(data$X1, data$j, data$k, sep=".")
  
  # cluster related effects
  # get cluster names
  clust_names <- unique(data$clust) 
  # generate cluster and cluster*time effects
  clust_effs <- rnorm(length(clust_names), 0, sd = sdb1) 
  clust_time_effs <- rnorm(length(clust_names), 0, sd = sdb2) 
  # assign names
  names(clust_effs) <- clust_names 
  names(clust_time_effs) <- clust_names
  
  # subject related effects
  # get subject names 
  subj_names <- unique(data$subj)
  # generate cluster and cluster*time effects
  subj_effs <- rnorm(length(subj_names), 0, sd = sdb3)
  subj_time_effs <- rnorm(length(subj_names), 0, sd = sdb4)
  # assign names
  names(subj_effs) <- subj_names
  names(subj_time_effs) <- subj_names
  
  # initialize random effect vectors
  data$b1 <- clust_effs[data$clust] # cluster random effect
  data$b2 <- clust_effs[data$clust] # cluster*time random effect
  data$b3 <- subj_effs[data$subj] # subject random effect
  data$b4 <- subj_effs[data$subj] # subject*time random effect
  
  data$e <- rnorm(nrow(data), 0, sde)
  
  if(ltf == 0){
    within(data = data, {y <- B0 + B1*X1 + B2*X2 + B3*X1*X2 + b1 + b2*X2 + b3 + b4*X2 + e})
  } else {
    # determine which subjects to keep at pre, post
    keep1 <- with(data, which(X2 == 0 & k <= (kmax - kloss)))
    keep2 <- with(data, which(X2 == 1 & !(k %in% ) < (kmax - kloss) | k > (kmax - kloss)) ))
    keep <- c(keep1, keep2)
    data <- data[keep,]
    rownames(data) <- NULL
    data
    # keep1 indicates the original k subjects
    # keep2 contains those of the original subject who stayed as well as replacements
    # mutate(keep1 = if_else(X2 == 0 & k <= (kmax - kdiff), 1, 0),
    #        keep2 = if_else(X2 == 1 & (k < (kmax - kdiff) | k > (kmax - kdiff)), 1, 0),
    #        keep = keep1 + keep2) %>%
    # filter(keep == 1) %>%
    # mutate(y = B0 + B1*X1 + B2*X2 + B3*X1*X2 + b1 + b2*X2 + b3 + b4*X2 + e)
  }
}

j <- 2
k <- 3
m <- 1
B0 <- 0
B1 <- 0
B2 <- 0
B3 <- 1
sdb1 <- 0.1
sdb2 <- 0.1
sdb3 <- 0.1
sdb4 <- 0.1
sde <- 0.1

ltf <- 1

test <- generate_did_data(j,k,m,B0,B1,B2,B3,sdb1,sdb2,sdb3,sdb4,sde,ltf)

mod <- lmer(y ~ X1*X2 + (X2|clust) + (X2|subj), data = test)

