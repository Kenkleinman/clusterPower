#' Updates the Shiny app options if the function arguments change
#'
#'
#' @author Alexandria C. Sakrejda (\email{acbro0@@umass.edu}
#'
#' @param fxnName String. The name of the function from which arguments ahould
#' be drawn.
#'
#' @return List of numericInput expressions.
#' @export



argMatch <- function(fxnName) {
  library(shiny)
  library("clusterPower")
  
  nsim <- numericInput("nsim", "Number of simulations", value = 1000)
  narms <- numericInput("narms", "Number of arms", value = 3)
  mu <- numericInput("mu", "Mean in arm 1", value = 1.4)
  mu2 <- numericInput("mu2", "Mean in arm 2", value = 3.2)
  sigma_sq <- numericInput("sigma_sq", "Within-cluster variance", value = 0.01)
  sigma_b_sq <-
    numericInput("sigma_b_sq", "Between-cluster variance", value = 0.1)
  sigma_sq2 <- numericInput("sigma_sq2", "Within-cluster variance (Arm 2)", value = 0.01)
  sigma_b_sq2 <-
    numericInput("sigma_b_sq2", "Between-cluster variance (Arm 2)", value = 0.1)  
  CV <- numericInput("CV", "Coefficient of variation (CV)", value = 0.01)
  d <- numericInput("d", "Means difference", value = 1)
  rho_c <- numericInput("rho_c", "Baseline and post-test cluster-level correlation", value = 0)
  rho_s <- numericInput("rho_s", "Baseline and post-test subject-level correlation", value = 0)
  ntimes <- numericInput("ntimes", "Number of measurement time points", value = 3)
  sigma_b_sq0 <- numericInput("sigma_b_sq0", "Pre-treatment between-cluster variance", value = 0)
  sigma_b_sq1 <- numericInput("sigma_b_sq1", "Pre-treatment between-cluster variance", value = 0)
  ICC <-
    numericInput("ICC", "Intracluster correlation coefficient (ICC)", value = 0.05)
  vart <-
    numericInput("vart", "Total variation of the outcome", value = 10)
  seed <- numericInput("seed", "Set the seed (for repeatability)", value = NA)
  ncontrols <- numericInput("ncontrols", "Number of control subjects", value = 10)
  varu <- numericInput("varu", "Intervention arm cluster random effect variance")
  varei <- numericInput("varei", "Intervention arm subject random error variance")
  varr <- numericInput("varr", "Control arm subject random error variance")
  vara <- numericInput("vara", "Between-arm variance")
  varc <- numericInput("varc", "Between-cluster variance")
  vare <- numericInput("vare", "Within-cluster variance")
  p1 <- numericInput("p1",  "Outcome proportion (Arm 1)", value = 0.1)
  p2 <- numericInput("p2",  "Outcome proportion (Arm 2)", value = 0.5)
  pooled <- checkboxInput("pooled", "Pooled standard error?", value = FALSE)
  p1inc <- checkboxInput("p1inc", "p1 > p2?", value = FALSE)
  multi.p.method <- numericInput("multi.p.method", "Multiple comparisons adjustment", value = "bonferroni")
  tdist <- checkboxInput("tdist", "Use t-distribution?", value = FALSE)
  means <- numericInput("means", "Expected absolute treatment effect for each arm", value = c(1.2, 4.3, 5.2))
  low.power.override <- checkboxInput("low.power.override", "Complete run even if power < 0.5", value = FALSE)
  steps <- numericInput("steps", "Number of crossover steps", value = 3)
  p.e <- numericInput("p.e",  "Outcome probability (Intervention arm)", value = 0.5)
  p.c <- numericInput("p.c",  "Outcome probability (Control arm)", value = 0.1)
  decrease <- checkboxInput("decrease", "Intervention probability < control probability?", value = FALSE)
  mu0 <- numericInput("mu0", "Baseline (arm 1) effect", value = 0.1)
  beta <- numericInput("beta", "Treatment (arm 2) effect", value = 0.4) 
  probs <- numericInput("probs", "treatment effect probabilities for each arm", value = c(0.15, 0.23, 0.22)) 
  r1 <- numericInput("r1",  "Mean event rate (Arm 1)", value = 0.1)
  r2 <- numericInput("r2",  "Mean event rate (Arm 2)", value = 0.2)
  r1inc <- checkboxInput("r1inc", "Intervention probability < control probability?", value = FALSE)
  CVB <- numericInput("CVB", "Between-cluster coefficient of variation (CV)", value = 0.01)
  lambda1 <- numericInput("lambda1", "Baseline rate for outcome of interest", value = 1.75)
  RR <- numericInput("RR", "Intervention relative risk", value = 0.9)
  c1 <- numericInput("c1",  "Expected outcome count (Arm 1)", value = 5)
  c2 <- numericInput("c2",  "Expected outcome count (Arm 2)", value = 7)
  counts <- numericInput("counts", "Mean event rates per unit time for each arm", value = c(30, 35, 70, 40))


  
  temp <- dplyr::intersect(objects(), names(formals(fxnName)))
  
  holder <- list()
  for (i in 1:length(temp)) {
    holder[[i]] <- eval(parse(text = temp[i]))
  }
  return(holder)
}