#' Updates the Shiny app options if the function arguments change
#'
#'
#' @author Alexandria C. Sakrejda (\email{acbro0@@umass.edu}
#'
#' @param fxnName String. The name of the function from which arguments should
#' be drawn.
#' @param justNames Logical. Default is FALSE. Should only argument names be returned
#'
#' @return List of numericInput expressions.
#' @export



argMatch <- function(fxnName, justNames = FALSE) {
  require("clusterPower")

  # The arguments from which to choose
  nsim <-
    numericInput(
      "nsim",
      "Number of simulations",
      value = 1000,
      max = 500000,
      min = 0
    )
  narms <-
    numericInput(
      "narms",
      "Number of arms",
      value = 3,
      min = 2,
      step = 1
    )
  mu <- numericInput("mu", "Mean in arm 1", value = 1.4)
  mu2 <- numericInput("mu2", "Mean in arm 2", value = 3.2)
  sigma_sq <-
    numericInput("sigma_sq",
                 "Within-cluster variance",
                 value = 0.01,
                 min = 0)
  sigma_b_sq <-
    numericInput("sigma_b_sq",
                 "Between-cluster variance",
                 value = 0.1,
                 min = 0)
  sigma_sq2 <-
    numericInput("sigma_sq2",
                 "Within-cluster variance (Arm 2)",
                 value = 0.01,
                 min = 0)
  sigma_b_sq2 <-
    numericInput(
      "sigma_b_sq2",
      "Between-cluster variance (Arm 2)",
      value = 0.1,
      min = 0
    )
  CV <-
    numericInput("CV", "Coefficient of variation (CV)", value = 0)
  d <- numericInput("d", "Means difference", value = 1)
  rho_c <-
    numericInput("rho_c",
                 "Baseline and post-test cluster-level correlation",
                 value = 0)
  rho_s <-
    numericInput("rho_s",
                 "Baseline and post-test subject-level correlation",
                 value = 0)
  ntimes <-
    numericInput(
      "ntimes",
      "Number of measurement time points",
      value = 3,
      min = 0,
      step = 1
    )
  sigma_b_sq0 <-
    numericInput(
      "sigma_b_sq0",
      "Pre-treatment between-cluster variance",
      value = 0,
      min = 0
    )
  sigma_b_sq1 <-
    numericInput(
      "sigma_b_sq1",
      "Pre-treatment between-cluster variance",
      value = 0,
      min = 0
    )
  if (fxnName == "cps.binary"){
    ICC <-
      numericInput(
        "ICC",
        "Intracluster correlation coefficient (ICC)",
        value = 0.05,
        min = 0,
        max = 1
      )
  } else {
  ICC <-
    numericInput(
      "ICC",
      "Intracluster correlation coefficient (ICC)",
      value = NA,
      min = 0,
      max = 1
    )
}
  vart <-
    numericInput("vart", "Total variation of the outcome", value = NA)
  seed <-
    numericInput("seed",
                 "Set the seed (for repeatability)",
                 value = NA,
                 step = 1)
  ncontrols <-
    numericInput(
      "ncontrols",
      "Number of control subjects",
      value = 10,
      step = 1,
      min = 0
    )
  varu <-
    numericInput(
      "varu",
      "Intervention arm cluster random effect variance",
      value = 0.1,
      min = 0
    )
  varei <-
    numericInput(
      "varei",
      "Intervention arm subject random error variance",
      value = 0.1,
      min = 0
    )
  varr <-
    numericInput(
      "varr",
      "Control arm subject random error variance",
      value = 0.1,
      min = 0
    )
  vara <-
    numericInput("vara",
                 "Between-arm variance",
                 value = 0.1,
                 min = 0)
  varc <-
    numericInput("varc",
                 "Between-cluster variance",
                 value = 0.1,
                 min = 0)
  vare <-
    numericInput("vare",
                 "Within-cluster variance",
                 value = 0.1,
                 min = 0)
  p1 <-
    numericInput(
      "p1",
      "Outcome proportion (Arm 1)",
      value = 0.1,
      min = 0,
      max = 1
    )
  p2 <-
    numericInput(
      "p2",
      "Outcome proportion (Arm 2)",
      value = 0.5,
      min = 0,
      max = 1
    )
  pooled <-
    checkboxInput("pooled", "Pooled standard error", value = FALSE)
  p1inc <- checkboxInput("p1inc", "p1 > p2", value = FALSE)
  multi_p_method <-
    selectInput(
      "multi.p.method",
      "Multiple comparisons adjustment",
      choices = c(
        "holm",
        "hochberg",
        "hommel",
        "bonferroni",
        "BH",
        "BY",
        "fdr",
        "none"
      ),
      selected = "bonferroni",
      multiple = FALSE
    )
  tdist <-
    checkboxInput("tdist", "Use t-distribution", value = FALSE)
  means <-
    textInput(
      "means",
      "Expected absolute treatment effect for each arm (comma delimited)",
      "1.2, 4.3, 5.2"
    )
  steps <-
    numericInput("steps", "Number of crossover steps", value = 3)
  pe <-
    numericInput(
      "pe",
      "Outcome probability (Intervention arm)",
      value = 0.5,
      min = 0,
      max = 1
    )
  pc <- numericInput(
    "pc",
    "Outcome probability (Control arm)",
    value = 0.1,
    min = 0,
    max = 1
  )
  decrease <- checkboxInput("decrease",
                            "Intervention probability < control probability",
                            value = FALSE)
  mu0 <- numericInput("mu0", "Baseline (arm 1) effect", value = 0.1)
  beta <-
    numericInput("beta", "Treatment (arm 2) effect", value = 0.4)
  probs <-
    textInput(
      "probs",
      "Treatment effect probabilities for each arm (comma delimited)",
      "0.15, 0.23, 0.22"
    )
  r1 <- numericInput("r1",  "Mean event rate (Arm 1)", value = 0.1)
  r2 <- numericInput("r2",  "Mean event rate (Arm 2)", value = 0.2)
  r1inc <-
    checkboxInput("r1inc",
                  "Intervention probability < control probability",
                  value = FALSE)
  CVB <-
    numericInput("CVB", "Between-cluster coefficient of variation (CV)", value = 0.01)
  lambda1 <-
    numericInput("lambda1", "Baseline rate for outcome of interest", value = 1.75)
  RR <- numericInput("RR", "Intervention relative risk", value = 0.9)
  c1 <-
    numericInput(
      "c1",
      "Expected outcome count (Arm 1)",
      value = 5,
      step = 1,
      min = 0
    )
  c2 <-
    numericInput(
      "c2",
      "Expected outcome count (Arm 2)",
      value = 7,
      step = 1,
      min = 0
    )
  counts <-
    textInput("counts",
              "Mean event per unit time for each arm (comma delimited)",
              "30, 35, 70")
  
  if (justNames == TRUE){
    return(objects())
  }
  
# compare the function arguments to the objects above  
  temp <- dplyr::intersect(objects(), names(formals(fxnName)))
  
# store the expressions in a list 
  holder <- list()
  for (i in 1:length(temp)) {
    if (!is.null(temp[i])) {
      holder[[i]] <- eval(parse(text = temp[i]))
    } else {
      next()
    }
  }
  return(holder)
}