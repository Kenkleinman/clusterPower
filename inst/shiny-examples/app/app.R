#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shinythemes)
library(shiny)
library(shinyBS)

get_vignette_link <- function(...) {
  x <- vignette(...)
  if (nzchar(out <- x$PDF)) {
    ext <- tools::file_ext(out)
    port <- if (tolower(ext) == "html")
      tools::startDynamicHelp(NA)
    else
      0L
    if (port > 0L) {
      out <- sprintf("http://127.0.0.1:%d/library/%s/doc/%s",
                     port,
                     basename(x$Dir),
                     out)
      return(out)
    }
  }
  stop("no html help found")
}


ui <- fluidPage(
  theme = shinytheme("united"),
  h1(
    id = "big-heading",
    "Power Estimation for Randomized Controlled Trials: clusterPower"
  ),
  tags$style(HTML("#big-heading{color: #337ab7;}")),
  shinyjs::useShinyjs(),
  sidebarLayout(
    sidebarPanel(
      selectInput(
        "type",
        "CRT Type",
        choices = c(
          "Parallel",
          "Multi-Arm",
          "Difference-in-Difference",
          "Stepped Wedge",
          "Individually-Randomized Group"
        )
      ),
      selectInput(
        "dist",
        "Outcome Distribution",
        choices = c("Normal", "Binary", "Count")
      ),
      selectInput("meth", "Method",
                  choices = c("Analytic", "Simulation")),
      # Below values can be reset to defaults using the restore defaults button
      div(
        id = "allValues",
        shinyjs::hidden(numericInput("power", "power", value = NA)),
        
        #input for cpa.normal
        conditionalPanel(
          "input.type == 'Parallel' & input.dist == 'Normal' & input.meth == 'Analytic'",
          numericInput("nclusters", "Number of Clusters", value = 10),
          numericInput("nsubjects", "Number of Observations (per cluster)", value = 20),
          numericInput("CVcpanormal", "Coefficient of variation (CV)", value = 0),
          numericInput("dcpanormal", "Means difference", value = 0.43),
          numericInput(
            "ICCcpanormal",
            "Intracluster correlation coefficient (ICC)",
            value = NA,
            step = 0.01,
            min = 0,
            max = 1
          ),
          numericInput(
            "sigma_sqcpanormal",
            "Within-cluster variance (Arm 1)",
            value = 0.01,
            step = 0.001,
            min = 0
          ),
          numericInput(
            "sigma_b_sqcpanormal",
            "Between-cluster variance (Arm 1)",
            step = 0.001,
            value = 0.1,
            min = 0
          ),
          numericInput("vartcpanormal", "Total variation of the outcome", value = NA)
        ),
        
        # input for cps.normal
        conditionalPanel(
          "input.type == 'Parallel' & input.dist == 'Normal' & input.meth == 'Simulation'",
          numericInput("nclusters", "Number of Clusters", value = 10),
          numericInput("nsubjects", "Number of Observations (per cluster)", value = 20),
          numericInput(
            "nsimcpsnormal",
            "Number of simulations",
            value = 100,
            max = 500000,
            min = 0
          ),
          numericInput(
            "ICCcpsnormal",
            "Intracluster correlation coefficient (ICC, Arm 1)",
            value = NA,
            step = 0.01,
            min = 0,
            max = 1
          ),
          numericInput(
            "ICC2cpsnormal",
            "Intracluster correlation coefficient (ICC, Arm 2)",
            value = NA,
            step = 0.01,
            min = 0,
            max = 1
          ),
          numericInput("mucpsnormal", "Mean in arm 1", value = 2.4),
          numericInput("mu2cpsnormal", "Mean in arm 2", value = 1.5),
          numericInput(
            "sigma_sqcpsnormal",
            "Within-cluster variance (Arm 1)",
            value = 0.2,
            step = 0.001,
            min = 0
          ),
          numericInput(
            "sigma_sq2cpsnormal",
            "Within-cluster variance (Arm 2)",
            value = 0.2,
            step = 0.001,
            min = 0
          ),
          numericInput(
            "sigma_b_sqcpsnormal",
            "Between-cluster variance (Arm 1)",
            value = 0.5,
            step = 0.001,
            min = 0
          ),
          numericInput(
            "sigma_b_sq2cpsnormal",
            "Between-cluster variance (Arm 2)",
            value = 0.5,
            step = 0.001,
            min = 0
          )
        ),
        
        # cpa.binary inputs start
        conditionalPanel(
          "input.type == 'Parallel' & input.dist == 'Binary' & input.meth == 'Analytic'",
          numericInput("nclusters", "Number of Clusters", value = 10),
          numericInput("nsubjects", "Number of Observations (per cluster)", value = 20),
          numericInput("CVcpabinary", "Coefficient of variation (CV)", value = 0),
          numericInput(
            "ICCcpabinary",
            "Intracluster correlation coefficient (ICC)",
            value = 0.05,
            step = 0.01,
            min = 0,
            max = 1
          ),
          numericInput(
            "p1cpabinary",
            "Outcome proportion (Arm 1)",
            value = 0.1,
            step = 0.001,
            min = 0,
            max = 1
          ),
          numericInput(
            "p2cpabinary",
            "Outcome proportion (Arm 2)",
            value = 0.24,
            step = 0.001,
            min = 0,
            max = 1
          ),
          checkboxInput("pooledcpabinary", "Pooled standard error", value = FALSE),
          checkboxInput("p1inccpabinary", "p1 > p2", value = FALSE),
          checkboxInput("tdistcpabinary", "Use t-distribution", value = FALSE)
        ),
        
        # cps.binary inputs start
        conditionalPanel(
          "input.type == 'Parallel' & input.dist == 'Binary' & input.meth == 'Simulation'",
          numericInput("nclusters", "Number of Clusters", value = 10),
          numericInput("nsubjects", "Number of Observations (per cluster)", value = 20),
          numericInput(
            "nsimcpsbinary",
            "Number of simulations",
            value = 100,
            max = 500000,
            min = 0
          ),
          numericInput(
            "p1cpsbinary",
            "Outcome proportion (Arm 1)",
            value = 0.8,
            step = 0.001,
            min = 0,
            max = 1
          ),
          numericInput(
            "p2cpsbinary",
            "Outcome proportion (Arm 2)",
            value = 0.5,
            step = 0.001,
            min = 0,
            max = 1
          ),
          numericInput(
            "sigma_b_sqcpsbinary",
            "Between-cluster variance (Arm 1)",
            value = 1,
            step = 0.001,
            min = 0
          ),
          numericInput(
            "sigma_b_sq2cpsbinary",
            "Between-cluster variance (Arm 2)",
            value = 1,
            step = 0.001,
            min = 0
          )
        ),
        
        # cpa.count input starts
        conditionalPanel(
          "input.type == 'Parallel' & input.dist == 'Count' & input.meth == 'Analytic'",
          numericInput("nclusters", "Number of Clusters", value = 10),
          numericInput("nsubjects", "Number of Observations (per cluster)", value = 20),
          numericInput(
            "CVBcpacount",
            "Between-cluster coefficient of variation (CV)",
            value = 0.01
          ),
          numericInput("r1cpacount",  "Mean event rate (Arm 1)", value = 0.2),
          numericInput("r2cpacount",  "Mean event rate (Arm 2)", value = 0.35),
          checkboxInput(
            "r1inccpacount",
            "Intervention probability < control probability",
            value = FALSE
          )
        ),
        
        # cps.count input starts
        conditionalPanel(
          "input.type == 'Parallel' & input.dist == 'Count' & input.meth == 'Simulation'",
          numericInput("nclusters", "Number of Clusters", value = 10),
          numericInput("nsubjects", "Number of Observations (per cluster)", value = 20),
          numericInput(
            "nsimcpscount",
            "Number of simulations",
            value = 100,
            max = 500000,
            min = 0
          ),
          numericInput(
            "sigma_b_sqcpscount",
            "Between-cluster variance (Arm 1)",
            value = 0.5,
            step = 0.001,
            min = 0
          ),
          numericInput(
            "sigma_b_sq2cpscount",
            "Between-cluster variance (Arm 2)",
            value = 0.5,
            step = 0.001,
            min = 0
          ),
          numericInput(
            "c1cpscount",
            "Expected outcome count (Arm 1)",
            value = 200,
            step = 1,
            min = 0
          ),
          numericInput(
            "c2cpscount",
            "Expected outcome count (Arm 2)",
            value = 80,
            step = 1,
            min = 0
          )
        ),
        
        # cpa.ma.normal input starts
        conditionalPanel(
          "input.type == 'Multi-Arm' & input.dist == 'Normal' & input.meth == 'Analytic'",
          numericInput("nclusterscpamanormal", "Number of Clusters", value = 10),
          numericInput("nsubjectscpamanormal", "Number of Observations (per cluster)", value = 20),
          numericInput(
            "narmscpamanormal",
            "Number of arms",
            value = 3,
            max = 50,
            min = 2,
            step = 1
          ),
          numericInput(
            "varacpamanormal",
            "Between-arm variance",
            value = 0.02,
            step = 0.001,
            min = 0
          ),
          numericInput(
            "varccpamanormal",
            "Between-cluster variance",
            value = 0.01,
            step = 0.001,
            min = 0
          ),
          numericInput(
            "varecpamanormal",
            "Within-cluster variance",
            value = 0.1,
            step = 0.001,
            min = 0
          )
        ),
        
        # cps.ma.normal input start
        conditionalPanel(
          "input.type == 'Multi-Arm' & input.dist == 'Normal' & input.meth == 'Simulation'",
          textInput("nclusterscpsmanormal", "Number of Clusters (comma delimited)", value = "10, 10, 10"),
          textInput(
            "nsubjectscpsmanormal",
            "Number of Observations (per cluster, comma delimited)",
            value = "20, 20, 20"
          ),
          numericInput(
            "nsimcpsmanormal",
            "Number of simulations",
            value = 100,
            max = 500000,
            min = 0
          ),
          
          sliderInput(
            "narmscpsmanormal",
            "Number of arms",
            value = 3,
            min = 2,
            max = 10,
            step = 1
          ),
          textInput(
            "meanscpsmanormal",
            "Expected absolute treatment effect for each arm (comma delimited)",
            "22.0, 21.0, 22.5"
          ),
          textInput("ICCcpsmanormal", "Intracluster correlation coefficient (ICC)", value = "NA, NA, NA"),
          textInput(
            "sigma_sqcpsmanormal",
            "Within-cluster variance (comma delimited)",
            value = "0.1, 0.1, 0.1"
          ),
          textInput(
            "sigma_b_sqcpsmanormal",
            "Between-cluster variance (comma delimited)",
            value = "0.01, 0.01, 0.01"
          ),
          selectInput(
            "multi_p_methodcpsmanormal",
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
          ),
          checkboxInput("tdistcpsmanormal", "Use t-distribution", value = FALSE)
        ),
        
        # cpa.ma.binary (no method)
        conditionalPanel(
          "input.type == 'Multi-Arm' & input.dist == 'Binary' & input.meth == 'Analytic'",
          HTML("No method exists. Use the simulation option instead.")
        ),
        
        # cps.ma.binary input start
        conditionalPanel(
          "input.type == 'Multi-Arm' & input.dist == 'Binary' & input.meth == 'Simulation'",
          numericInput("nclusters", "Number of Clusters", value = 10),
          numericInput("nsubjects", "Number of Observations (per cluster)", value = 20),
          numericInput(
            "nsimcpsmabinary",
            "Number of simulations",
            value = 100,
            max = 500000,
            min = 0
          ),
          sliderInput(
            "narmscpsmabinary",
            "Number of arms",
            value = 3,
            min = 2,
            max = 10,
            step = 1
          ),
          textInput(
            "probscpsmabinary",
            "Treatment effect probabilities for each arm (comma delimited)",
            "0.15, 0.23, 0.22"
          ),
          textInput(
            "sigma_b_sqcpsmabinary",
            "Between-cluster variance (comma delimited)",
            value = "0.01, 0.01, 0.01"
          ),
          selectInput(
            "multi_p_methodcpsmabinary",
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
          ),
          checkboxInput("tdistcpsmabinary", "Use t-distribution", value = FALSE)
        ),
        
        # cpa.ma.count (no method)
        conditionalPanel(
          "input.type == 'Multi-Arm' & input.dist == 'Count' & input.meth == 'Analytic'",
          HTML("No method exists. Use the simulation option instead.")
        ),
        
        # cps.ma.count input start
        conditionalPanel(
          "input.type == 'Multi-Arm' & input.dist == 'Count' & input.meth == 'Simulation'",
          numericInput("nclusters", "Number of Clusters", value = 10),
          numericInput("nsubjects", "Number of Observations (per cluster)", value = 20),
          numericInput(
            "nsimcpsmacount",
            "Number of simulations",
            value = 100,
            max = 500000,
            min = 0
          ),
          sliderInput(
            "narmscpsmacount",
            "Number of arms",
            value = 3,
            min = 2,
            max = 10,
            step = 1
          ),
          textInput(
            "countscpsmacount",
            "Mean event per unit time for each arm (comma delimited)",
            "30, 35, 70"
          ),
          textInput(
            "sigma_b_sqcpsmacount",
            "Between-cluster variance (comma delimited)",
            value = "0.01, 0.01, 0.01"
          ),
          selectInput(
            "multi_p_methodcpsmacount",
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
          ),
          checkboxInput("tdistcpsmacount", "Use t-distribution", value = FALSE)
        ),
        
        # cpa.did.normal input start
        conditionalPanel(
          "input.type == 'Difference-in-Difference' & input.dist == 'Normal' & input.meth == 'Analytic'",
          numericInput("nclusters", "Number of Clusters", value = 10),
          numericInput("nsubjects", "Number of Observations (per cluster)", value = 20),
          numericInput("dcpadidnormal", "Means difference", value = 1),
          numericInput(
            "ICCcpadidnormal",
            "Intracluster correlation coefficient (ICC)",
            value = NA,
            step = 0.01,
            min = 0,
            max = 1
          ),
          numericInput(
            "rho_ccpadidnormal",
            "Baseline and post-test cluster-level correlation",
            value = 0
          ),
          numericInput(
            "rho_scpadidnormal",
            "Baseline and post-test subject-level correlation",
            value = 0
          ),
          numericInput("vartcpadidnormal", "Total variation of the outcome", value = NA)
        ),
        
        # cps.did.normal input start
        conditionalPanel(
          "input.type == 'Difference-in-Difference' & input.dist == 'Normal' & input.meth == 'Simulation'",
          numericInput("nclusters", "Number of Clusters", value = 10),
          numericInput("nsubjects", "Number of Observations (per cluster)", value = 20),
          numericInput(
            "nsimcpsdidnormal",
            "Number of simulations",
            value = 100,
            max = 500000,
            min = 0
          ),
          numericInput("mucpsdidnormal", "Mean in arm 1", value = 1.4),
          numericInput("mu2cpsdidnormal", "Mean in arm 2", value = 3.2),
          numericInput(
            "sigma_sq",
            "Within-cluster variance (Arm 1)",
            step = 0.001,
            value = 0.01,
            min = 0
          ),
          numericInput(
            "sigma_b_sq0cpsdidnormal",
            "Pre-treatment between-cluster variance",
            value = 0,
            step = 0.001,
            min = 0
          ),
          numericInput(
            "sigma_b_sq1cpsdidnormal",
            "Post-treatment between-cluster variance",
            value = 0,
            step = 0.001,
            min = 0
          )
        ),
        
        # cpa.did.binary input start
        conditionalPanel(
          "input.type == 'Difference-in-Difference' & input.dist == 'Binary' & input.meth == 'Analytic'",
          numericInput("nclusters", "Number of Clusters", value = 10),
          numericInput("nsubjects", "Number of Observations (per cluster)", value = 20),
          numericInput("dcpadidbinary", "Means difference", value = 1),
          numericInput("pcpadidbinary", "Mean post-test expected proportion", value = 1),
          numericInput(
            "rho_ccpadidbinary",
            "Baseline and post-test cluster-level correlation",
            value = 0
          ),
          numericInput(
            "rho_scpadidbinary",
            "Baseline and post-test subject-level correlation",
            value = 0
          )
        ),
        
        # cps.did.binary input start
        conditionalPanel(
          "input.type == 'Difference-in-Difference' & input.dist == 'Binary' & input.meth == 'Simulation'",
          numericInput("nclusters", "Number of Clusters", value = 10),
          numericInput("nsubjects", "Number of Observations (per cluster)", value = 20),
          numericInput(
            "nsimcpsdidbinary",
            "Number of simulations",
            value = 100,
            max = 500000,
            min = 0
          ),
          numericInput(
            "p1cpsdidbinary",
            "Outcome proportion (Arm 1)",
            value = 0.1,
            step = 0.001,
            min = 0,
            max = 1
          ),
          numericInput(
            "p2cpsdidbinary",
            "Outcome proportion (Arm 2)",
            value = 0.5,
            step = 0.001,
            min = 0,
            max = 1
          ),
          numericInput(
            "sigma_b_sq0cpsdidbinary",
            "Pre-treatment between-cluster variance",
            value = 0,
            step = 0.001,
            min = 0
          ),
          numericInput(
            "sigma_b_sq1cpsdidbinary",
            "Post-treatment between-cluster variance",
            value = 0,
            step = 0.001,
            min = 0
          )
        ),
        
        # cpa.did.count (no method)
        conditionalPanel(
          "input.type == 'Difference-in-Difference' & input.dist == 'Count' & input.meth == 'Analytic'",
          HTML("No method exists. Use the simulation option instead.")
        ),
        
        # cps.did.count input start
        conditionalPanel(
          "input.type == 'Difference-in-Difference' & input.dist == 'Count' & input.meth == 'Simulation'",
          numericInput("nclusters", "Number of Clusters", value = 10),
          numericInput("nsubjects", "Number of Observations (per cluster)", value = 20),
          numericInput(
            "nsimcpsdidcount",
            "Number of simulations",
            value = 100,
            max = 500000,
            min = 0
          ),
          numericInput(
            "c1cpsdidcount",
            "Expected outcome count (Arm 1)",
            value = 5,
            step = 1,
            min = 0
          ),
          numericInput(
            "c2cpsdidcount",
            "Expected outcome count (Arm 2)",
            value = 7,
            step = 1,
            min = 0
          ),
          numericInput(
            "sigma_b_sq0cpsdidcount",
            "Pre-treatment between-cluster variance",
            value = 0,
            step = 0.001,
            min = 0
          ),
          numericInput(
            "sigma_b_sq1cpsdidcount",
            "Post-treatment between-cluster variance",
            value = 0,
            step = 0.001,
            min = 0
          )
        ),
        
        # cpa.sw.normal input start
        conditionalPanel(
          "input.type == 'Stepped Wedge' & input.dist == 'Normal' & input.meth == 'Analytic'",
          numericInput("nclusters", "Number of Clusters", value = 10),
          numericInput("nsubjects", "Number of Observations (per cluster)", value = 20),
          numericInput(
            "ntimescpaswnormal",
            "Number of measurement time points",
            value = 3,
            min = 0,
            step = 1
          ),
          numericInput("dcpaswnormal", "Means difference", value = 1),
          numericInput(
            "rho_ccpaswnormal",
            "Baseline and post-test cluster-level correlation",
            value = 0
          ),
          numericInput(
            "rho_scpaswnormal",
            "Baseline and post-test subject-level correlation",
            value = 0
          ),
          numericInput("vartcpaswnormal", "Total variation of the outcome", value = NA)
        ),
        
        # cps.sw.normal input start
        conditionalPanel(
          "input.type == 'Stepped Wedge' & input.dist == 'Normal' & input.meth == 'Simulation'",
          numericInput("nclusters", "Number of Clusters", value = 10),
          numericInput("nsubjects", "Number of Observations (per cluster)", value = 20),
          numericInput(
            "nsimcpsswnormal",
            "Number of simulations",
            value = 100,
            max = 500000,
            min = 0
          ),
          numericInput("stepscpsswnormal", "Number of crossover steps", value = 3),
          numericInput("mucpsswnormal", "Mean in arm 1", value = 1.4),
          numericInput("mu2cpsswnormal", "Mean in arm 2", value = 3.2),
          numericInput(
            "sigma_sqcpsswnormal",
            "Within-cluster variance (Arm 1)",
            value = 0.01,
            step = 0.001,
            min = 0
          ),
          numericInput(
            "sigma_b_sqcpsswnormal",
            "Between-cluster variance (Arm 1)",
            value = 0.1,
            step = 0.001,
            min = 0
          )
        ),
        
        # cpa.sw.binary input start
        conditionalPanel(
          "input.type == 'Stepped Wedge' & input.dist == 'Binary' & input.meth == 'Analytic'",
          numericInput("nclusters", "Number of Clusters", value = 10),
          numericInput("nsubjects", "Number of Observations (per cluster)", value = 20),
          numericInput("stepscpaswbinary", "Number of crossover steps", value = 3),
          numericInput("dcpaswbinary", "Means difference", value = 1),
          numericInput("mu0cpaswbinary", "Baseline (arm 1) effect", value = 0.1),
          numericInput("betacpaswbinary", "Treatment (arm 2) effect", value = 0.4),
          numericInput(
            "ICCcpaswbinary",
            "Intracluster correlation coefficient (ICC)",
            value = NA,
            step = 0.01,
            min = 0,
            max = 1
          )
        ),
        
        # cps.sw.binary input start
        conditionalPanel(
          "input.type == 'Stepped Wedge' & input.dist == 'Binary' & input.meth == 'Simulation'",
          numericInput("nclusters", "Number of Clusters", value = 10),
          numericInput("nsubjects", "Number of Observations (per cluster)", value = 20),
          numericInput(
            "nsimcpsswbinary",
            "Number of simulations",
            value = 100,
            max = 500000,
            min = 0
          ),
          numericInput("stepscpsswbinary", "Number of crossover steps", value = 3
          ),
          numericInput(
            "p1cpsswbinary",
            "Outcome proportion (Arm 1)",
            value = 0.1,
            step = 0.001,
            min = 0,
            max = 1
          ),
          numericInput(
            "p2cpsswbinary",
            "Outcome proportion (Arm 2)",
            value = 0.5,
            step = 0.001,
            min = 0,
            max = 1
          ),
          numericInput(
            "sigma_b_sqcpsswbinary",
            "Between-cluster variance (Arm 1)",
            value = 0.1,
            step = 0.001,
            min = 0
          )
        ),
        
        
        # cpa.sw.count input start
        conditionalPanel(
          "input.type == 'Stepped Wedge' & input.dist == 'Count' & input.meth == 'Analytic'",
          numericInput("nclusters", "Number of Clusters", value = 10),
          numericInput("nsubjects", "Number of Observations (per cluster)", value = 20),
          numericInput("stepscpaswcount", "Number of crossover steps", value = 3),
          numericInput(
            "ICCcpaswcount",
            "Intracluster correlation coefficient (ICC)",
            value = NA,
            step = 0.01,
            min = 0,
            max = 1
          ),
          numericInput(
            "lambda1cpaswcount",
            "Baseline rate for outcome of interest",
            value = 1.75
          ),
          numericInput("RRcpaswcount", "Intervention relative risk", value = 0.9)
        ),
        
        # cps.sw.count input start
        conditionalPanel(
          "input.type == 'Stepped Wedge' & input.dist == 'Count' & input.meth == 'Simulation'",
          numericInput("nclusters", "Number of Clusters", value = 10),
          numericInput("nsubjects", "Number of Observations (per cluster)", value = 20),
          numericInput(
            "nsimcpsmacount",
            "Number of simulations",
            value = 100,
            max = 500000,
            min = 0
          ),
          numericInput("stepscpsswnormal", "Number of crossover steps", value = 3),
          numericInput(
            "c1cpsdidcount",
            "Expected outcome count (Arm 1)",
            value = 5,
            step = 1,
            min = 0
          ),
          numericInput(
            "c2cpsdidcount",
            "Expected outcome count (Arm 2)",
            value = 7,
            step = 1,
            min = 0
          ),
          numericInput(
            "sigma_b_sqcpsdidcount",
            "Between-cluster variance (Arm 1)",
            value = 0.1,
            step = 0.001,
            min = 0
          )
        ),
        
        # cpa.irgtt.normal input start
        conditionalPanel(
          "input.type == 'Individually-Randomized Group' & input.dist == 'Normal' & input.meth == 'Analytic'",
          numericInput("nclusters", "Number of Clusters", value = 10),
          numericInput("nsubjects", "Number of Observations (per cluster)", value = 20),
          numericInput("dcpairgttnormal", "Means difference", value = 1),
          numericInput(
            "ncontrolscpairgttnormal",
            "Number of control subjects",
            value = 10,
            step = 1,
            min = 0
          ),
          numericInput(
            "vareicpairgttnormal",
            "Intervention arm subject random error variance",
            value = 0.1,
            step = 0.001,
            min = 0
          ),
          numericInput(
            "varrcpairgttnormal",
            "Control arm subject random error variance",
            value = 0.1,
            step = 0.001,
            min = 0
          ),
          numericInput(
            "varucpairgttnormal",
            "Intervention arm cluster random effect variance",
            value = 0.1,
            step = 0.001,
            min = 0
          )
        ),
        
        # cps.irgtt.normal input start
        conditionalPanel(
          "input.type == 'Individually-Randomized Group' & input.dist == 'Normal' & input.meth == 'Simulation'",
          numericInput("nclusters", "Number of Clusters", value = 10),
          numericInput("nsubjects", "Number of Observations (per cluster)", value = 20),
          numericInput(
            "nsimcpsirgttnormal",
            "Number of simulations",
            value = 100,
            max = 500000,
            min = 0
          ),
          numericInput("mucpsirgttnormal", "Mean in arm 1", value = 1.4),
          numericInput("mu2cpsirgttnormal", "Mean in arm 2", value = 3.2),
          numericInput(
            "ICCcpsirgttnormal",
            "Intracluster correlation coefficient (ICC)",
            value = NA,
            step = 0.01,
            min = 0,
            max = 1
          ),
          numericInput(
            "sigma_sqcpsirgttnormal",
            "Within-cluster variance (Arm 1)",
            value = 0.01,
            step = 0.001,
            min = 0
          ),
          numericInput(
            "sigma_sq2cpsirgttnormal",
            "Within-cluster variance (Arm 2)",
            value = 0.01,
            step = 0.001,
            min = 0
          ),
          numericInput(
            "sigma_b_sqcpsirgttnormal",
            "Between-cluster variance (Arm 1)",
            value = 0.1,
            step = 0.001,
            min = 0
          ),
          numericInput(
            "sigma_b_sq2cpsirgttnormal",
            "Between-cluster variance (Arm 2)",
            value = 0.1,
            step = 0.001,
            min = 0
          )
        ),
        
        # cpa.irgtt.binary input start
        conditionalPanel(
          "input.type == 'Individually-Randomized Group' & input.dist == 'Binary' & input.meth == 'Analytic'",
          numericInput("nclusters", "Number of Clusters", value = 10),
          numericInput("nsubjects", "Number of Observations (per cluster)", value = 20),
          numericInput(
            "p1cpairgttbinary",
            "Outcome proportion (Arm 1)",
            value = 0.1,
            step = 0.001,
            min = 0,
            max = 1
          ),
          numericInput(
            "p2cpairgttbinary",
            "Outcome proportion (Arm 2)",
            value = 0.5,
            step = 0.001,
            min = 0,
            max = 1
          ),
          numericInput(
            "ncontrolscpairgttbinary",
            "Number of control subjects",
            value = 10,
            step = 1,
            min = 0
          ),
          numericInput(
            "ICCcpairgttbinary",
            "Intracluster correlation coefficient (ICC)",
            value = NA,
            step = 0.01,
            min = 0,
            max = 1
          ),
          checkboxInput(
            "decreasecpairgttbinary",
            "Intervention probability < control probability",
            value = FALSE
          )
        ),
        
        # cps.irgtt.binary input start
        conditionalPanel(
          "input.type == 'Individually-Randomized Group' & input.dist == 'Binary' & input.meth == 'Simulation'",
          numericInput("nclusters", "Number of Clusters", value = 10),
          numericInput("nsubjects", "Number of Observations (per cluster)", value = 20),
          numericInput(
            "nsimcpsirgttbinary",
            "Number of simulations",
            value = 100,
            max = 500000,
            min = 0
          ),
          numericInput(
            "p1cpsirgttbinary",
            "Outcome proportion (Arm 1)",
            value = 0.1,
            step = 0.001,
            min = 0,
            max = 1
          ),
          numericInput(
            "p2cpsirgttbinary",
            "Outcome proportion (Arm 2)",
            value = 0.5,
            step = 0.001,
            min = 0,
            max = 1
          ),
          numericInput(
            "sigma_b_sqcpsirgttbinary",
            "Between-cluster variance (Arm 1)",
            value = 0.1,
            step = 0.001,
            min = 0
          ),
          numericInput(
            "sigma_b_sq2cpsirgttbinary",
            "Between-cluster variance (Arm 2)",
            value = 0.1,
            step = 0.001,
            min = 0
          )
        ),
        
        # cpa.irgtt.count (no method)
        conditionalPanel(
          "input.type == 'Individually-Randomized Group' & input.dist == 'Count' & input.meth == 'Analytic'",
          HTML("No method exists. Use the simulation option instead.")
        ),
        
        # cps.irgtt.count input start
        conditionalPanel(
          "input.type == 'Individually-Randomized Group' & input.dist == 'Count' & input.meth == 'Simulation'",
          numericInput("nclusters", "Number of Clusters", value = 10),
          numericInput("nsubjects", "Number of Observations (per cluster)", value = 20),
          numericInput(
            "nsimcpsirgttcount",
            "Number of simulations",
            value = 100,
            max = 500000,
            min = 0
          ),
          numericInput(
            "sigma_b_sqcpsirgttcount",
            "Between-cluster variance (Arm 1)",
            value = 0.1,
            step = 0.001,
            min = 0
          ),
          numericInput(
            "sigma_b_sq2cpsirgttcount",
            "Between-cluster variance (Arm 2)",
            value = 0.1,
            step = 0.001,
            min = 0
          ),
          numericInput(
            "c1cpsirgttcount",
            "Expected outcome count (Arm 1)",
            value = 5,
            step = 1,
            min = 0
          ),
          numericInput(
            "c2cpsirgttcount",
            "Expected outcome count (Arm 2)",
            value = 7,
            step = 1,
            min = 0
          )
        )
      ),
      
      #end of values that can be reset with the restore defaults button
      
      actionButton(
        "button",
        "Estimate Power",
        icon = icon("arrow-circle-right"),
        width = '100%',
        style = "color: #fff; background-color: #337ab7; border-color: #2e6da4"
      ),
      actionButton('cancel', 'Cancel'),
      checkboxInput("more", "Show advanced options", value = FALSE),
      conditionalPanel(
        "input.more == true",
        bookmarkButton("Save App State"),
        sliderInput(
          "alpha",
          "Significance level (alpha)",
          value = 0.05,
          min = 0.01,
          max = 0.1,
          step = 0.02
        )
      ),
      conditionalPanel(
        "input.more == true & input.meth == 'Simulation'",
        checkboxInput("timelimitOverride", "Allow unlimited calculation time", value = FALSE),
        checkboxInput("lowPowerOverride", "Allow completion when power is < 0.5", value = FALSE),
        checkboxInput(
          "poorFitOverride",
          "Allow completion when model fit is poor",
          value = FALSE
        ),
        textInput("optmethod", "Specify an optimization method", value = "NLOPT_LN_NELDERMEAD"),
        numericInput(
          "seed",
          "Set the seed (for repeatability)",
          value = NA,
          step = 1
        )
      ),
      conditionalPanel(
        "input.more == true",
        actionButton("restoreDefault", "Restore default parameters",
                     width = '100%',
                     style = "color: #fff; background-color: #337ab7; border-color: #2e6da4"),
        actionButton("reload", "Reset all", icon = icon("trash-alt"))
      )
    ),
 ########################   
    # Tabs start
 ########################   
    mainPanel(tabsetPanel(

      ####  DEBUG ACCESS PANEL START #####      
      
       tabPanel("DEBUG", tableOutput("show_inputs")),
      
      #### DEBUG ACCESS PANEL END #####
      
      tabPanel(
        "Results",
        shinycssloaders::withSpinner(verbatimTextOutput("CRTpower", placeholder = TRUE))
      ),
      tabPanel(
        "Parameters",
        dataTableOutput("tbl"),
        shinyjs::hidden(
          textInput("fxnName", "clusterPower function name", value = "cpa.normal")
        ),
        wellPanel(
          HTML(
            "<p>This table shows the values that the Shiny app passes
                   to the R functions.</p>"
          ),
          uiOutput("helpdetails"),
          HTML(
            "<p>Note: for more advanced features, see the clusterPower R package.</p>"
          )
        )),
      tabPanel(
        "Help",
        HTML(
          "<p>To use the calculator, select the trial type, outcome distribution, and calculation method.
        Then enter values for the quantities that appear below. When complete, select the ESTIMATE POWER button.</p>
          <p>Please contact ken.kleinman@gmail.com with any feedback.</p>"
        ),
        HTML("<h3>Getting started</h3>"),
        HTML(
          "<p>The clusterPower package is intended to perform power calculations for many of the most common
             randomized controlled trial (RCT) designs. This app does not allow the user to access all of the
             functions available in the clusterPower R package, such as calculating the numbers of clusters or
             subjects needed to obtain a specific power, returning the raw simulated datasets, or viewing the
             results of each model fit. For these functions, use clusterPower with the R console rather than
             from within the applet.</p>
             <p>The first and most important step for using this app is to choose the appropriate experimental
             design and outcome distribution for your RCT. For more on this topic, consult the "
        ),
        tags$a(
          "clusterPower vignette.",
          href = get_vignette_link("clusterpower", package = "clusterPower"),
          target = "_blank"
        ),
        HTML(
          "To return to this page, click the back or reload button at the top of your browser window.</p>
          <h4>Choosing a distribution</h4>
          <p>After you have selected the RCT type using the pulldown menu, select the outcome distribution using
             the next pulldown menu. The options here include normal, binary, and count distributions. Use normal
             distribution when your measurement of interest is a numeric value. This can include measurements like
             descriptive variables such as weights, lab results, density, etc. Choose binary distribution if your
             outcome is a yes/no type of response metric, which would be found in studies with outcomes
             represented by exactly two choices (sometimes qualitative), such as survived/deceased, uninfected/infected,
             or participated/withdrew. Choose count when the outcome has more than two possible outcomes, such as
             uninfected/infected/recovered/died.</p>
             <h4>Choosing a method</h4>
        <p>The user can choose the calculation method using the 'Method' dropdown menu. The choices are analytical
        or simulation. Analytical methods have the advantage that calculations are sometimes faster than
             simulated methods, but many make assumptions about variance or balance in design that may sacrifice
             the accuracy of the power estimation. Simulated methods take longer to run but are more flexible because
             they make fewer assumptions. Furthermore, analytical methods don't exist for all RCT types, meaning that in
             those cases the simulation approach may be the user's only option. However, analytical and simulated
             power estimation will likely produce similar results.</p>
             <h4>Parameters</h4>
        <p>Depending on the user's choices as outlined in the previous sections, different parameter entry
             options will appear in the left-hand panel. All of these include the number of observations (or
             subjects, depending on design requirements), and the number of clusters in each arm. Other options
             include variables representing the expected outcomes for each arm, and those representing measures
             of variablity among and within clusters. For simulated methods, the user can also supply a number
             of simulations they would like to use for calculation.</p>
             <p>Here the user may want to consult the 'Parameters' tab option, which displays the verbatim
             parameter values that the app passes to the clusterPower backend as they are entered by the user.
             This table also shows the internal argument names for each parameter.</p>
             <h4>Advanced options</h4>
        <p>The checkbox near the bottom of the left-hand panel opens the Advanced Options for the app. For
             all analysis types, these options include options to bookmark the app state, restore the default
             parameters, or reset (reload) the app. Users can also adjust the significance cutoff value alpha
             using the slider. Typically this value is set to 0.05, which is the default for this app.</p>
             <p>For the simulation methods, there are 3-4 additional controls. These include:</p>
             <p>1) Time limit override: Simulation methods sometimes take a long time to complete. Each
             operation will produce a time estimate based on a few model fits. If the estimate is longer than
             2 minutes, the function will stop and return the estimated time to produce the desired fits.
             Unchecking this option allows the app to run indefinitely, so users can remove the time limit
             constraint if they choose.</p>
             <p>2) Low power override: The simulated methods also have a built-in option to stop fitting
             models if the calculated power is lower than 0.5. Generally, 0.8 is the ideal power target
             for RCT power estimation, so estimations that are very low trigger the alogrithm to stop
             fitting so that the user doesn't have to wait for the run to finish. However, the low power
             error can be overridden by selecting the low power override option.</p>
             <p>3) Poor fit override: When more than 25% of models fail to converge, the default app state
             will stop the calculation with an error, again to prevent the user from waiting a long time
             for uninformative results. Lack of convergence generally is an indication that the data is
             not a good fit for the models. This is ideally addressed by adjusting the model parameters.
             However, if the user wants to allow the procedure to run despite lack of convergence, the
             poor fit override option will override this error. Use this option with caution, as models
             lacking convergence may produce unreliable estimates.</p>
             <p>4) Some simulated methods allow the user to specify an optimizer, which can sometimes
             address convergence issues as an alternative to overriding the poor fit checks or excluding
             non-convergent models.</p>
          <h4>Obtaining results</h4>
          <p>After selecting the desired parameters, submit the job by clicking the Estimate Power button
          near the bottom of the screen. When complete, results will appear on the Results tab. Please
          keep in mind that calculations may take up to 2 minutes to complete, unless the user has chosen
          to override the time limit, in which case the wait time may be longer."
        ),
      )
    )) # Tabs end
  )
)


######################################

#       SERVER

######################################
server <- function(input, output, session) {
  # Reload the app
  observeEvent(input$reload, {
    session$reload()
  })
  
  # Restore default values
  observeEvent(input$restoreDefault, {
    shinyjs::reset("allValues")
  })
  
  #change text input to numeric
  textToNum <- function(x) {
    result <- as.numeric(unlist(strsplit(x, split = ", ")))
    return(result)
  }
  
  #make some helpful fxns to extract arg names
  updateArgs <- function(fxnName) {
    argMatchResult <-
      c(
        clusterPower::argMatch(fxnName, justNames = TRUE),
        "lowPowerOverride",
        "poorFitOverride",
        "timelimitOverride",
        "power",
        "seed",
        "optmethod"
      )
    argNames <-
      c("nsubjects",
        "nclusters",
        "alpha",
        dplyr::intersect(argMatchResult,
                         names(formals(fxnName))))
    arghelper <- function(argname) {
      x <- paste0("input$", argname)
      x <- eval(parse(text = x))
      return(x)
    }
    holder <- list()
    for (i in 1:length(argNames)) {
      holder[[i]] <- arghelper(argNames[i])
    }
    names(holder) <- argNames
    return(holder)
  }
  
  printresult <- function(fxnName) {
    x <- rlang::exec(fxnName, !!!updateArgs(fxnName))
    return(x)
  }
  #end of fxns to extract argument names
  
  #which events to observe
  watchfor <- reactive({
    list(input$dist,input$meth, input$type, input$button)
  }) # end of which events to observe
  
  # update help documentation and params table when function is selected
  #observeEvent(watchfor(), {
  observe({
    if (input$type == 'Parallel' &&
        input$dist == 'Normal' && input$meth == 'Analytic') {
      updateTextInput(session, "fxnName", value = "cpa.normal")
    }
    if (input$type == 'Parallel' &&
        input$dist == 'Normal' && input$meth == 'Simulation') {
      updateTextInput(session, "fxnName", value = "cps.normal")
    }
    if (input$type == 'Parallel' &&
        input$dist == 'Binary' && input$meth == 'Analytic') {
      updateTextInput(session, "fxnName", value = "cpa.binary")
    }
    if (input$type == 'Parallel' &&
        input$dist == 'Binary' && input$meth == 'Simulation') {
      updateTextInput(session, "fxnName", value = "cps.binary")
    }
    if (input$type == 'Parallel' &&
        input$dist == 'Count' && input$meth == 'Analytic') {
      updateTextInput(session, "fxnName", value = "cpa.count")
    }
    if (input$type == 'Parallel' &&
        input$dist == 'Count' && input$meth == 'Simulation') {
      updateTextInput(session, "fxnName", value = "cps.count")
    }
    if (input$type == 'Multi-Arm' &&
        input$dist == 'Normal' && input$meth == 'Analytic') {
      updateTextInput(session, "fxnName", value = "cpa.ma.normal")
    }
    if (input$type == 'Multi-Arm' &&
        input$dist == 'Normal' && input$meth == 'Simulation') {
      updateTextInput(session, "fxnName", value = "cps.ma.normal")
    }
    if (input$type == 'Multi-Arm' &&
        input$dist == 'Binary' && input$meth == 'Analytic') {
      updateTextInput(session, "fxnName", value = "cpa.ma.binary")
    }
    if (input$type == 'Multi-Arm' &&
        input$dist == 'Binary' && input$meth == 'Simulation') {
      updateTextInput(session, "fxnName", value = "cps.ma.binary")
    }
    if (input$type == 'Multi-Arm' &&
        input$dist == 'Count' && input$meth == 'Analytic') {
      updateTextInput(session, "fxnName", value = "cpa.ma.count")
    }
    if (input$type == 'Multi-Arm' &&
        input$dist == 'Count' && input$meth == 'Simulation') {
      updateTextInput(session, "fxnName", value = "cps.ma.count")
    }
    if (input$type == 'Difference-in-Difference' &&
        input$dist == 'Normal' && input$meth == 'Analytic') {
      updateTextInput(session, "fxnName", value = "cpa.did.normal")
    }
    if (input$type == 'Difference-in-Difference' &&
        input$dist == 'Normal' && input$meth == 'Simulation') {
      updateTextInput(session, "fxnName", value = "cps.did.normal")
    }
    if (input$type == 'Difference-in-Difference' &&
        input$dist == 'Binary' && input$meth == 'Analytic') {
      updateTextInput(session, "fxnName", value = "cpa.did.binary")
    }
    if (input$type == 'Difference-in-Difference' &&
        input$dist == 'Binary' && input$meth == 'Simulation') {
      updateTextInput(session, "fxnName", value = "cps.did.binary")
    }
    if (input$type == 'Difference-in-Difference' &&
        input$dist == 'Count' && input$meth == 'Analytic') {
      updateTextInput(session, "fxnName", value = "cpa.did.count")
    }
    if (input$type == 'Difference-in-Difference' &&
        input$dist == 'Count' && input$meth == 'Simulation') {
      updateTextInput(session, "fxnName", value = "cps.did.count")
    }
    if (input$type == 'Stepped Wedge' &&
        input$dist == 'Normal' && input$meth == 'Analytic') {
      updateTextInput(session, "fxnName", value = "cpa.sw.normal")
    }
    if (input$type == 'Stepped Wedge' &&
        input$dist == 'Normal' && input$meth == 'Simulation') {
      updateTextInput(session, "fxnName", value = "cps.sw.normal")
    }
    if (input$type == 'Stepped Wedge' &&
        input$dist == 'Binary' && input$meth == 'Analytic') {
      updateTextInput(session, "fxnName", value = "cpa.sw.binary")
    }
    if (input$type == 'Stepped Wedge' &&
        input$dist == 'Binary' && input$meth == 'Simulation') {
      updateTextInput(session, "fxnName", value = "cps.sw.binary")
    }
    if (input$type == 'Stepped Wedge' &&
        input$dist == 'Count' && input$meth == 'Analytic') {
      updateTextInput(session, "fxnName", value = "cpa.sw.count")
    }
    if (input$type == 'Stepped Wedge' &&
        input$dist == 'Count' && input$meth == 'Simulation') {
      updateTextInput(session, "fxnName", value = "cps.sw.count")
    }
    if (input$type == 'Individually-Randomized Group' &&
        input$dist == 'Normal' && input$meth == 'Analytic') {
      updateTextInput(session, "fxnName", value = "cpa.irgtt.normal")
    }
    if (input$type == 'Individually-Randomized Group' &&
        input$dist == 'Normal' && input$meth == 'Simulation') {
      updateTextInput(session, "fxnName", value = "cps.irgtt.normal")
    }
    if (input$type == 'Individually-Randomized Group' &&
        input$dist == 'Binary' && input$meth == 'Analytic') {
      updateTextInput(session, "fxnName", value = "cpa.irgtt.binary")
    }
    if (input$type == 'Individually-Randomized Group' &&
        input$dist == 'Binary' && input$meth == 'Simulation') {
      updateTextInput(session, "fxnName", value = "cps.irgtt.binary")
    }
    if (input$type == 'Individually-Randomized Group' &&
        input$dist == 'Count' && input$meth == 'Analytic') {
      updateTextInput(session, "fxnName", value = "cpa.irgtt.count")
    }
    if (input$type == 'Individually-Randomized Group' &&
        input$dist == 'Count' && input$meth == 'Simulation') {
      updateTextInput(session, "fxnName", value = "cps.irgtt.count")
    }
  }) # end update help documentation and params table when function is selected
  
  # call the clusterPower functions
  answer <- eventReactive(input$button, {
    if (input$type == 'Parallel' &&
        input$dist == 'Normal' && input$meth == 'Analytic') {
      print(cpa.normal(alpha = input$alpha, power = input$power, 
                       nclusters = input$nclusters, 
                       nsubjects = input$nsubjects, 
                       sigma_sq = input$sigma_sqcpanormal, 
                       sigma_b_sq = input$sigma_b_sqcpanormal, 
                       CV = input$CVcpanormal, d = input$dcpanormal, 
                       ICC = input$ICCcpanormal, vart = input$vartcpanormal))
    }
    if (input$type == 'Parallel' &&
        input$dist == 'Normal' && input$meth == 'Simulation') {
      print(summary(cps.normal(nsim = input$nsimcpsnormal, 
                               nclusters = input$nclusters, 
                               nsubjects = input$nsubjects, 
                               mu = input$mucpsnormal, 
                               mu2 = input$mu2cpsnormal, 
                               ICC = input$ICCcpsnormal, 
                               sigma_sq = input$sigma_sqcpsnormal, 
                               sigma_b_sq = input$sigma_b_sqcpsnormal, 
                               ICC2 = input$ICC2cpsnormal, 
                               sigma_sq2 = input$sigma_sq2cpsnormal, 
                               sigma_b_sq2 = input$sigma_b_sq2cpsnormal, 
                               alpha = input$alpha, 
                               seed = input$seed, 
                               poorFitOverride = input$poorFitOverride, 
                               timelimitOverride = input$timelimitOverride, 
                               lowPowerOverride = input$lowPowerOverride)))
    }
    if (input$type == 'Parallel' &&
        input$dist == 'Binary' && input$meth == 'Analytic') {
      print(
        cpa.binary(
          alpha = input$alpha,
          power = input$power,
          nclusters = input$nclusters,
          nsubjects = input$nsubjects,
          CV = input$CVcpabinary,
          p1 = input$p1cpabinary,
          p2 = input$p2cpabinary,
          ICC = input$ICCcpabinary,
          pooled = input$pooledcpabinary,
          p1inc = input$p1inccpabinary,
          tdist = input$tdistcpabinary
        )
      )
    }
    if (input$type == 'Parallel' &&
        input$dist == 'Binary' && input$meth == 'Simulation') {
      print(summary(
        cps.binary(
          nsim = input$nsimcpsbinary,
          nsubjects = input$nsubjects,
          nclusters = input$nclusters,
          p1 = input$p1cpsbinary,
          p2 = input$p2cpsbinary,
          sigma_b_sq = input$sigma_b_sqcpsbinary,
          sigma_b_sq2 = input$sigma_b_sq2cpsbinary,
          alpha = input$alpha,
          seed = input$seed,
          poorFitOverride = input$poorFitOverride,
          lowPowerOverride = input$lowPowerOverride,
          timelimitOverride = input$timelimitOverride
        )
      ))
    }
    if (input$type == 'Parallel' &&
        input$dist == 'Count' && input$meth == 'Analytic') {
      print(
        cpa.count(
          alpha = input$alpha,
          power = input$power,
          nclusters = input$nclusters,
          nsubjects = input$nsubjects,
          r1 = input$r1cpacount,
          r2 = input$r2cpacount,
          CVB = input$CVBcpacount,
          r1inc = input$r1inccpacount
        )
      )
    }
    if (input$type == 'Parallel' &&
        input$dist == 'Count' && input$meth == 'Simulation') {
      print(summary(
        cps.count(
          nsim = input$nsimcpscount,
          nsubjects = input$nsubjects,
          nclusters = input$nclusters,
          c1 = input$c1cpscount,
          c2 = input$c2cpscount,
          sigma_b_sq = input$sigma_b_sqcpscount,
          sigma_b_sq2 = input$sigma_b_sq2cpscount,
          alpha = input$alpha,
          seed = input$seed,
          poorFitOverride = input$poorFitOverride,
          lowPowerOverride = input$lowPowerOverride,
          timelimitOverride = input$timelimitOverride
        )
      ))
    }
    if (input$type == 'Multi-Arm' &&
        input$dist == 'Normal' && input$meth == 'Analytic') {
      print(
        cpa.ma.normal(
          alpha = input$alpha,
          power = input$power,
          narms = input$narmscpamanormal,
          nclusters = input$nclusters,
          nsubjects = input$nsubjects,
          vara = input$varacpamanormal,
          varc = input$varccpamanormal,
          vare = input$varecpamanormal
        )
      )
    }
    if (input$type == 'Multi-Arm' &&
        input$dist == 'Normal' && input$meth == 'Simulation') {
      print(cps.ma.normal(          
        nsim = input$nsim,
        nsubjects = textToNum(input$nsubjectscpsmanormal),
        narms = input$narmscpsmanormal,
        nclusters = textToNum(input$nclusterscpsmanormal),
        means = textToNum(input$meanscpsmanormal),
        sigma_sq = textToNum(input$sigma_sqcpsmanormal),
        sigma_b_sq = textToNum(input$sigma_b_sqcpsmanormal),
        alpha = input$alpha,
        ICC = textToNum(input$ICCcpsmanormal),
        multi_p_method = input$multi_p_methodcpsmanormal,
        seed = input$seed,
        poorFitOverride = input$poorFitOverride,
        lowPowerOverride = input$lowPowerOverride,
        tdist = input$tdistcpsmanormal,
        optmethod = input$optmethod,
        timelimitOverride = input$timelimitOverride))
    }
    if (input$type == 'Multi-Arm' &&
        input$dist == 'Binary' && input$meth == 'Analytic') {
      print(cpa.ma.binary())
    }
    if (input$type == 'Multi-Arm' &&
        input$dist == 'Binary' && input$meth == 'Simulation') {
      print(
        cps.ma.binary(
          nsim = input$nsimcpsmabinary,
          nsubjects = input$nsubjects,
          narms = input$narmscpsmabinary,
          nclusters = input$nclusters,
          probs = input$probscpsmabinary,
          sigma_b_sq = input$sigma_b_sqcpsmabinary,
          alpha = input$alpha,
          multi_p_method = input$multi_p_methodcpsmabinary,
          seed = input$seed,
          tdist = input$tdistcpsmabinary,
          poorFitOverride = input$poorFitOverride,
          lowPowerOverride = input$lowPowerOverride,
          timelimitOverride = input$timelimitOverride
        ))
    }
    if (input$type == 'Multi-Arm' &&
        input$dist == 'Count' && input$meth == 'Analytic') {
      print(cpa.ma.count())
    }
    if (input$type == 'Multi-Arm' &&
        input$dist == 'Count' && input$meth == 'Simulation') {
      print(
        cps.ma.count(
          nsim = input$nsimcpsmacount,
          nsubjects = input$nsubjects,
          narms = input$narmscpsmacount,
          nclusters = input$nclusters,
          counts = input$countscpsmacount,
          sigma_b_sq = input$sigma_b_sqcpsmacount,
          alpha = input$alpha,
          multi_p_method = input$multi_p_methodcpsmacount,
          seed = input$seed,
          cores = input$corescpsmacount,
          tdist = input$tdistcpsmacount,
          poorFitOverride = input$poorFitOverride,
          lowPowerOverride = input$lowPowerOverride,
          timelimitOverride = input$timelimitOverride
        )
      )
    }
    if (input$type == 'Difference-in-Difference' &&
        input$dist == 'Normal' && input$meth == 'Analytic') {
      print(
        cpa.did.normal(
          alpha = input$alpha,
          power = input$power,
          nclusters = input$nclusters,
          nsubjects = input$nsubjects,
          d = input$dcpadidnormal,
          ICC = input$ICCcpadidnormal,
          rho_c = input$rho_ccpadidnormal,
          rho_s = input$rho_scpadidnormal,
          vart = input$vartcpadidnormal
        )
      )
    }
    if (input$type == 'Difference-in-Difference' &&
        input$dist == 'Normal' && input$meth == 'Simulation') {
      print(summary(
        cps.did.normal(
          nsim = input$nsimcpsdidnormal,
          nsubjects = input$nsubjects,
          nclusters = input$nclusters,
          mu = input$mucpsdidnormal,
          mu2 = input$mu2cpsdidnormal,
          sigma_sq = input$sigma_sqcpsdidnormal,
          sigma_b_sq0 = input$sigma_b_sq0cpsdidnormal,
          sigma_b_sq1 = input$sigma_b_sq1cpsdidnormal,
          alpha = input$alpha,
          poorFitOverride = input$poorFitOverride,
          lowPowerOverride = input$lowPowerOverride,
          timelimitOverride = input$timelimitOverride,
          seed = input$seed
        )
      ))
    }
    if (input$type == 'Difference-in-Difference' &&
        input$dist == 'Binary' && input$meth == 'Analytic') {
      print(
        cpa.did.binary(
          alpha = input$alpha,
          power = input$power,
          nclusters = input$nclusters,
          nsubjects = input$nsubjects,
          p = input$pcpadidbinary,
          d = input$dcpadidbinary,
          ICC = input$ICCcpadidbinary,
          rho_c = input$rho_ccpadidbinary,
          rho_s = input$rho_scpadidbinary
        )
      )
    }
    if (input$type == 'Difference-in-Difference' &&
        input$dist == 'Binary' && input$meth == 'Simulation') {
      print(summary(
        cps.did.binary(
          nsim = input$nsimcpsdidbinary,
          nsubjects = input$nsubjects,
          nclusters = input$nclusters,
          p1 = input$p1cpsdidbinary,
          p2 = input$p2cpsdidbinary,
          sigma_b_sq0 = input$sigma_b_sq0cpsdidbinary,
          sigma_b_sq1 = input$sigma_b_sq1cpsdidbinary,
          alpha = input$alpha,
          poorFitOverride = input$poorFitOverride,
          lowPowerOverride = input$lowPowerOverride,
          timelimitOverride = input$timelimitOverride,
          seed = input$seed
        )
      ))
    }
    if (input$type == 'Difference-in-Difference' &&
        input$dist == 'Count' && input$meth == 'Analytic') {
      print(cpa.did.count())
    }
    if (input$type == 'Difference-in-Difference' &&
        input$dist == 'Count' && input$meth == 'Simulation') {
      print(summary(
        cps.did.count(
          nsim = input$nsimcpsdidcount,
          nsubjects = input$nsubjects,
          nclusters = input$nclusters,
          c1 = input$c1cpsdidcount,
          c2 = input$c2cpsdidcount,
          sigma_b_sq0 = input$sigma_b_sq0cpsdidcount,
          sigma_b_sq1 = input$sigma_b_sq1cpsdidcount,
          alpha = input$alpha,
          poorFitOverride = input$poorFitOverride,
          lowPowerOverride = input$lowPowerOverride,
          timelimitOverride = input$timelimitOverride,
          seed = input$seed
        )
      ))
    }
    if (input$type == 'Stepped Wedge' &&
        input$dist == 'Normal' && input$meth == 'Analytic') {
      print(
        cpa.sw.normal(
          alpha = input$alpha,
          power = input$power,
          nclusters = input$nclusters,
          nsubjects = input$nsubjects,
          ntimes = input$ntimescpaswnormal,
          d = input$dcpaswnormal,
          ICC = input$ICCcpaswnormal,
          rho_c = input$rho_ccpaswnormal,
          rho_s = input$rho_scpaswnormal,
          vart = input$vartcpaswnormal
        )
      )
    }
    if (input$type == 'Stepped Wedge' &&
        input$dist == 'Normal' && input$meth == 'Simulation') {
      print(summary(
        cps.sw.normal(
          alpha = input$alpha,
          power = input$power,
          nclusters = input$nclusters,
          nsubjects = input$nsubjects,
          ntimes = input$ntimescpaswnormal,
          d = input$dcpaswnormal,
          ICC = input$ICCcpaswnormal,
          rho_c = input$rho_ccpaswnormal,
          rho_s = input$rho_scpaswnormal,
          vart = input$vartcpaswnormal
        )
      ))
    }
    if (input$type == 'Stepped Wedge' &&
        input$dist == 'Binary' && input$meth == 'Analytic') {
      print(
        cpa.sw.binary(
          nclusters = input$nclusters,
          steps = input$stepscpaswbinary,
          nsubjects = input$nsubjects,
          d = input$dcpaswbinary,
          ICC = input$ICCcpaswbinary,
          beta = input$betacpaswbinary,
          mu0 = input$mu0cpaswbinary
        )
      )
    }
    if (input$type == 'Stepped Wedge' &&
        input$dist == 'Binary' && input$meth == 'Simulation') {
      print(summary(
        cps.sw.binary(
          nsim = input$nsimcpsswbinary,
          nsubjects = input$nsubjects,
          nclusters = input$nclusters,
          p1 = input$p1cpsswbinary,
          p2 = input$p2cpsswbinary,
          steps = input$stepscpsswbinary,
          sigma_b_sq = input$sigma_b_sqcpsswbinary,
          alpha = input$alpha,
          poorFitOverride = input$poorFitOverride,
          lowPowerOverride = input$lowPowerOverride,
          timelimitOverride = input$timelimitOverride,
          seed = input$seed
        )
      ))
    }
    if (input$type == 'Stepped Wedge' &&
        input$dist == 'Count' && input$meth == 'Analytic') {
      print(
        cpa.sw.count(
          lambda1 = input$lambda1cpaswcount,
          RR = input$RRcpaswcount,
          nclusters = input$nclusters,
          steps = input$stepscpaswcount,
          nsubjects = input$nsubjects,
          ICC = input$ICCcpaswcount,
          alpha = input$alpha
        )
      )
    }
    if (input$type == 'Stepped Wedge' &&
        input$dist == 'Count' && input$meth == 'Simulation') {
      print(summary(
        cps.sw.count(
          nsim = input$nsimcpsswcount,
          nsubjects = input$nsubjects,
          nclusters = input$nclusters,
          c1 = input$c1cpsswcount,
          c2 = input$c2cpsswcount,
          steps = input$stepscpsswcount,
          sigma_b_sq = input$sigma_b_sqcpsswcount,
          alpha = input$alpha,
          poorFitOverride = input$poorFitOverride,
          lowPowerOverride = input$lowPowerOverride,
          timelimitOverride = input$timelimitOverride,
          seed = input$seed
        )
      ))
    }
    if (input$type == 'Individually-Randomized Group' &&
        input$dist == 'Normal' && input$meth == 'Analytic') {
      print(
        cpa.irgtt.normal(
          alpha = input$alpha,
          power = input$power,
          nclusters = input$nclusters,
          nsubjects = input$nsubjects,
          ncontrols = input$ncontrolscpairgttnormal,
          d = input$dcpairgttnormal,
          varu = input$varucpairgttnormal,
          varei = input$vareicpairgttnormal,
          varr = input$varrcpairgttnormal
        )
      )
    }
    if (input$type == 'Individually-Randomized Group' &&
        input$dist == 'Normal' && input$meth == 'Simulation') {
      print(summary(
        cps.irgtt.normal(
          nsim = input$nsimcpsirgttnormal,
          nsubjects = input$nsubjects,
          nclusters = input$nclusters,
          mu = input$mucpsirgttnormal,
          mu2 = input$mu2cpsirgttnormal,
          ICC = input$ICCcpsirgttnormal,
          sigma_sq = input$sigma_sqcpsirgttnormal,
          sigma_b_sq = input$sigma_b_sqcpsirgttnormal,
          ICC2 = input$ICC2cpsirgttnormal,
          sigma_sq2 = input$sigma_sq2cpsirgttnormal,
          sigma_b_sq2 = input$sigma_b_sq2cpsirgttnormal,
          alpha = input$alpha,
          seed = input$seed,
          poorFitOverride = input$poorFitOverride,
          lowPowerOverride = input$lowPowerOverride,
          timelimitOverride = input$timelimitOverride
        )
      ))
    }
    if (input$type == 'Individually-Randomized Group' &&
        input$dist == 'Binary' && input$meth == 'Analytic') {
      print(
        cpa.irgtt.binary(
          alpha = input$alpha,
          power = input$power,
          nclusters = input$nclusters,
          nsubjects = input$nsubjects,
          ncontrols = input$ncontrolscpairgttbinary,
          ICC = input$ICCcpairgttbinary,
          p2 = input$p2cpairgttbinary,
          p1 = input$p1cpairgttbinary,
          decrease = input$decreasecpairgttbinary
        )
      )
    }
    if (input$type == 'Individually-Randomized Group' &&
        input$dist == 'Binary' && input$meth == 'Simulation') {
      print(summary(
        cps.irgtt.binary(
          nsim = input$nsimcpsirgttbinary,
          nsubjects = input$nsubjects,
          nclusters = input$nclusters,
          p1 = input$p1cpsirgttbinary,
          p2 = input$p2cpsirgttbinary,
          sigma_b_sq = input$sigma_b_sqcpsirgttbinary,
          sigma_b_sq2 = input$sigma_b_sq2cpsirgttbinary,
          alpha = input$alpha,
          poorFitOverride = input$poorFitOverride,
          lowPowerOverride = input$lowPowerOverride,
          timelimitOverride = input$timelimitOverride,
          seed = input$seed
        )
      ))
    }
    if (input$type == 'Individually-Randomized Group' &&
        input$dist == 'Count' && input$meth == 'Analytic') {
      print(cpa.irgtt.count())
    }
    if (input$type == 'Individually-Randomized Group' &&
        input$dist == 'Count' && input$meth == 'Simulation') {
      print(summary(
        cps.irgtt.count(
          nsim = input$nsimcpsirgttcount,
          nsubjects = input$nsubjects,
          nclusters = input$nclusters,
          c1 = input$c1cpsirgttcount,
          c2 = input$c2cpsirgttcount,
          sigma_b_sq = input$sigma_b_sqcpsirgttcount,
          sigma_b_sq2 = input$sigma_b_sq2cpsirgttcount,
          alpha = input$alpha,
          poorFitOverride = input$poorFitOverride,
          lowPowerOverride = input$lowPowerOverride,
          timelimitOverride = input$timelimitOverride,
          seed = input$seed
        )
      ))
    }
  }) # end call the clusterPower functions
  ##################
  
  
#########################################
  #### DEBUG ACCESS PANEL #####
#########################################
  
  AllInputs <- reactive({
    x <- reactiveValuesToList(input)
    holder <- data.frame(
      names = names(x),
      values = unlist(x, use.names = FALSE),
      mode = unlist(lapply(x, mode))
    )
    dplyr::filter(holder, grepl(gsub("\\.", "", input$fxnName), names))
  })
  
  output$show_inputs <- renderTable({
    AllInputs()
  })
  
#########################################
  
  output$helpdetails <- renderUI({
    a(
      "Tell me more about this.",
      href = sprintf(
        "http://127.0.0.1:%d/library/clusterPower/html/%s",
        tools::startDynamicHelp(NA),
        paste0(input$fxnName, ".html")
      ),
      target = "_blank"
    )
  })
  
  # create reactive input data table
  args <- reactive({
    t(data.frame(unlist(updateArgs(input$fxnName))))
  })
  observe(args())
  output$tbl <- shiny::renderDataTable(args(),
                                       options = list(
                                         lengthChange = FALSE,
                                         searching = FALSE,
                                         paging = FALSE
                                       ))
  # end create reactive input data table
  
  output$CRTpower <- renderPrint({
    answer()
  })
  
} #end of server fxn



# Run the application
shinyApp(ui = ui, server = server)

