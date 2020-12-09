#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(clusterPower)
library(shinythemes)
library(shiny)
library(shinyBS)
library(promises)
library(future)
library(future.callr)
library(shinyjs)
library(V8)
library(ggplot2)

plan(callr)

# labels for arguments

analyticnsubjectstext <- "Mean observations per cluster (nsubjects)"
analyticnclusterstext <- "Mean clusters per arm (nclusters)"
analyticICCtext <- "Intracluster correlation coefficient (ICC)"
analyticsigma_sqtext <- "Within-cluster variance (sigma_sq)"
analyticsigma_b_sqtext <- "Between-cluster variance (sigma_b_sq)"
treatsigma_sqtext <-
  "Treatment arm within-cluster variance (sigma_sq)"
refsigma_b_sqtext <-
  "Reference arm between-cluster variance (sigma_b_sq)"
treatsigma_b_sqtext <-
  "Treatment arm between-cluster variance (sigma_b_sq)"
simnsimtext <- "Number of simulations (nsim)"
simnsubjectstext <- "Observations per cluster (nsubjects)"
simnclusterstext <- "Clusters per arm (nclusters)"
refmutext <- "Reference arm expected mean (mu)"
treatmutext <- "Treatment arm expected mean (mu2)"
refICCtext <- "Reference arm ICC (ICC)"
treatICCtext <- "Treatment arm ICC (ICC)"
refp1text <- "Reference arm estimated proportion (p1)"
treatp2text <- "Treatment arm estimated proportion (p2)"
refc1text <- "Reference arm expected count (c1)"
treatc2text <- "Treatment arm expected count (c2)"
narmstext <- "Number of arms (narms)"
delimtext <- "Note: comma delimited"
dtext <- "Expected difference in arm means (d)"
didsigma_b_sqposttext <-
  "Post-treatment between-cluster variance (sigma_b_sq1)"
didsigma_b_sqpretext <-
  "Pre-treatment between-cluster variance (sigma_b_sq0)"
ntimestext <- "Number of measurement time points (ntimes)"
varttext <- "Total variation of the outcome (vart)"
rho_ctext <-
  "Baseline and post-test cluster-level correlation (rho_c)"
rho_stext <-
  "Baseline and post-test subject-level correlation (rho_s)"
stepstext <- "Number of crossover steps (steps)"
mu0text <- "Estimated baseline mean (mu0)"
betatext <- "Estimated post-treatment effect (beta)"
mu1text <- "Estimated post-treatment effect (mu1)"
p0text <- "Estimated baseline proportion (p0)"
p1text <- "Estimated post-treatment effect (p1)"
c0text <- "Estimated baseline proportion (c0)"
c1text <- "Estimated post-treatment effect (c1)"
refsigma_sqtext <-
  "Reference arm within-cluster variance (sigma_sq)"
nsubjectsirgttunclusttext <-
  "Observations in unclustered arm (nsubjects)"
ncontrolsirgttunclusttext <-
  "Observations in unclustered arm (ncontrols)"
nclustersswtext <- "Number of clusters (nclusters)"
nsubjectsdelim <- 
  "Note: comma delimited. Length of nsubjects must be 1 (all equal) or equal to the sum of nclusters."



# returns the vignette for the help section link
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
        
        #input for cpa.normal  -------------------------------------------------------------
        conditionalPanel(
          "input.type == 'Parallel' && input.dist == 'Normal' && input.meth == 'Analytic'",
          
          # nclusters
          numericInput("nclusterscpanormal", analyticnclusterstext, value = 10),
          
          #nsubjects
          numericInput("nsubjectscpanormal",
                       analyticnsubjectstext,
                       value = 20),
          
          # CV
          numericInput("CVcpanormal", "Coefficient of variation (CV)", value = 0),
          bsTooltip(
            "CVcpanormal",
            "When CV equals 0, all clusters are the same size.",
            'right',
            options = list(container = "body")
          ),
          
          # d
          numericInput("dcpanormal", dtext, value = 0.43),
          
          # ICC
          numericInput(
            "ICCcpanormal",
            analyticICCtext,
            value = NA,
            step = 0.01,
            min = 0,
            max = 1
          ),
          
          # variance params
          numericInput(
            "sigma_sqcpanormal",
            analyticsigma_sqtext,
            value = 0.01,
            step = 0.001,
            min = 0
          ),
          numericInput(
            "sigma_b_sqcpanormal",
            analyticsigma_b_sqtext,
            step = 0.001,
            value = 0.1,
            min = 0
          ),
          numericInput("vartcpanormal", varttext, value = NA)
        ),
        
        # input for cps.normal
        conditionalPanel(
          "input.type == 'Parallel' && input.dist == 'Normal' && input.meth == 'Simulation'",
          
          # nsim
          numericInput(
            "nsimcpsnormal",
            simnsimtext,
            value = 100,
            max = 500000,
            min = 0
          ),
          
          ## REFERENCE VALUES
          h3(id = "armHead",
             "Reference Arm Parameters"),
          tags$style(HTML("#armHead{color: #337ab7;}")),
          
          # nclusters
          numericInput("nclusters1cpsnormal", simnclusterstext, value = 10),
          
          # nsubjects
          textInput("nsubjects1cpsnormal",
                    simnsubjectstext,
                    value = "20"),
          bsTooltip(
            "nsubjects1cpsnormal",
            nsubjectsdelim,
            'right',
            options = list(container = "body")
          ),
          
          # mu
          numericInput("mucpsnormal", refmutext, value = 2.4),
          
          # variance params
          numericInput(
            "ICCcpsnormal",
            refICCtext,
            value = NA,
            step = 0.01,
            min = 0,
            max = 1
          ),
          bsTooltip(
            "ICCcpsnormal",
            "Intracluster correlation coefficient",
            'right',
            options = list(container = "body")
          ),
          
          numericInput(
            "sigma_sqcpsnormal",
            refsigma_sqtext,
            value = 0.2,
            step = 0.001,
            min = 0
          ),
          numericInput(
            "sigma_b_sqcpsnormal",
            refsigma_b_sqtext,
            value = 0.5,
            step = 0.001,
            min = 0
          ),
          
          ### TREATMENT ARM
          h3(id = "armHead",
             "Treatment Arm Parameters"),
          
          # nclusters
          numericInput("nclusters2cpsnormal", simnclusterstext, value = 10),
          
          # nsubjects
          textInput("nsubjects2cpsnormal",
                    simnsubjectstext,
                    value = "20"),
          bsTooltip(
            "nsubjects2cpsnormal",
            nsubjectsdelim,
            'right',
            options = list(container = "body")
          ),
          
          # mu
          numericInput("mu2cpsnormal", treatmutext, value = 1.5),
          
          # variance params
          numericInput(
            "ICC2cpsnormal",
            treatICCtext,
            value = NA,
            step = 0.01,
            min = 0,
            max = 1
          ),
          bsTooltip(
            "ICC2cpsnormal",
            "Intracluster correlation coefficient",
            'right',
            options = list(container = "body")
          ),
          
          numericInput(
            "sigma_sq2cpsnormal",
            treatsigma_sqtext,
            value = 0.2,
            step = 0.001,
            min = 0
          ),
          
          numericInput(
            "sigma_b_sq2cpsnormal",
            treatsigma_b_sqtext,
            value = 0.5,
            step = 0.001,
            min = 0
          )
        ),
        
        # cpa.binary inputs start
        conditionalPanel(
          "input.type == 'Parallel' && input.dist == 'Binary' && input.meth == 'Analytic'",
          
          # nuclusters
          numericInput("nclusterscpabinary", analyticnclusterstext, value = 10),
          
          # nsubjects
          numericInput("nsubjectscpabinary",
                       analyticnsubjectstext,
                       value = 20),
          
          # CV
          numericInput("CVcpabinary", "Coefficient of variation (CV)", value = 0),
          
          # ICC
          numericInput(
            "ICCcpabinary",
            analyticICCtext,
            value = 0.05,
            step = 0.01,
            min = 0,
            max = 1
          ),
          
          # p1 & p2
          numericInput(
            "p1cpabinary",
            refp1text,
            value = 0.1,
            step = 0.001,
            min = 0,
            max = 1
          ),
          numericInput(
            "p2cpabinary",
            treatp2text,
            value = 0.24,
            step = 0.001,
            min = 0,
            max = 1
          ),
          
          # other selections
          checkboxInput("pooledcpabinary", "Pooled standard error", value = FALSE),
          checkboxInput("p1inccpabinary", "p1 > p2", value = FALSE),
          checkboxInput("tdistcpabinary", "Use t-distribution", value = FALSE)
        ),
        
        # cps.binary inputs start
        conditionalPanel(
          "input.type == 'Parallel' && input.dist == 'Binary' && input.meth == 'Simulation'",
          
          # nsim
          numericInput(
            "nsimcpsbinary",
            simnsimtext,
            value = 100,
            max = 500000,
            min = 0
          ),
          
          ## REFERENCE VALUES
          h3(id = "armHead",
             "Reference Arm Parameters"),
          
          # nclusters
          numericInput("nclusters1cpsbinary", simnclusterstext, value = 10),
          
          # nsubjects
          textInput("nsubjects1cpsbinary",
                    simnsubjectstext,
                    value = "20"),
          bsTooltip(
            "nsubjects1cpsbinary",
            nsubjectsdelim,
            'right',
            options = list(container = "body")
          ),
          
          # p1
          numericInput(
            "p1cpsbinary",
            refp1text,
            value = 0.8,
            step = 0.001,
            min = 0,
            max = 1
          ),
          
          # variance param
          numericInput(
            "sigma_b_sqcpsbinary",
            refsigma_b_sqtext,
            value = 1,
            step = 0.001,
            min = 0
          ),
          
          ## TREATMENT VALUES
          h3(id = "armHead",
             "Treatment Arm Parameters"),
          
          # nclusters
          numericInput("nclusters2cpsbinary", simnclusterstext, value = 10),
          
          # nsubjects
          textInput("nsubjects2cpsbinary",
                    simnsubjectstext,
                    value = "20"),
          bsTooltip(
            "nsubjects2cpsbinary",
            nsubjectsdelim,
            'right',
            options = list(container = "body")
          ),
          
          # p2
          numericInput(
            "p2cpsbinary",
            treatp2text,
            value = 0.5,
            step = 0.001,
            min = 0,
            max = 1
          ),
          
          #variance param
          numericInput(
            "sigma_b_sq2cpsbinary",
            treatsigma_b_sqtext,
            value = 1,
            step = 0.001,
            min = 0
          )
        ),
        
        # cpa.count input starts
        conditionalPanel(
          "input.type == 'Parallel' && input.dist == 'Count' && input.meth == 'Analytic'",
          
          # nclusters
          numericInput("nclusterscpacount", analyticnclusterstext, value = 10),
          
          # nsubjects
          numericInput("nsubjectscpacount",
                       analyticnsubjectstext,
                       value = 20),
          
          # variance param
          numericInput(
            "CVBcpacount",
            "Between-cluster coefficient of variation (CVB)",
            value = 0.01
          ),
          
          # expected outcome
          numericInput("r1cpacount",  "Reference arm mean event rate", value = 0.2),
          numericInput("r2cpacount",  "Treatment arm mean event rate", value = 0.35),
          
          # other params
          checkboxInput(
            "r1inccpacount",
            "Intervention probability < control probability",
            value = FALSE
          )
        ),
        
        # cps.count input starts
        conditionalPanel(
          "input.type == 'Parallel' && input.dist == 'Count' && input.meth == 'Simulation'",
          
          # nsim
          numericInput(
            "nsimcpscount",
            simnsimtext,
            value = 100,
            max = 500000,
            min = 0
          ),
          
          ## REFERENCE VALUES
          h3(id = "armHead",
             "Reference Arm Parameters"),
          
          # nclusters
          numericInput("nclusters1cpscount", simnclusterstext, value = 10),
          
          # nsubjects
          textInput("nsubjects1cpscount",
                    simnsubjectstext,
                    value = "20"),
          bsTooltip(
            "nsubjects1cpscount",
            nsubjectsdelim,
            'right',
            options = list(container = "body")
          ),
          
          # c1
          numericInput(
            "c1cpscount",
            refc1text,
            value = 200,
            step = 1,
            min = 0
          ),
          
          #variance param
          numericInput(
            "sigma_b_sqcpscount",
            refsigma_b_sqtext,
            value = 0.5,
            step = 0.001,
            min = 0
          ),
          
          ## Treatment VALUES
          h3(id = "armHead",
             "Treatment Arm Parameters"),
          
          # nclusters
          numericInput("nclusters2cpscount", simnclusterstext, value = 10),
          
          # nsubjects
          textInput("nsubjects2cpscount",
                    simnsubjectstext,
                    value = "20"),
          bsTooltip(
            "nsubjects2cpscount",
            nsubjectsdelim,
            'right',
            options = list(container = "body")
          ),
          
          # c2
          numericInput(
            "c2cpscount",
            treatc2text,
            value = 80,
            step = 1,
            min = 0
          ),
          
          #variance param
          numericInput(
            "sigma_b_sq2cpscount",
            treatsigma_b_sqtext,
            value = 0.5,
            step = 0.001,
            min = 0
          )
        ),
        
        # cpa.ma.normal input starts
        conditionalPanel(
          "input.type == 'Multi-Arm' && input.dist == 'Normal' && input.meth == 'Analytic'",
          
          #narms
          numericInput(
            "narmscpamanormal",
            narmstext,
            value = 3,
            max = 50,
            min = 2,
            step = 1
          ),
          
          # nclusters
          numericInput("nclusterscpamanormal", analyticnclusterstext, value = 10),
          
          # nsubjects
          numericInput("nsubjectscpamanormal",
                       analyticnsubjectstext,
                       value = 20),
          
          # variance params
          numericInput(
            "varacpamanormal",
            "Between-arm variance (vara)",
            value = 0.02,
            step = 0.001,
            min = 0
          ),
          numericInput(
            "varccpamanormal",
            "Between-cluster variance (varc)",
            value = 0.01,
            step = 0.001,
            min = 0
          ),
          numericInput(
            "varecpamanormal",
            "Within-cluster variance (vare)",
            value = 0.1,
            step = 0.001,
            min = 0
          )
        ),
        
        # cps.ma.normal input start
        conditionalPanel(
          "input.type == 'Multi-Arm' && input.dist == 'Normal' && input.meth == 'Simulation'",
          
          # nsim
          numericInput(
            "nsimcpsmanormal",
            simnsimtext,
            value = 100,
            max = 500000,
            min = 0
          ),
          
          # narms
          sliderInput(
            "narmscpsmanormal",
            narmstext,
            value = 3,
            min = 2,
            max = 10,
            step = 1
          ),
          
          # nclusters
          textInput("nclusterscpsmanormal",
                    simnclusterstext,
                    value = "10, 10, 10"),
          bsTooltip(
            "nclusterscpsmanormal",
            delimtext,
            'right',
            options = list(container = "body")
          ),
          
          # nsubjects
          textInput("nsubjectscpsmanormal",
                    simnsubjectstext,
                    value = "20, 20, 20"),
          bsTooltip(
            "nsubjectscpsmanormal",
            delimtext,
            'right',
            options = list(container = "body")
          ),
          
          # means
          textInput(
            "meanscpsmanormal",
            "Expected absolute treatment effect for each arm (means)",
            "22.0, 21.0, 22.5"
          ),
          bsTooltip(
            "meanscpsmanormal",
            delimtext,
            'right',
            options = list(container = "body")
          ),
          
          # ICC
          textInput("ICCcpsmanormal",
                    analyticICCtext,
                    value = NULL),
          bsTooltip(
            "ICCcpsmanormal",
            delimtext,
            'right',
            options = list(container = "body")
          ),
          
          # variance params
          textInput("sigma_sqcpsmanormal",
                    analyticsigma_sqtext,
                    value = "0.1, 0.1, 0.1"),
          bsTooltip(
            "sigma_sqcpsmanormal",
            delimtext,
            'right',
            options = list(container = "body")
          ),
          
          textInput("sigma_b_sqcpsmanormal",
                    analyticsigma_b_sqtext,
                    value = "0.1, 0.1, 0.1"),
          bsTooltip(
            "sigma_b_sqcpsmanormal",
            delimtext,
            'right',
            options = list(container = "body")
          ),
          
          #other choices
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
          "input.type == 'Multi-Arm' && input.dist == 'Binary' && input.meth == 'Analytic'",
          h4(id = "noMethod", "No method exists.  Use the simulation option instead."),
          tags$style(HTML("#noMethod{color: #d30000;}"))
        ),
        
        # cps.ma.binary input start
        conditionalPanel(
          "input.type == 'Multi-Arm' && input.dist == 'Binary' && input.meth == 'Simulation'",
          
          # nsim
          numericInput(
            "nsimcpsmabinary",
            simnsimtext,
            value = 100,
            max = 500000,
            min = 0
          ),
          
          # narms
          sliderInput(
            "narmscpsmabinary",
            narmstext,
            value = 3,
            min = 2,
            max = 10,
            step = 1
          ),
          
          # nclusters
          textInput("nclusterscpsmabinary",
                    simnclusterstext,
                    value = "10, 10, 10"),
          bsTooltip(
            "nclusterscpsmabinary",
            delimtext,
            'right',
            options = list(container = "body")
          ),
          
          # nsubjects
          textInput("nsubjectscpsmabinary",
                    simnsubjectstext,
                    value = "20, 20, 20"),
          bsTooltip(
            "nsubjectscpsmabinary",
            delimtext,
            'right',
            options = list(container = "body")
          ),
          
          # probs
          textInput(
            "probscpsmabinary",
            "Treatment effect probabilities for each arm (probs)",
            "0.30, 0.4, 0.5"
          ),
          bsTooltip(
            "probscpsmabinary",
            delimtext,
            'right',
            options = list(container = "body")
          ),
          
          #variance params
          textInput("sigma_b_sqcpsmabinary",
                    analyticsigma_b_sqtext,
                    value = "0.1, 0.1, 0.1"),
          bsTooltip(
            "sigma_b_sqcpsmabinary",
            delimtext,
            'right',
            options = list(container = "body")
          ),
          
          # other options
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
          "input.type == 'Multi-Arm' && input.dist == 'Count' && input.meth == 'Analytic'",
          h4(id = "noMethod", "No method exists.  Use the simulation option instead.")
        ),
        
        # cps.ma.count input start
        conditionalPanel(
          "input.type == 'Multi-Arm' && input.dist == 'Count' && input.meth == 'Simulation'",
          
          # nsim
          numericInput(
            "nsimcpsmacount",
            simnsimtext,
            value = 10,
            max = 500000,
            min = 0
          ),
          
          # narms
          sliderInput(
            "narmscpsmacount",
            narmstext,
            value = 3,
            min = 2,
            max = 10,
            step = 1
          ),
          
          # nclusters
          textInput("nclusterscpsmacount",
                    simnclusterstext,
                    value = "10, 10, 10"),
          bsTooltip(
            "nclusterscpsmacount",
            delimtext,
            'right',
            options = list(container = "body")
          ),
          
          # nsubjects
          textInput("nsubjectscpsmacount",
                    simnsubjectstext,
                    value = "20, 20, 20"),
          bsTooltip(
            "nsubjectscpsmacount",
            delimtext,
            'right',
            options = list(container = "body")
          ),
          
          # counts
          textInput(
            "countscpsmacount",
            "Expected count outcome for each arm (counts)",
            "30, 35, 70"
          ),
          bsTooltip(
            "countscpsmacount",
            "Mean event per unit time; comma delimited",
            'right',
            options = list(container = "body")
          ),
          
          # variance params
          textInput("sigma_b_sqcpsmacount",
                    analyticsigma_b_sqtext,
                    value = "1, 1.2, 1.9"),
          bsTooltip(
            "sigma_b_sqcpsmacount",
            delimtext,
            'right',
            options = list(container = "body")
          ),
          
          # other options
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
          "input.type == 'Difference-in-Difference' && input.dist == 'Normal' && input.meth == 'Analytic'",
          
          # nclusters
          numericInput("nclusterscpadidnormal", analyticnclusterstext, value = 10),
          
          # nsubjects
          numericInput("nsubjectscpadidnormal",
                       analyticnsubjectstext,
                       value = 20),
          
          # d
          numericInput("dcpadidnormal", dtext, value = 1.02),
          
          # ICC
          numericInput(
            "ICCcpadidnormal",
            analyticICCtext,
            value = 0.05,
            step = 0.01,
            min = 0,
            max = 1
          ),
          
          # variance params
          numericInput("rho_ccpadidnormal",
                       rho_ctext,
                       value = 0),
          numericInput("rho_scpadidnormal",
                       rho_stext,
                       value = 0),
          numericInput("vartcpadidnormal", varttext, value = 3)
        ),
        
        # cps.did.normal input start
        conditionalPanel(
          "input.type == 'Difference-in-Difference' && input.dist == 'Normal' && input.meth == 'Simulation'",
          
          # nsim
          numericInput(
            "nsimcpsdidnormal",
            simnsimtext,
            value = 100,
            max = 500000,
            min = 0
          ),
          
          ## REFERENCE VALUES
          h3(id = "armHead",
             "Reference Arm Parameters"),
          
          # nclusters
          numericInput("nclusters1cpsdidnormal", simnclusterstext, value = 10),
          
          # nsubjects
          textInput("nsubjects1cpsdidnormal",
                    simnsubjectstext,
                    value = "20"),
          bsTooltip(
            "nsubjects1cpsdidnormal",
            nsubjectsdelim,
            'right',
            options = list(container = "body")
          ),
          
          # mu
          numericInput("mucpsdidnormal", refmutext, value = 2.4),
          
          # variance params
          numericInput(
            "sigma_sqcpsdidnormal",
            refsigma_sqtext,
            value = 0.2,
            step = 0.001,
            min = 0
          ),
          numericInput(
            "sigma_b_sq01cpsdidnormal",
            didsigma_b_sqpretext,
            value = 0.1,
            step = 0.001,
            min = 0
          ),
          numericInput(
            "sigma_b_sq11cpsdidnormal",
            didsigma_b_sqposttext,
            value = 0.1,
            step = 0.001,
            min = 0
          ),
          
          ### TREATMENT ARM
          h3(id = "armHead",
             "Treatment Arm Parameters"),
          
          # nclusters
          numericInput("nclusters2cpsdidnormal", simnclusterstext, value = 10),

          # nsubjects
          textInput("nsubjects2cpsdidnormal",
                    simnsubjectstext,
                    value = "20"),
          bsTooltip(
            "nsubjects2cpsdidnormal",
            nsubjectsdelim,
            'right',
            options = list(container = "body")
          ),
          
          # mu
          numericInput("mu2cpsdidnormal", treatmutext, value = 1.5),
          
          # variance params
          numericInput(
            "sigma_sq2cpsdidnormal",
            treatsigma_sqtext,
            value = 0.2,
            step = 0.001,
            min = 0
          ),
          numericInput(
            "sigma_b_sq02cpsdidnormal",
            didsigma_b_sqpretext,
            value = 0.1,
            step = 0.001,
            min = 0
          ),
          numericInput(
            "sigma_b_sq12cpsdidnormal",
            didsigma_b_sqposttext,
            value = 0.1,
            step = 0.001,
            min = 0
          )
        ),
        
        # cpa.did.binary input start
        conditionalPanel(
          "input.type == 'Difference-in-Difference' && input.dist == 'Binary' && input.meth == 'Analytic'",
          
          # nclusters
          numericInput("nclusterscpadidbinary", analyticnclusterstext, value = 33),
          
          # nsubjects
          numericInput("nsubjectscpadidbinary",
                       analyticnsubjectstext,
                       value = 50),
          
          # d
          numericInput("dcpadidbinary", dtext, value = 0.1),
          
          # p
          numericInput(
            "pcpadidbinary",
            "Mean post-test expected proportion (p)",
            value = 0.5
          ),
          
          # variance parameters
          numericInput(
            "ICCcpadidbinary",
            analyticICCtext,
            value = 0.05,
            step = 0.01,
            min = 0,
            max = 1
          ),
          numericInput("rho_ccpadidbinary",
                       rho_ctext,
                       value = 0.3),
          numericInput("rho_scpadidbinary",
                       rho_stext,
                       value = 0.7)
        ),
        
        # cps.did.binary input start
        conditionalPanel(
          "input.type == 'Difference-in-Difference' && input.dist == 'Binary' && input.meth == 'Simulation'",
          
          # nsim
          numericInput(
            "nsimcpsdidbinary",
            simnsimtext,
            value = 100,
            max = 500000,
            min = 0
          ),
          
          ## REFERENCE VALUES
          h3(id = "armHead",
             "Reference Arm Parameters"),
          
          # nclusters
          numericInput("nclusters1cpsdidbinary", simnclusterstext, value = 10),
          
          # nsubjects
          textInput("nsubjects1cpsdidbinary",
                    simnsubjectstext,
                    value = "20"),
          bsTooltip(
            "nsubjects1cpsdidbinary",
            nsubjectsdelim,
            'right',
            options = list(container = "body")
          ),
          
          # p1
          numericInput(
            "p1cpsdidbinary",
            refp1text,
            value = 0.8,
            step = 0.001,
            min = 0,
            max = 1
          ),
          
          # variance param
          numericInput(
            "sigma_b_sq01cpsdidbinary",
            didsigma_b_sqpretext,
            value = 1,
            step = 0.001,
            min = 0
          ),
          numericInput(
            "sigma_b_sq11cpsdidbinary",
            didsigma_b_sqposttext,
            value = 0,
            step = 0.001,
            min = 0
          ),
          
          ## TREATMENT VALUES
          h3(id = "armHead",
             "Treatment Arm Parameters"),
          
          # nclusters
          numericInput("nclusters2cpsdidbinary", simnclusterstext, value = 10),
          
          # nsubjects
          textInput("nsubjects2cpsdidbinary",
                    simnsubjectstext,
                    value = "20"),
          bsTooltip(
            "nsubjects2cpsdidbinary",
            nsubjectsdelim,
            'right',
            options = list(container = "body")
          ),
          
          # p2
          numericInput(
            "p2cpsdidbinary",
            treatp2text,
            value = 0.5,
            step = 0.001,
            min = 0,
            max = 1
          ),
          
          #variance param
          numericInput(
            "sigma_b_sq02cpsdidbinary",
            didsigma_b_sqpretext,
            value = 1,
            step = 0.001,
            min = 0
          ),
          numericInput(
            "sigma_b_sq12cpsdidbinary",
            didsigma_b_sqposttext,
            value = 0,
            step = 0.001,
            min = 0
          )
        ),
        
        # cpa.did.count (no method)
        conditionalPanel(
          "input.type == 'Difference-in-Difference' && input.dist == 'Count' && input.meth == 'Analytic'",
          h4(id = "noMethod", "No method exists.  Use the simulation option instead.")
        ),
        
        # cps.did.count input start
        conditionalPanel(
          "input.type == 'Difference-in-Difference' && input.dist == 'Count' && input.meth == 'Simulation'",
          
          # nsim
          numericInput(
            "nsimcpscount",
            simnsimtext,
            value = 100,
            max = 500000,
            min = 0
          ),
          
          ## REFERENCE VALUES
          h3(id = "armHead",
             "Reference Arm Parameters"),
          
          # nclusters
          numericInput("nclusters1cpsdidcount", simnclusterstext, value = 10),
          
          # nsubjects
          textInput("nsubjects1cpsdidcount",
                    simnsubjectstext,
                    value = "20"),
          bsTooltip(
            "nsubjects1cpsdidcount",
            nsubjectsdelim,
            'right',
            options = list(container = "body")
          ),
          
          # c1
          numericInput(
            "c1cpsdidcount",
            refc1text,
            value = 200,
            step = 1,
            min = 0
          ),
          
          #variance param
          numericInput(
            "sigma_b_sq01cpsdidcount",
            didsigma_b_sqpretext,
            value = 1,
            step = 0.001,
            min = 0
          ),
          numericInput(
            "sigma_b_sq11cpsdidcount",
            didsigma_b_sqposttext,
            value = 0,
            step = 0.001,
            min = 0
          ),
          
          ## Treatment VALUES
          h3(id = "armHead",
             "Treatment Arm Parameters"),
          
          # nclusters
          numericInput("nclusters2cpsdidcount", simnclusterstext, value = 10),

          # nsubjects
          textInput("nsubjects2cpsdidcount",
                    simnsubjectstext,
                    value = "20"),
          bsTooltip(
            "nsubjects2cpsdidcount",
            nsubjectsdelim,
            'right',
            options = list(container = "body")
          ),
          
          # c2
          numericInput(
            "c2cpsdidcount",
            treatc2text,
            value = 80,
            step = 1,
            min = 0
          ),
          
          #variance param
          numericInput(
            "sigma_b_sq02cpsdidcount",
            didsigma_b_sqpretext,
            value = 1,
            step = 0.001,
            min = 0
          ),
          numericInput(
            "sigma_b_sq12cpsdidcount",
            didsigma_b_sqposttext,
            value = 0,
            step = 0.001,
            min = 0
          ),
        ),
        
        # cpa.sw.normal input start
        conditionalPanel(
          "input.type == 'Stepped Wedge' && input.dist == 'Normal' && input.meth == 'Analytic'",
          
          # nclusters
          numericInput("nclusterscpaswnormal", nclustersswtext, value = 8),
          
          # nsubjects
          numericInput("nsubjectscpaswnormal",
                       analyticnsubjectstext,
                       value = 8),
          
          # ntimes
          numericInput(
            "ntimescpaswnormal",
            ntimestext,
            value = 3,
            min = 0,
            step = 1
          ),
          
          # d
          numericInput("dcpaswnormal", dtext, value = 1.5),
          
          # variance params
          numericInput(
            "ICCcpaswnormal",
            analyticICCtext,
            value = 0.01,
            step = 0.01,
            min = 0,
            max = 1
          ),
          numericInput("rho_ccpaswnormal",
                       rho_ctext,
                       value = 0.8),
          numericInput("rho_scpaswnormal",
                       rho_stext,
                       value = 0),
          numericInput("vartcpaswnormal", varttext, value = 16)
        ),
        
        # cps.sw.normal input start
        conditionalPanel(
          "input.type == 'Stepped Wedge' && input.dist == 'Normal' && input.meth == 'Simulation'",
          
          # nsim
          numericInput(
            "nsimcpsswnormal",
            simnsimtext,
            value = 100,
            max = 500000,
            min = 0
          ),
          
          # steps
          numericInput("stepscpsswnormal", stepstext, value = 3),
          
          # nclusters
          numericInput("nclusterscpsswnormal", nclustersswtext, value = 12),

          # nsubjects
          textInput("nsubjectscpsswnormal",
                    simnsubjectstext,
                    value = "20"),
          bsTooltip(
            "nsubjectscpsswnormal",
            nsubjectsdelim,
            'right',
            options = list(container = "body")
          ),
          
          # expected outcomes
          numericInput("mu0cpsswnormal", mu0text, value = 1.4),
          numericInput("mu1cpsswnormal", mu1text, value = 3.2),
          
          # variance params
          numericInput(
            "sigma_sqcpsswnormal",
            analyticsigma_sqtext,
            value = 0.01,
            step = 0.001,
            min = 0
          ),
          numericInput(
            "sigma_b_sqcpsswnormal",
            analyticsigma_b_sqtext,
            value = 0.1,
            step = 0.001,
            min = 0
          )
        ),
        
        # cpa.sw.binary input start
        conditionalPanel(
          "input.type == 'Stepped Wedge' && input.dist == 'Binary' && input.meth == 'Analytic'",
          
          # steps
          numericInput("stepscpaswbinary", stepstext, value = 2),
          
          #nclusters
          numericInput("nclusterscpaswbinary", nclustersswtext, value = 50),
          
          # nsubjects
          numericInput("nsubjectscpaswbinary",
                       analyticnsubjectstext,
                       value = 100),
          
          # estimated outcomes
          numericInput("dcpaswbinary", dtext, value = -0.75),
          numericInput("mu0cpaswbinary", mu0text, value = 0.2),
          numericInput("betacpaswbinary", betatext, value = 0.4),
          
          # ICC
          numericInput(
            "ICCcpaswbinary",
            analyticICCtext,
            value = 0.01,
            step = 0.01,
            min = 0,
            max = 1
          )
        ),
        
        # cps.sw.binary input start
        conditionalPanel(
          "input.type == 'Stepped Wedge' && input.dist == 'Binary' && input.meth == 'Simulation'",
          
          # nsim
          numericInput(
            "nsimcpsswbinary",
            simnsimtext,
            value = 100,
            max = 500000,
            min = 0
          ),
          
          # steps
          numericInput("stepscpsswbinary", stepstext, value = 3),
          
          # nclusters
          numericInput("nclusterscpsswbinary", nclustersswtext, value = 12),

          # nsubjects
          textInput("nsubjectscpsswbinary",
                    simnsubjectstext,
                    value = "20"),
          bsTooltip(
            "nsubjectscpsswbinary",
            nsubjectsdelim,
            'right',
            options = list(container = "body")
          ),
          
          numericInput(
            "p0cpsswbinary",
            p0text,
            value = 0.1,
            step = 0.001,
            min = 0,
            max = 1
          ),
          numericInput(
            "p1cpsswbinary",
            p1text,
            value = 0.5,
            step = 0.001,
            min = 0,
            max = 1
          ),
          numericInput(
            "sigma_b_sqcpsswbinary",
            analyticsigma_b_sqtext,
            value = 0.1,
            step = 0.001,
            min = 0
          )
        ),
        
        # cpa.sw.count input start
        conditionalPanel(
          "input.type == 'Stepped Wedge' && input.dist == 'Count' && input.meth == 'Analytic'",
          
          #steps
          numericInput("stepscpaswcount", stepstext, value = 3),
          
          # nclusters
          numericInput("nclusterscpaswcount", nclustersswtext, value = 12),
          
          # nsubjects
          numericInput("nsubjectscpaswcount",
                       analyticnsubjectstext,
                       value = 20),
          
          # ICC
          numericInput(
            "ICCcpaswcount",
            analyticICCtext,
            value = 0.05,
            step = 0.01,
            min = 0,
            max = 1
          ),
          
          numericInput(
            "lambda1cpaswcount",
            "Baseline rate for outcome of interest (lambda1)",
            value = 0.65
          ),
          numericInput("RRcpaswcount", "Intervention relative risk (RR)", value = 0.6)
        ),
        
        # cps.sw.count input start
        conditionalPanel(
          "input.type == 'Stepped Wedge' && input.dist == 'Count' && input.meth == 'Simulation'",
          
          # nsim
          numericInput(
            "nsimcpsswcount",
            simnsimtext,
            value = 100,
            max = 500000,
            min = 0
          ),
          
          # steps
          numericInput("stepscpsswcount", stepstext, value = 3),
          
          # nclusters
          numericInput("nclusterscpsswcount", nclustersswtext, value = 12),
          
          # nsubjects
          textInput("nsubjectscpsswcount",
                    simnsubjectstext,
                    value = "20"),
          bsTooltip(
            "nsubjectscpsswcount",
            nsubjectsdelim,
            'right',
            options = list(container = "body")
          ),
          
          #est outcomes
          numericInput(
            "c0cpsswcount",
            c0text,
            value = 5,
            step = 1,
            min = 0
          ),
          numericInput(
            "c1cpsswcount",
            c1text,
            value = 7,
            step = 1,
            min = 0
          ),
          
          # variance params
          numericInput(
            "sigma_b_sqcpsswcount",
            analyticsigma_b_sqtext,
            value = 0.1,
            step = 0.001,
            min = 0
          )
        ),
        
        # cpa.irgtt.normal input start
        conditionalPanel(
          "input.type == 'Individually-Randomized Group' && input.dist == 'Normal' && input.meth == 'Analytic'",
          
          numericInput("dcpairgttnormal", dtext, value = 1.03),
          
          ## unclustered arm VALUES
          h3(id = "armHead",
             "Unclustered Arm Parameters"),
          
          # nsubjects
          numericInput(
            "ncontrolscpairgttnormal",
            ncontrolsirgttunclusttext,
            value = 40,
            step = 1,
            min = 0
          ),
          
          # varr
          numericInput(
            "varrcpairgttnormal",
            "Unclustered arm residual variance (varr)",
            value = 0.2,
            step = 0.001,
            min = 0
          ),
          
          ## Treatment VALUES
          h3(id = "armHead",
             "Clustered Arm Parameters"),
          
          # nclusters
          numericInput("nclusterscpairgttnormal",
                       nclustersswtext,
                       value = 8),
          
          # nsubjects
          numericInput("nsubjectscpairgttnormal",
                       analyticnsubjectstext,
                       value = 15),
          
          # variance params
          numericInput(
            "vareicpairgttnormal",
            "Subject-level random error variance (varei)",
            value = 0.5,
            step = 0.001,
            min = 0
          ),
          
          numericInput(
            "varucpairgttnormal",
            "Cluster-level random effect variance (varu)",
            value = 1,
            step = 0.001,
            min = 0
          )
        ),
        
        # cps.irgtt.normal input start
        conditionalPanel(
          "input.type == 'Individually-Randomized Group' && input.dist == 'Normal' && input.meth == 'Simulation'",
          
          # nsim
          numericInput(
            "nsimcpsirgttnormal",
            simnsimtext,
            value = 100,
            max = 500000,
            min = 0
          ),
          
          ## Unclustered VALUES
          h3(id = "armHead",
             "Unclustered Arm Parameters"),
          
          # nsubjects
          numericInput("nsubjectscpsirgttnormal",
                       nsubjectsirgttunclusttext,
                       value = 100),
          
          # mu
          numericInput("mucpsirgttnormal", "Unclustered arm expected mean (mu)", value = 1.1),
          
          # sigma_sq
          numericInput(
            "sigma_sqcpsirgttnormal",
            analyticsigma_sqtext,
            value = 0.1,
            step = 0.001,
            min = 0
          ),
          
          ## Clustered VALUES
          h3(id = "armHead",
             "Clustered Arm Parameters"),
          
          numericInput("nclusterscpsirgttnormal",
                       simnclusterstext,
                       value = 8),
          
          # nsubjects
          textInput("nsubjects2cpsirgttnormal",
                    simnsubjectstext,
                    value = "10"),
          bsTooltip(
            "nsubjects2cpsirgttnormal",
            nsubjectsdelim,
            'right',
            options = list(container = "body")
          ),
          
          # mu2
          numericInput("mu2cpsirgttnormal", "Clustered arm expected mean (mu2)", value = 1.5),
          
          # variance params
          numericInput(
            "ICC2cpsirgttnormal",
            analyticICCtext,
            value = NA,
            step = 0.01,
            min = 0,
            max = 1
          ),
          numericInput(
            "sigma_sq2cpsirgttnormal",
            analyticsigma_sqtext,
            value = 0.2,
            step = 0.001,
            min = 0
          ),
          numericInput(
            "sigma_b_sq2cpsirgttnormal",
            analyticsigma_b_sqtext,
            value = 0.1,
            step = 0.001,
            min = 0
          )
        ),
        
        # cpa.irgtt.binary input start
        conditionalPanel(
          "input.type == 'Individually-Randomized Group' && input.dist == 'Binary' && input.meth == 'Analytic'",
          
          ## unclustered VALUES
          h3(id = "armHead",
             "Unclustered Arm Parameters"),
          
          # ncontrols
          numericInput(
            "ncontrolscpairgttbinary",
            ncontrolsirgttunclusttext,
            value = 200,
            step = 1,
            min = 0
          ),
          bsTooltip(
            "ncontrolscpairgttbinary",
            "(in the unclustered arm)",
            'right',
            options = list(container = "body")
          ),
          
          # p1
          numericInput(
            "p1cpairgttbinary",
            "Unclustered arm estimated proportion (p1)",
            value = 0.1,
            step = 0.001,
            min = 0,
            max = 1
          ),
          
          ## Clustered VALUES
          h3(id = "armHead",
             "Clustered Arm Parameters"),
          
          # nclusters
          numericInput("nclusterscpairgttbinary", nclustersswtext, value = 10),
          
          # nsubjects
          numericInput("nsubjectscpairgttbinary",
                       analyticnsubjectstext,
                       value = 20),
          
          # p2
          numericInput(
            "p2cpairgttbinary",
            "Clustered arm estimated proportion (p2)",
            value = 0.21,
            step = 0.001,
            min = 0,
            max = 1
          ),
          
          numericInput(
            "ICCcpairgttbinary",
            analyticICCtext,
            value = 0.01,
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
          "input.type == 'Individually-Randomized Group' && input.dist == 'Binary' && input.meth == 'Simulation'",
          
          # nsim
          numericInput(
            "nsimcpsirgttbinary",
            simnsimtext,
            value = 100,
            max = 500000,
            min = 0
          ),
          
          ## Unclustered VALUES
          h3(id = "armHead",
             "Unclustered Arm Parameters"),
          
          # nsubjects
          numericInput("nsubjectscpsirgttbinary",
                       nsubjectsirgttunclusttext,
                       value = 200),
          
          # p1
          numericInput(
            "p1cpsirgttbinary",
            "Unclustered arm estimated proportion (p1)",
            value = 0.1,
            step = 0.001,
            min = 0,
            max = 1
          ),
          
          ## Clustered VALUES
          h3(id = "armHead",
             "Clustered Arm Parameters"),
          
          # nclusters
          numericInput("nclusterscpsirgttbinary",
                       simnclusterstext,
                       value = 10),
          
          # nsubjects
          textInput("nsubjects2cpsirgttbinary",
                    simnsubjectstext,
                    value = "20"),
          bsTooltip(
            "nsubjects2cpsirgttbinary",
            nsubjectsdelim,
            'right',
            options = list(container = "body")
          ),
          
          # p2
          numericInput(
            "p2cpsirgttbinary",
            "Clustered arm estimated proportion (p2)",
            value = 0.21,
            step = 0.001,
            min = 0,
            max = 1
          ),
          
          # variance params
          numericInput(
            "sigma_b_sq2cpsirgttbinary",
            analyticsigma_b_sqtext,
            value = 0.02,
            step = 0.001,
            min = 0
          )
        ),
        
        # cpa.irgtt.count (no method)
        conditionalPanel(
          "input.type == 'Individually-Randomized Group' && input.dist == 'Count' && input.meth == 'Analytic'",
          h4(id = "noMethod", "No method exists.  Use the simulation option instead.")
        ),
        
        # cps.irgtt.count input start
        conditionalPanel(
          "input.type == 'Individually-Randomized Group' && input.dist == 'Count' && input.meth == 'Simulation'",
          
          # nsim
          numericInput(
            "nsimcpsirgttcount",
            simnsimtext,
            value = 100,
            max = 500000,
            min = 0
          ),
          
          ## Unclustered VALUES
          h3(id = "armHead",
             "Unclustered Arm Parameters"),
          
          # nsubjects
          numericInput("nsubjectscpsirgttcount",
                       nsubjectsirgttunclusttext,
                       value = 200),
          
          # c1
          numericInput(
            "c1cpsirgttcount",
            "Unclustered arm estimated count (c1)",
            value = 5,
            step = 1,
            min = 0
          ),
          
          ## Clustered VALUES
          h3(id = "armHead",
             "Clustered Arm Parameters"),
          
          # nclusters
          numericInput("nclusterscpsirgttcount",
                       simnclusterstext,
                       value = 10),

          # nsubjects
          textInput("nsubjects2cpsirgttcount",
                    simnsubjectstext,
                    value = "20"),
          bsTooltip(
            "nsubjects2cpsirgttcount",
            nsubjectsdelim,
            'right',
            options = list(container = "body")
          ),
          
          # c2
          numericInput(
            "c2cpsirgttcount",
            "Clustered arm estimated count (c2)",
            value = 7,
            step = 1,
            min = 0
          ),
          
          # variance param
          numericInput(
            "sigma_b_sq2cpsirgttcount",
            analyticsigma_b_sqtext,
            value = 0.1,
            step = 0.001,
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
        ),
        checkboxInput("verbose", "Show verbose results", value = FALSE),
        checkboxInput("debug", "Show debug/diagnostics tab (advanced)", value = FALSE)
      ),
      conditionalPanel(
        "input.more == true && input.meth == 'Simulation'",
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
        actionButton(
          "restoreDefault",
          "Restore default parameters",
          width = '100%',
          style = "color: #fff; background-color: #337ab7; border-color: #2e6da4"
        ),
        actionButton("reload", "Reset all", icon = icon("trash-alt"))
      )
    ),
    ########################
    # Tabs start
    ########################
    mainPanel(tabsetPanel(
      tabPanel(
        "Results",
        conditionalPanel(
          "input.dismissMsg == false && input.dismissMsg2 == false && input.more == false",
          wellPanel(
            HTML(
              "<p>Note: If the estimated calculation time is longer than 3 minutes,
            the process will not run unless you override the time limit
            under advanced options.</p>"
            ),
            checkboxInput("dismissMsg", "dismiss this message", value = FALSE)
          )
        ),
        conditionalPanel(
          "input.type == 'Stepped Wedge' && input.dist == 'Binary' &&
          input.meth == 'Analytic' && input.stepscpaswbinary > 3 &&
          input.dismissMsgCrossover == false",
          wellPanel(
            HTML(
              "<p>Note: Crossover steps > 3 will substantially increase
              calculation time. </p>"
            ),
            checkboxInput("dismissMsgCrossover", "dismiss this message", value = FALSE)
          )
        ),
        verbatimTextOutput("CRTpower", placeholder = TRUE),
        
        ####  DEBUG ACCESS PANEL START #####
        conditionalPanel(
          "input.debug == true",
          actionButton("browser", "browser"),
          tableOutput("show_inputs")
        )
        
        #### DEBUG ACCESS PANEL END #####
      ),
      tabPanel(
        "Graphs",
        conditionalPanel(
          "input.dismissMsg == false && input.dismissMsg2 == false && input.more == false",
          wellPanel(
            HTML(
              "<p>Note: If the estimated calculation time is longer than 3 minutes,
            the process will not run unless you override the time limit
            under advanced options.</p>"
            ),
            checkboxInput("dismissMsg2", "dismiss this message", value = FALSE)
          )
        ),
        selectInput("axisname",
                    "Y-axis name",
                    choices = c("nclusters", "nsubjects")),
        plotOutput("graphic", click = "click"),
        tableOutput("dp"),
        actionButton("cleargraph", "Clear Data", icon = icon("trash-alt")),
        tags$style(type = 'text/css', "button#cleargraph { margin-top: 100px; }")
      ),
      tabPanel(
        "Parameters",
        h5(id = "helpnote", "Note: This tab shows help for the 
           clusterPower functions. For help with this app, see 
           the Help tab."),
        shinyjs::hidden(
          textInput("fxnName", "clusterPower function name", value = "cpa.normal")
        ),
        wellPanel(
          HTML(
            "<p>This table shows the values that the Shiny app passes
                   to the R functions.</p>
            <p>Note: for more advanced features, see the clusterPower R package.</p>"
          ),
          checkboxInput("showHelp", "Show help documentation", value = FALSE)
        ),
        fluidRow(
          column(width = 6, tableOutput("tracker")),
          column(width = 6),
          conditionalPanel("input.showHelp == true", htmlOutput("helpdetails"))
        ),
        downloadButton(
          "downloadData",
          "Download this table (.csv)",
          icon = icon("file-download"),
          style = "color: #fff; background-color: #337ab7; border-color: #2e6da4"
        ),
        actionButton("cleargraph2", "Clear Data", icon = icon("trash-alt")),
        tags$style(type = 'text/css', "button#cleargraph2 { margin-left: 250px; }")
      ),
      
      tabPanel(
        "Help",
        h5(id = "helpnote", "Note: This tab shows help for this app. 
           For help with the clusterPower functions, see the Parameters tab."),
        tags$style(HTML("#helpnote{color: #337ab7;}")),
        HTML(
          "<p>To use the calculator, select the trial type, outcome distribution, and calculation method.
        Then enter values for the quantities that appear below. When complete, select the ESTIMATE POWER button.</p>"
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
          href = get_vignette_link("clusterPower", package = "clusterPower"),
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
          <h3>Results</h3>
          <p>After selecting the desired parameters, submit the job by clicking the Estimate Power button
          near the bottom of the screen. When complete, results will appear on the Results tab. Please
          keep in mind that calculations may take up to 2 minutes to complete, unless the user has chosen
          to override the time limit, in which case the wait time may be longer. Wait times vary depending
          on the CRT design and complexity of the resulting model, whether the method is analytic or
          simulation, and the amount of RAM available in the host computer. The Results tab shows the power
          calculation resulting from the most recent clusterPower run, although consecutive runs are logged
          unless the cached data is manually cleared, or the CRT type, method, or distribution is
          changed by the user.</p>
          <h3>Graphics</h3>
          <p> As mentioned in the previous section, consecutive runs are logged
          unless the cached data is manually cleared, or the CRT type, method, or distribution is
          changed by the user. On the Graphics tab, the user can graph any of the user-selected
          parameters against the resulting power estimate using the drop-down menu. Exact coordinates
          for each point can be obtained by clicking on the point of interest. If the user would like
          to clear the graph manually, the Clear Data button at the bottom of the panel will clear the
          Results and Parameters tabs to their original state.</p>
          <h3>Parameters</h3>
          <p>On the parameters tab, consecutive clusterPower runs are logged until the data is cleared
          or the user selects a different CRT type, method, or distribution. Parameters are shown
          according to their argument names when passed to the clusterPower function. To learn more
          about each parameter, select the 'Tell Me More About This' link near the top of the panel
          to open a pop-up window showing the documentation for the clusterPower function in use.</p>
          <h3>Note</h3>
          <p>App created by Alexandria Sakrejda, Jon Moyer, and Ken Kleinman; support from NIGMS grant R01GM121370.
          Please contact ken.kleinman@gmail.com with any feedback.</p>"
        )
      )
    )) # Tabs end
  )
)

######################################
#       SERVER
######################################

server <- function(input, output, session) {
  disable("cancel")
  out1 <- reactiveValues(power = NULL)
  logargs <- reactiveValues(tab = NULL)
  
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
      #  x <- eval(parse(text = x))
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
  
  # update help documentation and params table when function is selected
  observe({
    if (input$type == 'Parallel' &&
        input$dist == 'Normal' && input$meth == 'Analytic') {
      updateTextInput(session, "fxnName", value = "cpa.normal")
      shinyjs::show(id = "button")
    }
    if (input$type == 'Parallel' &&
        input$dist == 'Normal' && input$meth == 'Simulation') {
      updateTextInput(session, "fxnName", value = "cps.normal")
      shinyjs::show(id = "button")
    }
    if (input$type == 'Parallel' &&
        input$dist == 'Binary' && input$meth == 'Analytic') {
      updateTextInput(session, "fxnName", value = "cpa.binary")
      shinyjs::show(id = "button")
    }
    if (input$type == 'Parallel' &&
        input$dist == 'Binary' && input$meth == 'Simulation') {
      updateTextInput(session, "fxnName", value = "cps.binary")
      shinyjs::show(id = "button")
    }
    if (input$type == 'Parallel' &&
        input$dist == 'Count' && input$meth == 'Analytic') {
      updateTextInput(session, "fxnName", value = "cpa.count")
      shinyjs::show(id = "button")
    }
    if (input$type == 'Parallel' &&
        input$dist == 'Count' && input$meth == 'Simulation') {
      updateTextInput(session, "fxnName", value = "cps.count")
      shinyjs::show(id = "button")
    }
    if (input$type == 'Multi-Arm' &&
        input$dist == 'Normal' && input$meth == 'Analytic') {
      updateTextInput(session, "fxnName", value = "cpa.ma.normal")
      shinyjs::show(id = "button")
    }
    if (input$type == 'Multi-Arm' &&
        input$dist == 'Normal' && input$meth == 'Simulation') {
      updateTextInput(session, "fxnName", value = "cps.ma.normal")
      shinyjs::show(id = "button")
    }
    if (input$type == 'Multi-Arm' &&
        input$dist == 'Binary' && input$meth == 'Analytic') {
      updateTextInput(session, "fxnName", value = "cpa.ma.binary")
      shinyjs::hide(id = "button")
    }
    if (input$type == 'Multi-Arm' &&
        input$dist == 'Binary' && input$meth == 'Simulation') {
      updateTextInput(session, "fxnName", value = "cps.ma.binary")
      shinyjs::show(id = "button")
    }
    if (input$type == 'Multi-Arm' &&
        input$dist == 'Count' && input$meth == 'Analytic') {
      updateTextInput(session, "fxnName", value = "cpa.ma.count")
      shinyjs::hide(id = "button")
    }
    if (input$type == 'Multi-Arm' &&
        input$dist == 'Count' && input$meth == 'Simulation') {
      updateTextInput(session, "fxnName", value = "cps.ma.count")
      shinyjs::show(id = "button")
    }
    if (input$type == 'Difference-in-Difference' &&
        input$dist == 'Normal' && input$meth == 'Analytic') {
      updateTextInput(session, "fxnName", value = "cpa.did.normal")
      shinyjs::show(id = "button")
    }
    if (input$type == 'Difference-in-Difference' &&
        input$dist == 'Normal' && input$meth == 'Simulation') {
      updateTextInput(session, "fxnName", value = "cps.did.normal")
      shinyjs::show(id = "button")
    }
    if (input$type == 'Difference-in-Difference' &&
        input$dist == 'Binary' && input$meth == 'Analytic') {
      updateTextInput(session, "fxnName", value = "cpa.did.binary")
      shinyjs::show(id = "button")
    }
    if (input$type == 'Difference-in-Difference' &&
        input$dist == 'Binary' && input$meth == 'Simulation') {
      updateTextInput(session, "fxnName", value = "cps.did.binary")
      shinyjs::show(id = "button")
    }
    if (input$type == 'Difference-in-Difference' &&
        input$dist == 'Count' && input$meth == 'Analytic') {
      updateTextInput(session, "fxnName", value = "cpa.did.count")
      shinyjs::hide(id = "button")
    }
    if (input$type == 'Difference-in-Difference' &&
        input$dist == 'Count' && input$meth == 'Simulation') {
      updateTextInput(session, "fxnName", value = "cps.did.count")
      shinyjs::show(id = "button")
    }
    if (input$type == 'Stepped Wedge' &&
        input$dist == 'Normal' && input$meth == 'Analytic') {
      updateTextInput(session, "fxnName", value = "cpa.sw.normal")
      shinyjs::show(id = "button")
    }
    if (input$type == 'Stepped Wedge' &&
        input$dist == 'Normal' && input$meth == 'Simulation') {
      updateTextInput(session, "fxnName", value = "cps.sw.normal")
      shinyjs::show(id = "button")
    }
    if (input$type == 'Stepped Wedge' &&
        input$dist == 'Binary' && input$meth == 'Analytic') {
      updateTextInput(session, "fxnName", value = "cpa.sw.binary")
      shinyjs::show(id = "button")
    }
    if (input$type == 'Stepped Wedge' &&
        input$dist == 'Binary' && input$meth == 'Simulation') {
      updateTextInput(session, "fxnName", value = "cps.sw.binary")
      shinyjs::show(id = "button")
    }
    if (input$type == 'Stepped Wedge' &&
        input$dist == 'Count' && input$meth == 'Analytic') {
      updateTextInput(session, "fxnName", value = "cpa.sw.count")
      shinyjs::show(id = "button")
    }
    if (input$type == 'Stepped Wedge' &&
        input$dist == 'Count' && input$meth == 'Simulation') {
      updateTextInput(session, "fxnName", value = "cps.sw.count")
      shinyjs::show(id = "button")
    }
    if (input$type == 'Individually-Randomized Group' &&
        input$dist == 'Normal' && input$meth == 'Analytic') {
      updateTextInput(session, "fxnName", value = "cpa.irgtt.normal")
      shinyjs::show(id = "button")
    }
    if (input$type == 'Individually-Randomized Group' &&
        input$dist == 'Normal' && input$meth == 'Simulation') {
      updateTextInput(session, "fxnName", value = "cps.irgtt.normal")
      shinyjs::show(id = "button")
    }
    if (input$type == 'Individually-Randomized Group' &&
        input$dist == 'Binary' && input$meth == 'Analytic') {
      updateTextInput(session, "fxnName", value = "cpa.irgtt.binary")
      shinyjs::show(id = "button")
    }
    if (input$type == 'Individually-Randomized Group' &&
        input$dist == 'Binary' && input$meth == 'Simulation') {
      updateTextInput(session, "fxnName", value = "cps.irgtt.binary")
      shinyjs::show(id = "button")
    }
    if (input$type == 'Individually-Randomized Group' &&
        input$dist == 'Count' && input$meth == 'Analytic') {
      updateTextInput(session, "fxnName", value = "cpa.irgtt.count")
      shinyjs::hide(id = "button")
    }
    if (input$type == 'Individually-Randomized Group' &&
        input$dist == 'Count' && input$meth == 'Simulation') {
      updateTextInput(session, "fxnName", value = "cps.irgtt.count")
      shinyjs::show(id = "button")
    }
  }) # end update help documentation and params table when function is selected
  
  # call the clusterPower functions
  observeEvent(input$button, {
    disable("button")
    enable("cancel")
    prog <- Progress$new(session)
    prog$set(message = "Analysis in progress",
             value = NULL)
    isolate({
      q <- reactiveValuesToList(input)
    })
    
    if (input$type == 'Parallel' &&
        input$dist == 'Normal' && input$meth == 'Analytic') {
      out1$power <- NULL
      out1$power$power <- "Calculating..."
      answer <<- future({
        val <-
          cpa.normal(
            alpha = q$alpha,
            power = q$power,
            nclusters = q$nclusterscpanormal,
            nsubjects = q$nsubjectscpanormal,
            sigma_sq = q$sigma_sqcpanormal,
            sigma_b_sq = q$sigma_b_sqcpanormal,
            CV = q$CVcpanormal,
            d = q$dcpanormal,
            ICC = q$ICCcpanormal,
            vart = q$vartcpanormal
          )
        return(val)
      }, seed = TRUE)
    }
    if (input$type == 'Parallel' &&
        input$dist == 'Normal' && input$meth == 'Simulation') {
      out1$power <- NULL
      out1$power$power <- "Calculating..."
      answer <<- future({
        val <- cps.normal(
          nsim = q$nsimcpsnormal,
          nclusters = c(q$nclusters1cpsnormal, q$nclusters2cpsnormal),
          nsubjects = textToNum(c(q$nsubjects1cpsnormal, q$nsubjects2cpsnormal)),
          mu = q$mucpsnormal,
          mu2 = q$mu2cpsnormal,
          ICC = q$ICCcpsnormal,
          sigma_sq = q$sigma_sqcpsnormal,
          sigma_b_sq = q$sigma_b_sqcpsnormal,
          ICC2 = q$ICC2cpsnormal,
          sigma_sq2 = q$sigma_sq2cpsnormal,
          sigma_b_sq2 = q$sigma_b_sq2cpsnormal,
          alpha = q$alpha,
          seed = q$seed,
          poorFitOverride = q$poorFitOverride,
          timelimitOverride = q$timelimitOverride,
          lowPowerOverride = q$lowPowerOverride
        )
        return(val)
      }, seed = TRUE)
    }
    if (input$type == 'Parallel' &&
        input$dist == 'Binary' && input$meth == 'Analytic') {
      out1$power <- NULL
      out1$power$power <- "Calculating..."
      answer <<- future({
        val <- cpa.binary(
          alpha = q$alpha,
          power = q$power,
          nclusters = q$nclusterscpabinary,
          nsubjects = q$nsubjectscpabinary,
          CV = q$CVcpabinary,
          p1 = q$p1cpabinary,
          p2 = q$p2cpabinary,
          ICC = q$ICCcpabinary,
          pooled = q$pooledcpabinary,
          p1inc = q$p1inccpabinary,
          tdist = q$tdistcpabinary
        )
        return(val)
      }, seed = TRUE)
    }
    if (input$type == 'Parallel' &&
        input$dist == 'Binary' && input$meth == 'Simulation') {
      out1$power <- NULL
      out1$power$power <- "Calculating..."
      answer <<- future({
        val <- cps.binary(
          nsim = q$nsimcpsbinary,
          nsubjects = textToNum(c(q$nsubjects1cpsbinary, q$nsubjects2cpsbinary)),
          nclusters = c(q$nclusters1cpsbinary, q$nclusters2cpsbinary),
          p1 = q$p1cpsbinary,
          p2 = q$p2cpsbinary,
          sigma_b_sq = q$sigma_b_sqcpsbinary,
          sigma_b_sq2 = q$sigma_b_sq2cpsbinary,
          alpha = q$alpha,
          seed = q$seed,
          poorFitOverride = q$poorFitOverride,
          lowPowerOverride = q$lowPowerOverride,
          timelimitOverride = q$timelimitOverride
        )
        return(val)
      }, seed = TRUE)
    }
    if (input$type == 'Parallel' &&
        input$dist == 'Count' && input$meth == 'Analytic') {
      out1$power <- NULL
      out1$power$power <- "Calculating..."
      answer <<- future({
        val <- cpa.count(
          alpha = q$alpha,
          power = q$power,
          nclusters = q$nclusterscpacount,
          nsubjects = q$nsubjectscpacount,
          r1 = q$r1cpacount,
          r2 = q$r2cpacount,
          CVB = q$CVBcpacount,
          r1inc = q$r1inccpacount
        )
        return(val)
      }, seed = TRUE)
    }
    if (input$type == 'Parallel' &&
        input$dist == 'Count' && input$meth == 'Simulation') {
      out1$power <- NULL
      out1$power$power <- "Calculating..."
      answer <<- future({
        val <- cps.count(
          nsim = q$nsimcpscount,
          nsubjects = textToNum(c(q$nsubjects1cpscount, q$nsubjects2cpscount)),
          nclusters = c(q$nclusters1cpscount, q$nclusters2cpscount),
          c1 = q$c1cpscount,
          c2 = q$c2cpscount,
          sigma_b_sq = q$sigma_b_sqcpscount,
          sigma_b_sq2 = q$sigma_b_sq2cpscount,
          alpha = q$alpha,
          seed = q$seed,
          poorFitOverride = q$poorFitOverride,
          lowPowerOverride = q$lowPowerOverride,
          timelimitOverride = q$timelimitOverride
        )
        return(val)
      }, seed = TRUE)
    }
    if (input$type == 'Multi-Arm' &&
        input$dist == 'Normal' && input$meth == 'Analytic') {
      out1$power <- NULL
      out1$power$power <- "Calculating..."
      answer <<- future({
        val <- cpa.ma.normal(
          alpha = q$alpha,
          power = q$power,
          narms = q$narmscpamanormal,
          nclusters = q$nclusterscpamanormal,
          nsubjects = q$nsubjectscpamanormal,
          vara = q$varacpamanormal,
          varc = q$varccpamanormal,
          vare = q$varecpamanormal
        )
        return(val)
      }, seed = TRUE)
    }
    if (input$type == 'Multi-Arm' &&
        input$dist == 'Normal' && input$meth == 'Simulation') {
      out1$power <- NULL
      out1$power$power <- "Calculating..."
      answer <<- future({
        val <- cps.ma.normal(
          nsim = q$nsimcpsmanormal,
          nsubjects = textToNum(q$nsubjectscpsmanormal),
          narms = q$narmscpsmanormal,
          nclusters = textToNum(q$nclusterscpsmanormal),
          means = textToNum(q$meanscpsmanormal),
          sigma_sq = textToNum(q$sigma_sqcpsmanormal),
          sigma_b_sq = textToNum(q$sigma_b_sqcpsmanormal),
          alpha = q$alpha,
          ICC = textToNum(q$ICCcpsmanormal),
          multi_p_method = q$multi_p_methodcpsmanormal,
          seed = q$seed,
          poorFitOverride = q$poorFitOverride,
          lowPowerOverride = q$lowPowerOverride,
          tdist = q$tdistcpsmanormal,
          optmethod = q$optmethod,
          timelimitOverride = q$timelimitOverride
        )
        return(val)
      }, seed = TRUE)
    }
    if (input$type == 'Multi-Arm' &&
        input$dist == 'Binary' && input$meth == 'Analytic') {
      out1$power <- NULL
      answer <<- future({
        val <- cpa.ma.binary()
        return(val)
      }, seed = TRUE)
    }
    if (input$type == 'Multi-Arm' &&
        input$dist == 'Binary' && input$meth == 'Simulation') {
      out1$power <- NULL
      out1$power$power <- "Calculating..."
      answer <<- future({
        val <- cps.ma.binary(
          nsim = q$nsimcpsmabinary,
          nsubjects = textToNum(q$nsubjectscpsmabinary),
          narms = q$narmscpsmabinary,
          nclusters = textToNum(q$nclusterscpsmabinary),
          probs = textToNum(q$probscpsmabinary),
          sigma_b_sq = textToNum(q$sigma_b_sqcpsmabinary),
          alpha = q$alpha,
          multi_p_method = q$multi_p_methodcpsmabinary,
          seed = q$seed,
          tdist = q$tdistcpsmabinary,
          poorFitOverride = q$poorFitOverride,
          lowPowerOverride = q$lowPowerOverride,
          timelimitOverride = q$timelimitOverride
        )
        return(val)
      }, seed = TRUE)
    }
    if (input$type == 'Multi-Arm' &&
        input$dist == 'Count' && input$meth == 'Analytic') {
      out1$power <- NULL
      answer <<- future({
        val <- cpa.ma.count()
        return(val)
      }, seed = TRUE)
    }
    if (input$type == 'Multi-Arm' &&
        input$dist == 'Count' && input$meth == 'Simulation') {
      out1$power <- NULL
      out1$power$power <- "Calculating..."
      answer <<- future({
        val <- cps.ma.count(
          nsim = q$nsimcpsmacount,
          nsubjects = textToNum(q$nsubjectscpsmacount),
          narms = q$narmscpsmacount,
          nclusters = textToNum(q$nclusterscpsmacount),
          counts = textToNum(q$countscpsmacount),
          sigma_b_sq = textToNum(q$sigma_b_sqcpsmacount),
          alpha = q$alpha,
          multi_p_method = q$multi_p_methodcpsmacount,
          seed = q$seed,
          tdist = q$tdistcpsmacount,
          poorFitOverride = q$poorFitOverride,
          lowPowerOverride = q$lowPowerOverride,
          timelimitOverride = q$timelimitOverride
        )
        return(val)
      }, seed = TRUE)
    }
    if (input$type == 'Difference-in-Difference' &&
        input$dist == 'Normal' && input$meth == 'Analytic') {
      out1$power <- NULL
      out1$power$power <- "Calculating..."
      answer <<- future({
        val <- cpa.did.normal(
          alpha = q$alpha,
          power = q$power,
          nclusters = q$nclusterscpadidnormal,
          nsubjects = q$nsubjectscpadidnormal,
          d = q$dcpadidnormal,
          ICC = q$ICCcpadidnormal,
          rho_c = q$rho_ccpadidnormal,
          rho_s = q$rho_scpadidnormal,
          vart = q$vartcpadidnormal
        )
        return(val)
      }, seed = TRUE)
    }
    if (input$type == 'Difference-in-Difference' &&
        input$dist == 'Normal' && input$meth == 'Simulation') {
      out1$power <- NULL
      out1$power$power <- "Calculating..."
      answer <<- future({
        val <- cps.did.normal(
          nsim = q$nsimcpsdidnormal,
          nsubjects = textToNum(c(q$nsubjects1cpsdidnormal, q$nsubjects2cpsdidnormal)),
          nclusters = c(q$nclusters1cpsdidnormal, q$nclusters2cpsdidnormal),
          mu = q$mucpsdidnormal,
          mu2 = q$mu2cpsdidnormal,
          sigma_sq = q$sigma_sqcpsdidnormal,
          sigma_b_sq0 = c(
            q$sigma_b_sq01cpsdidnormal,
            q$sigma_b_sq02cpsdidnormal
          ),
          sigma_b_sq1 = c(
            q$sigma_b_sq11cpsdidnormal,
            q$sigma_b_sq12cpsdidnormal
          ),
          alpha = q$alpha,
          poorFitOverride = q$poorFitOverride,
          lowPowerOverride = q$lowPowerOverride,
          timelimitOverride = q$timelimitOverride,
          seed = q$seed
        )
        return(val)
      }, seed = TRUE)
    }
    if (input$type == 'Difference-in-Difference' &&
        input$dist == 'Binary' && input$meth == 'Analytic') {
      out1$power <- NULL
      out1$power$power <- "Calculating..."
      answer <<- future({
        val <- cpa.did.binary(
          alpha = q$alpha,
          power = q$power,
          nclusters = q$nclusterscpadidbinary,
          nsubjects = q$nsubjectscpadidbinary,
          p = q$pcpadidbinary,
          d = q$dcpadidbinary,
          ICC = q$ICCcpadidbinary,
          rho_c = q$rho_ccpadidbinary,
          rho_s = q$rho_scpadidbinary
        )
        return(val)
      }, seed = TRUE)
    }
    if (input$type == 'Difference-in-Difference' &&
        input$dist == 'Binary' && input$meth == 'Simulation') {
      out1$power <- NULL
      out1$power$power <- "Calculating..."
      answer <<- future({
        val <- cps.did.binary(
          nsim = q$nsimcpsdidbinary,
          nsubjects = textToNum(c(q$nsubjects1cpsdidbinary, q$nsubjects2cpsdidbinary)),
          nclusters = c(q$nclusters1cpsdidbinary, q$nclusters2cpsdidbinary),
          p1 = q$p1cpsdidbinary,
          p2 = q$p2cpsdidbinary,
          sigma_b_sq0 = c(
            q$sigma_b_sq01cpsdidbinary,
            q$sigma_b_sq02cpsdidbinary
          ),
          sigma_b_sq1 = c(
            q$sigma_b_sq11cpsdidbinary,
            q$sigma_b_sq12cpsdidbinary
          ),
          alpha = q$alpha,
          poorFitOverride = q$poorFitOverride,
          lowPowerOverride = q$lowPowerOverride,
          timelimitOverride = q$timelimitOverride,
          seed = q$seed
        )
        return(val)
      }, seed = TRUE)
    }
    if (input$type == 'Difference-in-Difference' &&
        input$dist == 'Count' && input$meth == 'Analytic') {
      out1$power <- NULL
      answer <<- future({
        val <- cpa.did.count()
        return(val)
      }, seed = TRUE)
    }
    if (input$type == 'Difference-in-Difference' &&
        input$dist == 'Count' && input$meth == 'Simulation') {
      out1$power <- NULL
      out1$power$power <- "Calculating..."
      answer <<- future({
        val <- cps.did.count(
          nsim = q$nsimcpsdidcount,
          nsubjects = textToNum(c(q$nsubjects1cpsdidcount, q$nsubjects2cpsdidcount)),
          nclusters = c(q$nclusters1cpsdidcount, q$nclusters2cpsdidcount),
          c1 = q$c1cpsdidcount,
          c2 = q$c2cpsdidcount,
          sigma_b_sq0 = c(
            q$sigma_b_sq01cpsdidcount,
            q$sigma_b_sq02cpsdidcount
          ),
          sigma_b_sq1 = c(
            q$sigma_b_sq11cpsdidcount,
            q$sigma_b_sq12cpsdidcount
          ),
          alpha = q$alpha,
          poorFitOverride = q$poorFitOverride,
          lowPowerOverride = q$lowPowerOverride,
          timelimitOverride = q$timelimitOverride,
          seed = q$seed
        )
        return(val)
      }, seed = TRUE)
    }
    if (input$type == 'Stepped Wedge' &&
        input$dist == 'Normal' && input$meth == 'Analytic') {
      out1$power <- NULL
      out1$power$power <- "Calculating..."
      answer <<- future({
        val <- cpa.sw.normal(
          alpha = q$alpha,
          power = q$power,
          nclusters = q$nclusterscpaswnormal,
          nsubjects = q$nsubjectscpaswnormal,
          ntimes = q$ntimescpaswnormal,
          d = q$dcpaswnormal,
          ICC = q$ICCcpaswnormal,
          rho_c = q$rho_ccpaswnormal,
          rho_s = q$rho_scpaswnormal,
          vart = q$vartcpaswnormal
        )
        return(val)
      }, seed = TRUE)
    }
    if (input$type == 'Stepped Wedge' &&
        input$dist == 'Normal' && input$meth == 'Simulation') {
      out1$power <- NULL
      out1$power$power <- "Calculating..."
      answer <<- future({
        val <- cps.sw.normal(
          alpha = q$alpha,
          nsim = q$nsimcpsswnormal,
          nclusters = q$nclusterscpsswnormal,
          nsubjects = textToNum(q$nsubjectscpsswnormal),
          steps = q$stepscpsswnormal,
          mu0 = q$mu0cpsswnormal,
          mu1 = q$mu1cpsswnormal,
          sigma_sq = q$sigma_sqcpsswnormal,
          sigma_b_sq = q$sigma_b_sqcpsswnormal,
          poorFitOverride = q$poorFitOverride,
          lowPowerOverride = q$lowPowerOverride,
          timelimitOverride = q$timelimitOverride,
          seed = q$seed
        )
        return(val)
      }, seed = TRUE)
    }
    if (input$type == 'Stepped Wedge' &&
        input$dist == 'Binary' && input$meth == 'Analytic') {
      out1$power <- NULL
      out1$power$power <- "Calculating..."
      answer <<- future({
        val <- cpa.sw.binary(
          alpha = q$alpha,
          nclusters = q$nclusterscpaswbinary,
          steps = q$stepscpaswbinary,
          nsubjects = q$nsubjectscpaswbinary,
          d = q$dcpaswbinary,
          ICC = q$ICCcpaswbinary,
          beta = q$betacpaswbinary,
          mu0 = q$mu0cpaswbinary
        )
        return(val)
      }, seed = TRUE)
    }
    if (input$type == 'Stepped Wedge' &&
        input$dist == 'Binary' && input$meth == 'Simulation') {
      out1$power <- NULL
      out1$power$power <- "Calculating..."
      answer <<- future({
        val <- cps.sw.binary(
          nsim = q$nsimcpsswbinary,
          nsubjects = textToNum(q$nsubjectscpsswbinary),
          nclusters = q$nclusterscpsswbinary,
          p0 = q$p0cpsswbinary,
          p1 = q$p1cpsswbinary,
          steps = q$stepscpsswbinary,
          sigma_b_sq = q$sigma_b_sqcpsswbinary,
          alpha = q$alpha,
          poorFitOverride = q$poorFitOverride,
          lowPowerOverride = q$lowPowerOverride,
          timelimitOverride = q$timelimitOverride,
          seed = q$seed
        )
        return(val)
      }, seed = TRUE)
    }
    if (input$type == 'Stepped Wedge' &&
        input$dist == 'Count' && input$meth == 'Analytic') {
      out1$power <- NULL
      out1$power$power <- "Calculating..."
      answer <<- future({
        val <- cpa.sw.count(
          lambda1 = q$lambda1cpaswcount,
          RR = q$RRcpaswcount,
          nclusters = q$nclusterscpaswcount,
          steps = q$stepscpaswcount,
          nsubjects = q$nsubjectscpaswcount,
          ICC = q$ICCcpaswcount,
          alpha = q$alpha
        )
        return(val)
      }, seed = TRUE)
    }
    if (input$type == 'Stepped Wedge' &&
        input$dist == 'Count' && input$meth == 'Simulation') {
      out1$power <- NULL
      out1$power$power <- "Calculating..."
      answer <<- future({
        val <- cps.sw.count(
          nsim = q$nsimcpsswcount,
          nsubjects = textToNum(q$nsubjectscpsswcount),
          nclusters = q$nclusterscpsswcount,
          c0 = q$c0cpsswcount,
          c1 = q$c1cpsswcount,
          steps = q$stepscpsswcount,
          sigma_b_sq = q$sigma_b_sqcpsswcount,
          alpha = q$alpha,
          poorFitOverride = q$poorFitOverride,
          lowPowerOverride = q$lowPowerOverride,
          timelimitOverride = q$timelimitOverride,
          seed = q$seed
        )
        return(val)
      }, seed = TRUE)
    }
    if (input$type == 'Individually-Randomized Group' &&
        input$dist == 'Normal' && input$meth == 'Analytic') {
      out1$power <- NULL
      out1$power$power <- "Calculating..."
      answer <<- future({
        val <- cpa.irgtt.normal(
          alpha = q$alpha,
          power = q$power,
          nclusters = q$nclusterscpairgttnormal,
          nsubjects = q$nsubjectscpairgttnormal,
          ncontrols = q$ncontrolscpairgttnormal,
          d = q$dcpairgttnormal,
          varu = q$varucpairgttnormal,
          varei = q$vareicpairgttnormal,
          varr = q$varrcpairgttnormal
        )
        return(val)
      }, seed = TRUE)
    }
    if (input$type == 'Individually-Randomized Group' &&
        input$dist == 'Normal' && input$meth == 'Simulation') {
      out1$power <- NULL
      out1$power$power <- "Calculating..."
      answer <<- future({
        val <- cps.irgtt.normal(
          nsim = q$nsimcpsirgttnormal,
          nsubjects = textToNum(c(
            q$nsubjectscpsirgttnormal,
            q$nsubjects2cpsirgttnormal
          )),
          nclusters = q$nclusterscpsirgttnormal,
          mu = q$mucpsirgttnormal,
          mu2 = q$mu2cpsirgttnormal,
          sigma_sq = q$sigma_sqcpsirgttnormal,
          ICC2 = q$ICC2cpsirgttnormal,
          sigma_sq2 = q$sigma_sq2cpsirgttnormal,
          sigma_b_sq2 = q$sigma_b_sq2cpsirgttnormal,
          alpha = q$alpha,
          seed = q$seed,
          poorFitOverride = q$poorFitOverride,
          lowPowerOverride = q$lowPowerOverride,
          timelimitOverride = q$timelimitOverride
        )
        return(val)
      }, seed = TRUE)
    }
    if (input$type == 'Individually-Randomized Group' &&
        input$dist == 'Binary' && input$meth == 'Analytic') {
      out1$power <- NULL
      out1$power$power <- "Calculating..."
      answer <<- future({
        val <- cpa.irgtt.binary(
          alpha = q$alpha,
          power = q$power,
          nclusters = q$nclusterscpairgttbinary,
          nsubjects = q$nsubjectscpairgttbinary,
          ncontrols = q$ncontrolscpairgttbinary,
          ICC = q$ICCcpairgttbinary,
          p2 = q$p2cpairgttbinary,
          p1 = q$p1cpairgttbinary,
          decrease = q$decreasecpairgttbinary
        )
        return(val)
      }, seed = TRUE)
    }
    if (input$type == 'Individually-Randomized Group' &&
        input$dist == 'Binary' && input$meth == 'Simulation') {
      out1$power <- NULL
      out1$power$power <- "Calculating..."
      answer <<- future({
        val <- cps.irgtt.binary(
          nsim = q$nsimcpsirgttbinary,
          nsubjects = textToNum(c(
            q$nsubjectscpsirgttbinary,
            q$nsubjects2cpsirgttbinary
          )),
          nclusters = q$nclusterscpsirgttbinary,
          p1 = q$p1cpsirgttbinary,
          p2 = q$p2cpsirgttbinary,
          sigma_b_sq2 = q$sigma_b_sq2cpsirgttbinary,
          alpha = q$alpha,
          poorFitOverride = q$poorFitOverride,
          lowPowerOverride = q$lowPowerOverride,
          timelimitOverride = q$timelimitOverride,
          seed = q$seed
        )
        return(val)
      }, seed = TRUE)
    }
    if (input$type == 'Individually-Randomized Group' &&
        input$dist == 'Count' && input$meth == 'Analytic') {
      out1$power <- NULL
      answer <<- future({
        val <- cpa.irgtt.count()
        return(val)
      }, seed = TRUE)
    }
    if (input$type == 'Individually-Randomized Group' &&
        input$dist == 'Count' && input$meth == 'Simulation') {
      out1$power <- NULL
      out1$power$power <- "Calculating..."
      answer <<- future({
        val <- cps.irgtt.count(
          nsim = q$nsimcpsirgttcount,
          nsubjects = textToNum(c(q$nsubjectscpsirgttcount, q$nsubjects2cpsirgttcount)),
          nclusters = q$nclusterscpsirgttcount,
          c1 = q$c1cpsirgttcount,
          c2 = q$c2cpsirgttcount,
          sigma_b_sq2 = q$sigma_b_sq2cpsirgttcount,
          alpha = q$alpha,
          poorFitOverride = q$poorFitOverride,
          lowPowerOverride = q$lowPowerOverride,
          timelimitOverride = q$timelimitOverride,
          seed = q$seed
        )
        return(val)
      }, seed = TRUE)
    }
    
    then(
      answer,
      onFulfilled = function(value) {
        for (i in names(value)) {
          out1$power[[i]] <- value[[i]]
        }
        out1 <<- out1
      },
      onRejected = function(error) {
        out1$power <- paste0("ERROR: ", error$message)
      }
    )
    
    
    finally(answer, function() {
      prog$close()
      disable("cancel")
      enable("button")
    })
    
    return(NULL)
  }, ignoreInit = TRUE) # end call the clusterPower functions
  
  ##################
  
  # cancel button
  observeEvent(input$cancel, {
    async_pid <- answer$process$get_pid()
    print(paste("Killing PID:", async_pid))
    answer$process$kill()
    out1$power <- NULL
    disable("cancel")
    enable("run1")
  }, ignoreInit = TRUE)
  
  #########################################
  #### DEBUG ACCESS PANEL #####
  #########################################
  
  AllInputs <- reactive({
    x <- reactiveValuesToList(input)
    holder <- NULL
    if (sum(grepl("click", names(x))) == 1) {
      x$click <- NULL
    }
    holder <- data.frame(
      names = names(x),
      values = unlist(x, use.names = FALSE),
      mode = unlist(lapply(x, mode))
    )
  })
  
  output$show_inputs <- renderTable({
    AllInputs()
  })
  
  observeEvent(input$browser, {
    browser()
  })
  
  ######################################### END DEBUG table
  
  #embedded documentation
  
  helplink <- function(fxnName = input$fxnName) {
    help <- sprintf(
      "http://127.0.0.1:%d/library/clusterPower/html/%s",
      tools::startDynamicHelp(NA),
      paste0(input$fxnName, ".html")
    )
    return(help)
  }
  
  output$helpdetails <- renderUI({
    tags$iframe(src = helplink(),
                height = 600,
                width = 600)
  })
  
  #create the graphing table
  observeEvent(req(out1$power$power), {
    if (is.character(isolate(out1$power$power)) == FALSE) {
      x <- reactiveValuesToList(input)
      holder <- NULL
      if (sum(grepl("click", names(x))) == 1) {
        x$click <- NULL
      }
      holder <- data.frame(
        argument = names(isolate(x)),
        values = unlist(isolate(x), use.names = FALSE)
      )
      specialnames <-
        dplyr::filter(holder, grepl(gsub("\\.", "", input$fxnName), argument))
      specialnames$argument <-
        gsub(gsub("\\.", "", input$fxnName),
             "",
             specialnames$argument)
      specialnames <- dplyr::arrange(specialnames, argument)
      if (x$meth == "Analytic") {
        tab <-
          rbind(specialnames,
                c("alpha", isolate(input$alpha)),
                c("power", round(out1$power$power, 3)))
      }
      if (x$meth == "Simulation") {
        tab <-
          rbind(
            specialnames,
            c("alpha", isolate(input$alpha)),
            c("upper CI", round(out1$power$power$Upper.95.CI, 3)),
            c("power", round(out1$power$power$Power, 3)),
            c("lower CI", round(out1$power$power$Lower.95.CI, 3))
          )
      }
      if (is.null(logargs$tab)) {
        logargs$tab <- tab
      } else {
        tab <- dplyr::select(tab, values)
        tab <- cbind.data.frame(logargs$tab, tab)
        logargs$tab <- data.frame(tab, check.names = TRUE)
      }
    } else {
      # if logargs$tab is an error, ignore it
      if (!is.null(logargs$tab)) {
        tab <- cbind.data.frame(logargs$tab)
        logargs$tab <- data.frame(tab, check.names = TRUE)
      }
    }
  })   # END create the graphing table
  
  #clear the data log under certain circumstances
  observeEvent(input$cleargraph, {
    logargs$tab <- NULL
  })
  
  observeEvent(input$cleargraph2, {
    logargs$tab <- NULL
  })
  
  observeEvent(input$fxnName, {
    logargs$tab <- NULL
  })
  
  # END clear the data log under certain circumstances
  
  # make the graph
  #update the axis choices
  observeEvent(input$fxnName, {
    x <- reactiveValuesToList(input)
    holder <- names(x)
    specialnames <-
      grep(gsub("\\.", "", isolate(input$fxnName)), holder, value = TRUE)
    specialnames <-
      gsub(gsub("\\.", "", isolate(input$fxnName)), "", specialnames)
    specialnames <-
      specialnames[grepl("nsim", specialnames) == FALSE]
    args_ <- c("alpha", specialnames)
    updateSelectInput(session, "axisname",
                      choices = args_)
  })
  
  plot_this <-
    eventReactive(list(input$axisname, req(logargs$tab)), {
      data <- data.frame(isolate(logargs$tab), check.names = TRUE)
      data$values <- as.numeric(data$values)
      data$argument <- as.factor(data$argument)
      var <- isolate(input$axisname)
      data <-
        dplyr::filter(data, argument == !!var |
                        argument == "power")
      # first remember the names
      n <- data$argument
      # transpose all but the first column (name)
      data <- data.frame(t(data[, -1]))
      colnames(data) <- n
      data[, 1:2] <- apply(data[, 1:2], 2,
                           function(x)
                             as.numeric(as.character(x)))
      return(data)
    })
  
  dpfun <- function(x) {
    var <- isolate(input$axisname)
    data <- plot_this()
    fun <- function(x) {
      x <- enquo(x)
      sol <-
        ggplot(data, aes(x = !!x, y = power)) +
        geom_point(aes(colour = "fff"), size = 2.5) +
        theme_minimal() + theme(legend.position = "none")
      return(sol)
    }
    power_plot <- suppressWarnings(fun(get(var)))
    if (nrow(data) > 1) {
      power_plot <-
        power_plot + geom_line(aes(colour = "fff"), size = 1.25) + xlab(var)
    } else {
      power_plot <- power_plot + xlab(var)
    }
    return(power_plot)
  }
  
  output$graphic <- renderPlot({
    dpfun()
  }, res = 96)
  
  output$dp <- renderTable({
    q <- plot_this()
    nearPoints(q, isolate(input$click), xvar = isolate(input$axisname))
  })
  
  # create reactive input data table
  output$tracker <-
    renderTable(logargs$tab)
  
  # Downloadable csv of reactive input data table
  output$downloadData <- downloadHandler(
    filename = function() {
      paste("clusterPower_", isolate(input$fxnName), ".csv", sep = "")
    },
    content = function(file) {
      write.csv(logargs$tab, file, row.names = FALSE)
    }
  )
  
  resultdisplay <-  reactive({
    q <- reactiveValuesToList(out1)
    if (input$verbose == FALSE)
      return(q$power$power)
    else
      return(q$power)
  })
  
  # present the output verbose/not verbose
  output$CRTpower <- renderPrint(resultdisplay())
  
} #end of server fxn



# Run the application
shinyApp(ui = ui, server = server)
