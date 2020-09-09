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
  h1(id="big-heading", "Power Estimation for Randomized Controlled Trials: clusterPower"),
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
        numericInput("nclusters", "Number of Clusters", value = 10),
        numericInput("nsubjects", "Number of Observations (per cluster)", value = 20),
        shinyjs::hidden(numericInput("power", "power", value = NA)),
        conditionalPanel(
          "input.type == 'Parallel' & input.dist == 'Normal' & input.meth == 'Analytic'",
          clusterPower::argMatch("cpa.normal")
        ),
        conditionalPanel(
          "input.type == 'Parallel' & input.dist == 'Normal' & input.meth == 'Simulation'",
          clusterPower::argMatch("cps.normal")
        ),
        conditionalPanel(
          "input.type == 'Parallel' & input.dist == 'Binary' & input.meth == 'Analytic'",
          clusterPower::argMatch("cpa.binary")
        ),
        conditionalPanel(
          "input.type == 'Parallel' & input.dist == 'Binary' & input.meth == 'Simulation'",
          clusterPower::argMatch("cps.binary")
        ),
        conditionalPanel(
          "input.type == 'Parallel' & input.dist == 'Count' & input.meth == 'Analytic'",
          clusterPower::argMatch("cpa.count")
        ),
        conditionalPanel(
          "input.type == 'Parallel' & input.dist == 'Count' & input.meth == 'Simulation'",
          clusterPower::argMatch("cps.count")
        ),
        conditionalPanel(
          "input.type == 'Multi-Arm' & input.dist == 'Normal' & input.meth == 'Analytic'",
          clusterPower::argMatch("cpa.ma.normal")
        ),
        conditionalPanel(
          "input.type == 'Multi-Arm' & input.dist == 'Normal' & input.meth == 'Simulation'",
          clusterPower::argMatch("cps.ma.normal")
        ),
        conditionalPanel(
          "input.type == 'Multi-Arm' & input.dist == 'Binary' & input.meth == 'Analytic'",
          HTML("No method exists. Use the simulation option instead.")
        ),
        conditionalPanel(
          "input.type == 'Multi-Arm' & input.dist == 'Binary' & input.meth == 'Simulation'",
          clusterPower::argMatch("cps.ma.binary")
        ),
        conditionalPanel(
          "input.type == 'Multi-Arm' & input.dist == 'Count' & input.meth == 'Analytic'",
          HTML("No method exists. Use the simulation option instead.")
        ),
        conditionalPanel(
          "input.type == 'Multi-Arm' & input.dist == 'Count' & input.meth == 'Simulation'",
          clusterPower::argMatch("cps.ma.count")
        ),
        conditionalPanel(
          "input.type == 'Difference-in-Difference' & input.dist == 'Normal' & input.meth == 'Analytic'",
          clusterPower::argMatch("cpa.did.normal")
        ),
        conditionalPanel(
          "input.type == 'Difference-in-Difference' & input.dist == 'Normal' & input.meth == 'Simulation'",
          clusterPower::argMatch("cps.did.normal")
        ),
        conditionalPanel(
          "input.type == 'Difference-in-Difference' & input.dist == 'Binary' & input.meth == 'Analytic'",
          clusterPower::argMatch("cpa.did.binary")
        ),
        conditionalPanel(
          "input.type == 'Difference-in-Difference' & input.dist == 'Binary' & input.meth == 'Simulation'",
          clusterPower::argMatch("cps.did.binary")
        ),
        conditionalPanel(
          "input.type == 'Difference-in-Difference' & input.dist == 'Count' & input.meth == 'Analytic'",
          HTML("No method exists. Use the simulation option instead.")
        ),
        conditionalPanel(
          "input.type == 'Difference-in-Difference' & input.dist == 'Count' & input.meth == 'Simulation'",
          clusterPower::argMatch("cps.did.count")
        ),
        conditionalPanel(
          "input.type == 'Stepped Wedge' & input.dist == 'Normal' & input.meth == 'Analytic'",
          clusterPower::argMatch("cpa.sw.normal")
        ),
        conditionalPanel(
          "input.type == 'Stepped Wedge' & input.dist == 'Normal' & input.meth == 'Simulation'",
          clusterPower::argMatch("cps.sw.normal")
        ),
        conditionalPanel(
          "input.type == 'Stepped Wedge' & input.dist == 'Binary' & input.meth == 'Analytic'",
          clusterPower::argMatch("cpa.sw.binary")
        ),
        conditionalPanel(
          "input.type == 'Stepped Wedge' & input.dist == 'Binary' & input.meth == 'Simulation'",
          clusterPower::argMatch("cps.sw.binary")
        ),
        conditionalPanel(
          "input.type == 'Stepped Wedge' & input.dist == 'Count' & input.meth == 'Analytic'",
          clusterPower::argMatch("cpa.sw.count")
        ),
        conditionalPanel(
          "input.type == 'Stepped Wedge' & input.dist == 'Count' & input.meth == 'Simulation'",
          clusterPower::argMatch("cps.sw.count")
        ),
        conditionalPanel(
          "input.type == 'Individually-Randomized Group' & input.dist == 'Normal' & input.meth == 'Analytic'",
          clusterPower::argMatch("cpa.irgtt.normal")
        ),
        conditionalPanel(
          "input.type == 'Individually-Randomized Group' & input.dist == 'Normal' & input.meth == 'Simulation'",
          clusterPower::argMatch("cps.irgtt.normal")
        ),
        conditionalPanel(
          "input.type == 'Individually-Randomized Group' & input.dist == 'Binary' & input.meth == 'Analytic'",
          clusterPower::argMatch("cpa.irgtt.binary")
        ),
        conditionalPanel(
          "input.type == 'Individually-Randomized Group' & input.dist == 'Binary' & input.meth == 'Simulation'",
          clusterPower::argMatch("cps.irgtt.binary")
        ),
        conditionalPanel(
          "input.type == 'Individually-Randomized Group' & input.dist == 'Count' & input.meth == 'Analytic'",
          HTML("No method exists. Use the simulation option instead.")
        ),
        conditionalPanel(
          "input.type == 'Individually-Randomized Group' & input.dist == 'Count' & input.meth == 'Simulation'",
          clusterPower::argMatch("cps.irgtt.count")
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
    
    # Tabs start
    
    mainPanel(tabsetPanel(
      tabPanel("DEBUG", textOutput("powe")),
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
                   to the R functions based on user input. It is updated 
          when clusterPower returns a result.</p>"
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
  observeEvent(watchfor(), {
    if (input$type == 'Parallel' &&
        input$dist == 'Normal' && input$meth == 'Analytic') {
      updateTextInput(session, "fxnName", value = "cpa.normal")
    }
    if (input$type == 'Parallel' &&
        input$dist == 'Normal' && input$meth == 'Simulation') {
      updateTextInput(session, "fxnName", value = "cps.normal")
    }
    if (input$type == 'Parallel' &
        input$dist == 'Binary' & input$meth == 'Analytic') {
      updateTextInput(session, "fxnName", value = "cpa.binary")
    }
    if (input$type == 'Parallel' &
        input$dist == 'Binary' & input$meth == 'Simulation') {
      updateTextInput(session, "fxnName", value = "cps.binary")
    }
    if (input$type == 'Parallel' &
        input$dist == 'Count' & input$meth == 'Analytic') {
      updateTextInput(session, "fxnName", value = "cpa.count")
    }
    if (input$type == 'Parallel' &
        input$dist == 'Count' & input$meth == 'Simulation') {
      updateTextInput(session, "fxnName", value = "cps.count")
    }
    if (input$type == 'Multi-Arm' &
        input$dist == 'Normal' & input$meth == 'Analytic') {
      updateTextInput(session, "fxnName", value = "cpa.ma.normal")
    }
    if (input$type == 'Multi-Arm' &
        input$dist == 'Normal' & input$meth == 'Simulation') {
      updateTextInput(session, "fxnName", value = "cps.ma.normal")
    }
    if (input$type == 'Multi-Arm' &
        input$dist == 'Binary' & input$meth == 'Analytic') {
      updateTextInput(session, "fxnName", value = "cpa.ma.binary")
    }
    if (input$type == 'Multi-Arm' &
        input$dist == 'Binary' & input$meth == 'Simulation') {
      updateTextInput(session, "fxnName", value = "cps.ma.binary")
    }
    if (input$type == 'Multi-Arm' &
        input$dist == 'Count' & input$meth == 'Analytic') {
      updateTextInput(session, "fxnName", value = "cpa.ma.count")
    }
    if (input$type == 'Multi-Arm' &
        input$dist == 'Count' & input$meth == 'Simulation') {
      updateTextInput(session, "fxnName", value = "cps.ma.count")
    }
    if (input$type == 'Difference-in-Difference' &
        input$dist == 'Normal' & input$meth == 'Analytic') {
      updateTextInput(session, "fxnName", value = "cpa.did.normal")
    }
    if (input$type == 'Difference-in-Difference' &
        input$dist == 'Normal' & input$meth == 'Simulation') {
      updateTextInput(session, "fxnName", value = "cps.did.normal")
    }
    if (input$type == 'Difference-in-Difference' &
        input$dist == 'Binary' & input$meth == 'Analytic') {
      updateTextInput(session, "fxnName", value = "cpa.did.binary")
    }
    if (input$type == 'Difference-in-Difference' &
        input$dist == 'Binary' & input$meth == 'Simulation') {
      updateTextInput(session, "fxnName", value = "cps.did.binary")
    }
    if (input$type == 'Difference-in-Difference' &
        input$dist == 'Count' & input$meth == 'Analytic') {
      updateTextInput(session, "fxnName", value = "cpa.did.count")
    }
    if (input$type == 'Difference-in-Difference' &
        input$dist == 'Count' & input$meth == 'Simulation') {
      updateTextInput(session, "fxnName", value = "cps.did.count")
    }
    if (input$type == 'Stepped Wedge' &
        input$dist == 'Normal' & input$meth == 'Analytic') {
      updateTextInput(session, "fxnName", value = "cpa.sw.normal")
    }
    if (input$type == 'Stepped Wedge' &
        input$dist == 'Normal' & input$meth == 'Simulation') {
      updateTextInput(session, "fxnName", value = "cps.sw.normal")
    }
    if (input$type == 'Stepped Wedge' &
        input$dist == 'Binary' & input$meth == 'Analytic') {
      updateTextInput(session, "fxnName", value = "cpa.sw.binary")
    }
    if (input$type == 'Stepped Wedge' &
        input$dist == 'Binary' & input$meth == 'Simulation') {
      updateTextInput(session, "fxnName", value = "cps.sw.binary")
    }
    if (input$type == 'Stepped Wedge' &
        input$dist == 'Count' & input$meth == 'Analytic') {
      updateTextInput(session, "fxnName", value = "cpa.sw.count")
    }
    if (input$type == 'Stepped Wedge' &
        input$dist == 'Count' & input$meth == 'Simulation') {
      updateTextInput(session, "fxnName", value = "cps.sw.count")
    }
    if (input$type == 'Individually-Randomized Group' &
        input$dist == 'Normal' & input$meth == 'Analytic') {
      updateTextInput(session, "fxnName", value = "cpa.irgtt.normal")
    }
    if (input$type == 'Individually-Randomized Group' &
        input$dist == 'Normal' & input$meth == 'Simulation') {
      updateTextInput(session, "fxnName", value = "cps.irgtt.normal")
    }
    if (input$type == 'Individually-Randomized Group' &
        input$dist == 'Binary' & input$meth == 'Analytic') {
      updateTextInput(session, "fxnName", value = "cpa.irgtt.binary")
    }
    if (input$type == 'Individually-Randomized Group' &
        input$dist == 'Binary' & input$meth == 'Simulation') {
      updateTextInput(session, "fxnName", value = "cps.irgtt.binary")
    }
    if (input$type == 'Individually-Randomized Group' &
        input$dist == 'Count' & input$meth == 'Analytic') {
      updateTextInput(session, "fxnName", value = "cpa.irgtt.count")
    }
    if (input$type == 'Individually-Randomized Group' &
        input$dist == 'Count' & input$meth == 'Simulation') {
      updateTextInput(session, "fxnName", value = "cps.irgtt.count")
    }
  }) # end update help documentation and params table when function is selected
  
  # call the clusterPower functions
  answer <- eventReactive(input$button, {
    if (input$type == 'Parallel' &&
        input$dist == 'Normal' && input$meth == 'Analytic') {
      print(printresult("cpa.normal"))
    }
    if (input$type == 'Parallel' &&
        input$dist == 'Normal' && input$meth == 'Simulation') {
      print(summary(printresult("cps.normal")))
    }
    if (input$type == 'Parallel' &
        input$dist == 'Binary' & input$meth == 'Analytic') {
      print(printresult("cpa.binary"))
    }
    if (input$type == 'Parallel' &
        input$dist == 'Binary' & input$meth == 'Simulation') {
      print(summary(printresult("cps.binary")))
    }
    if (input$type == 'Parallel' &
        input$dist == 'Count' & input$meth == 'Analytic') {
      print(printresult("cpa.count"))
    }
    if (input$type == 'Parallel' &
        input$dist == 'Count' & input$meth == 'Simulation') {
      print(summary(printresult("cps.count")))
    }
    if (input$type == 'Multi-Arm' &
        input$dist == 'Normal' & input$meth == 'Analytic') {
      print(printresult("cpa.ma.normal"))
    }
    if (input$type == 'Multi-Arm' &
        input$dist == 'Normal' & input$meth == 'Simulation') {
      print(printresult("cps.ma.normal"))
    }
    if (input$type == 'Multi-Arm' &
        input$dist == 'Binary' & input$meth == 'Analytic') {
      print(printresult("cpa.ma.binary"))
    }
    if (input$type == 'Multi-Arm' &
        input$dist == 'Binary' & input$meth == 'Simulation') {
      print(printresult("cps.ma.binary"))
    }
    if (input$type == 'Multi-Arm' &
        input$dist == 'Count' & input$meth == 'Analytic') {
      print(printresult("cpa.ma.count"))
    }
    if (input$type == 'Multi-Arm' &
        input$dist == 'Count' & input$meth == 'Simulation') {
      print(printresult("cps.ma.count"))
    }
    if (input$type == 'Difference-in-Difference' &
        input$dist == 'Normal' & input$meth == 'Analytic') {
      print(printresult("cpa.did.normal"))
    }
    if (input$type == 'Difference-in-Difference' &
        input$dist == 'Normal' & input$meth == 'Simulation') {
      print(summary(printresult("cps.did.normal")))
    }
    if (input$type == 'Difference-in-Difference' &
        input$dist == 'Binary' & input$meth == 'Analytic') {
      print(printresult("cpa.did.binary"))
    }
    if (input$type == 'Difference-in-Difference' &
        input$dist == 'Binary' & input$meth == 'Simulation') {
      print(summary(printresult("cps.did.binary")))
    }
    if (input$type == 'Difference-in-Difference' &
        input$dist == 'Count' & input$meth == 'Analytic') {
      print(printresult("cpa.did.count"))
    }
    if (input$type == 'Difference-in-Difference' &
        input$dist == 'Count' & input$meth == 'Simulation') {
      print(summary(printresult("cps.did.count")))
    }
    if (input$type == 'Stepped Wedge' &
        input$dist == 'Normal' & input$meth == 'Analytic') {
      print(printresult("cpa.sw.normal"))
    }
    if (input$type == 'Stepped Wedge' &
        input$dist == 'Normal' & input$meth == 'Simulation') {
      print(summary(printresult("cps.sw.normal")))
    }
    if (input$type == 'Stepped Wedge' &
        input$dist == 'Binary' & input$meth == 'Analytic') {
      print(printresult("cpa.sw.binary"))
    }
    if (input$type == 'Stepped Wedge' &
        input$dist == 'Binary' & input$meth == 'Simulation') {
      print(summary(printresult("cps.sw.binary")))
    }
    if (input$type == 'Stepped Wedge' &
        input$dist == 'Count' & input$meth == 'Analytic') {
      print(printresult("cpa.sw.count"))
    }
    if (input$type == 'Stepped Wedge' &
        input$dist == 'Count' & input$meth == 'Simulation') {
      print(summary(printresult("cps.sw.count")))
    }
    if (input$type == 'Individually-Randomized Group' &
        input$dist == 'Normal' & input$meth == 'Analytic') {
      print(printresult("cpa.irgtt.normal"))
    }
    if (input$type == 'Individually-Randomized Group' &
        input$dist == 'Normal' & input$meth == 'Simulation') {
      print(summary(printresult("cps.irgtt.normal")))
    }
    if (input$type == 'Individually-Randomized Group' &
        input$dist == 'Binary' & input$meth == 'Analytic') {
      print(printresult("cpa.irgtt.binary"))
    }
    if (input$type == 'Individually-Randomized Group' &
        input$dist == 'Binary' & input$meth == 'Simulation') {
      print(summary(printresult("cps.irgtt.binary")))
    }
    if (input$type == 'Individually-Randomized Group' &
        input$dist == 'Count' & input$meth == 'Analytic') {
      print(printresult("cpa.irgtt.count"))
    }
    if (input$type == 'Individually-Randomized Group' &
        input$dist == 'Count' & input$meth == 'Simulation') {
      print(summary(printresult("cps.irgtt.count")))
    }
  }) # end call the clusterPower functions
  ##################
  
  output$powe <- renderText(input$ICC)
  
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
