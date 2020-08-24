#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#



ui <- fluidPage(
  shinyjs::useShinyjs(),
  HTML(
    "<h3>Estimate Power for a Randomized Controlled Trial with clusterPower</h3>
        <p>To use the calculator, select the trial type, outcome distribution, and calculation method.<p>
        <p>Then enter values for the quantities that appear below. When complete, select the ESTIMATE POWER button.</p>
        <p>You may specify more than one input quantity by separating numbers with commas.<p>"
  ),
  HTML(
    "This Beta has minimal documentation; please contact ken.kleinman@gmail.com with any feedback."
  ),
  column(12, bookmarkButton("Save App State")),
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
        clusterPower::argMatch("cpa.normal"),
        shinyjs::hidden(textInput("fxnName", "clusterPower function name", value = "cpa.normal"))
      ),
      conditionalPanel(
        "input.type == 'Parallel' & input.dist == 'Normal' & input.meth == 'Simulation'",
        clusterPower::argMatch("cps.normal"),
        shinyjs::hidden(textInput("fxnName", "clusterPower function name", value = "cps.normal"))
      ),
      conditionalPanel(
        "input.type == 'Parallel' & input.dist == 'Binary' & input.meth == 'Analytic'",
        clusterPower::argMatch("cpa.binary"),
        shinyjs::hidden(textInput("fxnName", "clusterPower function name", value = "cpa.binary"))
      ),
      conditionalPanel(
        "input.type == 'Parallel' & input.dist == 'Binary' & input.meth == 'Simulation'",
        clusterPower::argMatch("cps.binary"),
        shinyjs::hidden(textInput("fxnName", "clusterPower function name", value = "cps.binary"))
      ),
      conditionalPanel(
        "input.type == 'Parallel' & input.dist == 'Count' & input.meth == 'Analytic'",
        clusterPower::argMatch("cpa.count"),
        shinyjs::hidden(textInput("fxnName", "clusterPower function name", value = "cpa.count"))
      ),
      conditionalPanel(
        "input.type == 'Parallel' & input.dist == 'Count' & input.meth == 'Simulation'",
        clusterPower::argMatch("cps.count"),
        shinyjs::hidden(textInput("fxnName", "clusterPower function name", value = "cps.count"))
      ),
      conditionalPanel(
        "input.type == 'Multi-Arm' & input.dist == 'Normal' & input.meth == 'Analytic'",
        clusterPower::argMatch("cpa.ma.normal"),
        shinyjs::hidden(textInput("fxnName", "clusterPower function name", value = "cpa.ma.normal"))
      ),
      conditionalPanel(
        "input.type == 'Multi-Arm' & input.dist == 'Normal' & input.meth == 'Simulation'",
        clusterPower::argMatch("cps.ma.normal"),
        shinyjs::hidden(textInput("fxnName", "clusterPower function name", value = "cps.ma.normal"))
      ),
      conditionalPanel(
        "input.type == 'Multi-Arm' & input.dist == 'Binary' & input.meth == 'Analytic'",
        print("No method exists. Use the simulation option instead.")
      ),
      conditionalPanel(
        "input.type == 'Multi-Arm' & input.dist == 'Binary' & input.meth == 'Simulation'",
        clusterPower::argMatch("cps.ma.binary"),
        shinyjs::hidden(textInput("fxnName", "clusterPower function name", value = "cps.ma.binary"))
      ),
      conditionalPanel(
        "input.type == 'Multi-Arm' & input.dist == 'Count' & input.meth == 'Analytic'",
        print("No method exists. Use the simulation option instead.")
      ),
      conditionalPanel(
        "input.type == 'Multi-Arm' & input.dist == 'Count' & input.meth == 'Simulation'",
        clusterPower::argMatch("cps.ma.count"),
        shinyjs::hidden(textInput("fxnName", "clusterPower function name", value = "cps.ma.count"))
      ),
      conditionalPanel(
        "input.type == 'Difference-in-Difference' & input.dist == 'Normal' & input.meth == 'Analytic'",
        clusterPower::argMatch("cpa.did.normal"),
        shinyjs::hidden(textInput("fxnName", "clusterPower function name", value = "cpa.did.normal"))
      ),
      conditionalPanel(
        "input.type == 'Difference-in-Difference' & input.dist == 'Normal' & input.meth == 'Simulation'",
        clusterPower::argMatch("cps.did.normal"),
        shinyjs::hidden(textInput("fxnName", "clusterPower function name", value = "cps.did.normal"))
      ),
      conditionalPanel(
        "input.type == 'Difference-in-Difference' & input.dist == 'Binary' & input.meth == 'Analytic'",
        clusterPower::argMatch("cpa.did.binary"),
        shinyjs::hidden(textInput("fxnName", "clusterPower function name", value = "cpa.did.binary"))
      ),
      conditionalPanel(
        "input.type == 'Difference-in-Difference' & input.dist == 'Binary' & input.meth == 'Simulation'",
        clusterPower::argMatch("cps.did.binary"),
        shinyjs::hidden(textInput("fxnName", "clusterPower function name", value = "cps.did.binary"))
      ),
      conditionalPanel(
        "input.type == 'Difference-in-Difference' & input.dist == 'Count' & input.meth == 'Analytic'",
        print("No method exists. Use the simulation option instead.")
      ),
      conditionalPanel(
        "input.type == 'Difference-in-Difference' & input.dist == 'Count' & input.meth == 'Simulation'",
        clusterPower::argMatch("cps.did.count"),
        shinyjs::hidden(textInput("fxnName", "clusterPower function name", value = "cps.did.count"))
      ),
      conditionalPanel(
        "input.type == 'Stepped Wedge' & input.dist == 'Normal' & input.meth == 'Analytic'",
        clusterPower::argMatch("cpa.sw.normal"),
        shinyjs::hidden(textInput("fxnName", "clusterPower function name", value = "cpa.sw.normal"))
      ),
      conditionalPanel(
        "input.type == 'Stepped Wedge' & input.dist == 'Normal' & input.meth == 'Simulation'",
        clusterPower::argMatch("cps.sw.normal"),
        shinyjs::hidden(textInput("fxnName", "clusterPower function name", value = "cps.sw.normal"))
      ),
      conditionalPanel(
        "input.type == 'Stepped Wedge' & input.dist == 'Binary' & input.meth == 'Analytic'",
        clusterPower::argMatch("cpa.sw.binary"),
        shinyjs::hidden(textInput("fxnName", "clusterPower function name", value = "cpa.sw.binary"))
      ),
      conditionalPanel(
        "input.type == 'Stepped Wedge' & input.dist == 'Binary' & input.meth == 'Simulation'",
        clusterPower::argMatch("cps.sw.binary"),
        shinyjs::hidden(textInput("fxnName", "clusterPower function name", value = "cps.sw.binary"))
      ),
      conditionalPanel(
        "input.type == 'Stepped Wedge' & input.dist == 'Count' & input.meth == 'Analytic'",
        clusterPower::argMatch("cpa.sw.count"),
        shinyjs::hidden(textInput("fxnName", "clusterPower function name", value = "cpa.sw.count"))
      ),
      conditionalPanel(
        "input.type == 'Stepped Wedge' & input.dist == 'Count' & input.meth == 'Simulation'",
        clusterPower::argMatch("cps.sw.count"),
        shinyjs::hidden(textInput("fxnName", "clusterPower function name", value = "cps.sw.count"))
      ),
      conditionalPanel(
        "input.type == 'Individually-Randomized Group' & input.dist == 'Normal' & input.meth == 'Analytic'",
        clusterPower::argMatch("cpa.irgtt.normal"),
        shinyjs::hidden(textInput("fxnName", "clusterPower function name", value = "cpa.irgtt.normal"))
      ),
      conditionalPanel(
        "input.type == 'Individually-Randomized Group' & input.dist == 'Normal' & input.meth == 'Simulation'",
        clusterPower::argMatch("cps.irgtt.normal"),
        shinyjs::hidden(textInput("fxnName", "clusterPower function name", value = "cps.irgtt.normal"))
      ),
      conditionalPanel(
        "input.type == 'Individually-Randomized Group' & input.dist == 'Binary' & input.meth == 'Analytic'",
        clusterPower::argMatch("cpa.irgtt.binary"),
        shinyjs::hidden(textInput("fxnName", "clusterPower function name", value = "cpa.irgtt.binary"))
      ),
      conditionalPanel(
        "input.type == 'Individually-Randomized Group' & input.dist == 'Binary' & input.meth == 'Simulation'",
        clusterPower::argMatch("cps.irgtt.binary"),
        shinyjs::hidden(textInput("fxnName", "clusterPower function name", value = "cps.irgtt.binary"))
      ),
      conditionalPanel(
        "input.type == 'Individually-Randomized Group' & input.dist == 'Count' & input.meth == 'Analytic'",
        print("No method exists. Use the simulation option instead.")
      ),
      conditionalPanel(
        "input.type == 'Individually-Randomized Group' & input.dist == 'Count' & input.meth == 'Simulation'",
        clusterPower::argMatch("cps.irgtt.count"),
        shinyjs::hidden(textInput("fxnName", "clusterPower function name", value = "cps.irgtt.count"))
      )
      ), #end of values that can be reset with the restore defaults button
      
      actionButton(
        "button",
        "Estimate Power",
        icon = icon("arrow-circle-right"),
        width = '100%',
        style = "color: #fff; background-color: #337ab7; border-color: #2e6da4"
      ),
      actionButton('cancel', 'Cancel'),
      checkboxInput("more", "Show advanced options", value = FALSE), 
      conditionalPanel("input.more == true",
                       sliderInput(
                         "alpha",
                         "Significance level (alpha)",
                         value = 0.05,
                         min = 0.01,
                         max = 0.1,
                         step = 0.02
                       )
      ),
      conditionalPanel("input.more == true & input.meth == 'Simulation'",
        checkboxInput("timelimitOverride", "Allow unlimited calculation time", value = FALSE),
        checkboxInput("lowPowerOverride", "Allow completion when power is < 0.5", value = FALSE),
        checkboxInput("poorFitOverride", "Allow completion when model fit is poor", value = FALSE),
        textInput("optmethod", "Specify an optimization method", value = "NLOPT_LN_NELDERMEAD"),
        numericInput("seed",
                       "Set the seed (for repeatability)",
                       value = NA,
                       step = 1)
    ),
    conditionalPanel("input.more == true",
    actionButton("restoreDefault", "Restore default parameters"),
    actionButton("reload", "Reset all", icon = icon("trash-alt"))
    )),
    mainPanel(dataTableOutput("tbl"),
              shinycssloaders::withSpinner(verbatimTextOutput("CRTpower", placeholder = TRUE))
              )
  )
)


######################################

    #       SERVER

######################################
server <- function(input, output, session) {

  # Reload the app
  observeEvent(input$reload,{
    session$reload()
  })
  
  # Restore default values
  observeEvent(input$restoreDefault,{
    shinyjs::reset("allValues")
  })
  
  #make some helpful fxns to extract arg names
  updateArgs <- function(fxnName) {
    argMatchResult <- c(clusterPower::argMatch(fxnName, justNames = TRUE), 
                        "lowPowerOverride", "poorFitOverride", "timelimitOverride", 
                        "power", "seed", "optmethod")
    argNames <-
      c(
        "nsubjects",
        "nclusters",
        "alpha",
        dplyr::intersect(
          argMatchResult,
          names(formals(fxnName))
        )
      )
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
  
  answer <- eventReactive(input$button, {
    if (input$type == 'Parallel' &&
        input$dist == 'Normal' && input$meth == 'Analytic') {
      print(printresult("cpa.normal"))
      fxnName <- "cpa.normal"
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
      print("No method exists. Use the simulation option instead.")
    }
    if (input$type == 'Multi-Arm' &
        input$dist == 'Binary' & input$meth == 'Simulation') {
      print(printresult("cps.ma.binary"))
    }
    if (input$type == 'Multi-Arm' &
        input$dist == 'Count' & input$meth == 'Analytic') {
      print("No method exists. Use the simulation option instead.")
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
      print("No method exists. Use the simulation option instead.")
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
      print("No method exists. Use the simulation option instead.")
    }
    if (input$type == 'Individually-Randomized Group' &
        input$dist == 'Count' & input$meth == 'Simulation') {
      print(summary(printresult("cps.irgtt.count")))
    }
  })
  

  # create input data table
  args <- reactive({t(data.frame(unlist(updateArgs(input$fxnName))))
  })
  output$tbl <- shiny::renderDataTable(
    args(), options = list(lengthChange = FALSE, searching = FALSE, paging = FALSE)
  )
  output$CRTpower <- renderPrint({
    answer()
  })
} #end of server fxn



# Run the application
shinyApp(ui = ui, server = server)
