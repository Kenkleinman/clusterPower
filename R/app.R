library(shiny)
library(DT)
library(tidyverse)
library(stringr)

source("crtpwr.2mean.R")
source("crtpwr.2prop.R")
source("crtpwr.2rate.R")

names2mean <- names(as.list(args(crtpwr.2mean)))[1:8]
names2prop <- names(as.list(args(crtpwr.2prop)))[1:8]
names2rate <- names(as.list(args(crtpwr.2rate)))[1:7]

umass <- "font-family: 'Open Sans', Helvetica, Arial, sans-serif; font-weight: bold; color: #ffffff; background-color: #881c1c; border: 3px solid #000000;"

make_sequence <- function(x){
  # trim any spaces before and after the supplied string
  x <- str_trim(x)
  if(x == ""){
    # if x is an empty string, just return NA
    as.numeric(x)
  } else if(str_detect(x, "^[0-9.]")){
    # if the initial character of x is a digit or period, use the supplied
    # numbers to create a vector of possibilities
    temp <- as.numeric(str_split(x,"[^0-9.]",simplify=TRUE))
    temp[!is.na(temp)]
  } else {
    # if the initial character of x is not a digit or period, use the three
    # supplied numbers to create a sequence 'from x to y by z'
    temp <- as.numeric(str_split(x,"[^0-9.]",simplify=TRUE))
    temp <- temp[!is.na(temp)]
    seq(temp[1], temp[2], by = temp[3])
  }
}

get_graph_width <- function(x){
  if(is.na(x)){
    "auto"
  } else {
    x
  }
}

ui <- fluidPage(
  HTML("<h3>Simple Two-Arm Designs</h3>
        <p>To use the two-arm calculator, leave the desired quantity blank and enter values for the other quantities.</p>
        <p>You may specify more than one input quantity by separating numbers with spaces or commas.<p>
        <p>You may specify a sequence of values by typing 'from X to Y by Z', where 'X' is the starting value, 'Y' is the ending value, and 'Z' is the increment.</p>"),
  tabsetPanel(
    tabPanel("Continuous",
             column(2,
                    fluidRow(textInput("alpha2mean", HTML("Type I rate, &alpha; (alpha)"),
                                       value = "0.05", width = "100%")),
                    fluidRow(textInput("power2mean", "Power (power)",
                                       value = "", width = "100%")),
                    fluidRow(textInput("cv2mean", "Cluster size CV (cv)",
                                       value = "0", width = "100%")),
                    fluidRow(textInput("d2mean", "Difference (d)",
                                       value = ".2 .5", width = "100%")),
                    fluidRow(textInput("m2mean", "Clusters per arm (m)",
                                       value = "3 6 9 12", width = "100%")),
                    fluidRow(textInput("n2mean", "Cluster size (n)",
                                       value = "100 300 200", width = "100%")),
                    fluidRow(textInput("icc2mean", "ICC (icc)",
                                       value = ".001 .005 .01", width = "100%")),
                    fluidRow(textInput("varw2mean", "Within variance (varw)",
                                       value = "1", width = "100%")),
                    fluidRow(
                      column(6, style='padding:0px;', actionButton("default2mean", "Defaults", width = "100%")),
                      column(6, style='padding:0px;', actionButton("clear2mean", "Clear All", width = "100%"))
                    ),
                    fluidRow(
                      column(12, style='padding: 0px;', actionButton("calc2mean", "Calculate", width = "100%",
                                                                     style = umass))
                    ),
                    conditionalPanel(condition = "output.table2mean != null",
                                     fluidRow(downloadButton("dl2mean", "Download")))
             ),
             column(10,
                    tabsetPanel(
                      tabPanel("Data",
                               DT::dataTableOutput("table2mean")
                      ),
                      tabPanel("Graphs",
                               column(2,
                                      fluidRow(selectInput("y2mean", "Y",
                                                           choices = names2mean, selected = "power")),
                                      fluidRow(selectInput("x2mean", "X",
                                                           choices = names2mean, selected = "m")),
                                      fluidRow(selectInput("group2mean", "Group",
                                                           choices = c(None = ".", names2mean), selected = "n")),
                                      fluidRow(checkboxInput("color2mean", "Color by Group",value = TRUE)),
                                      fluidRow(selectInput("row2mean", "Facet Row",
                                                           choices = c(None = ".", names2mean))),
                                      fluidRow(selectInput("col2mean", "Facet Column",
                                                           choices = c(None = ".", names2mean))),
                                      fluidRow(numericInput("height2mean", "Plot Height", value = 400,
                                                            min = 100, max = 2000, step = 10)),
                                      fluidRow(numericInput("psize2mean", "Point Size", value = 3,
                                                            min = 0.5, max = 4, step = 0.25)),
                                      fluidRow(numericInput("lsize2mean", "Line Width", value = 1,
                                                            min = 0.5, max = 2, step = 0.25))
                               ),
                               column(10,
                                      plotOutput("graph2mean", height = "auto")
                               )
                      ) # end tabPanel("Graphs"...
                    ) # end tabsetPanel(...
             ) # end column(10,...
    ), # end tabPanel("Continuous ...
    tabPanel("Binary",
             column(2,
                    fluidRow(textInput("alpha2prop", HTML("Type I rate, &alpha; (alpha)"),
                                       value = "0.05", width = "100%")),
                    fluidRow(textInput("power2prop", "Power (power)",
                                       value = "0.80", width = "100%")),
                    fluidRow(textInput("cv2prop", "Cluster size CV (cv)",
                                       value = "0", width = "100%")),
                    fluidRow(textInput("icc2prop", "ICC (icc)",
                                       value = "", width = "100%")),
                    fluidRow(textInput("m2prop", "Clusters per arm (m)",
                                       value = "", width = "100%")),
                    fluidRow(textInput("n2prop", "Cluster size (n)",
                                       value = "", width = "100%")),
                    fluidRow(textInput("p12prop", "Trt Proportion (p1)",
                                       value = "", width = "100%")),
                    fluidRow(textInput("p22prop", "Ctrl Proportion (p2)",
                                       value = "", width = "100%")),
                    fluidRow(
                      column(6,checkboxInput("p1inc2prop", "p1 > p2")),
                      column(6,checkboxInput("pooled2prop", "Pooled"))
                    ),
                    fluidRow(
                      column(6, style='padding:0px;', actionButton("default2prop", "Defaults", width = "100%")),
                      column(6, style='padding:0px;', actionButton("clear2prop", "Clear All", width = "100%"))
                    ),
                    fluidRow(
                      column(12, style='padding:0px;', actionButton("calc2prop", "Calculate", width = "100%",
                                                                    style = umass))
                    ),
                    conditionalPanel(condition = "output.table2prop != null",
                                     fluidRow(downloadButton("dl2prop", "Download")))
             ), # end column(2, ..
             column(10,
                    tabsetPanel(
                      tabPanel("Data",
                               DT::dataTableOutput("table2prop")
                      ),
                      tabPanel("Graphs",
                               column(2,
                                      fluidRow(selectInput("y2prop", "Y",
                                                           choices = names2prop, selected = "power")),
                                      fluidRow(selectInput("x2prop", "X",
                                                           choices = names2prop, selected = "m")),
                                      fluidRow(selectInput("group2prop", "Group",
                                                           choices = c(None = ".", names2prop), selected = "n")),
                                      fluidRow(checkboxInput("color2prop", "Color by Group",value = TRUE)),
                                      fluidRow(selectInput("row2prop", "Facet Row",
                                                           choices = c(None = ".", names2prop))),
                                      fluidRow(selectInput("col2prop", "Facet Column",
                                                           choices = c(None = ".", names2prop))),
                                      fluidRow(numericInput("height2prop", "Plot Height", value = 400,
                                                            min = 100, max = 2000, step = 10)),
                                      fluidRow(numericInput("psize2mean", "Point Size", value = 3,
                                                            min = 0.5, max = 4, step = 0.25)),
                                      fluidRow(numericInput("lsize2mean", "Line Width", value = 1,
                                                            min = 0.5, max = 2, step = 0.25))
                               ),
                               column(10,
                                      plotOutput("graph2prop", height = "auto")
                               )
                      ) # end tabPanel("Graphs"...
                    ) # end tabsetPanel(...
             ) # end column(10,...
    ), # end tabPanel("Binary ...
    tabPanel("Count",
             column(2,
                    fluidRow(textInput("alpha2rate", HTML("Type I rate, &alpha; (alpha)"),
                                       value = "0.05", width = "100%")),
                    fluidRow(textInput("power2rate", "Power (power)",
                                       value = "0.80", width = "100%")),
                    fluidRow(textInput("m2rate", "Clusters per arm (m)",
                                       value = "", width = "100%")),
                    fluidRow(textInput("py2rate", "PY per cluster (py)",
                                       value = "", width = "100%")),
                    fluidRow(textInput("r12rate", "Trt Rate (r1)",
                                       value = "", width = "100%")),
                    fluidRow(textInput("r22rate", "Ctrl Rate (r2)",
                                       value = "", width = "100%")),
                    fluidRow(textInput("cvb2rate", "Btwn-cluster CV (cvb)",
                                       value = "", width = "100%")),
                    fluidRow(
                      column(6, style='padding:0px;', actionButton("default2rate", "Defaults", width = "100%")),
                      column(6, style='padding:0px;', actionButton("clear2rate", "Clear All", width = "100%"))),
                    fluidRow(
                      column(12, style='padding:0px;', actionButton("calc2rate", "Calculate", width = "100%",
                                                                    style = umass))
                    ),
                    conditionalPanel(condition = "output.table2rate != null",
                                     fluidRow(downloadButton("dl2rate", "Download")))
             ),
             column(10,
                    tabsetPanel(
                      tabPanel("Data",
                               DT::dataTableOutput("table2rate")
                      ),
                      tabPanel("Graphs",
                               column(2,
                                      fluidRow(selectInput("y2rate", "Y",
                                                           choices = names2rate, selected = "power")),
                                      fluidRow(selectInput("x2rate", "X",
                                                           choices = names2rate, selected = "m")),
                                      fluidRow(selectInput("group2rate", "Group",
                                                           choices = names2rate, selected = "py")),
                                      fluidRow(checkboxInput("color2rate", "Color by Group",value = TRUE)),
                                      fluidRow(selectInput("row2rate", "Facet Row",
                                                           choices = c(None = ".", names2rate))),
                                      fluidRow(selectInput("col2rate", "Facet Column",
                                                           choices = c(None = ".", names2rate))),
                                      fluidRow(numericInput("height2rate", "Plot Height", value = 400,
                                                            min = 100, max = 2000, step = 10)),
                                      fluidRow(numericInput("psize2rate", "Point Size", value = 3,
                                                            min = 0.5, max = 4, step = 0.25)),
                                      fluidRow(numericInput("lsize2rate", "Line Width", value = 1,
                                                            min = 0.5, max = 2, step = 0.25))
                               ),
                               column(10,
                                      plotOutput("graph2rate", height = "auto")
                               )
                      ) # end tabPanel("Graphs"...
                    ) # end tabsetPanel(...
             ) # end column(10,...
    ) # end tabPanel("Count ...
  )
) # end fluidPage

server <- function(input, output, session){
  
  #----------------------------------------------------------------------------
  # Two means
  #----------------------------------------------------------------------------
  
  # reset 2mean inputs to default values
  observeEvent(
    input$default2mean,
    {
      updateNumericInput(session, inputId = "alpha2mean", value = "0.05")
      updateNumericInput(session, inputId = "power2mean", value = "0.80")
      updateNumericInput(session, inputId = "cv2mean", value = "0")
      updateNumericInput(session, inputId = "m2mean", value = "")
      updateNumericInput(session, inputId = "n2mean", value = "")
      updateNumericInput(session, inputId = "d2mean", value = "")
      updateNumericInput(session, inputId = "icc2mean", value = "")
      updateNumericInput(session, inputId = "varw2mean", value = "")
    }
  ) # end observeEvent(input$default2mean ...
  
  
  # clear 2mean inputs 
  observeEvent(
    input$clear2mean,
    {
      updateTextInput(session, inputId = "alpha2mean", value = "")
      updateTextInput(session, inputId = "power2mean", value = "")
      updateTextInput(session, inputId = "cv2mean", value = "")
      updateTextInput(session, inputId = "m2mean", value = "")
      updateTextInput(session, inputId = "n2mean", value = "")
      updateTextInput(session, inputId = "d2mean", value = "")
      updateTextInput(session, inputId = "icc2mean", value = "")
      updateTextInput(session, inputId = "varw2mean", value = "")
    }
  ) # end observeEvent(input$clear2mean ...
  
  # create 2mean data
  res2mean <- eventReactive(
    input$calc2mean,
    {
      alpha <- make_sequence(isolate(input$alpha2mean))
      power <- make_sequence(isolate(input$power2mean))
      m <- make_sequence(isolate(input$m2mean))
      n <- make_sequence(isolate(input$n2mean))
      cv <- make_sequence(isolate(input$cv2mean))
      d <- make_sequence(isolate(input$d2mean))
      icc <- make_sequence(isolate(input$icc2mean))
      varw <- make_sequence(isolate(input$varw2mean))
      
      tab <- expand.grid(alpha,
                         power,
                         m,
                         n,
                         cv,
                         d,
                         icc,
                         varw,
                         stringsAsFactors = FALSE)
      
      needind <- which(is.na(tab[1,]))
      names(tab) <- names2mean
      target <- names2mean[needind]
      
      tab[[target]] <- signif(pmap_dbl(tab, crtpwr.2mean),4)
      mutate_if(tab, names(tab) != target, factor)
    })
  
  # create 2mean output table
  output$table2mean = DT::renderDataTable(
    res2mean(),
    filter = 'top',
    options = list(
      columnDefs = list(list(className = 'dt-center', targets = '_all')),
      pageLength = 10,
      autoWidth = TRUE
    )
    
  )
  
  # set up 2mean download handler
  output$dl2mean <- downloadHandler(
    filename = function() {
      paste('data-2mean-', Sys.time(), '.csv', sep='')
    },
    content = function(file) {
      write.csv(res2mean(), file, row.names = FALSE)
    }
  )
  
  # create 2mean graph
  output$graph2mean <- renderPlot({
    if(input$color2mean){
      p <- ggplot(res2mean(),
                  aes_string(x = input$x2mean,
                             y = input$y2mean,
                             group = input$group2mean,
                             color = input$group2mean))
    } else {
      p <- ggplot(res2mean(),
                  aes_string(x = input$x2mean,
                             y = input$y2mean,
                             group = input$group2mean))
    }
    p <- p + geom_line(size = input$lsize2mean) +
      geom_point(size = input$psize2mean) +
      theme_grey(base_size = 18) 
    
    facets <- paste(input$row2mean,"~",input$col2mean)
    if(facets != ". ~ .") p <- p + facet_grid(facets)
    
    p
  },
  height = reactive({input$height2mean}))
  
  
  #----------------------------------------------------------------------------
  # Two proportions
  #----------------------------------------------------------------------------
  
  # reset 2prop inputs to default values
  observeEvent(
    input$default2prop,
    {
      updateTextInput(session, inputId = "alpha2prop", value = "0.05")
      updateTextInput(session, inputId = "power2prop", value = "0.80")
      updateTextInput(session, inputId = "cv2prop", value = "0")
      updateTextInput(session, inputId = "m2prop", value = "")
      updateTextInput(session, inputId = "n2prop", value = "")
      updateTextInput(session, inputId = "icc2prop", value = "")
      updateTextInput(session, inputId = "p12prop", value = "")
      updateTextInput(session, inputId = "p22prop", value = "")
      updateTextInput(session, inputId = "pooled2prop", value = FALSE)
      updateTextInput(session, inputId = "p1inc2prop", value = FALSE)
    } # end observeEvent(input$default2prop ...
  )
  
  # clear 2prop inputs
  observeEvent(
    input$clear2prop,
    {
      updateTextInput(session, inputId = "alpha2prop", value = "")
      updateTextInput(session, inputId = "power2prop", value = "")
      updateTextInput(session, inputId = "cv2prop", value = "")
      updateTextInput(session, inputId = "m2prop", value = "")
      updateTextInput(session, inputId = "n2prop", value = "")
      updateTextInput(session, inputId = "icc2prop", value = "")
      updateTextInput(session, inputId = "p12prop", value = "")
      updateTextInput(session, inputId = "p22prop", value = "")
      updateTextInput(session, inputId = "pooled2prop", value = FALSE)
      updateTextInput(session, inputId = "p1inc2prop", value = FALSE)
    } 
  ) # end observeEvent(input$clear2prop ...
  
  # create 2prop data
  res2prop <- eventReactive(
    input$calc2prop,
    {
      alpha <- make_sequence(isolate(input$alpha2prop))
      power <- make_sequence(isolate(input$power2prop))
      m <- make_sequence(isolate(input$m2prop))
      n <- make_sequence(isolate(input$n2prop))
      cv <- make_sequence(isolate(input$cv2prop))
      p1 <- make_sequence(isolate(input$p12prop))
      p2 <- make_sequence(isolate(input$p22prop))
      icc <- make_sequence(isolate(input$icc2prop))
      pooled <- isolate(input$pooled2prop)
      p1inc <- isolate(input$p1inc2prop)
      
      tab <- expand.grid(alpha,
                         power,
                         m,
                         n,
                         cv,
                         p1,
                         p2,
                         icc,
                         stringsAsFactors = FALSE)
      
      needind <- which(is.na(tab[1,]))
      names(tab) <- names2prop
      target <- names2prop[needind]
      
      tab[[target]] <- signif(pmap_dbl(tab, crtpwr.2prop, pooled = pooled, p1inc = p1inc),4)
      mutate_if(tab, names(tab) != target, factor)
    })
  
  # create 2prop output table
  output$table2prop = DT::renderDataTable(
    res2prop(),
    filter = 'top',
    options = list(
      columnDefs = list(list(className = 'dt-center', targets = '_all')),
      pageLength = 10,
      autoWidth = TRUE
    )
  )
  
  # setup 2prop download handler
  output$dl2prop <- downloadHandler(
    filename = function() {
      paste('data-2prop-', Sys.time(), '.csv', sep='')
    },
    content = function(con) {
      write.csv(res2prop(), file, row.names = FALSE)
    }
  )
  
  # create 2prop graph
  output$graph2prop <- renderPlot({
    if(input$color2prop){
      p <- ggplot(res2prop(),
                  aes_string(x = input$x2prop,
                             y = input$y2prop,
                             group = input$group2prop,
                             color = input$group2prop))
    } else {
      p <- ggplot(res2prop(),
                  aes_string(x = input$x2prop,
                             y = input$y2prop,
                             group = input$group2prop))
    }
    p <- p + geom_line(size = input$lsize2prop) +
      geom_point(size = input$psize2prop) +
      theme_grey(base_size = 18) 
    
    facets <- paste(input$row2prop,"~",input$col2prop)
    if(facets != ". ~ .") p <- p + facet_grid(facets)
    
    p
  },
  height = reactive({input$height2prop}))
  
  
  #----------------------------------------------------------------------------
  # Two rates
  #----------------------------------------------------------------------------
  
  # reset 2rate inputs to default values
  observeEvent(
    input$default2rate,
    {
      updateTextInput(session, inputId = "alpha2rate", value = "0.05")
      updateTextInput(session, inputId = "power2rate", value = "0.80")
      updateTextInput(session, inputId = "m2rate", value = "")
      updateTextInput(session, inputId = "py2rate", value = "")
      updateTextInput(session, inputId = "r12rate", value = "")
      updateTextInput(session, inputId = "r22rate", value = "")
      updateTextInput(session, inputId = "cvb2rate", value = "")
    } # end observeEvent(input$default2rate ...
  )
  
  # clear 2rate inputs 
  observeEvent(
    input$clear2rate,
    {
      updateTextInput(session, inputId = "alpha2rate", value = "")
      updateTextInput(session, inputId = "power2rate", value = "")
      updateTextInput(session, inputId = "m2rate", value = "")
      updateTextInput(session, inputId = "py2rate", value = "")
      updateTextInput(session, inputId = "r12rate", value = "")
      updateTextInput(session, inputId = "r22rate", value = "")
      updateTextInput(session, inputId = "cvb2rate", value = "")
    } # end observeEvent(input$clear2rate ...
  )
  
  # create 2rate data
  res2rate <- eventReactive(
    input$calc2rate,
    {
      alpha <- make_sequence(isolate(input$alpha2rate))
      power <- make_sequence(isolate(input$power2rate))
      m <- make_sequence(isolate(input$m2rate))
      py <- make_sequence(isolate(input$py2rate))
      r1 <- make_sequence(isolate(input$r12rate))
      r2 <- make_sequence(isolate(input$r22rate))
      cvb <- make_sequence(isolate(input$cvb2rate))
      
      tab <- expand.grid(alpha,
                         power,
                         m,
                         py,
                         r1,
                         r2,
                         cvb,
                         stringsAsFactors = FALSE)
      
      needind <- which(is.na(tab[1,]))
      names(tab) <- names2rate
      target <- names2rate[needind]
      
      tab[[target]] <- signif(pmap_dbl(tab, crtpwr.2rate),4)
      mutate_if(tab, names(tab) != target, factor)
    })
  
  # create 2rate output table
  output$table2rate = DT::renderDataTable(
    res2rate(),
    filter = 'top',
    options = list(
      columnDefs = list(list(className = 'dt-center', targets = '_all')),
      pageLength = 10,
      autoWidth = TRUE
    )
  )
  
  # setup 2rate download handler
  output$dl2rate <- downloadHandler(
    filename = function() {
      paste('data-2rate-', Sys.time(), '.csv', sep='')
    },
    content = function(con) {
      write.csv(res2rate(), file, row.names = FALSE)
    }
  )
  
  # create 2rate graph
  output$graph2rate <- renderPlot({
    if(input$color2rate){
      p <- ggplot(res2rate(),
                  aes_string(x = input$x2rate,
                             y = input$y2rate,
                             group = input$group2rate,
                             color = input$group2rate))
    } else {
      p <- ggplot(res2rate(),
                  aes_string(x = input$x2rate,
                             y = input$y2rate,
                             group = input$group2rate))
    }
    p <- p + geom_line(size = input$lsize2rate) +
      geom_point(size = input$psize2rate) +
      theme_grey(base_size = 18) 
    
    facets <- paste(input$row2rate,"~",input$col2rate)
    if(facets != ". ~ .") p <- p + facet_grid(facets)
    
    p
  },
  height = reactive({input$height2rate}))
  
  
}


shinyApp(ui = ui, server = server)