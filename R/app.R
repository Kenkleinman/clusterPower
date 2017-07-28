library(shiny)
library(tidyverse)
library(stringr)

# source("crtpwr.2mean.R")
# source("crtpwr.2prop.R")
# source("crtpwr.2prop.R")

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

ui <- fluidPage(
  HTML("<h3>Simple Two-Arm Designs</h3>
        <p>To use the two-arm calculator, leave the desired quantity blank and entering values for the other quantities.</p>
        <p>You may specify more than one input quantity by separating numbers with spaces or commas.<p>
        <p>You may specify a sequence of values by typing 'from X to Y by Z', where 'X' is the starting value, 'Y' is the ending value, and 'Z' is the increment.</p>"),
  tabsetPanel(
    tabPanel("Continous",
             column(3,
                    fluidRow(textInput("alpha2mean", HTML("Type I rate, &alpha; (alpha)"),
                                          value = "0.05", width = "100%")),
                    fluidRow(textInput("power2mean", "Power (power)",
                                          value = "0.80", width = "100%")),
                    fluidRow(textInput("cv2mean", "Cluster size CV (cv)",
                                          value = "0", width = "100%")),
                    fluidRow(textInput("d2mean", "Difference (d)",
                                          value = "", width = "100%")),
                    fluidRow(textInput("m2mean", "Clusters per arm (m)",
                                          value = "", width = "100%")),
                    fluidRow(textInput("n2mean", "Cluster size (n)",
                                          value = "", width = "100%")),
                    fluidRow(textInput("icc2mean", "ICC (icc)",
                                          value = "", width = "100%")),
                    fluidRow(textInput("varw2mean", "Within variance (varw)",
                                          value = "", width = "100%")),
                    fluidRow(
                      column(6, actionButton("default2mean", "Defaults", width = "100%")),
                      column(6, actionButton("clear2mean", "Clear All", width = "100%"))
                      ),
                    fluidRow(
                      column(12, actionButton("calc2mean", "Calculate", width = "100%"))
                    )
             ),
             column(9,
                    fluidRow(
                      conditionalPanel(condition = "output.table2mean != null",
                                       downloadButton("dl2mean", "Download"))
                      
                    ),
                    dataTableOutput("table2mean")
             )
    ), # end tabPanel("Continuous ...
    tabPanel("Binary",
             column(3,
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
                      column(6,checkboxInput("p1inc2prop", "Expect p1 > p2?")),
                      column(6,checkboxInput("pooled2prop", "Pooled"))
                    ),
                    fluidRow(
                      column(6, actionButton("default2prop", "Defaults", width = "100%")),
                      column(6, actionButton("clear2prop", "Clear All", width = "100%"))
                    ),
                    fluidRow(
                      column(12, actionButton("calc2prop", "Calculate", width = "100%"))
                    )
             ),
             column(9,
                    fluidRow(
                      conditionalPanel(condition = "output.table2prop != null",
                                       downloadButton("dl2prop", "Download"))
                      
                    ),
                    dataTableOutput("table2prop")
             )
    ), # end tabPanel("Binary ...
    tabPanel("Count",
             column(3,
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
                      column(6, actionButton("default2rate", "Defaults", width = "100%")),
                      column(6, actionButton("clear2rate", "Clear All", width = "100%"))),
                    fluidRow(
                      column(12, actionButton("calc2rate", "Calculate", width = "100%"))
                    )
             ),
             column(9,
                    fluidRow(
                      conditionalPanel(condition = "output.table2rate != null",
                                       downloadButton("dl2rate", "Download"))
                      
                    ),
                    dataTableOutput("table2rate")
             )
    ) # end tabPanel("Count ...
  )
) # end fluidPage

server <- function(input, output, session){
  
  # observeEvent to reset 2mean inputs to default values
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
  
  # observeEvent to reset 2prop inputs to default values
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
  
  # observeEvent to reset 2rate inputs to default values
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
  
  # observeEvent to clear 2mean inputs 
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
  
  # observeEvent to clear 2prop inputs
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
  
  
  # observeEvent to clear 2rate inputs 
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
  
  #----------------------------------------------------------------------------
  # Two means
  #----------------------------------------------------------------------------
  
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
      names(tab) <- names(as.list(args(crtpwr.2mean)))[1:8]
      target <- names(tab)[needind]
      
      tab[[target]] <- signif(pmap_dbl(tab, crtpwr.2mean),4)
      as.data.frame(tab)
    })
  
  output$table2mean = renderDataTable({
    res2mean()
  })
  
  output$dl2mean <- downloadHandler(
    filename = function() {
      paste('data-2mean-', Sys.time(), '.csv', sep='')
    },
    content = function(file) {
      write.csv(res2mean(), file)
    }
  )
  
  
  #----------------------------------------------------------------------------
  # Two proportions
  #----------------------------------------------------------------------------
  
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
      names(tab) <- names(as.list(args(crtpwr.2prop)))[1:8]
      target <- names(tab)[needind]
      
      tab[[target]] <- signif(pmap_dbl(tab, crtpwr.2prop, pooled = pooled, p1inc = p1inc),4)
      as.data.frame(tab)[1:8]
    })
  
  output$table2prop = renderDataTable({
    res2prop()
  })
  
  output$dl2prop <- downloadHandler(
    filename = function() {
      paste('data-2prop-', Sys.time(), '.csv', sep='')
    },
    content = function(con) {
      write.csv(data, con)
    }
  )
  
  
  #----------------------------------------------------------------------------
  # Two rates
  #----------------------------------------------------------------------------
  
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
      names(tab) <- names(as.list(args(crtpwr.2rate)))[1:7]
      target <- names(tab)[needind]
      
      tab[[target]] <- signif(pmap_dbl(tab, crtpwr.2rate),4)
      as.data.frame(tab)
    })
  
  output$table2rate = renderDataTable({
    res2rate()
  })
  
  output$dl2rate <- downloadHandler(
    filename = function() {
      paste('data-2rate-', Sys.time(), '.csv', sep='')
    },
    content = function(con) {
      write.csv(data, con)
    }
  )
  
  
}


shinyApp(ui = ui, server = server)