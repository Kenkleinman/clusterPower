library(shiny)
library(tidyverse)
library(stringr)

# source("crtpwr.2mean.R")
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
  navlistPanel(
    tabPanel("Two-Arm",
             tabsetPanel(
               tabPanel("Continous",
                        column(12,
                               fluidRow(
                                 column(4, textInput("alpha2mean", HTML("Type I error rate, &alpha; (alpha)"),
                                                     value = "0.05")),
                                 column(4, textInput("power2mean", "Power (power)",
                                                     value = "0.80")),
                                 column(4, textInput("cv2mean", "Cluster size CV (cv)",
                                                     value = "0"))
                               ),
                               fluidRow(
                                 column(4, textInput("m2mean", "Clusters per arm (m)",
                                                     value = "")),
                                 column(4, textInput("n2mean", "Cluster size (n)",
                                                     value = "")),
                                 column(4, textInput("d2mean", "Difference (d)",
                                                     value = ""))
                               ),
                               fluidRow(
                                 column(4, textInput("icc2mean", "ICC (icc)",
                                                     value = "")),
                                 column(4, textInput("varw2mean", "Within-cluster variance (varw)",
                                                     value = ""))
                               ),
                               fluidRow(
                                 column(4, actionButton("default2mean", "Defaults")),
                                 column(4, actionButton("clear2mean", "Clear All")),
                                 column(4, actionButton("calc2mean", "Calculate"))),
                               fluidRow(
                                 dataTableOutput("table2mean")
                               )
                        ) # end column(12 ...
               ), # end tabPanel("Continuous ...
               tabPanel("Binary",
                        column(12,
                               fluidRow(
                                 column(4, textInput("alpha2prop", HTML("Type I error rate, &alpha; (alpha)"),
                                                     value = "0.05")),
                                 column(4, textInput("power2prop", "Power (power)",
                                                     value = "0.80")),
                                 column(4, textInput("cv2prop", "Cluster size CV (cv)",
                                                     value = "0"))
                               ),
                               fluidRow(
                                 column(4, textInput("m2prop", "Clusters per arm (m)",
                                                     value = "")),
                                 column(4, textInput("n2prop", "Cluster size (n)",
                                                     value = "")),
                                 column(4, textInput("icc2prop", "ICC (icc)",
                                                     value = ""))
                               ),
                               fluidRow(
                                 column(4, textInput("p12prop", "Treatment Proportion (p1)",
                                                     value = "")),
                                 column(4, textInput("p22prop", "Control Proportion (p2)",
                                                     value = "")),
                                 column(4, 
                                        fluidRow(
                                          column(12,
                                                 checkboxInput("pooled2prop", "Pooled")),
                                          column(12,
                                                 checkboxInput("p1inc2prop", "Expect increase? (p1 > p2)"))
                                        )
                                 )),
                               fluidRow(
                                 column(4, actionButton("default2prop", "Defaults")),
                                 column(4, actionButton("clear2prop", "Clear All")),
                                 column(4, actionButton("calc2prop", "Calculate"))),
                               fluidRow(
                                 dataTableOutput("table2prop")
                               )
                        ) # end column(12 ...
               ), # end tabPanel("Binary ...
               tabPanel("Count",
                        column(12,
                               fluidRow(
                                 column(4, textInput("alpha2rate", HTML("Type I error rate, &alpha; (alpha)"),
                                                     value = "0.05")),
                                 column(4, textInput("power2rate", "Power (power)",
                                                     value = "0.80")),
                                 column(4, textInput("m2rate", "Clusters per arm (m)",
                                                     value = ""))
                               ),
                               fluidRow(
                                 column(4, textInput("py2rate", "Person-years per cluster (py)",
                                                     value = "")),
                                 column(4, textInput("r12rate", "Treatment rate (r1)",
                                                     value = "")),
                                 column(4, textInput("r22rate", "Control rate (r2)",
                                                     value = ""))
                               ),
                               fluidRow(
                                 column(4, textInput("cvb2rate", "Between-cluster CV (cvb)",
                                                     value = ""))
                               ),
                               fluidRow(
                                 column(4, actionButton("default2rate", "Defaults")),
                                 column(4, actionButton("clear2rate", "Clear All")),
                                 column(4, actionButton("calc2rate", "Calculate"))),
                               fluidRow(
                                 dataTableOutput("table2rate")
                               )
                        ) # end column(12 ...
               ) # end tabPanel("Count ...
             )
    ),
    tabPanel("Other stuff",
             "Under construction ..."
    )
  ) # end navlistPanel
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
  
  
}


shinyApp(ui = ui, server = server)