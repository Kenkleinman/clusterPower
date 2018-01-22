library(shiny)
library(shinyBS)
library(DT)
library(tidyverse)
library(stringr)
library(clusterPower)

source("labels.R")
source("helpers.R")
source("builders.R")

# vectors of names, needed for graphs and target param selection
names2mean <- c("alpha","power","m","n","cv","d","icc","varw","method")
names2meanD <- c("alpha","power","m","n","d","icc","rho_c","rho_s","varw")
names2meanM <- c("alpha","power","m","n","d","icc","varw","rho_m")
names2prop <- c("alpha","power","m","n","cv","p1","p2","icc","pooled","p1inc")
names2propD <- c("alpha","power","m","n","p","d","icc","rho_c","rho_s","covdf","pvar_c","pvar_s")
names2propM <- c("alpha","power","m","n","p1","p2","cvm","p1inc")
names2rate <- c("alpha","power","m","py","r1","r2","cvb","r1inc")

umass <- "font-family: 'Open Sans', Helvetica, Arial, sans-serif; font-weight: bold; color: #ffffff; background-color: #881c1c; border: 3px solid #000000;"

ui <- fluidPage(
  HTML("<h3>Simple Two-Arm Designs</h3>
        <p>To use the two-arm calculator, leave the desired quantity blank and enter values for the other quantities.</p>
        <p>You may specify more than one input quantity by separating numbers with spaces or commas.<p>
        <p>You may specify a sequence of values by typing 'X to Y by Z', where 'X' is the starting value, 'Y' is the ending value, and 'Z' is the increment.</p>"),
  HTML("This Beta has minimal documentation; please contact ken.kleinman@gmail.com with any feedback."),
  tabsetPanel(
    #-----------------------------------------------------------------------------------------------------------
    tabPanel("Continuous",
             column(2,
                    #----------------------------------------------------------
                    fluidRow(textInput("alpha2mean", alphatext,
                                       value = "0.05", width = "100%")),
                    bsTooltip("alpha2mean", alphatooltip,
                              'right', options = list(container = "body")),
                    #----------------------------------------------------------
                    fluidRow(textInput("power2mean", powertext,
                                       value = "", width = "100%")),
                    bsTooltip("power2mean", powertooltip,
                              'right', options = list(container = "body")),
                    #----------------------------------------------------------
                    fluidRow(textInput("d2mean", dtext,
                                       value = "", width = "100%")),
                    bsTooltip("d2mean", dtooltip,
                              'right', options = list(container = "body")),
                    #----------------------------------------------------------
                    fluidRow(textInput("m2mean", mtext,
                                       value = "", width = "100%")),
                    bsTooltip("m2mean", mtooltip,
                              'right', options = list(container = "body")),
                    #----------------------------------------------------------
                    fluidRow(textInput("n2mean", ntext,
                                       value = "", width = "100%")),
                    bsTooltip("n2mean", ntooltip,
                              'right', options = list(container = "body")),
                    #----------------------------------------------------------
                    fluidRow(textInput("icc2mean", icctext,
                                       value = "", width = "100%")),
                    bsTooltip("icc2mean", icctooltip,
                              'right', options = list(container = "body")),
                    #----------------------------------------------------------
                    fluidRow(textInput("varw2mean", varwtext,
                                       value = "", width = "100%")),
                    bsTooltip("varw2mean", varwtooltip,
                              'right', options = list(container = "body")),
                    #----------------------------------------------------------
                    fluidRow(textInput("cv2mean", cvtext,
                                       value = "0", width = "100%")),
                    bsTooltip("cv2mean", cvtooltip,
                              'right', options = list(container = "body")),
                    #----------------------------------------------------------
                    fluidRow(checkboxGroupInput("method2mean", methodtext,
                                                choices = c(Taylor = "taylor", Weighted = "weighted"),
                                                selected = "taylor")),
                    bsTooltip("method2mean", methodtooltip,
                              'right', options = list(container = "body")),
                    #----------------------------------------------------------
                    fluidRow(
                      column(6, style='padding:0px;', actionButton("default2mean", defaulttext, width = "100%")),
                      column(6, style='padding:0px;', actionButton("clear2mean", clearalltext, width = "100%"))
                    ),
                    fluidRow(
                      column(12, style='padding: 0px;', actionButton("calc2mean", calctext, width = "100%",
                                                                     style = umass)),
                    fluidRow(column(12, credittext))

                    )
             ),
             column(10,
                    make_table_and_graph("2mean", names2mean)
             ) # end column(10,...
    ), # end tabPanel("Continuous ...
    #-----------------------------------------------------------------------------------------------------------
    tabPanel("Continuous DID",
             column(2,
                    #----------------------------------------------------------
                    fluidRow(textInput("alpha2meanD", alphatext,
                                       value = "0.05", width = "100%")),
                    bsTooltip("alpha2meanD", alphatooltip,
                              'right', options = list(container = "body")),
                    #----------------------------------------------------------
                    fluidRow(textInput("power2meanD", powertext,
                                       value = "", width = "100%")),
                    bsTooltip("power2meanD", powertooltip,
                              'right', options = list(container = "body")),
                    #----------------------------------------------------------
                    fluidRow(textInput("d2meanD", dDtext,
                                       value = "", width = "100%")),
                    bsTooltip("d2meanD", dDtooltip,
                              'right', options = list(container = "body")),
                    #----------------------------------------------------------
                    fluidRow(textInput("m2meanD", mtext,
                                       value = "", width = "100%")),
                    bsTooltip("m2meanD", mtooltip,
                              'right', options = list(container = "body")),
                    #----------------------------------------------------------
                    fluidRow(textInput("n2meanD", ntext,
                                       value = "", width = "100%")),
                    bsTooltip("n2meanD", ntooltip,
                              'right', options = list(container = "body")),
                    #----------------------------------------------------------
                    fluidRow(textInput("icc2meanD", icctext,
                                       value = "", width = "100%")),
                    bsTooltip("icc2meanD", icctooltip,
                              'right', options = list(container = "body")),
                    #----------------------------------------------------------
                    fluidRow(textInput("rho_c2meanD", rho_ctext,
                                       value = "", width = "100%")),
                    bsTooltip("rho_c2meanD", rho_ctooltip,
                              'right', options = list(container = "body")),
                    #----------------------------------------------------------
                    fluidRow(textInput("rho_s2meanD", rho_stext,
                                       value = "", width = "100%")),
                    bsTooltip("rho_s2meanD", rho_stooltip,
                              'right', options = list(container = "body")),
                    #----------------------------------------------------------
                    fluidRow(textInput("varw2meanD", varwtext,
                                       value = "", width = "100%")),
                    bsTooltip("varw2meanD", varwtooltip,
                              'right', options = list(container = "body")),
                    #----------------------------------------------------------
                    fluidRow(
                      column(6, style='padding:0px;', actionButton("default2meanD", defaulttext, width = "100%")),
                      column(6, style='padding:0px;', actionButton("clear2meanD", clearalltext, width = "100%"))
                    ),
                    fluidRow(
                      column(12, style='padding: 0px;', actionButton("calc2meanD", calctext, width = "100%",
                                                                     style = umass)),
                      fluidRow(column(12, credittext))
                    )
             ),
             column(10,
                    make_table_and_graph("2meanD", names2meanD)
             ) # end column(10,...
    ), # end tabPanel("Continuous DID ...
    #-----------------------------------------------------------------------------------------------------------
    tabPanel("Continuous Matched",
             column(2,
                    #----------------------------------------------------------
                    fluidRow(textInput("alpha2meanM", alphatext,
                                       value = "0.05", width = "100%")),
                    bsTooltip("alpha2meanM", alphatooltip,
                              'right', options = list(container = "body")),
                    #----------------------------------------------------------
                    fluidRow(textInput("power2meanM", powertext,
                                       value = "", width = "100%")),
                    bsTooltip("power2meanM", powertooltip,
                              'right', options = list(container = "body")),
                    #----------------------------------------------------------
                    fluidRow(textInput("d2meanM", dtext,
                                       value = "", width = "100%")),
                    bsTooltip("d2meanM", dtooltip,
                              'right', options = list(container = "body")),
                    #----------------------------------------------------------
                    fluidRow(textInput("m2meanM", mtext,
                                       value = "", width = "100%")),
                    bsTooltip("m2meanM", mtooltip,
                              'right', options = list(container = "body")),
                    #----------------------------------------------------------
                    fluidRow(textInput("n2meanM", ntext,
                                       value = "", width = "100%")),
                    bsTooltip("n2meanM", ntooltip,
                              'right', options = list(container = "body")),
                    #----------------------------------------------------------
                    fluidRow(textInput("icc2meanM", icctext,
                                       value = "", width = "100%")),
                    bsTooltip("icc2meanM", icctooltip,
                              'right', options = list(container = "body")),
                    #----------------------------------------------------------
                    fluidRow(textInput("rho_m2meanM", rho_mtext,
                                       value = "", width = "100%")),
                    bsTooltip("rho_m2meanM", rho_mtooltip,
                              'right', options = list(container = "body")),
                    #----------------------------------------------------------
                    fluidRow(textInput("varw2meanM", varwtext,
                                       value = "", width = "100%")),
                    bsTooltip("varw2meanM", varwtooltip,
                              'right', options = list(container = "body")),
                    #----------------------------------------------------------
                    fluidRow(
                      column(6, style='padding:0px;', actionButton("default2meanM", defaulttext, width = "100%")),
                      column(6, style='padding:0px;', actionButton("clear2meanM", clearalltext, width = "100%"))
                    ),
                    fluidRow(
                      column(12, style='padding: 0px;', actionButton("calc2meanM", calctext, width = "100%",
                                                                     style = umass)),
                      fluidRow(column(12, credittext))
                    )
             ),
             column(10,
                    make_table_and_graph("2meanM", names2meanM)
             ) # end column(10,...
    ), # end tabPanel("Continuous Matched ...
    #-----------------------------------------------------------------------------------------------------------
    tabPanel("Binary",
             column(2,
                    #----------------------------------------------------------
                    fluidRow(textInput("alpha2prop", alphatext,
                                       value = "0.05", width = "100%")),
                    bsTooltip("alpha2prop", alphatooltip,
                              'right', options = list(container = "body")),
                    #----------------------------------------------------------
                    fluidRow(textInput("power2prop", powertext,
                                       value = "0.80", width = "100%")),
                    bsTooltip("power2prop", powertooltip,
                              'right', options = list(container = "body")),
                    #----------------------------------------------------------
                    fluidRow(textInput("p12prop", p1text,
                                       value = "", width = "100%")),
                    bsTooltip("p12prop", p1tooltip,
                              'right', options = list(container = "body")),
                    #----------------------------------------------------------
                    fluidRow(textInput("p22prop", p2text,
                                       value = "", width = "100%")),
                    bsTooltip("p22prop", p2tooltip,
                              'right', options = list(container = "body")),
                    #----------------------------------------------------------
                    fluidRow(checkboxInput("p1inc2prop", p1inctext)),
                    bsTooltip("p1inc2prop", p1inctooltip,
                              'top', options = list(container = "body")),
                    #----------------------------------------------------------
                    fluidRow(checkboxInput("pooled2prop", pooledtext)),
                    bsTooltip("pooled2prop", pooledtooltip,
                              'top', options = list(container = "body")),
                    #----------------------------------------------------------
                    fluidRow(textInput("m2prop", mtext,
                                       value = "", width = "100%")),
                    bsTooltip("m2prop", mtooltip,
                              'right', options = list(container = "body")),
                    #----------------------------------------------------------
                    fluidRow(textInput("n2prop", ntext,
                                       value = "", width = "100%")),
                    bsTooltip("n2prop", ntooltip,
                              'right', options = list(container = "body")),
                    #----------------------------------------------------------
                    fluidRow(textInput("icc2prop", icctext,
                                       value = "", width = "100%")),
                    bsTooltip("icc2prop", icctooltip,
                              'right', options = list(container = "body")),
                    #----------------------------------------------------------
                    fluidRow(textInput("cv2prop", cvtext,
                                       value = "0", width = "100%")),
                    bsTooltip("cv2prop", cvtooltip,
                              'right', options = list(container = "body")),
                    #----------------------------------------------------------
                    fluidRow(
                      column(6, style='padding:0px;', actionButton("default2prop", defaulttext, width = "100%")),
                      column(6, style='padding:0px;', actionButton("clear2prop", clearalltext, width = "100%"))
                    ),
                    fluidRow(
                      column(12, style='padding:0px;', actionButton("calc2prop", calctext, width = "100%",
                                                                    style = umass)),
                    fluidRow(column(12, credittext))
                      
                    )
             ), # end column(2, ..
             column(10,
                    make_table_and_graph("2prop", names2prop)
             ) # end column(10,...
    ), # end tabPanel("Binary ...
    #-----------------------------------------------------------------------------------------------------------
    tabPanel("Binary DID",
             column(2,
                    #----------------------------------------------------------
                    fluidRow(textInput("alpha2propD", alphatext,
                                       value = "0.05", width = "100%")),
                    bsTooltip("alpha2propD", alphatooltip,
                              'right', options = list(container = "body")),
                    #----------------------------------------------------------
                    fluidRow(textInput("power2propD", powertext,
                                       value = "", width = "100%")),
                    bsTooltip("power2propD", powertooltip,
                              'right', options = list(container = "body")),
                    #----------------------------------------------------------
                    fluidRow(textInput("p2propD", ptext,
                                       value = "", width = "100%")),
                    bsTooltip("p2propD", ptooltip,
                              'right', options = list(container = "body")),
                    #----------------------------------------------------------
                    fluidRow(textInput("d2propD", dDtext,
                                       value = "", width = "100%")),
                    bsTooltip("d2propD", dDtooltip,
                              'right', options = list(container = "body")),
                    #----------------------------------------------------------
                    fluidRow(textInput("m2propD", mtext,
                                       value = "", width = "100%")),
                    bsTooltip("m2propD", mtooltip,
                              'right', options = list(container = "body")),
                    #----------------------------------------------------------
                    fluidRow(textInput("n2propD", ntext,
                                       value = "", width = "100%")),
                    bsTooltip("n2propD", ntooltip,
                              'right', options = list(container = "body")),
                    #----------------------------------------------------------
                    fluidRow(textInput("icc2propD", icctext,
                                       value = "", width = "100%")),
                    bsTooltip("icc2propD", icctooltip,
                              'right', options = list(container = "body")),
                    #----------------------------------------------------------
                    fluidRow(textInput("rho_c2propD", rho_ctext,
                                       value = "", width = "100%")),
                    bsTooltip("rho_c2propD", rho_ctooltip,
                              'right', options = list(container = "body")),
                    #----------------------------------------------------------
                    fluidRow(textInput("rho_s2propD", rho_stext,
                                       value = "", width = "100%")),
                    bsTooltip("rho_s2propD", rho_stooltip,
                              'right', options = list(container = "body")),
                    #----------------------------------------------------------
                    fluidRow(textInput("covdf2propD", covdftext,
                                       value = "0", width = "100%")),
                    bsTooltip("covdf2propD", covdftooltip,
                              'right', options = list(container = "body")),
                    #----------------------------------------------------------
                    fluidRow(textInput("pvar_c2propD", pvar_ctext,
                                       value = "0", width = "100%")),
                    bsTooltip("pvar_c2propD", pvar_ctooltip,
                              'right', options = list(container = "body")),
                    #----------------------------------------------------------
                    fluidRow(textInput("pvar_s2propD", pvar_stext,
                                       value = "0", width = "100%")),
                    bsTooltip("pvar_s2propD", pvar_stooltip,
                              'right', options = list(container = "body")),
                    #----------------------------------------------------------
                    fluidRow(
                      column(6, style='padding:0px;', actionButton("default2propD", defaulttext, width = "100%")),
                      column(6, style='padding:0px;', actionButton("clear2propD", clearalltext, width = "100%"))
                    ),
                    fluidRow(
                      column(12, style='padding: 0px;', actionButton("calc2propD", calctext, width = "100%",
                                                                     style = umass)),
                      fluidRow(column(12, credittext))
                    )
             ),
             column(10,
                    make_table_and_graph("2propD", names2propD)
             ) # end column(10,...
    ), # end tabPanel("Binary DID ...
    #-----------------------------------------------------------------------------------------------------------
    tabPanel("Binary Matched",
             column(2,
                    #----------------------------------------------------------
                    fluidRow(textInput("alpha2propM", alphatext,
                                       value = "0.05", width = "100%")),
                    bsTooltip("alpha2propM", alphatooltip,
                              'right', options = list(container = "body")),
                    #----------------------------------------------------------
                    fluidRow(textInput("power2propM", powertext,
                                       value = "0.80", width = "100%")),
                    bsTooltip("power2propM", powertooltip,
                              'right', options = list(container = "body")),
                    #----------------------------------------------------------
                    fluidRow(textInput("p12propM", p1text,
                                       value = "", width = "100%")),
                    bsTooltip("p12propM", p1tooltip,
                              'right', options = list(container = "body")),
                    #----------------------------------------------------------
                    fluidRow(textInput("p22propM", p2text,
                                       value = "", width = "100%")),
                    bsTooltip("p22propM", p2tooltip,
                              'right', options = list(container = "body")),
                    #----------------------------------------------------------
                    fluidRow(checkboxInput("p1inc2propM", p1inctext)),
                    bsTooltip("p1inc2propM", p1inctooltip,
                              'top', options = list(container = "body")),
                    #----------------------------------------------------------
                    fluidRow(textInput("m2propM", mtext,
                                       value = "", width = "100%")),
                    bsTooltip("m2propM", mtooltip,
                              'right', options = list(container = "body")),
                    #----------------------------------------------------------
                    fluidRow(textInput("n2propM", ntext,
                                       value = "", width = "100%")),
                    bsTooltip("n2propM", ntooltip,
                              'right', options = list(container = "body")),
                    #----------------------------------------------------------
                    fluidRow(textInput("cvm2propM", cvmtext,
                                       value = "", width = "100%")),
                    bsTooltip("cvm2propM", cvmtooltip,
                              'right', options = list(container = "body")),
                    #----------------------------------------------------------
                    fluidRow(
                      column(6, style='padding:0px;', actionButton("default2propM", defaulttext, width = "100%")),
                      column(6, style='padding:0px;', actionButton("clear2propM", clearalltext, width = "100%"))
                    ),
                    fluidRow(
                      column(12, style='padding:0px;', actionButton("calc2propM", calctext, width = "100%",
                                                                    style = umass)),
                      fluidRow(column(12, credittext))
                      
                    )
             ), # end column(2, ..
             column(10,
                    make_table_and_graph("2propM", names2propM)
             ) # end column(10,...
    ), # end tabPanel("Binary Matched ...
    #-----------------------------------------------------------------------------------------------------------
    tabPanel("Count",
             column(2,
                    #----------------------------------------------------------
                    fluidRow(textInput("alpha2rate", alphatext,
                                       value = "0.05", width = "100%")),
                    bsTooltip("alpha2rate", alphatooltip,
                              'right', options = list(container = "body")),
                    #----------------------------------------------------------
                    fluidRow(textInput("power2rate", powertext,
                                       value = "0.80", width = "100%")),
                    bsTooltip("power2rate", powertooltip,
                              'right', options = list(container = "body")),
                    #----------------------------------------------------------
                    fluidRow(textInput("r12rate", r1text,
                                       value = "", width = "100%")),
                    bsTooltip("r12rate", r1tooltip,
                              'right', options = list(container = "body")),
                    #----------------------------------------------------------
                    fluidRow(textInput("r22rate", r2text,
                                       value = "", width = "100%")),
                    bsTooltip("r22rate", r2tooltip,
                              'right', options = list(container = "body")),
                    #----------------------------------------------------------
                    fluidRow(checkboxInput("r1inc2rate", r1inctext)),
                    bsTooltip("r1inc2rate", r1inctooltip,
                              'top', options = list(container = "body")),
                    #----------------------------------------------------------
                    fluidRow(textInput("m2rate", mtext,
                                       value = "", width = "100%")),
                    bsTooltip("m2rate", mtooltip,
                              'right', options = list(container = "body")),
                    #----------------------------------------------------------
                    fluidRow(textInput("py2rate", pytext,
                                       value = "", width = "100%")),
                    bsTooltip("py2rate", pytooltip,
                              'right', options = list(container = "body")),
                    #----------------------------------------------------------
                    fluidRow(textInput("cvb2rate", cvbtext,
                                       value = "", width = "100%")),
                    bsTooltip("cvb2rate", cvbtooltip,
                              'right', options = list(container = "body")),
                    #----------------------------------------------------------
                    fluidRow(
                      column(6, style='padding:0px;', actionButton("default2rate", defaulttext, width = "100%")),
                      column(6, style='padding:0px;', actionButton("clear2rate", clearalltext, width = "100%"))),
                    fluidRow(
                      column(12, style='padding:0px;', actionButton("calc2rate", calctext, width = "100%",
                                                                    style = umass)),
                    fluidRow(column(12,credittext))
                    )
             ),
             column(10,
                    make_table_and_graph("2rate", names2rate)
             ) # end column(10,...
    ) # end tabPanel("Count ...
    #-----------------------------------------------------------------------------------------------------------
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
      updateTextInput(session, inputId = "alpha2mean", value = "0.05")
      updateTextInput(session, inputId = "power2mean", value = "")
      updateTextInput(session, inputId = "cv2mean", value = "0")
      updateCheckboxGroupInput(session, inputId = "method2mean", selected = "taylor")
      updateTextInput(session, inputId = "m2mean", value = "")
      updateTextInput(session, inputId = "n2mean", value = "")
      updateTextInput(session, inputId = "d2mean", value = "")
      updateTextInput(session, inputId = "icc2mean", value = "")
      updateTextInput(session, inputId = "varw2mean", value = "")
    }
  ) # end observeEvent(input$default2mean ...
  
  
  # clear 2mean inputs 
  observeEvent(
    input$clear2mean,
    {
      updateTextInput(session, inputId = "alpha2mean", value = "")
      updateTextInput(session, inputId = "power2mean", value = "")
      updateTextInput(session, inputId = "cv2mean", value = "")
      updateCheckboxGroupInput(session, inputId = "method2mean",
                               choices = c(Taylor = "taylor", Weighted = "weighted"))
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
      # convert inputs to numeric vectors
      alpha <- make_sequence(isolate(input$alpha2mean))
      power <- make_sequence(isolate(input$power2mean))
      m <- make_sequence(isolate(input$m2mean))
      n <- make_sequence(isolate(input$n2mean))
      cv <- make_sequence(isolate(input$cv2mean))
      d <- make_sequence(isolate(input$d2mean))
      icc <- make_sequence(isolate(input$icc2mean))
      varw <- make_sequence(isolate(input$varw2mean))
      method <- na.omit(isolate(input$method2mean))
      
      if(!is.na(power)){
        validate(
          need(power >= 0 & power <= 1,
               powervalidmsg)
        )
      }

      if(!is.na(alpha)){
        validate(
          need(alpha >= 0 & alpha <= 1,
               alphavalidmsg)
        )
      }

      
      # create a table of input values
      tab <- expand.grid(alpha,
                         power,
                         m,
                         n,
                         cv,
                         d,
                         icc,
                         varw,
                         method,
                         stringsAsFactors = FALSE)
      
      # record the column index of the target parameter
      needind <- which(is.na(tab[1,]))
      # validate that only one input is blank
      validate(
        need(length(needind) == 1,
             "Exactly one of 'alpha', 'power', 'd', 'm', 'n', 'icc', 'varw', and 'cv' must be left blank."
        )
      )
      names(tab) <- names2mean
      target <- names2mean[needind]
      
      # apply function over table of input values
      temp <-pmap_df(tab, crtpwr.2mean.safe)
      
      tab[[target]] <- signif(temp$result, 4)
      tab$error <- map_chr(temp$error, shorten_error, target = target)
      
      # make a column to store target variable for use in graphing
      tab$target <- target
      
      # convert all input values to factors for ggplot
      mutate_if(tab, !(names(tab) %in% c(target,"error","target")), factor)
    })
  
  # create 2mean output table
  output$table2mean <- DT::renderDataTable(
    res2mean()[,1:10],
    server = FALSE,
    extensions = 'Buttons',
    filter = 'top',
    options = list(
      # create the button
      dom = 'fBrtlip',
      buttons = list(list(extend = 'csv', filename = paste('data-2mean-', Sys.time(), sep=''), text = 'Download')),
      autoWidth = TRUE,
      columnDefs = list(list(className = 'dt-center', targets = '_all'),
                        list(width = '500px', targets = 1:9)),
      lengthMenu = list(c(10, 25, 100, -1), list('10','25','100','All')),
      pageLength = 10
    )
  )

  # update graph UI
  observeEvent(res2mean(),
               {
                 update_graph_ui(session, res2mean(), "2mean", names2mean)
               })
  
  # create 2mean graph
  output$graph2mean <- renderPlot({
    create_graph(res2mean(), input$x2mean, input$y2mean, input$group2mean,
                 input$lsize2mean, input$psize2mean, input$row2mean, input$col2mean)
  },
  height = reactive({input$height2mean})
  )
  
  #----------------------------------------------------------------------------
  # Two means: Difference in difference
  #----------------------------------------------------------------------------
  
  # reset 2meanD inputs to default values
  observeEvent(
    input$default2meanD,
    {
      updateTextInput(session, inputId = "alpha2meanD", value = "0.05")
      updateTextInput(session, inputId = "power2meanD", value = "")
      updateTextInput(session, inputId = "rho_c2meanD", value = "")
      updateTextInput(session, inputId = "rho_s2meanD", value = "")
      updateTextInput(session, inputId = "m2meanD", value = "")
      updateTextInput(session, inputId = "n2meanD", value = "")
      updateTextInput(session, inputId = "d2meanD", value = "")
      updateTextInput(session, inputId = "icc2meanD", value = "")
      updateTextInput(session, inputId = "varw2meanD", value = "")
    }
  ) # end observeEvent(input$default2meanD ...
  
  
  # clear 2meanD inputs 
  observeEvent(
    input$clear2meanD,
    {
      updateTextInput(session, inputId = "alpha2meanD", value = "")
      updateTextInput(session, inputId = "power2meanD", value = "")
      updateTextInput(session, inputId = "rho_c2meanD", value = "")
      updateTextInput(session, inputId = "rho_s2meanD", value = "")
      updateTextInput(session, inputId = "m2meanD", value = "")
      updateTextInput(session, inputId = "n2meanD", value = "")
      updateTextInput(session, inputId = "d2meanD", value = "")
      updateTextInput(session, inputId = "icc2meanD", value = "")
      updateTextInput(session, inputId = "varw2meanD", value = "")
    }
  ) # end observeEvent(input$clear2meanD ...
  
  # create 2meanD data
  res2meanD <- eventReactive(
    input$calc2meanD,
    {
      # convert inputs to numeric vectors
      alpha <- make_sequence(isolate(input$alpha2meanD))
      power <- make_sequence(isolate(input$power2meanD))
      m <- make_sequence(isolate(input$m2meanD))
      n <- make_sequence(isolate(input$n2meanD))
      rho_c <- make_sequence(isolate(input$rho_c2meanD))
      rho_s <- make_sequence(isolate(input$rho_s2meanD))
      d <- make_sequence(isolate(input$d2meanD))
      icc <- make_sequence(isolate(input$icc2meanD))
      varw <- make_sequence(isolate(input$varw2meanD))
      
      if(!is.na(power)){
        validate(
          need(power >= 0 & power <= 1,
               powervalidmsg)
        )
      }
      
      if(!is.na(alpha)){
        validate(
          need(alpha >= 0 & alpha <= 1,
               alphavalidmsg)
        )
      }
      
      
      # create a table of input values
      tab <- expand.grid(alpha,
                         power,
                         m,
                         n,
                         d,
                         icc,
                         rho_c,
                         rho_s,
                         varw,
                         stringsAsFactors = FALSE)
      
      # record the column index of the target parameter
      needind <- which(is.na(tab[1,]))
      # validate that only one input is blank
      validate(
        need(length(needind) == 1,
             "Exactly one of 'alpha', 'power', 'd', 'm', 'n', 'icc', 'rho_c', 'rho_s', and 'varw' must be left blank."
        )
      )
      names(tab) <- names2meanD
      target <- names2meanD[needind]
      
      # apply function over table of input values
      temp <-pmap_df(tab, crtpwr.2meanD.safe)
      
      tab[[target]] <- signif(temp$result, 4)
      tab$error <- map_chr(temp$error, shorten_error, target = target)
      
      # make a column to store target variable for use in graphing
      tab$target <- target
      
      # check to see if there are errors, if so, set maxcol to 10 to show error in datatable
      # if not, set colmax to 9 so that error column not displayed
      #tab$colmax <- ifelse(sum(!is.na(res2meanD()$error) != 0, 10, 9))
      
      # convert all input values to factors for ggplot
      mutate_if(tab, !(names(tab) %in% c(target,"error","target")), factor)
    })
  
  # create 2meanD output table
  output$table2meanD <- DT::renderDataTable(
    res2meanD()[, 1:10],
    server = FALSE,
    extensions = 'Buttons',
    filter = 'top',
    options = list(
      # create the button
      dom = 'fBrtlip',
      buttons = list(list(extend = 'csv', filename = paste('data-2meanD-', Sys.time(), sep=''), text = 'Download')),
      autoWidth = TRUE,
      columnDefs = list(list(className = 'dt-center', targets = '_all'),
                        list(width = '500px', targets = 10)),
      pageLength = 10,
      lengthMenu =  list(c(10, 25, 100, -1), list('10', '25', '100', 'All')) 
    )
  )
  
  # update graph UI
  observeEvent(res2meanD(),
               {
                 update_graph_ui(session, res2meanD(), "2meanD", names2meanD)
               })
  
  # create 2meanD graph
  output$graph2meanD <- renderPlot({
    create_graph(res2meanD(), input$x2meanD, input$y2meanD, input$group2meanD,
                 input$lsize2meanD, input$psize2meanD, input$row2meanD, input$col2meanD)
  },
  height = reactive({input$height2meanD})
  )
  
  
  #----------------------------------------------------------------------------
  # Two means: Matching
  #----------------------------------------------------------------------------
  
  # reset 2meanM inputs to default values
  observeEvent(
    input$default2meanM,
    {
      updateTextInput(session, inputId = "alpha2meanM", value = "0.05")
      updateTextInput(session, inputId = "power2meanM", value = "")
      updateTextInput(session, inputId = "rho_m2meanM", value = "")
      updateTextInput(session, inputId = "m2meanM", value = "")
      updateTextInput(session, inputId = "n2meanM", value = "")
      updateTextInput(session, inputId = "d2meanM", value = "")
      updateTextInput(session, inputId = "icc2meanM", value = "")
      updateTextInput(session, inputId = "varw2meanM", value = "")
    }
  ) # end observeEvent(input$default2meanM ...
  
  
  # clear 2meanM inputs 
  observeEvent(
    input$clear2meanM,
    {
      updateTextInput(session, inputId = "alpha2meanM", value = "")
      updateTextInput(session, inputId = "power2meanM", value = "")
      updateTextInput(session, inputId = "rho_m2meanM", value = "")
      updateTextInput(session, inputId = "m2meanM", value = "")
      updateTextInput(session, inputId = "n2meanM", value = "")
      updateTextInput(session, inputId = "d2meanM", value = "")
      updateTextInput(session, inputId = "icc2meanM", value = "")
      updateTextInput(session, inputId = "varw2meanM", value = "")
    }
  ) # end observeEvent(input$clear2meanM ...
  
  # create 2meanM data
  res2meanM <- eventReactive(
    input$calc2meanM,
    {
      # convert inputs to numeric vectors
      alpha <- make_sequence(isolate(input$alpha2meanM))
      power <- make_sequence(isolate(input$power2meanM))
      m <- make_sequence(isolate(input$m2meanM))
      n <- make_sequence(isolate(input$n2meanM))
      rho_m <- make_sequence(isolate(input$rho_m2meanM))
      d <- make_sequence(isolate(input$d2meanM))
      icc <- make_sequence(isolate(input$icc2meanM))
      varw <- make_sequence(isolate(input$varw2meanM))
      
      if(!is.na(power)){
        validate(
          need(power >= 0 & power <= 1,
               powervalidmsg)
        )
      }
      
      if(!is.na(alpha)){
        validate(
          need(alpha >= 0 & alpha <= 1,
               alphavalidmsg)
        )
      }
      
      
      # create a table of input values
      tab <- expand.grid(alpha,
                         power,
                         m,
                         n,
                         d,
                         icc,
                         varw,
                         rho_m,
                         stringsAsFactors = FALSE)
      
      # record the column index of the target parameter
      needind <- which(is.na(tab[1,]))
      # validate that only one input is blank
      validate(
        need(length(needind) == 1,
             "Exactly one of 'alpha', 'power', 'd', 'm', 'n', 'icc', 'rho_m', and 'varw' must be left blank."
        )
      )
      names(tab) <- names2meanM
      target <- names2meanM[needind]
      
      # apply function over table of input values
      temp <-pmap_df(tab, crtpwr.2meanM.safe)
      
      tab[[target]] <- signif(temp$result, 4)
      tab$error <- map_chr(temp$error, shorten_error, target = target)
      
      # make a column to store target variable for use in graphing
      tab$target <- target
      
      # check to see if there are errors, if so, set maxcol to 10 to show error in datatable
      # if not, set colmax to 9 so that error column not displayed
      #tab$colmax <- ifelse(sum(!is.na(res2meanM()$error) != 0, 10, 9))
      
      # convert all input values to factors for ggplot
      mutate_if(tab, !(names(tab) %in% c(target,"error","target")), factor)
    })
  
  # create 2meanM output table
  output$table2meanM <- DT::renderDataTable(
    res2meanM()[, 1:10],
    server = FALSE,
    extensions = 'Buttons',
    filter = 'top',
    options = list(
      # create the button
      dom = 'fBrtlip',
      buttons = list(list(extend = 'csv', filename = paste('data-2meanM-', Sys.time(), sep=''), text = 'Download')),
      autoWidth = TRUE,
      columnDefs = list(list(className = 'dt-center', targets = '_all'),
                        list(width = '500px', targets = 10)),
      pageLength = 10,
      lengthMenu =  list(c(10, 25, 100, -1), list('10', '25', '100', 'All')) 
    )
  )
  
  # update graph UI
  observeEvent(res2meanM(),
               {
                 update_graph_ui(session, res2meanM(), "2meanM", names2meanM)
               })
  
  # create 2meanM graph
  output$graph2meanM <- renderPlot({
    create_graph(res2meanM(), input$x2meanM, input$y2meanM, input$group2meanM,
                 input$lsize2meanM, input$psize2meanM, input$row2meanM, input$col2meanM)
  },
  height = reactive({input$height2meanM})
  )
  
  
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
                         pooled,
                         p1inc,
                         stringsAsFactors = FALSE)
      
      if(!is.na(power)){
        validate(
          need(power >= 0 & power <= 1,
               powervalidmsg)
        )
      }
      
      if(!is.na(alpha)){
        validate(
          need(alpha >= 0 & alpha <= 1,
               alphavalidmsg)
        )
      }
      
      # record the column index of the target parameter
      needind <- which(is.na(tab[1,]))
      # validate that only one input is blank
      validate(
        need(length(needind) == 1,
             "Exactly one of 'alpha', 'power', 'p1', 'p2', 'm', 'n', 'icc' or 'cv' must be left blank."
        )
      )
      names(tab) <- names2prop
      target <- names2prop[needind]
      
      # apply function over table of input values
      temp <-pmap_df(tab, crtpwr.2prop.safe)
      
      tab[[target]] <- signif(temp$result, 4)
      tab$error <- map_chr(temp$error, shorten_error, target = target)
      
      # make a column to store target variable for use in graphing
      tab$target <- target
      
      # convert all input values to factors for ggplot
      mutate_if(tab, !(names(tab) %in% c(target,"error","target")), factor)
    })
  
  # create 2prop output table
  output$table2prop <- DT::renderDataTable(
    res2prop()[,c(1:8,11)],
    server = FALSE,
    extensions = 'Buttons',
    filter = 'top',
    options = list(
      # create the button
      dom = 'fBrtlip',
      buttons = list(list(extend = 'csv', filename = paste('data-2prop-', Sys.time(), sep=''), text = 'Download')),
      autoWidth = TRUE,
      columnDefs = list(list(className = 'dt-center', targets = '_all'),
                        list(width = '500px', targets = 9)),
      pageLength = 10,
      lengthMenu =  list(c(10, 25, 100, -1), list('10', '25', '100', 'All')) 
    )
  )

  # update graph UI
  observeEvent(res2prop(),
               {
                 update_graph_ui(session, res2prop(), "2prop", names2prop)
               })
  
  # create 2prop graph
  output$graph2prop <- renderPlot({
    create_graph(res2prop(), input$x2prop, input$y2prop, input$group2prop,
                 input$lsize2prop, input$psize2prop, input$row2prop, input$col2prop)
  },
  height = reactive({input$height2prop})
  )
  
  
  #----------------------------------------------------------------------------
  # Two proportions, difference in difference
  #----------------------------------------------------------------------------
  
  # reset 2propD inputs to default values
  observeEvent(
    input$default2propD,
    {
      updateTextInput(session, inputId = "alpha2propD", value = "0.05")
      updateTextInput(session, inputId = "power2propD", value = "")
      updateTextInput(session, inputId = "rho_c2propD", value = "")
      updateTextInput(session, inputId = "rho_s2propD", value = "")
      updateTextInput(session, inputId = "m2propD", value = "")
      updateTextInput(session, inputId = "n2propD", value = "")
      updateTextInput(session, inputId = "p2propD", value = "")
      updateTextInput(session, inputId = "d2propD", value = "")
      updateTextInput(session, inputId = "icc2propD", value = "")
      updateTextInput(session, inputId = "covdf2propD", value = "0")
      updateTextInput(session, inputId = "pvar_c2propD", value = "0")
      updateTextInput(session, inputId = "pvar_s2propD", value = "0")
    }
  ) # end observeEvent(input$default2propD ...
  
  
  # clear 2propD inputs 
  observeEvent(
    input$clear2propD,
    {
      updateTextInput(session, inputId = "alpha2propD", value = "")
      updateTextInput(session, inputId = "power2propD", value = "")
      updateTextInput(session, inputId = "rho_c2propD", value = "")
      updateTextInput(session, inputId = "rho_s2propD", value = "")
      updateTextInput(session, inputId = "m2propD", value = "")
      updateTextInput(session, inputId = "n2propD", value = "")
      updateTextInput(session, inputId = "p2propD", value = "")
      updateTextInput(session, inputId = "d2propD", value = "")
      updateTextInput(session, inputId = "icc2propD", value = "")
      updateTextInput(session, inputId = "covdf2propD", value = "")
      updateTextInput(session, inputId = "pvar_c2propD", value = "")
      updateTextInput(session, inputId = "pvar_s2propD", value = "")
    }
  ) # end observeEvent(input$clear2propD ...
  
  # create 2propD data
  res2propD <- eventReactive(
    input$calc2propD,
    {
      # convert inputs to numeric vectors
      alpha <- make_sequence(isolate(input$alpha2propD))
      power <- make_sequence(isolate(input$power2propD))
      m <- make_sequence(isolate(input$m2propD))
      n <- make_sequence(isolate(input$n2propD))
      rho_c <- make_sequence(isolate(input$rho_c2propD))
      rho_s <- make_sequence(isolate(input$rho_s2propD))
      p <- make_sequence(isolate(input$p2propD))
      d <- make_sequence(isolate(input$d2propD))
      icc <- make_sequence(isolate(input$icc2propD))
      covdf <- make_sequence(isolate(input$covdf2propD))
      pvar_c <- make_sequence(isolate(input$pvar_c2propD))
      pvar_s <- make_sequence(isolate(input$pvar_s2propD))

      if(!is.na(power)){
        validate(
          need(power >= 0 & power <= 1,
               powervalidmsg)
        )
      }
      
      if(!is.na(alpha)){
        validate(
          need(alpha >= 0 & alpha <= 1,
               alphavalidmsg)
        )
      }
      
      
      # create a table of input values
      tab <- expand.grid(alpha,
                         power,
                         m,
                         n,
                         p,
                         d,
                         icc,
                         rho_c,
                         rho_s,
                         covdf,
                         pvar_c,
                         pvar_s,
                         stringsAsFactors = FALSE)
      
      # record the column index of the target parameter
      needind <- which(is.na(tab[1,]))
      # validate that only one input is blank
      validate(
        need(length(needind) == 1,
             "Exactly one of 'alpha', 'power', 'p', 'd', 'm', 'n', 'icc', 'rho_c', 'rho_s', 'covdf', 'pvar_c', and 'pvar_s' must be left blank."
        )
      )
      names(tab) <- names2propD
      target <- names2propD[needind]
      
      # apply function over table of input values
      temp <-pmap_df(tab, crtpwr.2propD.safe)
      
      tab[[target]] <- signif(temp$result, 4)
      tab$error <- map_chr(temp$error, shorten_error, target = target)
      
      # make a column to store target variable for use in graphing
      tab$target <- target
      
      # check to see if there are errors, if so, set maxcol to 10 to show error in datatable
      # if not, set colmax to 9 so that error column not displayed
      #tab$colmax <- ifelse(sum(!is.na(res2propD()$error) != 0, 10, 9))
      
      # convert all input values to factors for ggplot
      mutate_if(tab, !(names(tab) %in% c(target,"error","target")), factor)
    })
  
  # create 2propD output table
  output$table2propD <- DT::renderDataTable(
    res2propD()[, c(1:13)],
    server = FALSE,
    extensions = 'Buttons',
    filter = 'top',
    options = list(
      # create the button
      dom = 'fBrtlip',
      buttons = list(list(extend = 'csv', filename = paste('data-2propD-', Sys.time(), sep=''), text = 'Download')),
      autoWidth = TRUE,
      columnDefs = list(list(className = 'dt-center', targets = '_all'),
                        list(width = '500px', targets = 13)),
      pageLength = 10,
      lengthMenu =  list(c(10, 25, 100, -1), list('10', '25', '100', 'All')) 
    )
  )
  
  # update graph UI
  observeEvent(res2propD(),
               {
                 update_graph_ui(session, res2propD(), "2propD", names2propD)
               })
  
  # create 2propD graph
  output$graph2propD <- renderPlot({
    create_graph(res2propD(), input$x2propD, input$y2propD, input$group2propD,
                 input$lsize2propD, input$psize2propD, input$row2propD, input$col2propD)
  },
  height = reactive({input$height2propD})
  )
  
  #----------------------------------------------------------------------------
  # Two proportions, matched
  #----------------------------------------------------------------------------
  
  # reset 2propM inputs to default values
  observeEvent(
    input$default2propM,
    {
      updateTextInput(session, inputId = "alpha2propM", value = "0.05")
      updateTextInput(session, inputId = "power2propM", value = "0.80")
      updateTextInput(session, inputId = "cvm2propM", value = "")
      updateTextInput(session, inputId = "m2propM", value = "")
      updateTextInput(session, inputId = "n2propM", value = "")
      updateTextInput(session, inputId = "p12propM", value = "")
      updateTextInput(session, inputId = "p22propM", value = "")
      updateTextInput(session, inputId = "p1inc2propM", value = FALSE)
    } # end observeEvent(input$default2propM ...
  )
  
  # clear 2propM inputs
  observeEvent(
    input$clear2propM,
    {
      updateTextInput(session, inputId = "alpha2propM", value = "")
      updateTextInput(session, inputId = "power2propM", value = "")
      updateTextInput(session, inputId = "cvm2propM", value = "")
      updateTextInput(session, inputId = "m2propM", value = "")
      updateTextInput(session, inputId = "n2propM", value = "")
      updateTextInput(session, inputId = "p12propM", value = "")
      updateTextInput(session, inputId = "p22propM", value = "")
      updateTextInput(session, inputId = "p1inc2propM", value = FALSE)
    } 
  ) # end observeEvent(input$clear2propM ...
  
  # create 2propM data
  res2propM <- eventReactive(
    input$calc2propM,
    {
      alpha <- make_sequence(isolate(input$alpha2propM))
      power <- make_sequence(isolate(input$power2propM))
      m <- make_sequence(isolate(input$m2propM))
      n <- make_sequence(isolate(input$n2propM))
      p1 <- make_sequence(isolate(input$p12propM))
      p2 <- make_sequence(isolate(input$p22propM))
      cvm <- make_sequence(isolate(input$cvm2propM))
      p1inc <- isolate(input$p1inc2propM)
      
      tab <- expand.grid(alpha,
                         power,
                         m,
                         n,
                         p1,
                         p2,
                         cvm,
                         p1inc,
                         stringsAsFactors = FALSE)
      
      if(!is.na(power)){
        validate(
          need(power >= 0 & power <= 1,
               powervalidmsg)
        )
      }
      
      if(!is.na(alpha)){
        validate(
          need(alpha >= 0 & alpha <= 1,
               alphavalidmsg)
        )
      }
      
      # record the column index of the target parameter
      needind <- which(is.na(tab[1,]))
      # validate that only one input is blank
      validate(
        need(length(needind) == 1,
             "Exactly one of 'alpha', 'power', 'p1', 'p2', 'm', 'n', or 'cvm' must be left blank."
        )
      )
      names(tab) <- names2propM
      target <- names2propM[needind]
      
      # apply function over table of input values
      temp <-pmap_df(tab, crtpwr.2propM.safe)
      
      tab[[target]] <- signif(temp$result, 4)
      tab$error <- map_chr(temp$error, shorten_error, target = target)
      
      # make a column to store target variable for use in graphing
      tab$target <- target
      
      # convert all input values to factors for ggplot
      mutate_if(tab, !(names(tab) %in% c(target,"error","target")), factor)
    })
  
  # create 2propM output table
  output$table2propM <- DT::renderDataTable(
    res2propM()[,c(1:7,9)],
    server = FALSE,
    extensions = 'Buttons',
    filter = 'top',
    options = list(
      # create the button
      dom = 'fBrtlip',
      buttons = list(list(extend = 'csv', filename = paste('data-2propM-', Sys.time(), sep=''), text = 'Download')),
      autoWidth = TRUE,
      columnDefs = list(list(className = 'dt-center', targets = '_all'),
                        list(width = '500px', targets = 8)),
      pageLength = 10,
      lengthMenu =  list(c(10, 25, 100, -1), list('10', '25', '100', 'All')) 
    )
  )
  
  # update graph UI
  observeEvent(res2propM(),
               {
                 update_graph_ui(session, res2propM(), "2propM", names2propM)
               })
  
  # create 2propM graph
  output$graph2propM <- renderPlot({
    create_graph(res2propM(), input$x2propM, input$y2propM, input$group2propM,
                 input$lsize2propM, input$psize2propM, input$row2propM, input$col2propM)
  },
  height = reactive({input$height2propM})
  )
  
  
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
      updateTextInput(session, inputId = "r1inc2rate", value = FALSE)
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
      updateTextInput(session, inputId = "r1inc2rate", value = FALSE)
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
      r1inc <- isolate(input$r1inc2rate)
      
      tab <- expand.grid(alpha,
                         power,
                         m,
                         py,
                         r1,
                         r2,
                         cvb,
                         r1inc,
                         stringsAsFactors = FALSE)
      
      if(!is.na(power)){
        validate(
          need(power >= 0 & power <= 1,
               powervalidmsg)
        )
      }
      
      if(!is.na(alpha)){
        validate(
          need(alpha >= 0 & alpha <= 1,
               alphavalidmsg)
        )
      }
      
      # record column index of target parameter
      needind <- which(is.na(tab[1,]))
      # validate that only one input is blank
      validate(
        need(length(needind) == 1,
             "Exactly one of 'alpha', 'power', 'r1', 'r2', 'm', 'py', or 'cvb' must be left blank."
        )
      )
      names(tab) <- names2rate
      target <- names2rate[needind]
      
      # apply function over table of input values
      temp <-pmap_df(tab, crtpwr.2rate.safe)
      
      tab[[target]] <- signif(temp$result, 4)
      tab$error <- map_chr(temp$error, shorten_error, target = target)
      
      # make a column to store target variable for use in graphing
      tab$target <- target
      
      # convert all input values to factors for ggplot
      mutate_if(tab, !(names(tab) %in% c(target,"error","target")), factor)
    })
  
  # create 2rate output table
  output$table2rate <- DT::renderDataTable(
    res2rate()[,c(1:7,9)],
    server = FALSE,
    extensions = 'Buttons',
    filter = 'top',
    options = list(
      # create the button
      dom = 'fBrtlip',
      buttons = list(list(extend = 'csv', filename = paste('data-2rate-', Sys.time(), sep=''), text = 'Download')),
      autoWidth = TRUE,
      columnDefs = list(list(className = 'dt-center', targets = '_all'),
                        list(width = '500px', targets = 8)),
      pageLength = 10,
      lengthMenu =  list(c(10, 25, 100, -1), list('10', '25', '100', 'All')) 
    )
  )

  # update graph UI
  observeEvent(res2rate(),
               {
                 update_graph_ui(session, res2rate(), "2rate", names2rate)
               })
  
  # create 2rate graph
  output$graph2rate <- renderPlot({
    create_graph(res2rate(), input$x2rate, input$y2rate, input$group2rate,
                 input$lsize2rate, input$psize2rate, input$row2rate, input$col2rate)
  },
  height = reactive({input$height2rate})
  )
  
}


shinyApp(ui = ui, server = server)