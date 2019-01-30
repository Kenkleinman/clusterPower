library(shiny)
library(shinyBS)
library(shinycssloaders)
library(DT)
library(tidyverse)
library(stringr)
library(clusterPower)

source("labels.R")
source("helpers.R")
source("builders.R")

# vectors of names, needed for graphs and target param selection
names2mean <- c("alpha","power","nclusters","nsubjects","cv","d","icc","vart","method")
names2meanD <- c("alpha","power","nclusters","nsubjects","d","icc","rho_c","rho_s","vart")
names2meanM <- c("alpha","power","nclusters","nsubjects","d","icc","vart","rho_m")
namesnmean <- c("alpha","power","narms","nclusters","nsubjects","vara","varc","vare")
names2prop <- c("alpha","power","nclusters","nsubjects","cv","p1","p2","icc","pooled","p1inc")
names2propD <- c("alpha","power","nclusters","nsubjects","p","d","icc","rho_c","rho_s")
names2propM <- c("alpha","power","nclusters","nsubjects","p1","p2","cvm","p1inc")
names2rate <- c("alpha","power","nclusters","py","r1","r2","cvb","r1inc")

umass <- "font-family: 'Open Sans', Helvetica, Arial, sans-serif; font-weight: bold; color: #ffffff; background-color: #881c1c; border: 3px solid #000000;"

ui <- function(request){
  fluidPage(
  HTML("<h3>Simple Two-Arm Designs</h3>
        <p>To use the two-arm calculator, leave the desired quantity blank and enter values for the other quantities.</p>
        <p>You may specify more than one input quantity by separating numbers with spaces or commas.<p>
        <p>You may specify a sequence of values by typing 'X to Y by Z', where 'X' is the starting value, 'Y' is the ending value, and 'Z' is the increment.</p>"),
  HTML("This Beta has minimal documentation; please contact ken.kleinman@gmail.com with any feedback."),
  column(12, bookmarkButton("Save App State")),
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
                    fluidRow(textInput("nclusters2mean", nclusterstext,
                                       value = "", width = "100%")),
                    bsTooltip("nclusters2mean", nclusterstooltip,
                              'right', options = list(container = "body")),
                    #----------------------------------------------------------
                    fluidRow(textInput("n2mean", nsubjectstext,
                                       value = "", width = "100%")),
                    bsTooltip("n2mean", nsubjectstooltip,
                              'right', options = list(container = "body")),
                    #----------------------------------------------------------
                    fluidRow(textInput("icc2mean", icctext,
                                       value = "", width = "100%")),
                    bsTooltip("icc2mean", icctooltip,
                              'right', options = list(container = "body")),
                    #----------------------------------------------------------
                    fluidRow(textInput("vart2mean", varttext,
                                       value = "", width = "100%")),
                    bsTooltip("vart2mean", varttooltip,
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
                    fluidRow(textInput("nclusters2meanD", nclusterstext,
                                       value = "", width = "100%")),
                    bsTooltip("nclusters2meanD", nclusterstooltip,
                              'right', options = list(container = "body")),
                    #----------------------------------------------------------
                    fluidRow(textInput("n2meanD", nsubjectstext,
                                       value = "", width = "100%")),
                    bsTooltip("n2meanD", nsubjectstooltip,
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
                    fluidRow(textInput("vart2meanD", varttext,
                                       value = "", width = "100%")),
                    bsTooltip("vart2meanD", varttooltip,
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
                    fluidRow(textInput("nclusters2meanM", nclusterstext,
                                       value = "", width = "100%")),
                    bsTooltip("nclusters2meanM", nclusterstooltip,
                              'right', options = list(container = "body")),
                    #----------------------------------------------------------
                    fluidRow(textInput("n2meanM", nsubjectstext,
                                       value = "", width = "100%")),
                    bsTooltip("n2meanM", nsubjectstooltip,
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
                    fluidRow(textInput("vart2meanM", varttext,
                                       value = "", width = "100%")),
                    bsTooltip("vart2meanM", varttooltip,
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
    tabPanel("Continuous Multiarm",
             column(2,
                    #----------------------------------------------------------
                    fluidRow(textInput("alphanmean", alphatext,
                                       value = "0.05", width = "100%")),
                    bsTooltip("alphanmean", alphatooltip,
                              'right', options = list(container = "body")),
                    #----------------------------------------------------------
                    fluidRow(textInput("powernmean", powertext,
                                       value = "", width = "100%")),
                    bsTooltip("powernmean", powertooltip,
                              'right', options = list(container = "body")),
                    #----------------------------------------------------------
                    fluidRow(textInput("narmsnmean", narmstext,
                                       value = "", width = "100%")),
                    bsTooltip("narmsnmean", narmstooltip,
                              'right', options = list(container = "body")),
                    #----------------------------------------------------------
                    fluidRow(textInput("nclustersnmean", nclusterstext,
                                       value = "", width = "100%")),
                    bsTooltip("nclustersnmean", nclusterstooltip,
                              'right', options = list(container = "body")),
                    #----------------------------------------------------------
                    fluidRow(textInput("nnmean", nsubjectstext,
                                       value = "", width = "100%")),
                    bsTooltip("nnmean", nsubjectstooltip,
                              'right', options = list(container = "body")),
                    #----------------------------------------------------------
                    fluidRow(textInput("varanmean", varatext,
                                       value = "", width = "100%")),
                    bsTooltip("varanmean", varatooltip,
                              'right', options = list(container = "body")),
                    #----------------------------------------------------------
                    fluidRow(textInput("varcnmean", varctext,
                                       value = "", width = "100%")),
                    bsTooltip("varcnmean", varctooltip,
                              'right', options = list(container = "body")),
                    #----------------------------------------------------------
                    fluidRow(textInput("varenmean", varetext,
                                       value = "", width = "100%")),
                    bsTooltip("varenmean", varetooltip,
                              'right', options = list(container = "body")),
                    #----------------------------------------------------------
                    fluidRow(
                      column(6, style='padding:0px;', actionButton("defaultnmean", defaulttext, width = "100%")),
                      column(6, style='padding:0px;', actionButton("clearnmean", clearalltext, width = "100%"))
                    ),
                    fluidRow(
                      column(12, style='padding: 0px;', actionButton("calcnmean", calctext, width = "100%",
                                                                     style = umass)),
                      fluidRow(column(12, credittext))
                      
                    )
             ),
             column(10,
                    make_table_and_graph("nmean", namesnmean)
             ) # end column(10,...
    ), # end tabPanel("Continuous ...
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
                    fluidRow(textInput("nclusters2prop", nclusterstext,
                                       value = "", width = "100%")),
                    bsTooltip("nclusters2prop", nclusterstooltip,
                              'right', options = list(container = "body")),
                    #----------------------------------------------------------
                    fluidRow(textInput("n2prop", nsubjectstext,
                                       value = "", width = "100%")),
                    bsTooltip("n2prop", nsubjectstooltip,
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
                    fluidRow(textInput("nclusters2propD", nclusterstext,
                                       value = "", width = "100%")),
                    bsTooltip("nclusters2propD", nclusterstooltip,
                              'right', options = list(container = "body")),
                    #----------------------------------------------------------
                    fluidRow(textInput("n2propD", nsubjectstext,
                                       value = "", width = "100%")),
                    bsTooltip("n2propD", nsubjectstooltip,
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
                    fluidRow(textInput("nclusters2propM", nclusterstext,
                                       value = "", width = "100%")),
                    bsTooltip("nclusters2propM", nclusterstooltip,
                              'right', options = list(container = "body")),
                    #----------------------------------------------------------
                    fluidRow(textInput("n2propM", nsubjectstext,
                                       value = "", width = "100%")),
                    bsTooltip("n2propM", nsubjectstooltip,
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
                    fluidRow(textInput("nclusters2rate", nclusterstext,
                                       value = "", width = "100%")),
                    bsTooltip("nclusters2rate", nclusterstooltip,
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
    ) # end tabsetPanel
)} # end fluidPage

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
      updateTextInput(session, inputId = "nclusters2mean", value = "")
      updateTextInput(session, inputId = "n2mean", value = "")
      updateTextInput(session, inputId = "d2mean", value = "")
      updateTextInput(session, inputId = "icc2mean", value = "")
      updateTextInput(session, inputId = "vart2mean", value = "")
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
      updateTextInput(session, inputId = "nclusters2mean", value = "")
      updateTextInput(session, inputId = "n2mean", value = "")
      updateTextInput(session, inputId = "d2mean", value = "")
      updateTextInput(session, inputId = "icc2mean", value = "")
      updateTextInput(session, inputId = "vart2mean", value = "")
    }
  ) # end observeEvent(input$clear2mean ...
  
  # create 2mean data
  res2mean <- eventReactive(
    input$calc2mean,
    {
      # convert inputs to numeric vectors
      alpha <- make_sequence(isolate(input$alpha2mean))
      power <- make_sequence(isolate(input$power2mean))
      nclusters <- make_sequence(isolate(input$nclusters2mean))
      nsubjects <- make_sequence(isolate(input$n2mean))
      cv <- make_sequence(isolate(input$cv2mean))
      d <- make_sequence(isolate(input$d2mean))
      icc <- make_sequence(isolate(input$icc2mean))
      vart <- make_sequence(isolate(input$vart2mean))
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
                         nclusters,
                         nsubjects,
                         cv,
                         d,
                         icc,
                         vart,
                         method,
                         stringsAsFactors = FALSE)
      
      # record the column index of the target parameter
      needind <- which(is.na(tab[1,]))
      # validate that only one input is blank
      validate(
        need(length(needind) == 1,
             "Exactly one of 'alpha', 'power', 'd', 'nclusters', 'nsubjects', 'icc', 'vart', and 'cv' must be left blank."
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
      #autoWidth = TRUE,
      columnDefs = list(list(className = 'dt-center', targets = '_all'),
                        list(width = '700px', targets = 10)),
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
  width = reactive({input$width2mean}),
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
      updateTextInput(session, inputId = "nclusters2meanD", value = "")
      updateTextInput(session, inputId = "n2meanD", value = "")
      updateTextInput(session, inputId = "d2meanD", value = "")
      updateTextInput(session, inputId = "icc2meanD", value = "")
      updateTextInput(session, inputId = "vart2meanD", value = "")
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
      updateTextInput(session, inputId = "nclusters2meanD", value = "")
      updateTextInput(session, inputId = "n2meanD", value = "")
      updateTextInput(session, inputId = "d2meanD", value = "")
      updateTextInput(session, inputId = "icc2meanD", value = "")
      updateTextInput(session, inputId = "vart2meanD", value = "")
    }
  ) # end observeEvent(input$clear2meanD ...
  
  # create 2meanD data
  res2meanD <- eventReactive(
    input$calc2meanD,
    {
      # convert inputs to numeric vectors
      alpha <- make_sequence(isolate(input$alpha2meanD))
      power <- make_sequence(isolate(input$power2meanD))
      nclusters <- make_sequence(isolate(input$nclusters2meanD))
      nsubjects <- make_sequence(isolate(input$n2meanD))
      rho_c <- make_sequence(isolate(input$rho_c2meanD))
      rho_s <- make_sequence(isolate(input$rho_s2meanD))
      d <- make_sequence(isolate(input$d2meanD))
      icc <- make_sequence(isolate(input$icc2meanD))
      vart <- make_sequence(isolate(input$vart2meanD))
      
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
                         nclusters,
                         nsubjects,
                         d,
                         icc,
                         rho_c,
                         rho_s,
                         vart,
                         stringsAsFactors = FALSE)
      
      # record the column index of the target parameter
      needind <- which(is.na(tab[1,]))
      # validate that only one input is blank
      validate(
        need(length(needind) == 1,
             "Exactly one of 'alpha', 'power', 'd', 'nclusters', 'nsubjects', 'icc', 'rho_c', 'rho_s', and 'vart' must be left blank."
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
      #autoWidth = TRUE,
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
  width = reactive({input$width2meanD}),
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
      updateTextInput(session, inputId = "nclusters2meanM", value = "")
      updateTextInput(session, inputId = "n2meanM", value = "")
      updateTextInput(session, inputId = "d2meanM", value = "")
      updateTextInput(session, inputId = "icc2meanM", value = "")
      updateTextInput(session, inputId = "vart2meanM", value = "")
    }
  ) # end observeEvent(input$default2meanM ...
  
  
  # clear 2meanM inputs 
  observeEvent(
    input$clear2meanM,
    {
      updateTextInput(session, inputId = "alpha2meanM", value = "")
      updateTextInput(session, inputId = "power2meanM", value = "")
      updateTextInput(session, inputId = "rho_m2meanM", value = "")
      updateTextInput(session, inputId = "nclusters2meanM", value = "")
      updateTextInput(session, inputId = "n2meanM", value = "")
      updateTextInput(session, inputId = "d2meanM", value = "")
      updateTextInput(session, inputId = "icc2meanM", value = "")
      updateTextInput(session, inputId = "vart2meanM", value = "")
    }
  ) # end observeEvent(input$clear2meanM ...
  
  # create 2meanM data
  res2meanM <- eventReactive(
    input$calc2meanM,
    {
      # convert inputs to numeric vectors
      alpha <- make_sequence(isolate(input$alpha2meanM))
      power <- make_sequence(isolate(input$power2meanM))
      nclusters <- make_sequence(isolate(input$nclusters2meanM))
      nsubjects <- make_sequence(isolate(input$n2meanM))
      rho_m <- make_sequence(isolate(input$rho_m2meanM))
      d <- make_sequence(isolate(input$d2meanM))
      icc <- make_sequence(isolate(input$icc2meanM))
      vart <- make_sequence(isolate(input$vart2meanM))
      
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
                         nclusters,
                         nsubjects,
                         d,
                         icc,
                         vart,
                         rho_m,
                         stringsAsFactors = FALSE)
      
      # record the column index of the target parameter
      needind <- which(is.na(tab[1,]))
      # validate that only one input is blank
      validate(
        need(length(needind) == 1,
             "Exactly one of 'alpha', 'power', 'd', 'nclusters', 'nsubjects', 'icc', 'rho_m', and 'vart' must be left blank."
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
      #autoWidth = TRUE,
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
  width = reactive({input$width2meanM}),
  height = reactive({input$height2meanM})
  )
  
  #----------------------------------------------------------------------------
  # Multi-Arm (nmean)
  #----------------------------------------------------------------------------
  
  # reset nmean inputs to default values
  observeEvent(
    input$defaultnmean,
    {
      updateTextInput(session, inputId = "alphanmean", value = "0.05")
      updateTextInput(session, inputId = "powernmean", value = "")
      updateTextInput(session, inputId = "narmnmean", value = "")
      updateTextInput(session, inputId = "nclustersnmean", value = "")
      updateTextInput(session, inputId = "nnmean", value = "")
      updateTextInput(session, inputId = "varanmean", value = "")
      updateTextInput(session, inputId = "varcnmean", value = "")
      updateTextInput(session, inputId = "varenmean", value = "")
    }
  ) # end observeEvent(input$defaultnmean ...
  
  
  # clear nmean inputs 
  observeEvent(
    input$clearnmean,
    {
      updateTextInput(session, inputId = "alphanmean", value = "")
      updateTextInput(session, inputId = "powernmean", value = "")
      updateTextInput(session, inputId = "narmnmean", value = "")
      updateTextInput(session, inputId = "nclustersnmean", value = "")
      updateTextInput(session, inputId = "nnmean", value = "")
      updateTextInput(session, inputId = "varanmean", value = "")
      updateTextInput(session, inputId = "varcnmean", value = "")
      updateTextInput(session, inputId = "varenmean", value = "")
    }
  ) # end observeEvent(input$clearnmean ...
  
  # create nmean data
  resnmean <- eventReactive(
    input$calcnmean,
    {
      # convert inputs to numeric vectors
      alpha <- make_sequence(isolate(input$alphanmean))
      power <- make_sequence(isolate(input$powernmean))
      narms <- make_sequence(isolate(input$narmsnmean))
      nclusters <- make_sequence(isolate(input$nclustersnmean))
      nsubjects <- make_sequence(isolate(input$nnmean))
      vara <- make_sequence(isolate(input$varanmean))
      varc <- make_sequence(isolate(input$varcnmean))
      vare <- make_sequence(isolate(input$varenmean))
      
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
                         narms,
                         nclusters,
                         nsubjects,
                         vara,
                         varc,
                         vare,
                         stringsAsFactors = FALSE)
      
      # record the column index of the target parameter
      needind <- which(is.na(tab[1,]))
      # validate that only one input is blank
      validate(
        need(length(needind) == 1,
             "Exactly one of 'alpha', 'power', 'narms', 'nclusters', 'nsubjects', 'vara', 'varc', and 'vare' must be left blank."
        )
      )
      names(tab) <- namesnmean
      target <- namesnmean[needind]
      
      # apply function over table of input values
      temp <-pmap_df(tab, crtpwr.nmean.safe)
      
      tab[[target]] <- signif(temp$result, 4)
      tab$error <- map_chr(temp$error, shorten_error, target = target)
      
      # make a column to store target variable for use in graphing
      tab$target <- target
      
      # convert all input values to factors for ggplot
      mutate_if(tab, !(names(tab) %in% c(target,"error","target")), factor)
    })
  
  # create nmean output table
  output$tablenmean <- DT::renderDataTable(
    resnmean()[,1:10],
    server = FALSE,
    extensions = 'Buttons',
    filter = 'top',
    options = list(
      # create the button
      dom = 'fBrtlip',
      buttons = list(list(extend = 'csv', filename = paste('data-nmean-', Sys.time(), sep=''), text = 'Download')),
      #autoWidth = TRUE,
      columnDefs = list(list(className = 'dt-center', targets = '_all'),
                        list(width = '700px', targets = 10)),
      lengthMenu = list(c(10, 25, 100, -1), list('10','25','100','All')),
      pageLength = 10
    )
  )
  
  # update graph UI
  observeEvent(resnmean(),
               {
                 update_graph_ui(session, resnmean(), "nmean", namesnmean)
               })
  
  # create nmean graph
  output$graphnmean <- renderPlot({
    create_graph(resnmean(), input$xnmean, input$ynmean, input$groupnmean,
                 input$lsizenmean, input$psizenmean, input$rownmean, input$colnmean)
  },
  width = reactive({input$widthnmean}),
  height = reactive({input$heightnmean})
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
      updateTextInput(session, inputId = "nclusters2prop", value = "")
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
      updateTextInput(session, inputId = "nclusters2prop", value = "")
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
      nclusters <- make_sequence(isolate(input$nclusters2prop))
      nsubjects <- make_sequence(isolate(input$n2prop))
      cv <- make_sequence(isolate(input$cv2prop))
      p1 <- make_sequence(isolate(input$p12prop))
      p2 <- make_sequence(isolate(input$p22prop))
      icc <- make_sequence(isolate(input$icc2prop))
      pooled <- isolate(input$pooled2prop)
      p1inc <- isolate(input$p1inc2prop)
      
      tab <- expand.grid(alpha,
                         power,
                         nclusters,
                         nsubjects,
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
             "Exactly one of 'alpha', 'power', 'p1', 'p2', 'nclusters', 'nsubjects', 'icc' or 'cv' must be left blank."
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
      #autoWidth = TRUE,
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
  width = reactive({input$width2prop}),
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
      updateTextInput(session, inputId = "nclusters2propD", value = "")
      updateTextInput(session, inputId = "n2propD", value = "")
      updateTextInput(session, inputId = "p2propD", value = "")
      updateTextInput(session, inputId = "d2propD", value = "")
      updateTextInput(session, inputId = "icc2propD", value = "")
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
      updateTextInput(session, inputId = "nclusters2propD", value = "")
      updateTextInput(session, inputId = "n2propD", value = "")
      updateTextInput(session, inputId = "p2propD", value = "")
      updateTextInput(session, inputId = "d2propD", value = "")
      updateTextInput(session, inputId = "icc2propD", value = "")
    }
  ) # end observeEvent(input$clear2propD ...
  
  # create 2propD data
  res2propD <- eventReactive(
    input$calc2propD,
    {
      # convert inputs to numeric vectors
      alpha <- make_sequence(isolate(input$alpha2propD))
      power <- make_sequence(isolate(input$power2propD))
      nclusters <- make_sequence(isolate(input$nclusters2propD))
      nsubjects <- make_sequence(isolate(input$n2propD))
      rho_c <- make_sequence(isolate(input$rho_c2propD))
      rho_s <- make_sequence(isolate(input$rho_s2propD))
      p <- make_sequence(isolate(input$p2propD))
      d <- make_sequence(isolate(input$d2propD))
      icc <- make_sequence(isolate(input$icc2propD))

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
                         nclusters,
                         nsubjects,
                         p,
                         d,
                         icc,
                         rho_c,
                         rho_s,
                         stringsAsFactors = FALSE)
      
      # record the column index of the target parameter
      needind <- which(is.na(tab[1,]))
      # validate that only one input is blank
      validate(
        need(length(needind) == 1,
             "Exactly one of 'alpha', 'power', 'p', 'd', 'nclusters', 'nsubjects', 'icc', 'rho_c', and 'rho_s' must be left blank."
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
    res2propD()[, c(1:10)],
    server = FALSE,
    extensions = 'Buttons',
    filter = 'top',
    options = list(
      # create the button
      dom = 'fBrtlip',
      buttons = list(list(extend = 'csv', filename = paste('data-2propD-', Sys.time(), sep=''), text = 'Download')),
      #autoWidth = TRUE,
      columnDefs = list(list(className = 'dt-center', targets = '_all'),
                        list(width = '500px', targets = 10)),
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
  width = reactive({input$width2propD}),
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
      updateTextInput(session, inputId = "nclusters2propM", value = "")
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
      updateTextInput(session, inputId = "nclusters2propM", value = "")
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
      nclusters <- make_sequence(isolate(input$nclusters2propM))
      nsubjects <- make_sequence(isolate(input$n2propM))
      p1 <- make_sequence(isolate(input$p12propM))
      p2 <- make_sequence(isolate(input$p22propM))
      cvm <- make_sequence(isolate(input$cvm2propM))
      p1inc <- isolate(input$p1inc2propM)
      
      tab <- expand.grid(alpha,
                         power,
                         nclusters,
                         nsubjects,
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
             "Exactly one of 'alpha', 'power', 'p1', 'p2', 'nclusters', 'nsubjects', or 'cvm' must be left blank."
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
      #autoWidth = TRUE,
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
  width = reactive({input$width2propM}),
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
      updateTextInput(session, inputId = "nclusters2rate", value = "")
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
      updateTextInput(session, inputId = "nclusters2rate", value = "")
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
      nclusters <- make_sequence(isolate(input$nclusters2rate))
      py <- make_sequence(isolate(input$py2rate))
      r1 <- make_sequence(isolate(input$r12rate))
      r2 <- make_sequence(isolate(input$r22rate))
      cvb <- make_sequence(isolate(input$cvb2rate))
      r1inc <- isolate(input$r1inc2rate)
      
      tab <- expand.grid(alpha,
                         power,
                         nclusters,
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
             "Exactly one of 'alpha', 'power', 'r1', 'r2', 'nclusters', 'py', or 'cvb' must be left blank."
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
      #autoWidth = TRUE,
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
  width = reactive({input$width2rate}),
  height = reactive({input$height2rate})
  )
  
}


shinyApp(ui = ui, server = server, enableBookmarking = "server")