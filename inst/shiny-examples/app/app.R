library(shiny)
library(shinyBS)
library(DT)
library(tidyverse)
library(stringr)
library(clusterPower)

source("labels.R")

# vectors of names for graphin purposes
names2mean <- c("alpha","power","m","n","cv","d","icc","varw","method")
names2prop <- c("alpha","power","m","n","cv","p1","p2","icc")
names2rate <- c("alpha","power","m","py","r1","r2","cvb")

# "safe" versions of functions to catch errors
crtpwr.2mean.safe <- function(alpha,power,m,n,cv,d,icc,varw,method){
  # make safe version
  fun <- safely(crtpwr.2mean, otherwise = NA)
  # store result
  res <- fun(alpha,power,m,n,cv,d,icc,varw,method)
  # if res$error NULL, set to NA, otherwise set to message
  if(is.null(res$error)){
    res$error = NA
  } else {
    res$error <- res$error$message
  }
  res
}

crtpwr.2prop.safe <- function(alpha,power,m,n,cv,p1,p2,icc){
  # make safe version
  fun <- safely(crtpwr.2prop, otherwise = NA)
  # store result
  res <- fun(alpha,power,m,n,cv,p1,p2,icc)
  # if res$error NULL, set to NA, otherwise set to message
  if(is.null(res$error)){
    res$error = NA
  } else {
    res$error <- res$error$message
  }
  res
}

crtpwr.2rate.safe <- function(alpha,power,m,py,r1,r2,cvb){
  # make safe version
  fun <- safely(crtpwr.2rate, otherwise = NA)
  # store result
  res <- fun(alpha,power,m,py,r1,r2,cvb)
  # if res$error NULL, set to NA, otherwise set to message
  if(is.null(res$error)){
    res$error = NA
  } else {
    res$error <- res$error$message
  }
  res
}

# function to shorten error messages
shorten_error <- function(err, target){
  ifelse(err == "did not succeed extending the interval endpoints for f(lower) * f(upper) <= 0",
         paste0("No '", target, "' found for given inputs."),
         err)
}
  
umass <- "font-family: 'Open Sans', Helvetica, Arial, sans-serif; font-weight: bold; color: #ffffff; background-color: #881c1c; border: 3px solid #000000;"

# function to convert app input into numeric vectors
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
        <p>To use the two-arm calculator, leave the desired quantity blank and enter values for the other quantities.</p>
        <p>You may specify more than one input quantity by separating numbers with spaces or commas.<p>
        <p>You may specify a sequence of values by typing 'from X to Y by Z', where 'X' is the starting value, 'Y' is the ending value, and 'Z' is the increment.</p>"),
  HTML("This Beta has minimal documentation; please contact ken.kleinman@gmail.com with any feedback."),
  tabsetPanel(
    #-----------------------------------------------------------------------------------------------------------
    tabPanel("Continuous",
             column(2,
                    #----------------------------------------------------------
                    fluidRow(textInput("alpha2mean", HTML(alphatext),
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
                      
                    ),
                    conditionalPanel(condition = "output.table2mean != null",
                                     fluidRow(downloadButton("dl2mean", dltext)))
             ),
             column(10,
                    tabsetPanel(
                      tabPanel("Data",
                               DT::dataTableOutput("table2mean")
                      ),
                      tabPanel("Graphs",
                               column(2,
                                      fluidRow(selectInput("y2mean", ylab,
                                                           choices = c(None = ".", names2mean), selected = ".")),
                                      fluidRow(selectInput("x2mean", xlab,
                                                           choices = c(None = ".", names2mean), selected = ".")),
                                      fluidRow(selectInput("group2mean", grouplab,
                                                           choices = c(None = ".", names2mean), selected = ".")),
                                      fluidRow(checkboxInput("color2mean", colorlab,value = TRUE)),
                                      fluidRow(selectInput("row2mean", rowlab,
                                                           choices = c(None = ".", names2mean))),
                                      fluidRow(selectInput("col2mean", collab,
                                                           choices = c(None = ".", names2mean))),
                                      fluidRow(numericInput("height2mean", heightlab, value = 400,
                                                            min = 100, max = 2000, step = 10)),
                                      fluidRow(numericInput("psize2mean", psizelab, value = 3,
                                                            min = 0.5, max = 4, step = 0.25)),
                                      fluidRow(numericInput("lsize2mean", lsizelab, value = 1,
                                                            min = 0.5, max = 2, step = 0.25))
                               ),
                               column(10,
                                      plotOutput("graph2mean", height = "auto")
                               )
                      ) # end tabPanel("Graphs"...
                    ) # end tabsetPanel(...
             ) # end column(10,...
    ), # end tabPanel("Continuous ...
    #-----------------------------------------------------------------------------------------------------------
    tabPanel("Binary",
             column(2,
                    #----------------------------------------------------------
                    fluidRow(textInput("alpha2prop", HTML(alphatext),
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
                      
                    ),
                    conditionalPanel(condition = "output.table2prop != null",
                                     fluidRow(downloadButton("dl2prop", dltext)))
             ), # end column(2, ..
             column(10,
                    tabsetPanel(
                      tabPanel("Data",
                               DT::dataTableOutput("table2prop")
                      ),
                      tabPanel("Graphs",
                               column(2,
                                      fluidRow(selectInput("y2prop", ylab,
                                                           choices = names2prop, selected = "power")),
                                      fluidRow(selectInput("x2prop", xlab,
                                                           choices = names2prop, selected = "m")),
                                      fluidRow(selectInput("group2prop", grouplab,
                                                           choices = c(None = ".", names2prop), selected = "n")),
                                      fluidRow(checkboxInput("color2prop", colorlab,value = TRUE)),
                                      fluidRow(selectInput("row2prop", rowlab,
                                                           choices = c(None = ".", names2prop))),
                                      fluidRow(selectInput("col2prop", collab,
                                                           choices = c(None = ".", names2prop))),
                                      fluidRow(numericInput("height2prop", heightlab, value = 400,
                                                            min = 100, max = 2000, step = 10)),
                                      fluidRow(numericInput("psize2prop", psizelab, value = 3,
                                                            min = 0.5, max = 4, step = 0.25)),
                                      fluidRow(numericInput("lsize2prop", lsizelab, value = 1,
                                                            min = 0.5, max = 2, step = 0.25))
                               ),
                               column(10,
                                      plotOutput("graph2prop", height = "auto")
                               )
                      ) # end tabPanel("Graphs"...
                    ) # end tabsetPanel(...
             ) # end column(10,...
    ), # end tabPanel("Binary ...
    #-----------------------------------------------------------------------------------------------------------
    tabPanel("Count",
             column(2,
                    #----------------------------------------------------------
                    fluidRow(textInput("alpha2rate", HTML(alphatext),
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
                    ),
                    conditionalPanel(condition = "output.table2rate != null",
                                     fluidRow(downloadButton("dl2rate", dltext)))
             ),
             column(10,
                    tabsetPanel(
                      tabPanel("Data",
                               DT::dataTableOutput("table2rate")
                      ),
                      tabPanel("Graphs",
                               column(2,
                                      fluidRow(selectInput("y2rate", ylab,
                                                           choices = names2rate, selected = "power")),
                                      fluidRow(selectInput("x2rate", xlab,
                                                           choices = names2rate, selected = "m")),
                                      fluidRow(selectInput("group2rate", grouplab,
                                                           choices = names2rate, selected = "py")),
                                      fluidRow(checkboxInput("color2rate", colorlab,value = TRUE)),
                                      fluidRow(selectInput("row2rate", rowlab,
                                                           choices = c(None = ".", names2rate))),
                                      fluidRow(selectInput("col2rate", collab,
                                                           choices = c(None = ".", names2rate))),
                                      fluidRow(numericInput("height2rate", heightlab, value = 400,
                                                            min = 100, max = 2000, step = 10)),
                                      fluidRow(numericInput("psize2rate", psizelab, value = 3,
                                                            min = 0.5, max = 4, step = 0.25)),
                                      fluidRow(numericInput("lsize2rate", lsizelab, value = 1,
                                                            min = 0.5, max = 2, step = 0.25))
                               ),
                               column(10,
                                      plotOutput("graph2rate", height = "auto")
                               )
                      ) # end tabPanel("Graphs"...
                    ) # end tabsetPanel(...
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
      updateNumericInput(session, inputId = "alpha2mean", value = "0.05")
      updateNumericInput(session, inputId = "power2mean", value = "")
      updateNumericInput(session, inputId = "cv2mean", value = "0")
      updateCheckboxGroupInput(session, inputId = "method2mean", selected = "taylor")
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
      
      validate(
        need(power >= 0 & power <= 1,
             powervalidmsg)
      )
      
      validate(
        need(alpha >= 0 & alpha <= 1,
             alphavalidmsg)
      )
      
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
    filter = 'top',
    options = list(
      autoWidth = TRUE,
      columnDefs = list(list(className = 'dt-center', targets = '_all'),
                        list(width = '500px', targets = 10)),
      pageLength = 10
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
  
  # update graph UI
  observeEvent(res2mean(),
               {
                 # store target
                 target <- unique(res2mean()$target)
                 # update y-axis to target variable
                 updateSelectInput(session, "y2mean", label = "Y",
                                   choices = c(None = ".", names2mean), selected = target)
                 # if the target is m, set default x-axis to n, otherwise set default axis to m
                 if(target == "m"){
                   updateSelectInput(session, "x2mean", label = "X",
                                     choices = c(None = ".", names2mean), selected = "n")
                   updateSelectInput(session, "group2mean", label = "X",
                                     choices = c(None = ".", names2mean), selected = "n")
                 } else {
                   updateSelectInput(session, "x2mean", label = "X",
                                     choices = c(None = ".", names2mean), selected = "m")
                   updateSelectInput(session, "group2mean", label = "X",
                                     choices = c(None = ".", names2mean), selected = "n")
                 }
               })
  
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
      
      validate(
        need(power >= 0 & power <= 1,
             powervalidmsg)
      )
      
      validate(
        need(alpha >= 0 & alpha <= 1,
             alphavalidmsg)
      )
      
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
    res2prop()[,1:9],
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
  
  # update graph UI
  observeEvent(res2prop(),
               {
                 # store target
                 target <- unique(res2prop()$target)
                 # update y-axis to target variable
                 updateSelectInput(session, "y2prop", label = "Y",
                                   choices = c(None = ".", names2prop), selected = target)
                 # if the target is m, set default x-axis to n, otherwise set default axis to m
                 if(target == "m"){
                   updateSelectInput(session, "x2prop", label = "X",
                                     choices = c(None = ".", names2prop), selected = "n")
                   updateSelectInput(session, "group2prop", label = "X",
                                     choices = c(None = ".", names2prop), selected = "n")
                 } else {
                   updateSelectInput(session, "x2prop", label = "X",
                                     choices = c(None = ".", names2prop), selected = "m")
                   updateSelectInput(session, "group2prop", label = "X",
                                     choices = c(None = ".", names2prop), selected = "n")
                 }
               })
  
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
      
      validate(
        need(power >= 0 & power <= 1,
             powervalidmsg)
      )
      
      validate(
        need(alpha >= 0 & alpha <= 1,
             alphavalidmsg)
      )
      
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
    res2rate()[,1:8],
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
  
  # update graph UI
  observeEvent(res2rate(),
               {
                 # store target
                 target <- unique(res2rate()$target)
                 # update y-axis to target variable
                 updateSelectInput(session, "y2rate", label = "Y",
                                   choices = c(None = ".", names2rate), selected = target)
                 # if the target is m, set default x-axis to n, otherwise set default axis to m
                 if(target == "m"){
                   updateSelectInput(session, "x2rate", label = "X",
                                     choices = c(None = ".", names2rate), selected = "n")
                   updateSelectInput(session, "group2rate", label = "X",
                                     choices = c(None = ".", names2rate), selected = "n")
                 } else {
                   updateSelectInput(session, "x2rate", label = "X",
                                     choices = c(None = ".", names2rate), selected = "m")
                   updateSelectInput(session, "group2rate", label = "X",
                                     choices = c(None = ".", names2rate), selected = "n")
                 }
               })
  
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