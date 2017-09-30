# Helper functions for the app

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


# function to generate graph output
make_table_and_graph <- function(outcome, paramnames){
  tabsetPanel(
    tabPanel("Data",
             DT::dataTableOutput(paste0("table2",outcome))
    ),
    tabPanel("Graphs",
             column(2,
                    fluidRow(selectInput(paste0("y2",outcome), ylab,
                                         choices = c(None = ".", paramnames), selected = ".")),
                    fluidRow(selectInput(paste0("x2",outcome), xlab,
                                         choices = c(None = ".", paramnames), selected = ".")),
                    fluidRow(selectInput(paste0("group2",outcome), grouplab,
                                         choices = c(None = ".", paramnames), selected = ".")),
                    fluidRow(checkboxInput(paste0("color2",outcome), colorlab,value = TRUE)),
                    fluidRow(selectInput(paste0("row2",outcome), rowlab,
                                         choices = c(None = ".", paramnames))),
                    fluidRow(selectInput(paste0("col2",outcome), collab,
                                         choices = c(None = ".", paramnames))),
                    fluidRow(numericInput(paste0("height2",outcome), heightlab, value = 400,
                                          min = 100, max = 2000, step = 10)),
                    fluidRow(numericInput(paste0("psize2",outcome), psizelab, value = 3,
                                          min = 0.5, max = 4, step = 0.25)),
                    fluidRow(numericInput(paste0("lsize2",outcome), lsizelab, value = 1,
                                          min = 0.5, max = 2, step = 0.25))
             ),
             column(10,
                    plotOutput(paste0("graph2",outcome), height = "auto")
             )
    ) # end tabPanel("Graphs"...
  ) # end tabsetPanel(...
}