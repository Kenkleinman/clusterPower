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
                    #fluidRow(checkboxInput(paste0("color2",outcome), colorlab,value = TRUE)),
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


# function to update graph UI
update_graph_ui <- function(sessionname, dataset, outcome, paramnames){
  # update graph UI
  observeEvent(dataset,
               {
                 # store target
                 target <- unique(dataset$target)
                 # update y-axis to target variable
                 updateSelectInput(sessionname, paste0("y2",outcome), label = ylab,
                                   choices = c(None = ".", paramnames), selected = target)
                 # get the levels of quantities in data
                 levs <- map(dataset,levels)
                 # sort levels 
                 sortlevs <- sort(map_int(levs,length), decreasing = TRUE)
                 # put the quantity with most levels, levels > 1, on x
                 if(sortlevs[1] > 1){
                   updateSelectInput(sessionname, paste0("x2",outcome), label = xlab,
                                     choices = c(None = ".", paramnames), selected = names(sortlevs[1]))
                 }
                 # put the quantity with the second most levels, levels > 1, on group
                 if(sortlevs[2] > 1){
                   updateSelectInput(sessionname, paste0("group2",outcome), label = grouplab,
                                     choices = c(None = ".", paramnames), selected = names(sortlevs[2]))
                 }
                 # put the quantity with the third most levels, levels > 1, on row
                 if(sortlevs[3] > 1){
                   updateSelectInput(sessionname, paste0("row2",outcome), label = rowlab,
                                     choices = c(None = ".", paramnames), selected = names(sortlevs[3]))
                 }
                 # put the quantity with the fourth most levels, levels > 1, on column
                 if(sortlevs[4] > 1){
                   updateSelectInput(sessionname, paste0("col2",outcome), label = collab,
                                     choices = c(None = ".", paramnames), selected = names(sortlevs[4]))
                 }
                 # if no quantity has more than one levels, just target on x-axis
                 if(max(sortlevs) == 1){
                   updateSelectInput(sessionname, paste0("x2",outcome), label = xlab,
                                     choices = c(None = ".", paramnames), selected = target)
                 }
               })
}

# function to create graph
create_graph <- function(dataset, xvar, yvar, groupvar, lsizevar, psizevar, rowvar, colvar){
  # if either x or y are ".", or if no more than one level of a quantity, render a blank plot
  if(xvar == "." | yvar == "."){
    p <- ggplot(dataset) + geom_blank()
  } else {
    
    if(groupvar != "."){
      p <- ggplot(dataset,
                  aes_string(x = xvar,
                             y = yvar,
                             group = groupvar,
                             color = groupvar))
    } else {
      p <- ggplot(dataset,
                  aes_string(x = xvar,
                             y = yvar,
                             group = 1))
    }
    p <- p + geom_line(size = lsizevar) +
      geom_point(size = psizevar) +
      theme_grey(base_size = 18) 
    
    facets <- paste(rowvar,"~",colvar)
    if(facets != ". ~ .") p <- p + facet_grid(facets)
  }
  p
}
