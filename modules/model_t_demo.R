label_gene <- c("S100P","PDIA2","PRB3","KLK14","SCGB3A2","ZACN") 
load("data/example")
ui_t_demo <- function(id) {
  ns <- NS(id)
  tagList(
    fluidRow(
      column(12,
             box(width = 4, height = 1000,
                 tableOutput(ns("example_eset")),
                 tableOutput(ns("example_group")),
                 tableOutput(ns("example_annot")),
                 textOutput(ns('example_gene'))
             ),
             box(width = 4, height = 1000,
             ),
             box(width = 4, height = 1000,
             ),
      )  
    )
  )
}



server_t_demo <- function(id) {
  moduleServer(
    id,
    function(input, output, session) {

    }
  ))
}