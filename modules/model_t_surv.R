###
list.select <- c("median", "mean", "auto cut point", "upper75 vs low25", 
                 "upper25 vs low75", "upper25 vs low25")
list.dataset <- c("LUAD", "ESCC", "ESAD")


ui_t_surv <- function(id) {
  ns <- NS(id)
  tagList(
    fluidRow(
      column(12,
             box(width = 3, height = 1000,
                 #fileInput(ns("file_eset"), "Upload Eset File (CSV)", accept = c(".csv"),
                 #         multiple = F, buttonLabel = "Load"),
                 selectInput(ns("surv_dataset"), "DataSet", list.dataset, selected="LUAD"),
                 textInput(ns("surv_gene"), "Please type Gene name", value = "EGFR", placeholder = "EGFR"),# separated by commas!!!,BRAF,TP53
                 selectInput(ns("surv_method"), "Method", list.select, selected="median"),
                 actionButton(ns("surv_btn"), "Analyze")
             ),
             box(width = 9, height = 1000,
                 #plotOutput(ns("surv_plot"))
                plotOutput(ns("surv_plot"))
             ),
      )  
    )
  )
}

server_t_surv <- function(id, eset_os=eset_clin) {
  moduleServer(
    id,
    function(input, output, session) {
      #options(shiny.maxRequestSize=1000*1024^2)
      options(shiny.useragg = FALSE)
      #c("S100P","PDIA2","PRB3","KLK14","SCGB3A2","ZACN") 
      input_value <- reactiveValues(method="median", gene="EGFR")
      
      #readannot_csv <- reactive({
      #  #read.csv(input$file_annot$datapath) 
      #})

      loaddataset <- reactive({
        # list.dataset <- c("LUAD", "ESCC", "ESAD")
        if (input$surv_dataset == "LUAD") {
          eset_os <- shiny.readeset("data/tcga/DataClean_LUAD.tpm&os.csv.gz", gene1.col = 4)
        } else if (input$surv_dataset == "ESCC") {
          eset_os <- shiny.readeset("data/tcga/DataClean_ESCC.tpm&clic.csv.gz", gene1.col = 13)
        } else if (input$surv_dataset == "ESAD") {
          eset_os <- shiny.readeset("data/tcga/DataClean_ESAD.tpm&clic.csv.gz", gene1.col = 13)
        }
        return(eset_os)
      })

      #observeEvent(input$surv_btn, {
      #  input_value$method = input$surv_method# %>% as.numeric()
      #  input_value$gene = input$surv_gene # %>% str_split(., pattern = ",", simplify = T) %>% .[1,]
      #  # 生存图
      #  eset_os.cut <- xena.surv.cut(eset_os.df = loaddataset(), tagetGene=input_value$gene, 
      #                              method=input_value$method)
      #  ps <- xena.surv.getPvale(rt = eset_os.cut, num.tran = 365, main.text = input_value$gene)
      #  #
      #  output$surv_plot <- renderPlot(width = 600, height=500, res = 100,
      #        ps$pic %>% print() #ggplotGrob() %>% cowplot::plot_grid() 
      #  )
      #  #output$surv_plot <- renderPlot(width = 600, height=500, res = 100,
      #  #  ps$pic %>% print() #ggplotGrob() %>% cowplot::plot_grid() 
      #  #)
      #  #output$surv_plot <- renderUI({
      #  #  if (calcStatus() == "complete") {
      #  #    # If the calculation is complete, show the plot
      #  #      plotOutput("my_surv_plot")
      #  #  } else {
      #  #    # If the calculation is in progress, show the spinner
      #  #    withSpinner(p("Plotting in progress..."))
      #  #  }
      #  #})
      #})
      observeEvent(input$surv_btn, {
        id <- showNotification('通知', 
                               tags$div(
                                       "正在计算中，请勿关闭网页...",
                                       class = "notification-center"
                                       ),
                               duration = NULL, closeButton = FALSE)
        on.exit(removeNotification(id), add = TRUE)
        input_value$method = input$surv_method# %>% as.numeric()
        input_value$gene = input$surv_gene #%>% as.character() #%>% str_split(., pattern = ",", simplify = T) %>% .[1,]
        # 生存图
        eset_os.cut <- xena.surv.cut(eset_os.df = loaddataset(), tagetGene=input_value$gene, 
                                    method=input_value$method)
        ps <- xena.surv.getPvale(rt = eset_os.cut, num.tran = 365, main.text = input_value$gene)
        #
        output$surv_plot <- renderPlot(width = 600, height=500, res = 100,
              ps$pic %>% print() #ggplotGrob() %>% cowplot::plot_grid() 
        )
      })

    }
  )
}