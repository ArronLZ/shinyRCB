eset_clin <- fread("data/tcga/DataClean_LUAD.tpm&os.csv.gz", data.table = F)
rindex <- str_sub(eset_clin[,1],14,16) == "01A"
eset_clin <- eset_clin[rindex, ]
eset_clin[,1] <- str_sub(eset_clin[,1], 1, 12)
eset_clin <- eset_clin %>% distinct(., sample_id, .keep_all = T)
eset_clin[, 4:ncol(eset_clin)] <- log2(eset_clin[, 4:ncol(eset_clin)] + 1)
###
list.select <- c("median", "mean", "auto cut point", "upper75 vs low25", 
                 "upper25 vs low75", "upper25 vs low25")


ui_t_surv <- function(id) {
  ns <- NS(id)
  tagList(
    fluidRow(
      column(12,
             box(width = 3, height = 1000,
                 #fileInput(ns("file_eset"), "Upload Eset File (CSV)", accept = c(".csv"),
                 #         multiple = F, buttonLabel = "Load"),
                 textInput(ns("surv_gene"), "Please type Gene name", placeholder = "EGFR"),# separated by commas!!!,BRAF,TP53
                 selectInput(ns("surv_method"), "Method", list.select, selected="median"),
                 actionButton(ns("surv_btn"), "Analyze")
             ),
             box(width = 9, height = 1000,
                 plotOutput(ns("surv_plot")),
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
      #  read.csv(input$file_annot$datapath) 
      #})
      observeEvent(input$surv_btn, {
        input_value$method = input$surv_method# %>% as.numeric()
        input_value$gene = input$surv_gene # %>% str_split(., pattern = ",", simplify = T) %>% .[1,]
        # 生存图
        eset_os.cut <- xena.surv.cut(eset_os.df = eset_os, tagetGene=input_value$gene, 
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
