p.range <- c(0.0001, 0.001, 0.01, 0.05, seq(0.1, 1, 0.1))

ui_t_diff <- function(id) {
  ns <- NS(id)
  tagList(
    fluidRow(
      column(12,
             box(width = 2, height = 1000,
                 fileInput(ns("file_eset"), "Upload Eset File (CSV)", accept = c(".csv"),
                           multiple = F, buttonLabel = "Load"),
                 fileInput(ns("file_group"), "Upload Group File (CSV)", accept = c(".csv"),
                           multiple = F, buttonLabel = "Load"),
                 fileInput(ns("file_annot"), "Upload Annot File (CSV)", accept = c(".csv"),
                           multiple = F, buttonLabel = "Load"),
                 selectInput(ns("pval"), "pvalue", p.range, selected=0.05),
                 selectInput(ns("fdr"), "FDR", p.range, selected=0.1),
                 sliderInput(ns("logfc"), "log2FC", value = 1, min = 0, max = 5, step = 0.1),
                 textInput(ns("marker"), "Please type Gene name separated by commas!!!", placeholder = "EGFR,BRAF,TP53"),
                 actionButton(ns("btn"), "Analyze")
             ),
             box(width = 4, height = 1000,
                 tableOutput(ns("table_eset")),
                 tableOutput(ns("table_group")),
                 tableOutput(ns("table_annot")),
                 textOutput(ns('marker_gene'))
             ),
             box(width = 6, height = 1000,
                 plotOutput(ns("volcano_plot")),
                 br(),
                 textOutput(ns('deg_nrow')) # tableOutput("table")
             ),
      )  
    )
  )
}

server_t_diff <- function(id) {
  moduleServer(
    id,
    function(input, output, session) {
      options(shiny.maxRequestSize=1000*1024^2)
      options(shiny.useragg = FALSE)
      label_gene <- c("S100P","PDIA2","PRB3","KLK14","SCGB3A2","ZACN") 
      # S100P,PDIA2,PRB3,KLK14,SCGB3A2,ZACN
      input_value <- reactiveValues(pval=0.05, fdr=0.1,logfc=1,marker=label_gene)
      
      readannot_csv <- reactive({
        read.csv(input$file_annot$datapath) 
      })
      readeset_csv <- reactive({
        fread(input$file_eset$datapath, data.table = F)
      })
      readgroup_csv <- reactive({
        read.csv(input$file_group$datapath, row.names = 1, stringsAsFactors = T)
      })
      
      diffan_edger <- reactive({
        eset <- Xena.process(df = readeset_csv(), annot = readannot_csv())
        eset.group <- list(eset=eset, group=readgroup_csv(), f_mark="luad.mRNA")
        diffan <- DEG.edgeR(eset.group, pval = 0.05, fdr = 0.1, logfc = 1)
        # save(diffan, file = "diffan.RData", compress = T)
        return(diffan)
      })
      
      output$table_eset <- renderTable({
        if(is.null(input$file_eset))
          return(NULL)
        else{
          df <- readeset_csv()
          df_sub <- df[1:3, 1:3]
          return(df_sub)
        }
      })
      output$table_group <- renderTable({
        if(is.null(input$file_group))
          return(NULL)
        else{
          df <- readgroup_csv()
          df_sub <- df[1:3, 1:2]
          return(df_sub)
        }
      })
      output$table_annot <- renderTable({
        if(is.null(input$file_annot))
          return(NULL)
        else{
          df <- readannot_csv()
          df_sub <- df[1:3, 1:2]
          return(df_sub)
        }
      })
      
      observeEvent(input$btn, {
        input_value$pval = input$pval %>% as.numeric()
        input_value$fdr = input$fdr %>% as.numeric()
        input_value$logfc = input$logfc %>% as.numeric()
        input_value$marker = input$marker %>% str_split(., pattern = ",", simplify = T) %>% .[1,]
        
        # 火山图
        resdf <- diffan_edger()
        df_valcano <- resdf$resdf %>% distinct(Gene, .keep_all = T) %>% 
          dplyr::select(Gene, log2FC, PValue, FDR) %>% 
          dplyr::rename(log2FoldChange=log2FC, padj=FDR) %>% 
          na.omit()
        output$volcano_plot <- renderPlot(res = 100,
          DEGplot.volcano(result = df_valcano, logFC = input_value$logfc, 
                          adj_P = input_value$pval, 
                          label_geneset = intersect(input_value$marker, df_valcano$Gene)) %>% 
            print() #ggplotGrob() %>% cowplot::plot_grid() 
        )
        
        # 差异基因
        newdeg <- resdf$resdf %>% filter(abs(log2FC) > input_value$logfc & 
                                           PValue < input_value$pval & 
                                           FDR < input_value$fdr)
        output$deg_nrow <- renderText({
          paste0("pval= ",  input_value$pval, "\n",
                 "FDR= ",  input_value$fdr, "\n",
                 "Log2FC= ",  input_value$logfc, "\n",
                 "DEG numbers: ", nrow(newdeg)
          )
        })
        
        output$marker_gene<- renderText({
          input_value$marker
        })
      })
    }
  )
}
