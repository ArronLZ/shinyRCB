suppressMessages(library(tidyverse))
suppressMessages(library(shiny))
suppressMessages(library(shinydashboard))
suppressMessages(library(data.table))
suppressMessages(library(scales))
suppressMessages(library(ggplot2))
suppressMessages(library(edgeR))
suppressMessages(library(DESeq2))
suppressMessages(library(limma))
suppressMessages(library(ggsci))

source("src/predata_diffanalysis.R")
source("src/LZ_plot.valcano.R")
## a <- read.csv("data/example/example.annot.csv") 
## b <- fread("data/example/example.RNAseq.csv", data.table = F)
## c <- read.csv("data\\example\\example.group.csv", row.names = 1,
##  stringsAsFactors = T)
## eset <- Xena.process(df = b, annot = a)
## eset.group <- list(eset=eset, group=c, f_mark="luad.mRNA")
## diffan <- DEG.edgeR(eset.group, pval = 0.05, fdr = 0.1, logfc = 1)
## resdf = diffan$resdf
## DEGplot.volcano(df_valcano)

ui <- dashboardPage(
  # skin = "black",
  dashboardHeader(title = "Refined Cell Biotech SMU"),
  dashboardSidebar(
    sidebarMenu(
      menuItem("差异分析", tabName ="t1" , icon = icon("th")),#icon("dashboard")
      menuItem("化疗计算", tabName ="t2" , icon = icon("th")),
      menuItem("帮助", tabName = "t3", icon = icon("th"))
    )
  ),
  dashboardBody(
    tags$head(
      tags$link(rel = "stylesheet", type = "text/css", href = "custom.css")
    ),
    tabItems(
      # First tab content
      tabItem(tabName = "t1",
              fluidRow(
                column(12,
                box(width = 3,
                    fileInput("file_eset", "Upload Eset File (CSV)", accept = c(".csv"), multiple = F),
                    fileInput("file_group", "Upload Group File (CSV)", accept = c(".csv"), multiple = F),
                    fileInput("file_annot", "Upload Annot File (CSV)", accept = c(".csv"), multiple = F),
                    selectInput("pval", "pvalue", c(0.0001, 0.001, 0.01, 0.05,0.1,0.3,0.5,0.7, 1), selected=0.05),
                    selectInput("fdr", "FDR", c(0.0001, 0.001, 0.01, 0.05, 0.1, 0.2,0.3,0.5,0.7, 1), selected=0.1),
                    sliderInput("logfc", "log2FC", 0, 5, 0.1),
                    actionButton("btn", "Analyze")
                ),
                box(width = 9,
                  plotOutput("volcano_plot"),
                  br(),
                  textOutput('deg_nrow') # tableOutput("table")
                )
              )),
              fluidRow(
                column(9,
                box(width = 3,
                  tableOutput("table_eset")
                ),
                box(width = 3,
                  tableOutput("table_group")
                ),
                box(width = 3,
                  tableOutput("table_annot")
                ),
              )),
      ),
    # Second tab content
      tabItem(tabName = "t2",
        fluidRow(
                column(10,
                box(
                    selectInput('g','性别', c("男","女"), selected="男"),
                    sliderInput("h", "身高(cm):", 80, 220, 160),
                    sliderInput("w", "体重(kg):", 30, 120, 60),
                ),
                box(
                  textOutput('bsa'),
                  br(),
                  textOutput('bz')
                )
              )),
              fluidRow(
                column(10,
                box(
                  sliderInput("y", "年龄(岁):", 18, 90, 50),
                  sliderInput("cr", "肌酐(umol/L):", 0, 900, 30),
                  sliderInput("auc", "AUC(mg/ml/min):", 5, 10, 5)
                ),
                box(
                  textOutput('kb1'),
                  textOutput('kb2')
                )
              ))
      ),
      tabItem(tabName = "t3",
               h2("敬请期待")
      )
    )
  )
)

server <- function(input, output) {
  options(shiny.maxRequestSize=1000*1024^2)
  label_gene <- c("S100P","PDIA2","PRB3"," KLK14","SCGB3A2","ZACN") 
  input_value <- reactiveValues(pval=0.05, fdr=0.1,logfc=1)
  
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
    if(is.null(input$file_eset)&is.null(input$file_annot))
      return(NULL)
    else{
      df <- readeset_csv()
      df_sub <- df[1:3, 1:3]
      return(df_sub)
    }
  })
  output$table_group <- renderTable({
    if(is.null(input$file_eset)&is.null(input$file_annot))
      return(NULL)
    else{
      df <- readgroup_csv()
      df_sub <- df[1:3, 1:2]
      return(df_sub)
    }
  })
  output$table_annot <- renderTable({
    if(is.null(input$file_eset)&is.null(input$file_annot))
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
    
    # 火山图
    resdf <- diffan_edger()
    df_valcano <- resdf$resdf %>% distinct(Gene, .keep_all = T) %>% 
                     dplyr::select(Gene, log2FC, PValue, FDR) %>% 
                     dplyr::rename(log2FoldChange=log2FC, padj=FDR) %>% 
                     na.omit()
    output$volcano_plot <- renderPlot({
        DEGplot.volcano(result = df_valcano, logFC = input_value$logfc, 
                       adj_P = input_value$pval, 
                       label_geneset = intersect(label_gene, df_valcano$Gene)) %>% 
                       ggplotGrob() %>% cowplot::plot_grid() 
    })

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
    
    })
  
  ###
  c_bsa <- reactive(
    if (input$g == '男') {
      bsa <- round(0.00607*input$h + 0.0127*input$w - 0.0698, digits = 2)
    } else {
      bsa <- round(0.00586*input$h + 0.0126*input$w - 0.0461, digits = 2)
    }
  )
  
  c_crqcl <- reactive(
    if (input$g == '男') {
      crq <- 1.23*input$w*(140-input$y)/input$cr
    } else {
      crq <- 1.03*input$w*(140-input$y)/input$cr
    }
  )

  ###out
  #体表面积
  output$bsa <- renderText({
    paste0("体表面积(", input$g, "):",c_bsa())
  })
  
  # 白紫
  output$bz <- renderText({
    paste0("白蛋白紫杉醇(", input$g, ")(260):", 
           round(260*c_bsa(), digits = 2), 
           "mg")
  })
  
  #卡铂剂量（AUC方法，体表面积法）
  output$kb1 <- renderText({
    paste0("卡铂(", input$g, ")AUC:",
           round(input$auc*(c_crqcl() + 25), digits = 2), 
           "mg")
  })
  output$kb2 <- renderText({
    paste0("卡铂(", input$g, ")m2:    ",
           round(400*c_bsa(), digits = 2), 
           "mg")
  })
  
}

#options(shiny.port="0.0.0.0", shiny.port = 9012)
cat("按住Ctrl并单击下面网址即可打开应用\n")
browseURL("http://127.0.0.1:9012")
shinyApp(ui = ui, server = server)
