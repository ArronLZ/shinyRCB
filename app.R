#library(tidyverse)
#library(shiny)
#library(shinydashboard)
#library(DESeq2)

load("E:/OneDrive/Desktop/DM.R1/v5/.diffan.RData")

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
                column(10,
                box(
                    fileInput("file_eset", "Upload Eset File (CSV)", accept = c(".csv"), multiple = F),
                    fileInput("file_group", "Upload Group File (CSV)", accept = c(".csv"), multiple = F),
                    selectInput("pval", "pvalue", c(0.0001, 0.001, 0.01, 0.05, 1), selected=0.05),
                    selectInput("fdr", "FDR", c(0.0001, 0.001, 0.01, 0.05, 0.1, 0.2, 1), selected=0.1),
                    sliderInput("logfc", "log2FC", 0, 5, 1),
                    actionButton("btn", "Analyze")
                ),
                box(
                  plotOutput("volcano_plot"),
                  br(),
                  textOutput('deg_nrow') # tableOutput("table")
                )
              )),
              fluidRow(
                column(10,
                box(
                ),
                box(
                  tableOutput("table_eset")
                )
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
  readfile_csv <- reactive({
    read.csv(input$file_eset$datapath, header = TRUE, stringsAsFactors = FALSE)
  })

  input_value <- reactiveValues(pval=1, fdr=1,logfc=1)

  output$table_eset <- renderTable({
    if(is.null(input$file_eset))
      return(NULL)
    else{
      df <- readfile_csv()
      df_sub <- df[1:5, 1:5]
      return(df_sub)
    }
  })

  observeEvent(input$btn, {
    input_value$pval = input$pval %>% as.numeric()
    input_value$fdr = input$fdr %>% as.numeric()
    input_value$logfc = input$logfc %>% as.numeric()
    
    res_df <- diffan$resdf
    deg <- res_df %>% filter(pvalue < input_value$pval & padj < input_value$fdr & 
                               abs(log2FoldChange) > input_value$logfc)
    
    output$volcano_plot <- renderPlot({
        ggplot(res_df, aes(x = log2FoldChange, y = -log10(pvalue))) + 
          geom_point() + 
          xlab("Log2 Fold Change") + ylab("-log10(p-value)") + 
          ggtitle("Volcano Plot") +
          geom_vline(xintercept = c(-input_value$logfc, input_value$logfc), 
                     linetype = "dashed") + 
          geom_hline(yintercept = -log10(input_value$pval), linetype = "dashed")
    })

    output$deg_nrow <- renderText({
      # st()
      paste0("pval=",  input_value$pval, 
             "FDR=",  input_value$fdr, 
             "Log2FC=",  input_value$logfc, 
             ": ", nrow(deg) )
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
