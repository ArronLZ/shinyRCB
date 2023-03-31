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
suppressMessages(library(survival))
suppressMessages(library(survivalROC))
suppressMessages(library(survminer))
suppressMessages(library(waiter))

source("src/predata_diffanalysis.R")
source("src/LZ_plot.valcano.R")
source("src/LZ_plot.surv.R")
source("modules/model_t_diff.R")
source("modules/model_t_hualiao.R")
source("modules/model_t_demo.R")
source("modules/model_t_surv.R")

ui <- dashboardPage(
  # skin = "black",
  dashboardHeader(title = "Refined Cell Biotech SMU"),
  dashboardSidebar(
    sidebarMenu(
      menuItem("操作说明", tabName = "t_demo", icon = icon("th")),
      menuItem("差异分析", tabName ="t_diff" , icon = icon("th")),
      menuItem("富集分析", tabName ="t_pathway" , icon = icon("th")),
      menuItem("生存分析", tabName ="t_surv" , icon = icon("th")),#icon("dashboard")
      menuItem("上游分析", tabName ="t0" , icon = icon("th")),
      menuItem("化疗计算", tabName ="t_hualiao" , icon = icon("th")),
      menuItem("帮助", tabName = "t_help", icon = icon("th"))
    )
  ),
  dashboardBody(
    tags$head(
      tags$link(rel = "stylesheet", type = "text/css", href = "custom.css"),
      tags$style(
       ".notification-center {
           position: fixed;
           top: 50%;
           left: 50%;
           transform: translate(-50%, -50%);
           background-color: rgba(128, 128, 128, 0.3);
           color: darkred;
           font-weight: bold;
           font-size: 24px;
       }"
      )
    ),
    tabItems(
      tabItem(tabName = "t_demo",
              uiOutput(outputId = "ui_demo")
      ),
      tabItem(tabName = "t0",
              h2("仅管理员可用")
      ),
      tabItem(tabName = "t_diff",
              uiOutput(outputId = "ui_diff")
      ),
      tabItem(tabName = "t_pathway",
              h2("即将上线")
      ),
      tabItem(tabName = "t_surv",
              uiOutput(outputId = "ui_surv")
      ),
      # Second tab content
      tabItem(tabName = "t_hualiao",
              uiOutput(outputId = "ui_hualiao")
      ),
      tabItem(tabName = "t_help",
              h2("敬请期待")
      )
    )
  )
)

server <- function(input, output, session) {
  # t_demo
  output$ui_demo <- renderUI({
    ui_t_demo("demo")
  })
  server_t_demo("demo")
  # t_diff
  output$ui_diff <- renderUI({
    ui_t_diff("diff")
  })
  server_t_diff("diff")
  # t_surv
  output$ui_surv <- renderUI({
    ui_t_surv("surv")
  })
  server_t_surv("surv")
  # t_hualiao
  output$ui_hualiao <- renderUI({
    ui_t_hualiao("hualiao")
  })
  server_t_hualiao("hualiao")
}

#options(shiny.port="0.0.0.0", shiny.port = 9012)
cat("按住Ctrl并单击下面网址即可打开应用\n")
browseURL("http://127.0.0.1:9012")
shinyApp(ui = ui, server = server)
