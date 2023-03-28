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
source("modules/model_t_diff.R")
source("modules/model_t_hualiao.R")

ui <- dashboardPage(
  # skin = "black",
  dashboardHeader(title = "Refined Cell Biotech SMU"),
  dashboardSidebar(
    sidebarMenu(
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
      tags$link(rel = "stylesheet", type = "text/css", href = "custom.css")
    ),
    tabItems(
      tabItem(tabName = "t0",
              h2("仅管理员可用")
      ),
      # First tab content
      tabItem(tabName = "t_diff",
              uiOutput(outputId = "ui_t_diff")
      ),
      tabItem(tabName = "t_pathway",
              h2("即将上线")
      ),
      tabItem(tabName = "t_surv",
              h2("正在维护")
      ),
      # Second tab content
      tabItem(tabName = "t_hualiao",
              uiOutput(outputId = "ui_t_hualiao")
      ),
      tabItem(tabName = "t_help",
              h2("敬请期待")
      )
    )
  )
)

server <- function(input, output, session) {
  output$ui_t_diff <- renderUI({
    ui_t_diff("diff")
  })
  server_t_diff("diff")
  ###
  output$ui_t_hualiao <- renderUI({
    ui_t_hualiao("hualiao")
  })
  server_t_hualiao("hualiao")
}

#options(shiny.port="0.0.0.0", shiny.port = 9012)
cat("按住Ctrl并单击下面网址即可打开应用\n")
browseURL("http://127.0.0.1:9012")
shinyApp(ui = ui, server = server)
