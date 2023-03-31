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


header <- dashboardHeader(title = "Refined Cell Biotech SMU")
sidebar <- dashboardSidebar(
    sidebarMenu(
      menuItem("登录账号", tabName = "t_demo", icon = icon("th")),
      menuItem("差异分析", tabName ="t_diff" , icon = icon("th")),
      menuItem("富集分析", tabName ="t_pathway" , icon = icon("th")),
      menuItem("生存分析", tabName ="t_surv" , icon = icon("th")),#icon("dashboard")
      menuItem("上游分析", tabName ="t0" , icon = icon("th")),
      menuItem("化疗计算", tabName ="t_hualiao" , icon = icon("th")),
      menuItem("敬请期待", tabName = "t_help", icon = icon("th"))
    )
  )
body <- dashboardBody(
    tags$head(
      tags$link(rel = "stylesheet", type = "text/css", href = "custom.css"),
      tags$style(
       ".notification-center {
           position: fixed;
           top: 50%;
           left: 50%;
           transform: translate(-50%, -50%);
           background-color: rgba(128, 128, 128, 0.1);
           color: rgba(70,130,180,1);
           font-weight: bold;
           font-size: 24px;
           position: fixed;
           padding: 8px;
           border-radius: 4px;
           box-shadow: 0 2px 4px rgba(0, 0, 0, 0.2);
           z-index: 9999;
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
              h2("敬请期待...")
      )
    )
  )
ui <- dashboardPage(header, sidebar, body)

##password 
login_details <- data.frame(user = c("jun", "rcb", "guest"),
                            pswd = c("litchi123", "xwk123", "xwk"))

login <- box(
  title = "Welcome to RCB",
  textInput("userName", "Please Enter your UserName:"),
  passwordInput("passwd", "Please Enter your PassWord:"),
  br(),
  actionButton("Login", "Log in")
)

server <- function(input, output, session) {
    # To logout back to login page
  login.page = paste(
    isolate(session$clientData$url_protocol),
    "//",
    isolate(session$clientData$url_hostname),
    ":",
    isolate(session$clientData$url_port),
    sep = ""
  )
  
  USER <- reactiveValues(Logged = F)
  observe({
    if (USER$Logged == FALSE) {
      if (!is.null(input$Login)) {
        if (input$Login > 0) {
          Username <- isolate(input$userName)
          Password <- isolate(input$passwd)
          Id.username <- which(login_details$user == Username)
          Id.password <- which(login_details$pswd == Password)
          if (length(Id.username) > 0 & length(Id.password) > 0){
            if (Id.username == Id.password) {
              USER$Logged <- TRUE
            }
          }
        }
      }
    }
  })
  
  # t_demo
  output$ui_demo <- renderUI({
    if (USER$Logged == TRUE) {
      ui_t_demo("demo")
    } else {
      dashboardBody(login)
    }
  })
  server_t_demo("demo")

  # t_diff
  output$ui_diff <- renderUI({
    if (USER$Logged == TRUE) {
    ui_t_diff("diff")
    } else {
       h3("请先登录账号")
    }
  })
  server_t_diff("diff")
  # t_surv
  output$ui_surv <- renderUI({
    if (USER$Logged == TRUE) {
    ui_t_surv("surv") 
    } else {
       h3("请先登录账号")
    }
  })
  server_t_surv("surv")
  # t_hualiao
  output$ui_hualiao <- renderUI({
    if (USER$Logged == TRUE) {
    ui_t_hualiao("hualiao") 
    } else {
       h3("请先登录账号")
    }
  })
  server_t_hualiao("hualiao")
}

#options(shiny.port="0.0.0.0", shiny.port = 9012)
cat("按住Ctrl并单击下面网址即可打开应用\n")
browseURL("http://127.0.0.1:9012")
shinyApp(ui = ui, server = server)
