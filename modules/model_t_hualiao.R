counterButton2 <- function(id) {
  ns <- NS(id)
  tagList(
    fluidRow(
      column(10,
             box(
               selectInput(ns('g'),'性别', c("男","女"), selected="男"),
               sliderInput(ns("h"), "身高(cm):", 80, 220, 160),
               sliderInput(ns("w"), "体重(kg):", 30, 120, 60),
             ),
             box(
               textOutput(ns('bsa')),
               br(),
               textOutput(ns('bz'))
             )
      )),
    fluidRow(
      column(10,
             box(
               sliderInput(ns("y"), "年龄(岁):", 18, 90, 50),
               sliderInput(ns("cr"), "肌酐(umol/L):", 0, 900, 30),
               sliderInput(ns("auc"), "AUC(mg/ml/min):", 5, 10, 5)
             ),
             box(
               textOutput(ns('kb1')),
               textOutput(ns('kb2'))
             )
      ))
  )
}

counterServer2 <- function(id) {
  moduleServer(
    id,
    function(input, output, session) {
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
  )
}