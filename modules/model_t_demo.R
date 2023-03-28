label_gene <- c("S100P","PDIA2","PRB3","KLK14","SCGB3A2","ZACN") 
load("data/example/example.diff3set.env.RData") # env_data

ui_t_demo <- function(id) {
  ns <- NS(id)
  tagList(
    fluidRow(
      column(12,
             box(width = 6, height = 1000,
                 h3("差异分析需要上传三个csv文件："),
                 textOutput(ns("my_text2")),
                 br(),
                 "1. 基因表达文件格式要求：第一列为gene_id名称，第二至N列为样本名称",
                 tableOutput(ns("example_eset")),
                 "2. 样本分组文件格式要求：第一列为样本名称(对应表达文件的样本名)，后续列为分组信息，目前暂仅支持type进行差异分析，请将需要分析的表型信息列设置为type",
                 tableOutput(ns("example_group")),
                 "3. 基因注释文件格式要求：第一列为gene id(对应表达文件的第一列),第二列为要转换成的目标gene symbol",
                 tableOutput(ns("example_annot")),
                 "4. 标记基因：请以英文输入法状态下的逗号间隔每个基因",
                 textOutput(ns('example_gene'))
             ),
             box(width = 6, height = 1000,
             ),
      )  
    )
  )
}


server_t_demo <- function(id, exdata=env_exdata) {
  moduleServer(
    id,
    function(input, output, session) {
      ##
      F_renderTable <- function(data) {
        renderTable({
          data[1:3,1:3]
        })
      }
      ##
      output$example_eset <- F_renderTable(exdata$eset)
      output$example_group <- F_renderTable(exdata$group)
      output$example_annot <- F_renderTable(exdata$annot)
      output$example_gene <- renderText({
        paste(label_gene, collapse = ',')
        })
      ##
      output$my_text2 <- renderText({
        "  csv格式文件可在excel中打开表格另存为csv格式即可"
      })
    }
  )
}