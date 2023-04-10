library(shiny)
library(shinycssloaders)

ui <- fluidPage(
  actionButton("plotButton", "Plot"),
  uiOutput("plotOutput")
)

server <- function(input, output) {
  
  # Define a reactive variable to store the status of the calculation
  calcStatus <- reactiveVal("not started")
  
  # Define the plotting function
  plotFunction <- function() {
    # Set the calculation status to "in progress"
    calcStatus("in progress")
    
    # Do some time-consuming calculation
    Sys.sleep(5)
    
    # Create the plot
    plot(1:10)
    
    # Set the calculation status to "complete"
    calcStatus("complete")
  }
  
  # Define the action to take when the button is clicked
  observeEvent(input$plotButton, {
    # Run the plotting function
    plotFunction()
   
  })
  # Define the UI output to display the plot or the "in progress" message
  output$plotOutput <- renderUI({
    if (calcStatus() == "complete") {
      # If the calculation is complete, show the plot
      plotOutput("myPlot")
    } else {
      # If the calculation is in progress, show the spinner
      withSpinner(p("Plotting in progress..."))
    }
  })
  
  # Define the plot output
  output$myPlot <- renderPlot({
    plot(1:10)
  })
  
}
ui <- fluidPage(
  waiter::use_waiter(),
  actionButton("go", "go"),
  plotOutput("plot"),
)

server <- function(input, output, session) {
  data <- eventReactive(input$go, {
    waiter::Waiter$new(id = "plot")$show()
    
    Sys.sleep(3)
    data.frame(x = runif(50), y = runif(50))
  })
  
  output$plot <- renderPlot(plot(data()), res = 96)
}

shinyApp(ui, server)

install.packages('waiter')
?showNotification()
eset_clin[,20] %>% fivenum()
eset_clin %>% names() %>% .[20]

mytheme_prism <- function() {
  # max, by=5 # breaks = seq(0, max, by = by) ##limits = c(0, 40),
  list(
    scale_y_continuous(expand = c(0,0)),
    scale_x_continuous(expand = c(0,0)),
    #ggprism::theme_prism()
    theme_classic()
  )
}
LZplot.dens <- function(data, col) {
  five <- fivenum(data[, col])
  ggplot(data, aes(x = data[, col])) +
    geom_density(fill = "gray", alpha = 0.5) +
    geom_vline(xintercept = five, col = "blue", size=1) +
    geom_vline(xintercept = median(data[, col]), col = "red", size=1) +
    geom_vline(xintercept = mean(data[, col]), col = "yellow", size=1, lty=2) +
    annotate("text", x = five, y = 0, 
             label = c(paste0("Min:", round(five[1],2)), 
                       paste0("Q1:", round(five[2],2)), 
                       paste0("Median:", round(five[3],2)), 
                       paste0("Q3:", round(five[4],2)), 
                       paste0("Max:", round(five[5],2)) ), 
             color = c("black", "black", "red", "black", "black"),
             hjust = 1, vjust = c(-1,-2.5,-4.5,-6,-8)) +
    annotate("text", x = mean(data[, col]), y = 0, 
             label = paste0("Mean:", round(mean(data[, col]))), color="yellow",
             hjust = 1, vjust = -1) +
    mytheme_prism() + xlab(label = col) +ylab(label = NULL)
}
#ggplot(eset_clin, aes(x=PIGV)) + geom_density(fill="grey")
LZplot.dens(eset_clin, "PIGV")

library(survminer)
?ggsurvplot
