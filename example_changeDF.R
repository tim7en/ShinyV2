summarized <- data.frame(id = 1:20, group = letters[1:4], TY_COMP = runif(20), LY_COMP = runif(20))

library(shiny)

ui <- fluidPage(
  verbatimTextOutput("text"),
  actionButton("btn", "Add the 'cac' column to summarized")
)

server <- function(input, output){
  rv <- reactiveValues(summarized = summarized)
  
  output$text <- renderPrint(rv$summarized)
  
  observeEvent(input$btn,{
    print (ncol (rv$summarized))
  })
  
  observeEvent(input$btn, {
    rv$summarized$cac <- (summarized$TY_COMP / summarized$LY_COMP - 1)*(rnorm (1,2))
  })
  
  summarized_mod <- reactive({
    summarized()$TY_COMP / summarized()$LY_COMP-1
  })
}

shinyApp(ui, server)