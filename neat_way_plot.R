library(shiny)

ui <- fluidPage(
  
  sidebarLayout(
    sidebarPanel(
      actionButton("update", "Update")
    ),
    
    mainPanel(
      column(4, tableOutput('mytable')),
      column(4, tableOutput('mytable2')),
      column(4, plotOutput("autoupdate_plot"))
    )
    
    
  )
)

server <- function(input, output) {
  
  values <- reactiveValues(df = RenderMyTable())
  
  observe({
    invalidateLater(1000)
    values$df <- RenderMyTable()  # This does not update after 1 sec
  })
  
  observeEvent(input$update, {
    values$df <- RenderMyTable()  # This does update upon clicking
  })
  
  output$mytable  <- renderTable(values$df)  # Depends on reactiveValues
  
  autoInvalidate <- reactiveTimer(1000)
  
  output$mytable2 <- renderTable({
    autoInvalidate()
    RenderMyTable()  # >This does update after 1 sec
  })
  
  output$autoupdate_plot <- renderPlot({
    
    invalidateLater(2000)
    hist(rnorm(isolate(RenderMyTable())), main = "Auto Hist")
    
  })
}

time1 <- Sys.time()  # Start time
df <- data.frame(a = 1:1000)  # Some data

RenderMyTable <- function(){
  # Seconds since start time
  time2 <- as.integer(difftime(Sys.time(), time1, units="secs"))
  
  df.now <- df[1:time2,]  # Updates each second
  
  df.now
}

shinyApp(ui = ui, server = server)