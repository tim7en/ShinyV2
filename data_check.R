library(shiny)
library(DT)
library(plyr)
library(shinydashboard)
library (zCompositions)

# checks for negatives
is.negative <- function(x) {
  options(warn = -1)
  x <- try(as.numeric(x))
  if ((is.numeric(x) || is.integer(x))) {
    return(length(which(na.omit(x) < 0)))
  } else {
    return(NULL)
  }
  options(warn = 1)
}

# checks for zeros
is.zero <- function(x) {
  options(warn = -1)
  x <- try(as.numeric(x))
  if ((is.numeric(x) || is.integer(x))) {
    return(length(which(x == 0)))
  } else {
    return(NULL)
  }
  options(warn = 1)
}

# read in the data table
x <- read.csv("SourceData.csv")

# function to output summary table
dfSummary <- function(x) {
  tabSummary <- NULL
  tabSummary <- rbind(tabSummary, sapply(x, function(x) class(x)))
  tabSummary <- rbind(tabSummary, sapply(x, function(x) sum(is.na(x))))
  tabSummary <- rbind(tabSummary, sapply(x, function(x) nlevels(x)))
  tabSummary <- rbind(tabSummary, apply(x, 2, is.negative))
  tabSummary <- rbind(tabSummary, apply(x, 2, is.zero))
  tabSummary
}


# summary table
dat <- dfSummary(x)
dat <- as.data.frame(dat)
rownames(dat) <- c("Type", "NAs", "Unique", "Negative", "Zeros")

# 3 columns, first column is selection available of columns with factors
# factors available in the selected column
# select column - rename with the name in the 3rd column, use apply and then reset
# column where 0 found, if there is a zero, warn
# column where negative found, if there is a negative, warn
# if it is missing, ask replace with mean ? , median ? or keep it and remove with outliers

# match source and target column names

# define UI for app that draws a histogram ----
ui <- dashboardPage(
  dashboardHeader(),
  dashboardSidebar(
    tabsetPanel(
      tabPanel(
        "Factors",
        uiOutput("dat_facts"),
        uiOutput("dat_levels"),
        uiOutput("new_levels"),
        uiOutput("accept_revalue")
      ),
      tabPanel(
        "Negatives",
        uiOutput ('dat_negatives'),
        uiOutput ('dat_neg_text')
      ),
      tabPanel(
        "Zeros",
        uiOutput ('dat_zeros'),
        uiOutput ('dat_zero_text')
      ),
      tabPanel(
        'Missing',
        uiOutput ('dat_missing'),
        uiOutput ('dat_missing_sub'),
        uiOutput ('dat_missing_adj')
      )
    )
  ),
  dashboardBody(
    tabsetPanel(
      tabPanel(
        "Summary",
        box(
          title = "Summary table", status = "danger", height = "630", width = 12, solidHeader = T,
          DTOutput(outputId = "distTable"),
          style = "height:'auto'; overflow-y: scroll;overflow-x: scroll"
        )
      ),
      tabPanel(
        "Data",
        box(
          width = 12,
          title = "Data table", solidHeader = TRUE, status = "primary",
          DTOutput("dat_levels_dt"),
          style = "height:'auto'; overflow-y: scroll;overflow-x: scroll;"
        )
      )
    )
  )
)

# Define server logic required to draw a histogram ----
server <- function(input, output) {
  dat_rv <- reactiveValues(dtab = x)
  dat_rv_sum <- reactiveValues (dtab = dat)
  
  output$distTable <- renderDT({
    dat_rv_sum$dtab
  })

  output$dat_facts <- renderUI({
    selectInput("facts", "Factors", choices = colnames(dat)[which(as.numeric(as.matrix(dat[3, ])) > 0)], selected = NULL)
  })

  output$dat_negatives <- renderUI ({
    selectInput("negts", "Negatives", choices = colnames(dat)[which(as.numeric(as.matrix(dat[4, ])) > 0)], selected = NULL, multiple = TRUE)
  })
  
  output$dat_neg_text <- renderUI ({
    textInput('neg_text', 'Input equation')
  })
  
  output$dat_zeros <- renderUI ({
    selectInput("zeros", "Zeros", choices = colnames(dat)[which(as.numeric(as.matrix(dat[5, ])) > 0)], selected = NULL, multiple = TRUE)
  })
  
  output$dat_missing <- renderUI ({
    selectInput("missing", "Missing", choices = colnames(dat)[which(as.numeric(as.matrix(dat[2, ])) > 0)], selected = NULL, multiple = TRUE)
  })
  
  output$dat_missing_sub <- renderUI ({
    col_indx <- which(names(dat_rv$dtab) == input$facts)
    selectInput("missing_lvl", "Levels", choices = levels(dat_rv$dtab[, col_indx]), multiple = TRUE, selected = FALSE)
  })
  
  output$dat_missing_adj <- renderUI ({
    selectInput("missing_adj", "Fill", choices = c('mean', 'median', 'single impute', 'multiple impute'), selected = NULL, multiple = TRUE)
  })
  
  output$dat_zero_text <- renderUI ({
    textInput('zer_text', 'Input equation')
  })

  output$dat_levels <- renderUI({
    col_indx <- which(names(dat_rv$dtab) == input$facts)
    selectInput("flevels", "Levels", choices = levels(dat_rv$dtab[, col_indx]), multiple = TRUE, selected = FALSE)
  })

  output$new_levels <- renderUI({
    textInput("new_names", "Names")
  })

  output$accept_revalue <- renderUI({
    actionButton("accept_revalue1", "Accept")
  })

  dat_levels_rename <- eventReactive(input$accept_revalue1, {
    col_indx <- which(names(dat_rv$dtab) == input$facts)
    new_names <- try(unlist(strsplit(input$new_names, ",")))
    if ((length(input$flevels) == 1) && (length(new_names) == 1)) {
      dat_rv$dtab[, col_indx] <- mapvalues(dat_rv$dtab[, col_indx], from = c(input$flevels), to = c(new_names))
    } else if ((length(input$flevels) > 1) && (length(new_names) == 1)) {
      new_names <- rep(new_names, length(input$flevels))
      dat_rv$dtab[, col_indx] <- mapvalues(dat_rv$dtab[, col_indx], from = c(input$flevels), to = c(new_names))
    } else if ((length(input$flevels) > 1) && (length(new_names) > 1)) {
      if (length(input$flevels) < length(new_names)) {
        NULL
      } else {
        dat_rv$dtab[, col_indx] <- mapvalues(dat_rv$dtab[, col_indx], from = c(input$flevels), to = c(new_names))
      }
    }

    dat_rv$dtab
  })

  observeEvent(input$accept_revalue1, {
    output$dat_levels_dt <- renderDT({
      dat_levels_rename()
      dat_rv$dtab
    })
  })
}

shinyApp(ui, server)
