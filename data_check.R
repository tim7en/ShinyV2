library(shiny)
library(DT)
library(plyr)
library(shinydashboard)
library(zCompositions)

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
  rownames (tabSummary) <- c("Type", "NAs", "Unique", "Negative", "Zeros")
  return (tabSummary)
}


# summary table
#dat <- dfSummary(x)
#dat <- as.data.frame(dat)
#rownames(dat) <- c("Type", "NAs", "Unique", "Negative", "Zeros")

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
  dashboardSidebar(),
  dashboardBody(
    sidebarLayout(
      sidebarPanel(
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
            uiOutput("dat_negatives"),
            uiOutput("dat_neg_text"),
            uiOutput("accept_negconv")
          ),
          tabPanel(
            "Zeros",
            uiOutput("dat_zeros"),
            uiOutput("dat_zero_text"),
            uiOutput("accept_zeroconv")
          ),
          tabPanel(
            "Missing",
            uiOutput("dat_missing"),
            uiOutput("dat_det_lim"),
            uiOutput("dat_facts_mis"),
            uiOutput("dat_missing_sub"),
            uiOutput("dat_missing_adj"),
            uiOutput("dat_miss_neg_sel"),
            uiOutput("imp_rb"),
            uiOutput("accept_missconv")
          )
        )
      ),
      mainPanel(
        tabsetPanel(
          tabPanel(
            "Data",
            box(
              width = 12,
              title = "Data table", solidHeader = TRUE, status = "primary",
              DTOutput("dat_levels_dt"),
              style = "height:'auto'; overflow-y: scroll;overflow-x: scroll;"
            ),
            uiOutput("dt_lim_ui")
          ),
          tabPanel(
            "Summary",
            box(
              title = "Summary table", status = "danger", height = "630", width = 12, solidHeader = T,
              DTOutput(outputId = "distTable"),
              style = "height:'auto'; overflow-y: scroll;overflow-x: scroll"
            )
          )
        )
      )
    )
  )
)

# Define server logic required to draw a histogram ----
server <- function(input, output) {
  dat_rv <- reactiveValues(dtab = x)

  output$distTable <- renderDT({
    data.frame(dfSummary (dat_rv$dtab))
  })
  
  dat_sum <- reactive ({
    dfSummary (dat_rv$dtab)
  })

  output$dat_facts <- renderUI({
    dat <- dat_sum ()
    selectInput("facts", "Factors", choices = colnames(dat)[which(as.numeric(as.matrix(dat[3, ])) > 0)], selected = NULL)
  })

  output$dat_facts_mis <- renderUI({
    dat <- dat_sum ()
    selectInput("facts2", "GROUPING", choices = colnames(dat)[which(as.numeric(as.matrix(dat[3, ])) > 0)], selected = NULL, multiple = T)
  })

  output$dat_negatives <- renderUI({
    dat <- dat_sum ()
    selectInput("negts", "Negatives", choices = colnames(dat)[which(as.numeric(as.matrix(dat[4, ])) > 0)], selected = NULL, multiple = TRUE)
  })

  output$dat_neg_text <- renderUI({
    textInput("neg_text", "Input equation")
  })

  output$dat_zeros <- renderUI({
    dat <- dat_sum ()
    selectInput("zeros", "ZEROS", choices = colnames(dat)[which(as.numeric(as.matrix(dat[5, ])) > 0)], selected = NULL, multiple = TRUE)
  })

  output$dat_missing <- renderUI({
    dat <- dat_sum ()
    selectInput("missing", "MISSING", choices = colnames(dat)[which(as.numeric(as.matrix(dat[2, ])) > 0)], selected = NULL, multiple = TRUE)
  })

  output$dat_missing_sub <- renderUI({
    col_indx <- which(names(dat_rv$dtab) == input$facts2)
    sel_names <- rownames_missing()[, col_indx]
    selectInput("missing_lvl", "LEVELS", choices = c(levels(dat_rv$dtab[, col_indx]), "JOIN"), multiple = TRUE, selected = sel_names) # one that are missing values
  })

  output$dat_missing_adj <- renderUI({
    selectInput("missing_adj", "FILL", choices = c("MEAN", "MEDIAN", "IMPUTE"), selected = NULL, multiple = FALSE)
  })

  # could be replaced to textOutput instead of selectInput
  output$dat_miss_neg_sel <- renderUI({
    req(input$missing_adj)
    dat <- dat_sum ()
    cl_remove <- unique(c(colnames(dat)[which(as.numeric(as.matrix(dat[4, ])) > 0)], colnames(dat)[which(as.numeric(as.matrix(dat[5, ])) > 0)]))
    if (input$missing_adj == "IMPUTE") {
      selectInput("missing_neg", "EXCLUDED FROM IMPUTATION (NEGATIVES, ZEROS)", choices = colnames(dat)[-c(1, 2)], selected = cl_remove, multiple = TRUE)
    } else {
      NULL
    }
  })

  output$dat_det_lim <- renderUI({
    req(input$missing_adj)
    req(input$impute_selected)

    if (input$impute_selected == "def") {
      if (input$missing_adj == "IMPUTE") {
        textInput("det_lim", "TYPE DETECTION LIMIT", value = '1')
      } else {
        NULL
      }
    } else { # This function is repsonsible for loading in the selected file
      fileInput("datafile", "Choose CSV file",
        accept = c("text/csv", "text/comma-separated-values,text/plain")
      )
    }
  })

  # This function is repsonsible for loading in the selected file
  filedata <- reactive({
    infile <- input$datafile
    if (is.null(infile)) {
      return(NULL)
    }
    read.csv(infile$datapath)
  })

  output$dt_lim <- renderDT({
    req(input$impute_selected)
    if (input$impute_selected == "def") {
      return(NULL)
    }
    else {
      filedata()
    }
  })

  output$dt_lim_ui <- renderUI({
    req(input$impute_selected != "def")
    box(
      width = 12,
      title = "DETECTION LIMIT", solidHeader = TRUE, status = "primary",
      DTOutput("dt_lim"),
      style = "height:'auto'; overflow-y: scroll;overflow-x: scroll;"
    )
  })

  output$dat_zero_text <- renderUI({
    textInput("zer_text", "Input equation")
  })

  output$dat_levels <- renderUI({
    col_indx <- which(names(dat_rv$dtab) == input$facts)
    selectInput("flevels", "Levels", choices = levels(dat_rv$dtab[, col_indx]), multiple = TRUE, selected = FALSE)
  })

  output$new_levels <- renderUI({
    textInput("new_names", "Names")
  })

  output$accept_revalue <- renderUI({
    actionButton("accept_revalue_bt", "Accept")
  })

  output$accept_negconv <- renderUI({
    actionButton("accept_negconv_bt", "Accept")
  })

  output$accept_zeroconv <- renderUI({
    actionButton("accept_zeroconv_bt", "Accept")
  })

  output$accept_missconv <- renderUI({
    actionButton("accept_missconv_bt", "Accept")
  })

  output$imp_rb <- renderUI({
    req(input$missing_adj)
    if (input$missing_adj == "IMPUTE") {
      radioButtons(
        "impute_selected", "Impute method :",
        c(
          "User type" = "def",
          "Upload matrix" = "sel_mat",
          "Upload table" = "sel_tab"
        )
      )
    } else {
      NULL
    }
  })

  # main function to rename factor levels
  dat_levels_rename <- eventReactive(input$accept_revalue_bt, {
    col_indx <- which(names(dat_rv$dtab) == input$facts) # find columns selected
    new_names <- try(unlist(strsplit(input$new_names, ","))) # unlist user input and use coma as separator
    if ((length(input$flevels) == 1) && (length(new_names) == 1)) { # if only one factor selected to rename and only one provided
      dat_rv$dtab[, col_indx] <- mapvalues(dat_rv$dtab[, col_indx], from = c(input$flevels), to = c(new_names))
    } else if ((length(input$flevels) > 1) && (length(new_names) == 1)) { # if more then one factor selected to rename and only one provided
      new_names <- rep(new_names, length(input$flevels))
      dat_rv$dtab[, col_indx] <- mapvalues(dat_rv$dtab[, col_indx], from = c(input$flevels), to = c(new_names))
    } else if ((length(input$flevels) > 1) && (length(new_names) > 1)) { # if more then one factor and more then one name provided
      if (length(input$flevels) != length(new_names)) { # if the length don't match
        NULL
      } else {
        dat_rv$dtab[, col_indx] <- mapvalues(dat_rv$dtab[, col_indx], from = c(input$flevels), to = c(new_names))
      }
    }
    dat_rv$dtab
  })

  # main function to convert negatives into positive
  dat_neg_conv <- eventReactive(input$accept_negconv_bt, {

  })

  # main function to accept or convert 0
  dat_zero_conv <- eventReactive(input$accept_zeroconv_bt, {

  })



  # main function to deal with missing data , impute nondetects using detection limits or other available options (mean, median)
  dat_miss_conv <- eventReactive(input$accept_missconv_bt, {
    col_ind <- which(names(dat_rv$dtab) %in% input$facts2)
    dats <- dat_rv$dtab[which(as.character(dat_rv$dtab[, col_ind]) %in% as.character(input$missing_lvl)), ]
    join_flag <- grep("JOIN", input$missing_lvl)

    d <- NULL
    if (length(join_flag) < 1) { # get the mean or median for a subset factors
      for (i in seq(1, length(input$missing_lvl))) {
        dats_sub <- dats[which(as.character(dats[, col_ind]) %in% as.character(input$missing_lvl[i])), ]
        if (nrow(dats_sub) == 1) {
          return(NULL)
        }
        for (j in seq(1, length(input$missing))) {
          if (input$missing_adj == "MEAN") {
            dats_sub[is.na(dats_sub[, input$missing[j]]), input$missing[j]] <- mean(dats_sub[, input$missing[j]], na.rm = TRUE)
          } else if (input$missing_adj == "MEDIAN") {
            dats_sub[is.na(dats_sub[, input$missing[j]]), input$missing[j]] <- median(dats_sub[, input$missing[j]], na.rm = TRUE)
          }
        }
        d <- rbind(d, dats_sub)
      }
      dat_rv$dtab[which(as.character(dat_rv$dtab[, col_ind]) %in% as.character(input$missing_lvl)), ] <- d
    
    } else { # if the join flag tagged, use all of them to get the mean or median from the entire data set
      if (nrow(dats) == length(input$missing_lvl) - 1) {
        return(NULL)
      }
      for (j in seq(1, length(input$missing))) {
        if (input$missing_adj == "MEAN") {
          dats[is.na(dats[, input$missing[j]]), input$missing[j]] <- mean(dats[, input$missing[j]], na.rm = TRUE)
        }
        if (input$missing_adj == "MEDIAN") {
          dats[is.na(dats[, input$missing[j]]), input$missing[j]] <- median(dats[, input$missing[j]], na.rm = TRUE)
        }
      }
      dat_rv$dtab[which(as.character(dat_rv$dtab[, col_ind]) %in% as.character(input$missing_lvl)), ] <- dats
    }
    
    
    ####CURRENTLY WORKING ON####
    if (input$missing_adj == 'IMPUTE' && input$impute_selected == 'def'){
      det_lim_vals <- as.numeric(unlist(strsplit(input$det_lim, ",")))
      dat <- (dat_rv$dtab)[,-c(1,2)] #remove first two columns (source and organic content)
      dat <- dat[,which(!names(dat) %in% input$missing_neg)]
      dl = rep (0, length(names(dat)))
      ind = which(names(dat) %in% input$missing)
      dl[ind] <- det_lim_vals
      dat[is.na(dat)]<-0
      dats <- lrEM(dat,label=0,dl=dl,ini.cov="multRepl")
      dat_rv$dtab[,c(input$missing)]<- dats[,c(input$missing)]
    } else if (input$missing_adj == 'IMPUTE' && input$impute_selected == 'sel_mat') {
      req (filedata())
      dat <- filedata ()
      

    } else if (input$missing_adj == 'IMPUTE' && input$impute_selected == 'sel_tab') {
        req (filedata())
        dats <- filedata ()
        det_lim_vals<-dats[,2][(match(input$missing, dats[,1]))]
        dat <- (dat_rv$dtab)[,-c(1,2)] #remove first two columns (source and organic content)
        dat <- dat[,which(!names(dat) %in% input$missing_neg)]
        dl = rep (0, length(names(dat)))
        ind = which(names(dat) %in% input$missing)
        dl[ind] <- det_lim_vals
        dat[is.na(dat)]<-0
        dats <- lrEM(dat,label=0,dl=dl,ini.cov="multRepl")
        dat_rv$dtab[,c(input$missing)]<- dats[,c(input$missing)]
      
      }
    # find column selected as missing
    # find levels selected as missing
    # based on the fill selected
    # for each level find a mean and fill in the data
    # for each level find a median and fill in the data
    # if all selected, use the mean from all of them
    # if 3 selected, subset to each of them
    #

    # subset
    # df[,c("A","B","E")]
  })

  observeEvent(input$accept_revalue_bt, {
    output$dat_levels_dt <- renderDT({
      dat_levels_rename()
      dat_rv$dtab
    })
  })

  observeEvent(input$accept_negconv_bt, {

  })

  observeEvent(input$accept_zeroconv_bt, {

  })

  observeEvent(input$accept_missconv_bt, {
    dat_miss_conv()
  })

  rownames_missing <- reactive({
    dats <- dat_rv$dtab[, c(input$missing)]
    dat_rv$dtab[as.vector(attributes(na.omit(dats))$na.action), ]
  })
}

shinyApp(ui, server)
