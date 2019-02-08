
#* ------------------------------------------------- Outliers Removal -----
#  |  Function: Type ShinyDashboard
#
#  |  Purpose:  Current code will run a shiny daashboard UI and perform outlier
#                           removal of source data for sediments fingerprinting.
#               Requires csv type input file.
#   | Author: Timur Sabitov
#|
#* -------------------------------------------------------------------*/
rm(list = ls())
source("wd.R")
source("pcg.R")
source("func.R")
source("Regressions.R")
source("function_cor.R")
source("modules_navtab1.R")
source("modules_navtab2.R")
# source ('correction_test.R')
source("input_mod.R")
library("shiny")
library("DT")
# source ('correction_test.R')


server <- function(input, output, session) {
  options(shiny.maxRequestSize = 70 * 1024^2) # Max csv data limit set to 60 mb

  src.Data <- callModule(csvFile, "file1", stringsAsFactors = FALSE)
  dats <- callModule(inputMod, "dat1", src.Data())

  trg.Data <- callModule(csvFile, "file2")
  trgs <- callModule(inputMod, "dat2", trg.Data())

  output$DT <- renderDT({
    DTout()
  })

  DTout <- reactive({
    dats()
  })

  output$TD <- renderDT({
    TTout()
  })

  output$TD_p2 <- renderDT({
    TTout()
  })

  TTout <- reactive({
    trgs()[, which(names(trgs()) %in% names(DT_p2()))]
  })

  vars <- callModule(checkD, "dat1", DTout)

  DT_p2 <- reactive({
    vars()
  })

  output$DT_p2 <- renderDT({
    DT_p2()
  })

  output$selcCor <- renderUI({
    selectInput("slcC", "Select columns with size, toc or together", choices = names(DT_p2())[-c(1, 2)], multiple = TRUE, width = "100%")
  })

  selcRem <- reactive({
    input$slcC
  })

  # option to pick columns that will not be used for correction
  output$selcRem <- renderUI({
    selectInput("slcR", "Select columns to remove from correction", choices = names(DT_p2())[-c(1, 2)], multiple = TRUE, width = "100%")
  })

  # select p value threshold of normality for residuals
  output$shapiroP <- renderUI({
    sliderInput("shapiroP", "Shapiro-Wilk Univariate Normality Test p-value:", value = 0.05, min = 0.001, max = 1, step = 0.01)
  })

  # sliderinput for corrplot R threshold
  output$corR <- renderUI({
    sliderInput("corR", "R", value = 0.6, min = 0, max = 0.99, step = 0.1, animate = F)
  })

  selcCor <- reactive({
    input$slcR
  })

  output$applyCor <- renderUI({
    actionButton("applyCor", "Apply corrections")
  })

  backframeDt <- eventReactive(input$applyCor, {
    datas <- DT_p2()
    datas <- corect.func(as.data.frame(datas), input$corR, input$shapiroP, input$slcC, input$slcR)
    datas <- datas[, -which(names(datas) %in% c("formula"))]
    datas$id <- seq(1, nrow(datas))
    if (!is.null(datas)) {
      datas
    } else {

    }
  })

  convOptions <- reactive({
    datas <- backframeDt()
  })

  output$convOptions <- renderDT({
    datas <- convOptions()
    colnames(datas)[5] <- "Formula"
    datas[, 5] <- (gsub("I(datas$", replacement = "(", datas[, 5], fixed = T))
    datas[, 5] <- (gsub("datas$", replacement = "", datas[, 5], fixed = T))
    datas[, 5] <- (gsub("I", replacement = "", datas[, 5], fixed = T))
    datas[, 5] <- (gsub("logI", replacement = "log", datas[, 5], fixed = T))
    # datas$id <- seq (1, nrow(datas))
    if (!is.null(datas)) {
      datas
    } else {

    }
  })

  fDt <- reactive({
    dat <- convOptions()
    valTable <- dat
    varr <- valTable %>% group_by(source, element)
    f <- varr %>% arrange(-desc(residualsSD), desc(Cooks), .by_group = TRUE)
    f <- f %>%
      group_by(source, element) %>%
      mutate(rank = rank(-desc(residualsSD), ties.method = "first"))
    f <- f[, -which(names(f) %in% c("formula", "grade", "Cooks", "residualsSD", "p-eq", "p-resi", "p.eq", "p.resi"))]
    f <- as.data.frame(f)
    f[, 3] <- (gsub("I(datas$", replacement = "(", f[, 3], fixed = T))
    f[, 3] <- (gsub("datas$", replacement = "", f[, 3], fixed = T))
    f[, 3] <- (gsub("I", replacement = "", f[, 3], fixed = T))
    f[, 3] <- (gsub("logI", replacement = "log", f[, 3], fixed = T))
    f[, 1] <- as.factor(f[, 1])
    f[, 3] <- as.factor(f[, 3])
    f <- f[!duplicated(f[, 3]), ]
    return(f)
  })

  convFDt <- reactive({
    dat1 <- convOptions()
    dat2 <- finalDtrg()
    dat <- dat1[which(dat1$id %in% dat2$id), ]
    dat <- cbind(dat[, 1], 0, dat[, 2:ncol(dat)])
    dat
  })

  output$convFDt <- renderDT({
    convFDt()
  })

  outputCorrected <- reactive({
    req(DT_p2()) # source data
    req(convFDt()) # formulas format
    req(TTout()) # Target data
    req(input$trgS)

    datas <- DT_p2()
    datas <- datas[, which(!colnames(datas) %in% input$slcR)]
    target <- as.data.frame(TTout())
    target <- target[, which(!colnames(target) %in% input$slcR)]
    y <- as.data.frame(convFDt())
    x <- as.data.frame(datas)
    output <- list()
    slopes <- list()
    drops <- list()

    # print (x) #source
    # print (target) #target
    # print (input$slcC) #formula

    for (i in seq(1, dim(target)[1])) {
      # print (x)
      # print (target[i,])
      # print (y)
      # print (input$slcC)
      dat <- correct(x, target[i, ], y, input$slcC)
      # print (dat)
      # print (dat)
      output[[i]] <- dat[[1]]
      slopes[[i]] <- dat[[2]]
      if (length(dat[[3]]) > 0) {
        drops [[i]] <- dat[[3]]
      } else {
        drops[[i]] <- NA
      }
    }


    y <- slopes
    slopes.DT <- data.frame(matrix(NA, nrow = 8, ncol = 3))

    for (i in seq(1, length(y[[1]]))) {
      if (length(y[[1]][[i]]) < 2) {
        y[[1]][[i]] <- c(y[[1]][[i]], NA)
        y[[1]][[i]] <- c(y[[1]][[i]], NA)
      } else if (length(y[[1]][[i]]) < 3) {
        y[[1]][[i]] <- c(y[[1]][[i]], NA)
      }
      slopes.DT[i, ] <- y[[1]][[i]]
    }

    colnames(slopes.DT) <- c("var1", "var2", "var1*var2")

    cleanT <- function(x) {
      x <- (gsub("I(datas$", replacement = "(", x, fixed = T))
      x <- (gsub("datas$", replacement = "", x, fixed = T))
      x <- (gsub("I", replacement = "", x, fixed = T))
      x <- (gsub("logI", replacement = "log", x, fixed = T))
    }
    drops <- lapply(drops, cleanT)

    drops <- plyr::ldply(drops, rbind)

    list(output, slopes.DT, drops)
  })

  output$trgS <- renderUI({
    selectInput("trgS", "Target ID", choices = seq(1, nrow(TTout())), selected = 1)
  })

  tabList <- reactive({
    outputCorrected() [[1]][[as.numeric(input$trgS)]]
  })
  
  output$tabList <- renderDT ({
    tabList ()
  })

  output$tabLslopes <- renderDT({
    outputCorrected() [[2]]
  })

  output$tabLdrops <- renderDT({
    outputCorrected() [[3]]
  })
  
  output$fDt <- renderDT({
    if (!is.null(fDt())) {
      dat <- fDt ()[,c(3,2,1,6,4,5)]
      datatable(dat, filter = "top", options = list(
        pageLength = 5, autoWidth = F
      ))
    } else {}
  })

  selectPick <- reactive({
    fDt()[input$fDt_rows_selected, ]
  })

  # observe rows selected
  output$selDat <- renderUI({
    radioButtons(
      "selDat", "Select data :",
      c(
        "Top rank" = "def",
        "User choice" = "sel"
      )
    )
  })

  output$finalDtrg <- renderDT({
    datatable(finalDtrg(), filter = "top", options = list(
      pageLength = 5, autoWidth = F
    ))
  })

  finalDtrg <- reactive({
    req (input$selDat)
    if (input$selDat == "def") {
      datas <- topPick()
    } else {
      datas <- selectPick()
    }
    datas
  })


  output$topPick <- renderDT({
    datatable(topPick(), filter = "top", options = list(
      pageLength = 5, autoWidth = F
    ))
  })

  topPick <- reactive({
    dat <- fDt()
    datas <- dat[which(dat$rank == 1), ]
  })

  output$selVarplot <- renderUI({
    # print (fDt()[which(fDt()$rank == 1),])
    if (!is.null(input$selDat)) {
      if (input$selDat == "def") {
        selectInput("selVarplot", label = "Select Element", choices = fDt()[which(fDt()$rank == 1), 3])
      } else {
        selectInput("selVarplot", label = "Select Element", choices = fDt()[input$fDt_rows_selected, 3])
      }
    } else {}
  })
  # callModule (corrAnalysis, 'dat1', DT_p2(),TTout(), convOptions ())

  # regressions plot function
  regPlotfunc <- reactive({
    req(topPick())
    # req (input$selVarplot)
    datas <- DT_p2()

    # print (input$selDat)

    if (!is.null(input$selDat)) {
      if (input$selDat == "def") {
        if (!is.null(input$selVarplot)) {
          sclass <- topPick()[which(topPick()[, 3] == input$selVarplot), 4]
          sclass <- as.character(sclass)
          datas <- datas[which(datas[, 2] == sclass), ]
        } else {}
      } else {
        if (!is.null(input$selVarplot)) {
          sclass <- selectPick()[which(selectPick()[, 3] == input$selVarplot), 4]
          sclass <- as.character(sclass)
          datas <- datas[which(datas[, 2] == sclass), ]
          datas
        } else {}
      }
    }

    # print (datas)
    datas
  })

  # regressions plot output
  output$regPlot <- renderPlot({
    req(input$selDat)
    req(input$selVarplot)
    req(selectPick())
    req(input$slcC)

    dats <- NULL
    dats <- selectPick()
    slcC <- input$slcC

    if (dim(dats)[1] < 1 & input$selDat == "sel") {} else {
      datas <- regPlotfunc()

      if (length(slcC) > 1) {
        var1 <- as.numeric(datas[, slcC[1]])
        var2 <- as.numeric(datas[, slcC[2]])
      } else {
        var1 <- as.numeric(datas[, slcC])
        var2 <- NULL
      }

      datas[, 3:dim(datas)[2]] <- apply(datas[, 3:dim(datas)[2]], 2, as.numeric)

      f <- input$selVarplot

      f <- sub("^.", "I(datas$", f)
      f <- sub("(var1", "I(var1", f, fixed = T)
      f <- sub("(var2", "I(var2", f, fixed = T)
      f <- sub("logI", "log", f, fixed = T)
      f <- sub("I(datas$og(", "log(datas$", f, fixed = T)
      fit <- lm(eval(parse(text = f)))
      par(mfrow = c(1, 3))


      cooksd <- cooks.distance(fit)
      plot(cooksd, cex = 2, main = "Influential Obs by Cooks distance")
      lines(lowess(cooksd), col = "blue")

      plot(fitted(fit), residuals(fit))
      lines(lowess(fitted(fit), residuals(fit)))
      title("Residual vs Fit. value ", col = "red")

      acf(residuals(fit), main = "")
      title("Residual Autocorrelation Plot")
    }
  })
  
  output$xvar <- renderUI ({
    selectInput ('xvar', 'Observed, Xvar', choices = names (tabList ())[-c(1,2)])
  })
  
  output$yvar <- renderUI ({
    selectInput ('yvar', 'Corrected, Yvar', choices = names (tabList ())[-c(1,2)])
  })
  
  output$scatterplot1 <- renderScatterD3({
    req (input$xvar)
    req (input$yvar)
    #scatterD3(data = data, x=mpg, y=carb,
    mtdf <- tabList ()
    dtorg <- DT_p2 ()
    
    x <- dtorg[[input$xvar]]
    Sclass <- dtorg [,2]
    y <- mtdf[[input$yvar]]
    
    scatterD3(x=x,y=y,
              labels_size= 9, point_opacity = 1, col_var = Sclass, 
              lines = data.frame(slope = 1, intercept = Inf),
              #col_var=cyl, symbol_var= data$Assay,
              #lab= paste(mpg, carb, sep="|") , lasso=TRUE,
              #xlab= "IFN-γ", ylab= "IL-10",
              #click_callback = "function(id, index) {
              #  alert('scatterplot ID: ' + id + ' - Point index: ' + index) 
              #  }", 
              transitions= T)
  })
  
  output$brack_range <- renderUI ({
    sliderInput ('brack_rng', 'Bracket range parameter', value = 0.1, min = 0, max = 3, step = 0.05)
  })
  
  tarBrackets <- reactive({
    req (outputCorrected ())
    req (TTout())
    req (input$brack_rng)
    x <- outputCorrected ()[[1]] #list of sources corrected
    datas <- DT_p2()
    datas <- datas[, which(!names(datas) %in% names(x[[1]]))]
    y <- trgs() #list of targets corrected
    
    dat <- NULL
    trg <- NULL
    for (i in seq (1,nrow(y))){
      x[[i]]<- cbind (x[[i]], datas)
      l <- bracketT (x[[i]],y[i,], input$brack_rng, input$slcR)
      dat <- rbind (dat, l[[1]])
      trg <- c (trg, l[[2]])
    }
    
    levels(dat[,1]) <- as.character (dat[,1])
    
    return (list(dat, trg))
    
  })
  
  output$trgBrkts <- renderDT ({
    req (tarBrackets ())
    dat <- tarBrackets ()
    selection <- dat[[2]]
    trgDat <<- dat
    
    DT::datatable(dat[[1]]) %>% formatStyle(
      c(colnames(dat[[1]])),
      backgroundColor = styleEqual(selection, rep("lightsalmon", length(selection)))
    )
  })
  
  
}

# User interface side of the user input
ui <- dashboardPage(
  skin = "blue",
  dashboardHeader(title = "SedSat_ShinyV2"),
  ## Sidebar content
  dashboardSidebar(
    sidebarMenu(
      menuItem("Input", tabName = "dataInput", icon = icon("upload")),
      menuItem("Size & TOC Correction", tabName = "regressions", icon = icon("random")),
      menuItem("Discriminant Function Analysis", tabName = "DFA", icon = icon("table")),
      menuItem("Mixing Model", tabName = "mixmod", icon = icon("cubes")) # ,
      # menuItem("ML Model", tabName = "mlmod", icon = icon("circle"))
    )
  ),
  dashboardBody( # Body content
    tabItems(
      tabItem( # First tab content
        tabName = "dataInput",
        tabsetPanel(
          tabPanel(
            "Source",
            fluidPage(
              fluidRow(
                sidebarLayout(
                  sidebarPanel(
                    fluidRow(
                      column(
                        width = 12,
                        csvFileInput("file1", "User data (.csv format)"),
                        columnChooserUI("dat1")
                      )
                    )
                  ),
                  mainPanel(
                    fluidRow(
                      tabsetPanel(
                        tabPanel(
                          "Data",
                          column(
                            width = 12,
                            DTOutput("DT"),
                            style = "height:'auto'; overflow-y: scroll;overflow-x: scroll;"
                          )
                        ),
                        tabPanel(
                          "Plots",
                          fluidRow(
                            column(
                              width = 9,
                              br(),
                              srcCor("dat1")
                            ),
                            column(
                              width = 3,
                              mat_par("dat1")
                            ),
                            column(
                              width = 12,
                              br(),
                              srcDCor("dat1")
                            )
                          )
                        )
                      )
                    )
                  )
                )
              ),
              tabsetPanel(
                tabPanel(
                  "Original Data p-values",
                  box(
                    title = textOutput("origTitle"), status = "success", height = "630", width = 12, solidHeader = T,
                    column(
                      width = 12,
                      srcSP("dat1"),
                      style = "height: 550px; overflow-y: scroll; overflow-x: scroll;"
                    )
                  )
                ),
                tabPanel(
                  "Transformations Advanced",
                  tabsetPanel(
                    tabPanel(
                      "Methods",
                      box(
                        title = "Methods p-values ", status = "danger", height = "630", width = 12, solidHeader = T,
                        column(
                          width = 12,
                          getspMethods("dat1"),
                          style = "height:550px; overflow-y: scroll;overflow-x: scroll;"
                        )
                      )
                    ),
                    tabPanel(
                      "Transformed p-values",
                      box(
                        title = "Methods p-values ", status = "danger", height = "630", width = 12, solidHeader = T,
                        column(
                          width = 12,
                          getspPval("dat1"),
                          style = "height: 550px; overflow-y: scroll; overflow-x: scroll;"
                        )
                      )
                    ),
                    tabPanel(
                      "QQ Plots",
                      box(
                        title = "QQ plot of Original Data", status = "success", height = "1250", width = 6, solidHeader = T,
                        column(
                          width = 10,
                          style = "height:100px;",
                          spPlotpick("dat1")
                        ),
                        getorigQQval("dat1"), style = "height:1200px;overflow-y: scroll;overflow-x: scroll;",
                        column(
                          width = 6
                        )
                      ),
                      box(
                        title = "QQ plot of Transformed Data", status = "success", height = "1250", width = 6, solidHeader = T,
                        column(
                          width = 10,
                          style = "height:100px;",
                          shapiroP("dat1")
                        ),
                        getspQQval("dat1"), style = "height:1200px;overflow-y: scroll;overflow-x: scroll;", # 570
                        column(
                          width = 6
                        )
                      )
                    )
                  )
                ),
                tabPanel(
                  "Outliers",
                  fluidPage(
                    fluidRow(
                      box(
                        title = "Data: Original ", status = "success", height =
                          "595", width = "12", solidHeader = T,
                        column(
                          width = 12,
                          downloadButton("downloadData", "Download"),
                          br(),
                          tags$hr(),
                          outliersTab1("dat1"), style = "height:500px; overflow-y: scroll;overflow-x: scroll;"
                        )
                      )
                    )
                  )
                ),
                tabPanel(
                  "Outliers Advanced",
                  tabsetPanel(
                    tabPanel(
                      "Transformed Data & Outliers",
                      box(
                        title = "Data: Transformed ", status = "warning", height =
                          "595", width = 12, solidHeader = T,
                        column(
                          width = 12,
                          outliersADtab1("dat1"), style = "height:500px; overflow-y: scroll;overflow-x: scroll;"
                        )
                      )
                    ),
                    tabPanel(
                      "Edit Table To Keep/Remove Outliers",
                      box(
                        # tags$style(HTML(".js-irs-1 .irs-single, .js-irs-1 .irs-bar-edge, .js-irs-1 .irs-bar {background: green;}")),
                        # tags$style(HTML(".js-irs-1 .irs-single, .js-irs-1 .irs-bar-edge, .js-irs-1 .irs-grid-text { font-size: 10pt;}")),
                        nStd("dat1"),
                        title = "Select Or Deselect Rows ", status = "success", height =
                          "695", width = 12, solidHeader = T,
                        column(
                          width = 12,
                          outliersADtab2("dat1"), style = "height:500px; overflow-y: scroll;overflow-x: scroll;"
                        )
                      ),
                      box(
                        title = "Standard Normal Deviate", status = "success", height = "600", width = 12, solideHeader = T,
                        column(
                          width = 12,
                          stdTab1("dat1"), style = "height:500px; overflow-y: scroll;overflow-x: scroll;"
                        )
                      )
                    ),
                    tabPanel(
                      "Outliers Removed",
                      box(
                        title = "These Rows Will Be Excluded From Final Output ", status = "danger", height =
                          "595", width = 12, solidHeader = T,
                        column(
                          width = 12,
                          outliersADtab3("dat1"), style = "height:500px; overflow-y: scroll;overflow-x: scroll;"
                        )
                      )
                    ),
                    tabPanel(
                      "Final Output",
                      box(
                        title = "Final Output Table", status = "success", height =
                          "595", width = 12, solidHeader = T,
                        column(
                          width = 12,
                          br(),
                          downloadButton("downloadData2", "Download"),
                          tags$hr(),
                          outliersADtab4("dat1"), style = "height:500px; overflow-y: scroll;overflow-x: scroll;"
                        )
                      )
                    )
                  )
                )
              )
            )
          ),
          tabPanel(
            "Target",
            fluidPage(
              fluidRow(
                sidebarLayout(
                  sidebarPanel(
                    fluidRow(
                      column(
                        width = 12,
                        csvFileInput("file2", "User data (.csv format)"),
                        columnChooserUI("dat2")
                      )
                    )
                  ),
                  mainPanel(
                    fluidRow(
                      tabsetPanel(
                        tabPanel(
                          "Data",
                          column(
                            width = 12,
                            DTOutput("TD"),
                            style = "height:'auto'; overflow-y: scroll;overflow-x: scroll;"
                          )
                        ),
                        tabPanel(
                          "Plots",
                          fluidRow(
                            column(
                              width = 9,
                              br(),
                              srcCor("dat2")
                            ),
                            column(
                              width = 3,
                              mat_par("dat2")
                            )
                          ),
                          column(
                            width = 12,
                            br(),
                            srcDCor("dat2")
                          )
                        )
                      )
                    )
                  )
                )
              )
            )
          )
        )
      ),
      tabItem( # First tab content
        tabName = "regressions",
        tabsetPanel(
          tabPanel(
            "Sources",
            column(
              width = 12,
              DTOutput("DT_p2"), style = "height:'auto'; overflow-y: scroll;overflow-x: scroll;"
            )
          ),
          tabPanel(
            "Target",
            column(
              width = 12,
              DTOutput("TD_p2"), style = "height:'auto'; overflow-y: scroll;overflow-x: scroll;"
            )
          ),
          tabPanel(
            "Corrections",
            sidebarLayout(
              sidebarPanel(
                
                uiOutput("selcCor"), # ,
                uiOutput("selcRem"),
                uiOutput("shapiroP"),
                uiOutput("corR"),
                uiOutput("applyCor")
              ),
              mainPanel(
                box(
                  title = "Available options", status = "success", height =
                    "auto", solidHeader = T, width = "auto",
                  withSpinner(DTOutput("convOptions")), style = "height:'auto'; overflow-y: scroll;overflow-x: scroll;"
                )
              )
              
            )
          ),
          tabPanel(
            "Target corrected",
            fluidRow(
              column(
                width = 12,
                box(
                  title = "Formulas", status = "success", height =
                    "auto", solidHeader = T,
                  DTOutput("fDt"), style = "height:'auto'; overflow-y: scroll;overflow-x: scroll;"
                ),
                box(
                  title = "Selected", status = "primary", height =
                    "auto", solidHeader = T,
                  DTOutput("finalDtrg"), style = "height:'auto'; overflow-y: scroll;overflow-x: scroll;"
                )
              ),
              column(
                width = 12,
                box(
                  title = "Top pick", status = "success", height =
                    "auto", solidHeader = T,
                  DTOutput("topPick"), style = "height:'auto'; overflow-y: scroll;overflow-x: scroll;"
                ),
                box(
                  title = "Plots", status = "success", height =
                    "auto", solidHeader = T,
                  uiOutput("selDat"),
                  uiOutput("selVarplot"),
                  plotOutput("regPlot")
                )
              )
            )
          ),
          tabPanel(
            "Corrected Data",
            tabsetPanel(
              tabPanel (
                'Data',
                fluidRow (
                  column(
                    width = 12,
                      uiOutput("trgS"),
                      DTOutput("tabList"), style = "height:'auto'; overflow-y: scroll;overflow-x: scroll;"
                    ),
                      column (
                        width = 6,
                        DTOutput("tabLslopes"), style = "height:'auto'; overflow-y: scroll;overflow-x: scroll;"
                      ),
                      column (
                        width = 6,
                        DTOutput("tabLdrops"), style = "height:'auto'; overflow-y: scroll;overflow-x: scroll;"
                    )
                  
                )
              ),
              tabPanel (
                'Plot',
                uiOutput ('xvar'),
                uiOutput ('yvar'),
                scatterD3Output("scatterplot1")
              )
            )
          ),
          tabPanel(
            "Bracket test",
            uiOutput ('brack_range'),
            DTOutput ('trgBrkts'), style = "height:'auto'; overflow-y: scroll;overflow-x: scroll;"
          )
        )
      )
    ),
    tabItems(
      tabItem( # First tab content
        tabName = "DFA",
        tabsetPanel()
      )
    ),
    tabItems(
      tabItem( # First tab content
        tabName = "mixmod",
        tabsetPanel()
      )
    )
  )
)



# Run the app ----
shinyApp(ui, server)
