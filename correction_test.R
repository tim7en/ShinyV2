library("shiny")
library("dplyr")
library("DT")
#source ('func.R')
# x <- read.csv("datasGlobal.csv")
# x <- x[, -1]
# x <- x[!duplicated(x[, 6]), ]
# x$id <- seq(1, nrow(x))
# y <- read.csv("targetData.csv")
# v <- read.csv("sourceData.csv")
# v <- na.omit(v)

outputCorrected <- function (id) {
  ns <- NS(id)
  DTOutput (ns('outputCorrected'))
}

selectClass <- function(id) {
  ns <- NS(id)
  uiOutput("selectClass")
}

selectVar <- function(id) {
  ns <- NS(id)
  uiOutput("selectVar")
}

fDt <- function(id) {
  ns <- NS(id)
  DTOutput(ns("fDt"))
}

topPick <- function(id) {
  ns <- NS(id)
  DTOutput(ns("topPick"))
}

selectPick <- function(id) {
  ns <- NS(id)
  DTOutput(ns("selectPick"))
}

selDat <- function(id) {
  ns <- NS(id)
  uiOutput(ns("selDat"))
}

selVarplot <- function(id) {
  ns <- NS(id)
  uiOutput(ns("selVarplot"))
}

regPlot <- function(id) {
  ns <- NS(id)
  plotOutput(ns("regPlot"))
}

outputCorrected <- function (id) {
  ns <- NS(id)
  DTOutput (ns('outputCorrected'))
}

corDrops <- function (id) {
  ns <- NS(id)
  DTOutput(ns('corDrops'))
}

corSlopes <- function (id) {
  ns <- NS(id)
  DTOutput(ns('corSlopes'))
}

# module server side function
corrAnalysis <- function(input, output, session, x, y,dat) {
  
  valTable <- dat
  varr <- valTable %>% group_by(source, element)
  f <- varr %>% arrange(-desc(residualsSD), desc(Cooks), .by_group = TRUE)
  f <- f %>%
    group_by(source, element) %>%
    mutate(rank = rank(-desc(residualsSD), ties.method = "first"))
  f <- f[, -which(names(f) %in% c("formula", "grade", "Cooks", "residualsSD", "p.eq", "p.resi"))]
  colnames(f)[3] <- "formula"
  f <- as.data.frame(f)
  f[, 3] <- (gsub("I(datas$", replacement = "(", f[, 3], fixed = T))
  f[, 3] <- (gsub("datas$", replacement = "", f[, 3], fixed = T))
  f[, 3] <- (gsub("I", replacement = "", f[, 3], fixed = T))
  f[, 3] <- (gsub("logI", replacement = "log", f[, 3], fixed = T))
  f[, 3] <- as.factor(f[, 3])


  output$fDt <- renderDT({
    datatable(f, filter = "top", options = list(
      pageLength = 10, autoWidth = TRUE
    ))
  })

  output$topPick <- renderDT({
    datatable(topPick(), filter = "top", options = list(
      pageLength = 5, autoWidth = F
    ))
  })

  topPick <- reactive({
    datas <- f[which(f$rank == 1), ]
  })

  output$selectPick <- renderDT({
    datatable(selectPick(), filter = "top", options = list(
      pageLength = 5, autoWidth = F
    ))
  })

  selectPick <- reactive({
    f[input$"fDt_rows_selected", ]
  })

  # observe rows selected
  an_observe_func <- observe(suspended = T, {
    input$"fDt_rows_selected"
    isolate({
      # do stuff here
      # print(input$"fDt_rows_selected")
    })
  })

  # start the observer, without "suspended=T" the observer
  #  will start on init instead of when needed
  an_observe_func$resume()

  output$selDat <- renderUI({
    ns <- session$ns
    radioButtons(
      ns("selDat"), "Select data :",
      c(
        "Top rank" = "def",
        "User choice" = "sel"
      )
    )
  })

  output$selVarplot <- renderUI({
    ns <- session$ns
    if (!is.null(input$selDat)) {
      if (input$selDat == "def") {
        selectInput(ns("selVarplot"), label = "Select Element", choices = f[which(f$rank == 1), 3])
      } else {
        selectInput(ns("selVarplot"), label = "Select Element", choices = f[input$"fDt_rows_selected", 3])
      }
    } else {}
  })

  # regressions plot function
  regPlotfunc <- reactive({
    req(topPick())
    # req (input$selVarplot)
    datas <- v

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

    datas
  })

  # regressions plot output
  output$regPlot <- renderPlot({
    req(input$selDat)
    req(input$selVarplot)
    req(regPlotfunc())
    req(slcC)

    if (dim(selectPick())[1] < 1 & input$selDat == "sel") {} else {
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
  
  
  #
  finalDtrg <- reactive ({
    if (input$selDat == "def"){
      datas <- topPick()
    } else {
      datas <- selectPick()
    }
    #print (datas)
    datas
  })
  
  output$finalDtrg <- renderDT ({
    finalDtrg ()
  })
  
  srcDT <- reactive ({
    x
  })
  
  trgDT <- reactive ({
    y
  })
  
  outputCorrected <- reactive({
    req(srcDT())
    req(finalDtrg())
    req(trgDT())
    
    datas <- srcDT()
    datas <-  datas[, which(!colnames(datas) %in% slcR)]
    target <- as.data.frame(trgDT())
    target <- target[,which(!colnames(target)%in% slcR)]
    
    y <- as.data.frame(valTable[which(valTable[,ncol(valTable)] %in% finalDtrg()[,(ncol(finalDtrg())-1)]),])
    x <- as.data.frame(datas)
    output <- list ()
    slopes <- list ()
    drops <- list ()
    
    y[,6] <- as.character(y[,6])
    
    for (i in seq(1, dim(target)[1])) {
      
      dat <- correct(x, target[i, ], y, slcC)
      #print (dat)
      output[[i]] <- dat[[1]]
      slopes[[i]] <- dat[[2]]
      if (length(dat[[3]])>0){
        drops [[i]] <- dat[[3]]
      } else {
        drops[[i]] <-NA
      }
    }
    
    
    y <- slopes
    slopes.DT <- data.frame(matrix(NA, nrow = 8, ncol = 3))
    
    for (i in seq (1,length(y[[1]]))){
      if (length(y[[1]][[i]])<2 ){
        y[[1]][[i]] <- c(y[[1]][[i]], NA)
        y[[1]][[i]] <- c(y[[1]][[i]], NA)
      } else if (length(y[[1]][[i]])<3 ){
        y[[1]][[i]] <- c(y[[1]][[i]], NA)
      }
      slopes.DT[i,] <- y[[1]][[i]]
    }
    
    colnames(slopes.DT) <- c('var1','var2','var1*var2')
    
    cleanT <- function (x) {
      x <- (gsub("I(datas$", replacement = "(", x, fixed = T))
      x <- (gsub("datas$", replacement = "", x, fixed = T))
      x <- (gsub("I", replacement = "", x, fixed = T))
      x <- (gsub("logI", replacement = "log", x, fixed = T))
    }
    drops <- lapply (drops, cleanT)
    
    drops <- plyr::ldply(drops, rbind)
    
    list(output, slopes.DT, drops)
  })
  
  output$outputCorrected <- renderDT({
    req(outputCorrected())
    targetD <- as.data.frame(trgDT())
    indCor <- which(targetD[, 1] == 'target1')
    outputCorrected()[[1]][[indCor]]
  })
  
  output$corSlopes <- renderDT ({
    req(outputCorrected())
    outputCorrected()[[2]]
  })
  
  output$corDrops <- renderDT ({
    req(outputCorrected())
    outputCorrected()[[3]]
  })
  
  return (outputCorrected)
}
# 
# ui <- tabsetPanel(
#   tabPanel(
#     "Data plot",
#     fluidPage(
#       fluidRow(
#         column(
#           width = 6,
#           fDt("dat1"), style = "height:'auto'; overflow-y: scroll;overflow-x: scroll;"
#         ),
#         column(
#           width = 6,
#           topPick("dat1"), style = "height:'auto'; overflow-y: scroll;overflow-x: scroll;",
#           br(),
#           br(),
#           selectPick("dat1"), style = "height:'auto'; overflow-y: scroll;overflow-x: scroll;"
#         )
#       ),
#       sidebarLayout(
#         sidebarPanel(
#           selDat("dat1"),
#           selVarplot("dat1")
#         ),
#         mainPanel(
#           regPlot("dat1")
#         )
#       )
#     )
#   ),
#   tabPanel(
#     "Model stats & conversion",
#     fluidPage(
#       fluidRow(
#         column(
#           width = 12,
#           outputCorrected ('dat1'),
#           corSlopes ('dat1'),
#           corDrops ('dat1')
#         )
#       )
#     )
#   )
# )
# 
# 
# server <- function(input, output, session) {
#   mod1 <- callModule(correctSrc, "dat1", v, y, x, c("Size", "OrganicContent"),c('Delta15N', 'Delta13COrganic') )
# }
# 
# shinyApp(ui, server)
