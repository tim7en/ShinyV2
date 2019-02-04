#module ui side function

#option to choose columns from dropdown names to correct for
selcRem <- function (id) {
  ns <- NS(id)
  uiOutput (ns('selcRem'))
}

#option to choose columns from dropdown names to remove from correction
selcCor <- function (id) {
  ns <- NS (id)
  uiOutput (ns('selcCor'))
}

#DT output, src
srcDT <- function (id) {
  ns <- NS (id)
  DTOutput(ns('srcDT'))
}

#DT output, trg
trgDT <- function (id) {
  ns <- NS (id)
  DTOutput(ns('trgDT'))
}

#apply corrections button
ApplyCor <- function (id) {
  ns <- NS (id)
  uiOutput (ns('ApplyCor'))
}

#parameters for user to choose for correlation relationship and residuals p value
corPar<- function(id) {
  ns <- NS(id)
  tagList(
    fluidRow (
      column (
        width = 6,
        br (),
        uiOutput(ns("corR"))
      ), 
      column (
        width = 6,
        br (),
        uiOutput(ns("shapiroP")
      )
    )
  )
)}

listDTcor <- function (id) {
  ns <- NS (id)
  DTOutput(ns('listDTcor'))
}

srcClass <- function (id) {
  ns <- NS (id)
  uiOutput (ns('srcClass'))
}

srcVar <- function (id) {
  ns <- NS(id)
  uiOutput (ns('srcVar'))
}

eqDatoutput <- function (id) {
  ns <- NS(id)
  DTOutput(ns('eqDatoutput'))
}

selEq <- function (id){
  ns <- NS(id)
  uiOutput(ns('selEq'))
}

defaultBut <- function (id) {
  ns <- NS(id)
  uiOutput(ns('defaultBut'))
}

resetBut <- function (id) {
  ns <- NS(id)
  uiOutput(ns('resetBut'))
}

selectBut <- function (id) {
  ns <- NS(id)
  uiOutput(ns('selectBut'))
}

dtDefout <- function(id) {
  ns <- NS(id)
  DTOutput(ns('dtDefout'))
}

eqPick <- function(id) {
  ns <- NS(id)
  DTOutput(ns('eqPick'))
}

regPlot <- function (id) {
  ns <- NS(id)
  plotOutput(ns('regPlot'))
}

selTrg <- function (id) {
  ns <- NS(id)
  uiOutput (ns('selTrg'))
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
#output$dtDefout #output$eqPick



#module server side function
correctSrc <- function(input, output, session, src, trg) {
  
  #by default, size and organic content are the 3rd and 4th columns of a data frame
  #usr has an option to skip correction, or select (size only, toc only, both)
  output$selcCor <- renderUI({
    ns <- session$ns
    selectInput(ns("slcC"), "Select columns with size, toc or together", choices = names(src[[2]]())[-c(1,2)], multiple = TRUE, width = '100%')
  })
  
  #option to pick columns that will not be used for correction
  output$selcRem <- renderUI({
    ns <- session$ns
    selectInput(ns("slcR"), "Select columns to remove from correction", choices = names(src[[2]]())[-c(1,2)], multiple = TRUE, width = '100%')
  })
  
  #target data input
  trgDT<- reactive ({
    trg <- trg[,colnames(srcDT())]
  })
  
  #source data input
  srcDT <- reactive ({
    src[[2]] ()
  })
  
  #outputs for source dt
  output$srcDT <- renderDT ({
    srcDT ()
  }) 
  
  #outputs for trg data
  output$trgDT <- renderDT({
    trgDT()
  })

  #ui slider input for corr apply
  output$ApplyCor <- renderUI ({
    ns <- session$ns
    actionButton(ns('ApplyCor'), 'Apply corrections')
  })
  
  #select p value threshold of normality for residuals
  output$shapiroP <- renderUI({
    ns <- session$ns
    sliderInput(ns("shapiroP"), "Shapiro-Wilk Univariate Normality Test p-value:", value = 0.05, min = 0.001, max = 1, step = 0.01)
  })
  
  # sliderinput for corrplot R threshold
  output$corR <- renderUI({
    ns <- session$ns
    sliderInput(ns("corR"), "R2", value = 0.6, min = 0, max = 0.99, step = 0.1, animate = F)
  })
  
  #x13
  #apply correction function over data frame 
  listDTcor <- eventReactive(input$ApplyCor, {
    datas <- srcDT ()
    datas <- corect.func(as.data.frame(datas), input$corR, input$shapiroP,input$slcC, input$slcR )
    datas 
  })
  
  #output$x13 #correct
  output$listDTcor <- renderDT({
    datasGlobal <<- listDTcor()
    datas <- listDTcor()
    datas <- datas[, -which(names(datas) %in% c("formula"))]
    colnames(datas)[5] <- "Formula"
    datas[, 5] <- (gsub("I(datas$", replacement = "(", datas[, 5], fixed = T))
    datas[, 5] <- (gsub("datas$", replacement = "", datas[, 5], fixed = T))
    datas[, 5] <- (gsub("I", replacement = "", datas[, 5], fixed = T))
    datas[, 5] <- (gsub("logI", replacement = "log", datas[, 5], fixed = T))
    datas
    if (!is.null(datas)) {
      datas
    } else {
      
    }
  })
  
  #option to pick class #correct #xE4
  output$srcClass <- renderUI ({
    ns <- session$ns
    selectInput(ns("srcClass"), "Select represented classes", choices = unique(listDTcor()[,8]), multiple = F, selected = unique(listDTcor()[,8])[1], width = '60%')
  })
  
  #option to pick unique element #correct
  output$srcVar <- renderUI ({
    ns <- session$ns
    selectInput(ns("srcVar"), "Select represented classes", choices = unique(eqClassDT()[,1]), multiple = F, width = '60%')
  })
  
  #output table of equations, subset by source class #x14
  eqClassDT <- reactive ({
    datas <- listDTcor()
    datas <- datas[which(datas[, 8] == as.character(input$srcClass)), ]
    datas <- datas[!duplicated(datas[, 6]), ]
    datas
  })
  
  #output table of equations, subset by name of the variable #x15 ()
  eqDatoutput <- reactive ({
    datas <- eqClassDT()
    if (!input$srcVar %in% datas[,1]){
      return ()
    } else {
      datas <- datas[which(datas[, 1] == input$srcVar), ]
      datas <- datas[with(datas, order(datas[, 10], -datas[, 9])), ] # ncol(datas)
      datas$rank <- seq(1, nrow(datas))
      datas
    }
  })
  
  #data frame output of equations subset by class and var #output$x15 ()
  output$eqDatoutput <- renderDT({
    req (input$srcVar)
    req (input$srcClass)
    req (eqDatoutput())
    datas <- eqDatoutput()
    s <- datas[which(datas$rank == 1), ]
    datas <- datas[, -which(names(datas) %in% c("formula", "grade", "source"))]
    colnames(datas)[5] <- "Formula"
    datas[, 5] <- (gsub("I(datas$", replacement = "(", datas[, 5], fixed = T))
    datas[, 5] <- (gsub("datas$", replacement = "", datas[, 5], fixed = T))
    datas[, 5] <- (gsub("I", replacement = "", datas[, 5], fixed = T))
    datas[, 5] <- (gsub("logI", replacement = "log", datas[, 5], fixed = T))
    
    DT::datatable(datas, selection = list(mode = "multiple", selected = s, taget = "row")) %>% formatStyle(
      columns = "rank",
      target = "row",
      backgroundColor = styleEqual(1, "cyan"), fontWeight = "bold"
    )
    
  })
  
  #reactive function select unique equation
  selEq <- reactive ({ #Elements7
    req (input$srcClass)
    req (input$srcVar)
    ns <- session$ns
    datas <- eqDatoutput ()
    f <- as.character(unique(datas[, 6]))
    f <- (gsub("I(datas$", replacement = "(", f, fixed = T))
    f <- (gsub("datas$", replacement = "", f, fixed = T))
    f <- (gsub("I", replacement = "", f, fixed = T))
    f <- (gsub("logI", replacement = "log", f, fixed = T))
    f
  })
  
  #select equation
  output$selEq <- renderUI ({
    ns <- session$ns
    selectInput(ns("slEq"), label = "Select Formula", choices = selEq(), width = '100%')
  })
  
  #correction back log (used for display of tables (hidden))
  corBlog <- reactive({ # If none selected, use default selection and apply default formulas
    req(listDTcor())
    datas <- listDTcor()
    datasU <- unique(datas[, 1])
    output <- NULL
    for (i in seq(1, length(datasU))) {
      x <- datas[which(datas[, 1] == datasU[i]), ]
      x <- x[with(x, order(x[, 10], -x[, 9])), ]
      x$rank <- seq(1, nrow(x))
      output <- rbind(output, x)
    }
    output <- output[which(output[, ncol(output)] == 1), ]
    output
  })
  
  #default button
  output$defaultBut <- renderUI ({
    ns <- session$ns
    actionButton(ns("dfBut"), label = "Accept Default")
  })
  
  #reset button
  output$resetBut <- renderUI ({
    ns <- session$ns
    actionButton(ns("reBut"), label = "reset")
  })
  
  #select button
  output$selectBut <- renderUI ({
    ns <- session$ns
    actionButton(ns("slBut"), label = "select")
  })
  
  #output default table showing possible options for correction
  #observe button default
  observeEvent(input$dfBut, {
    print ('def button clicked')
     dtDefout1<- reactive({
      datas <- corBlog()
      datas <- datas[, -which(names(datas) %in% c("formula", "grade"))]
      colnames(datas)[5] <- "Formula"
      datas[, 5] <- (gsub("I(datas$", replacement = "(", datas[, 5], fixed = T))
      datas[, 5] <- (gsub("datas$", replacement = "", datas[, 5], fixed = T))
      datas[, 5] <- (gsub("I", replacement = "", datas[, 5], fixed = T))
      datas[, 5] <- (gsub("logI", replacement = "log", datas[, 5], fixed = T))
      datas
    })
     
     output$dtDefout <- renderDT ({
       dtDefout1 ()
     })
  }) #output$dtDefout
  
  #observe button reset
  observeEvent(input$reBut, {
    datasGlobal <<- listDTcor()
    output$eqPick <- renderDT({
      datasGlobal <- datasGlobal[!duplicated(datasGlobal[, 6]), ]
      datasGlobal <- datasGlobal[,which(!colnames(datasGlobal)%in% c('formula', 'grade', 'Cooks'))]
      colnames(datasGlobal)[5] <- 'formula'
      datasGlobal
    })
  }) #output$eqPick
  
  #observe button select
  observeEvent(input$slBut, {
    f <- datasGlobal[,6]
    f <- (gsub("I(datas$", replacement = "(", f, fixed = T))
    f <- (gsub("datas$", replacement = "", f, fixed = T))
    f <- (gsub("I", replacement = "", f, fixed = T))
    f <- (gsub("logI", replacement = "log", f, fixed = T))
    datasGlobal[which(f == input$slEq & datasGlobal$source == input$srcClass), 11] <<- "SELECTED"
    datasGlobal <<- as.data.frame(datasGlobal)
    
    output$eqPick <- renderDT({
      datasGlobal <- datasGlobal[!duplicated(datasGlobal[, 6]), ]
      #datasGlobal <- subset(datasGlobal, -c(formula, grade, Cooks))
      datasGlobal <- datasGlobal[,which(!colnames(datasGlobal)%in% c('formula', 'grade', 'Cooks'))]
      colnames(datasGlobal)[5] <- 'formula'
      colnames(datasGlobal)[ncol(datasGlobal)] <- 'Selection'
      na.omit(datasGlobal)
    })
    
    dtDefout1 <- reactive({
      datasGlobal <- datasGlobal[, -which(names(datasGlobal) %in% c("formula", "grade"))]
      colnames(datasGlobal)[5] <- "Formula"
      datasGlobal[, 5] <- (gsub("I(datas$", replacement = "(", datasGlobal[, 5], fixed = T))
      datasGlobal[, 5] <- (gsub("datas$", replacement = "", datasGlobal[, 5], fixed = T))
      datasGlobal[, 5] <- (gsub("I", replacement = "", datasGlobal[, 5], fixed = T))
      datasGlobal[, 5] <- (gsub("logI", replacement = "log", datasGlobal[, 5], fixed = T))
      datasGlobal[which(datasGlobal[, ncol(datasGlobal)] == 'SELECTED' & datasGlobal[,6] == input$srcClass), ]
    })
    
    output$dtDefout <- renderDT ({
      dtDefout1 ()
    })
    
  })
  
  #regressions plot function
  regPlotfunc <- reactive({
    req (input$srcClass)
    req (input$srcVar)
    datas <- src [[2]]()
    datas <- datas[which(datas[, 2] == input$srcClass), ]
    datas <- as.data.frame(datas)
  })
  
  #regressions plot output
  output$regPlot <- renderPlot({
    req (input$slEq)
    datas <- regPlotfunc()
    if (length(input$slcC) > 1){
      var1 <- as.numeric(datas[, input$slcC[1]])
      var2 <- as.numeric(datas[, input$slcC[2]])
    } else {
      var1 <- as.numeric(datas[, input$slcC])
      var2 <- NULL
    }

    datas[, 3:dim(datas)[2]] <- apply(datas[, 3:dim(datas)[2]], 2, as.numeric)
    
    f <- input$slEq
    #print (f)
    f <- sub("^.", "I(datas$", f)
    f <- sub("(var1", "I(var1", f, fixed = T)
    f <- sub("(var2", "I(var2", f, fixed = T)
    f <- sub("logI", "log", f, fixed = T)
    f <- sub ('I(datas$og(', 'log(datas$', f, fixed = T)
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
  })
  
  outputCorrected <- reactive({
    req(srcDT())
    req(corBlog())
    req(trgDT())
    
    datas <- srcDT()
    datas <-  datas[, which(!colnames(datas) %in% input$slcR)]
    target <- as.data.frame(trgDT())
    target <- target[,which(!colnames(target)%in% input$slcR)]
    y <- as.data.frame(corBlog())
    x <- as.data.frame(datas)
    output <- list ()
    slopes <- list ()
    drops <- list ()
    
    print (x)
    print (target)
    print (y)
    for (i in seq(1, dim(target)[1])) {
      #print (target[i,])
      dat <- correct(x, target[i, ], y, input$slcC)
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
  
  selTrg <- reactive({
    datas <- trgDT()
    as.character(unique(datas[, 1]))
  })
  
  output$selTrg <- renderUI ({
    ns <-session$ns
    selectInput(ns('sltrg'),'Select target', choices = selTrg())
  })
  
  output$outputCorrected <- renderDT({
    req(outputCorrected())
    req(input$sltrg)
    targetD <- as.data.frame(trgDT())
    indCor <- which(targetD[, 1] == input$sltrg)
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
  
  #return (dtDefout1)
}


# #option to choose type of correction to be used
# radBselection <- function (id) {
#   ns <- NS(id)
#   uiOutput ('radBselection')
# }
# #radio buttons to choose type of correction (can be avaoided with just if statement and checking
# #the length of the input$selection. If null, use default, no correction, if only 1, use it for correction
# #if 2, first one is size, second one is toc)
# output$radBselection <- renderUI({
#   ns <- session$ns
#   radioButtons(ns('cortype'), 'Select type of correction',
#                c('Single' = 'sing',
#                  'Multiple' = 'multiple',
#                  'No correction' = 'default'
#                ))
# })