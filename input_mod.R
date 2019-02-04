library(shiny)
library(DT)

# Module UI function
csvFileInput <- function(id, label = "CSV file") {
  # Create a namespace function using the provided id
  ns <- NS(id)
  tagList(
    fileInput(ns("file"), label,
              multiple = F,
              accept = c(
                "text/csv",
                "text/comma-separated-values,text/plain",
                ".csv"
              )
    )
  )
}


# Module server function
csvFile <- function(input, output, session, stringsAsFactors) {
  # The selected file, if any
  userFile <- reactive({
    #cat ('usr inp')
    # If no file is selected, don't do anything
    validate(need(input$file, message = FALSE))

    input$file
  })
  # The user's data, parsed into a data frame
  dataframe <- reactive({
    options(warn = -1)
    #cat ('user inp read')
    x <- read.csv(userFile()$datapath)
    options(warn = 1)
    x
  })
  # We can run observers in here if we want to
  observe({
    msg <- sprintf("File %s was uploaded", userFile()$name)
  })
  # Return the reactive that yields the data frame
  return(dataframe)
}


# Used to select column of data.table or a data frame
# Module UI function
columnChooserUI <- function(id) {
  ns <- NS(id)
  print (ns('before error'))
  uiOutput(ns("columnChooserUI"))
}


DT_tab <- function(id) {
  ns <- NS(id)
  withSpinner(DTOutput(ns("DT_tab")))
}


mat_par <- function(id) {
  ns <- NS(id)
  tagList(
    uiOutput(ns("corR")),
    uiOutput(ns("corMethod")),
    uiOutput(ns("corType")),
    uiOutput(ns("tl.cex"))
  )
}


srcCor <- function(id) { #simple corrplot
  ns <- NS(id)
  withSpinner(plotOutput(ns("srcCor")))
}


srcDCor <- function(id) { #correlation plot with scatterplots
  ns <- NS(id)
  withSpinner(plotOutput(ns("srcDCor"), height = "1000px"))
}


# Module server function
inputMod <- function(input, output, session, data) {
  
  output$corR <- renderUI({
    ns <- session$ns
    req (data)
    sliderInput(ns("corR"), "R2", value = 0.6, min = 0, max = 0.99, step = 0.1, animate = T)
  })
  
  output$corMethod <- renderUI({
    ns <- session$ns
    req (data)
    selectInput(ns("corMethod"), "Method", c(
      "pie", "circle", "square", "ellipse", "number", "shade",
      "color"
    ))
  })
  
  output$corType <- renderUI({
    ns <- session$ns
    req (data)
    selectInput(ns("corType"), "Type", c( "lower","full", "upper"))
  })
  
  output$tl.cex <- renderUI({
    ns <- session$ns
    req (data)
    sliderInput(ns("tl.cex"), "FontSize", value = 1, min = 0.1, max = 1.5, step = 0.1)
  })
  

  clmns <- reactive({
    ns <- session$ns
    req (data)
    checkboxGroupInput(ns("col"), "Remove: ", names(data), selected = names(data))
  })
  
  # select columns for correlation plot
  output$columnChooserUI <- renderUI({
    clmns()
  })
  
  #cat ('running module')
  DT_tab <- reactive({
    req (input$col)
    if (length(input$col)<1){
      data
    } else {
      data[, input$col]
    }
  })
  
  # table output
  output$DT_tab <- renderDT({
    req(data)
    DT_tab ()
  })
  
  # Correlationp plot
  output$srcCor <- renderPlot({
    req(data)
    req (input$corR)
    dat <- na.omit(data)
    

      M <- cor(dat[, input$col[-c(1,2)]])
      M[M < input$corR & M > -input$corR] <- 0
      corrplot(M, method = input$corMethod, order = "hclust", input$corType, tl.cex = input$tl.cex, diag = FALSE)
    
  }, height = "auto")
  
  # corrplot of src with distributions
  output$srcDCor <- renderPlot({
    req(data)
    if (length(input$col) < 2) {
    } else {
      dat <- na.omit(data)
      dat <- dat[, input$col]
      pairs.panels(dat,
                   method = "pearson", # correlation method
                   hist.col = "#00AFBB",
                   density = TRUE, # show density plots
                   ellipses = TRUE # show correlation ellipses
      )
    }
  })
  
  return (DT_tab)
}









# data frame of shapiro-wilk test p-values
srcSP <- function(id) {
  ns <- NS(id)
  withSpinner(DTOutput(ns("srcSP")))
}

shapiroP <- function(id) {
  ns <- NS(id)
  tagList(
    uiOutput(ns("shapiroP"))
  )
}

getspMethods <- function(id) {
  ns <- NS(id)
  tagList(
    DTOutput(ns("getspMethods"))
  )
}

getspPval <- function(id) {
  ns <- NS(id)
  tagList(
    DTOutput(ns("getspPval"))
  )
}

spPlotpick <- function(id) {
  ns <- NS(id)
  uiOutput(ns("spPlotpick"))
}

getspQQval <- function(id) {
  ns <- NS(id)
  plotOutput(ns("getspQQval"))
}

getorigQQval <- function(id) {
  ns <- NS(id)
  plotOutput(ns("getorigQQval"))
}

# number of standard diviates from normal
nStd <- function(id) {
  ns <- NS(id)
  uiOutput(ns("nStd"))
}

outliersTab1 <- function(id) {
  ns <- NS(id)
  DTOutput(ns("outliersTab1"))
}

# number of standard diviates from normal table output
stdTab1 <- function(id) {
  ns <- NS(id)
  DTOutput(ns("stdTab1"))
}

outliersADtab1 <- function(id) {
  ns <- NS(id)
  DTOutput(ns("outliersADtab1"))
}

outliersADtab2 <- function(id) {
  ns <- NS(id)
  DTOutput(ns("outliersADtab2"))
}

outliersADtab3 <- function(id) {
  ns <- NS(id)
  DTOutput(ns("outliersADtab3"))
}

outliersADtab4 <- function(id) {
  ns <- NS(id)
  DTOutput(ns("outliersADtab4"))
} 
# Module server function
checkD <- function(input, output, session, datas) {

  options (warn = -1)
  ns <- session$ns
  
  #read in data input 
  DT_tab <- reactive ({
    datas ()
    
  })

  # compute shapiro wilk test p value  & output it as data frame
  compSP <- reactive({
    datas <- rawShapiro(DT_tab())
  })

  # output as data frame shapiro wilk table
  output$srcSP <- renderDT({
    # req (DT_tab())
    req(is.factor(DT_tab()[, 2]))
    datas <- compSP()

    if (!is.null(input$shapiroP)) {
      cut <- input$shapiroP
    } else {
      cut <- 0.05
    }
    DT::datatable(datas) %>% formatStyle(
      c(colnames(datas)),
      backgroundColor = styleInterval(c(cut), c("lightsalmon", "white")), fontWeight = "bold"
    )
  })

  # shapiro-wilk test p value slider
  output$shapiroP <- renderUI({
    ns <- session$ns
    sliderInput(ns("shapiroP"), "Shapiro-Wilk Univariate Normality Test p-value:", value = 0.05, min = 0.001, max = 1, step = 0.01)
  })

  # get best transformation methods applied to normalize data
  getspMethods <- reactive({
    req(DT_tab())
    req(compSP())
    req(input$col)
    datas <- DT_tab()

    if (!is.null(input$shapiroP)) {
      cut <- input$shapiroP
    } else {
      cut <- 0.05
    }

    res <- transform(DT_tab(), compSP(), cut)[[2]]
  })

  # output data frame of best methods applied to normalize each class
  output$getspMethods <- renderDT({
    req(is.factor(DT_tab()[, 2]))
    DT::datatable(getspMethods()) %>% formatStyle(
      c(colnames(getspMethods())),
      backgroundColor = styleEqual("None", "lightblue"), fontWeight = "bold"
    )
  })

  # get p values of methods applied after shapiro wilk normalization test
  getspPval <- reactive({
    req(DT_tab())
    req(compSP())
    req(input$col)

    datas <- DT_tab()

    if (!is.null(input$shapiroP)) {
      cut <- input$shapiroP
    } else {
      cut <- 0.05
    }

    res <- transform(datas, compSP(), cut)[[1]]
  })

  # output data frame of best methods applied to normalize each class
  output$getspPval <- renderDT({
    req(is.factor(DT_tab()[, 2]))
    if (!is.null(input$shapiroP)) {
      cut <- input$shapiroP
    } else {
      cut <- 0.05
    }
    DT::datatable(getspPval()) %>% formatStyle(
      c(colnames(getspPval())),
      backgroundColor = styleInterval(c(cut), c("lightsalmon", "white")), fontWeight = "bold"
    )
  })

  # ui to pick column name for qqplots
  output$spPlotpick <- renderUI({
    req(DT_tab())
    req(input$col)
    ns <- session$ns
    selectInput(ns("spPlotpick"), label = "Select Element", choices = colnames(DT_tab())[-c(1, 2)])
  })

  # get shapiro-wilk test applied methods
  getspQQval <- reactive({
    req(DT_tab())
    req(compSP())
    req(input$col)

    datas <- DT_tab()

    if (!is.null(input$shapiroP)) {
      cut <- input$shapiroP
    } else {
      cut <- 0.05
    }

    res <- transform(datas, compSP(), cut)
    y <- res[[4]]
    z <- res[[3]]
    datas <- replacevals(datas, y, z)
  })

  # plotOutput of shapiro wilk transformations, before and after
  output$getspQQval <- renderPlot({
    req(is.factor(DT_tab()[, 2]))

    req(getspQQval())
    datas <- getspQQval()
    tryCatch({
      dat <- datas
      dat <- dat[order(dat[, 2]), ]
      suppressMessages(attach(dat))
      suppressWarnings(assign("val", get(input$spPlotpick)))
      colnames(dat)[2] <- "Classes"
      ggplot(dat, aes(sample = val, colour = Classes)) +
        stat_qq() + facet_wrap(~Classes, ncol = 2, scales = "free") +
        stat_qq_line()
    }, warning = function(cond) {}, error = function(cond) {})
  })

  output$getorigQQval <- renderPlot({
    req(is.factor(DT_tab()[, 2]))

    tryCatch({
      dat <- DT_tab()
      dat <- dat[order(dat[, 2]), ]
      suppressMessages(attach(dat))
      suppressWarnings(assign("val", get(input$spPlotpick)))
      colnames(dat)[2] <- "Classes"
      ggplot(dat, aes(sample = val, colour = Classes)) +
        stat_qq() + facet_wrap(~Classes, ncol = 2, scales = "free") +
        stat_qq_line()
    }, warning = function(cond) {}, error = function(cond) {})
  })

  # renderUI number of standard diviates to be considered as an outlier
  output$nStd <- renderUI({
    ns <- session$ns
    sliderInput(ns("nStd"), "Deviates From Standard Normal Mean For Outliers Detection:", value = 2.576, min = 0, max = 6, step = 0.01)
  })

  # showOutliers reactive function output
  # Output transformed data and values with with adjusent asterix if outlier, without order
  shwOutliers <- reactive({
    dat <- DT_tab()

    if (!is.null(input$nStd)) {
      cut <- input$nStd
    } else {
      cut <- 2.58
    }

    dfoutput <- detectoutliers(getspQQval(), cut)
    datas <- dfoutput[[1]]
    l <- datas[dfoutput[[2]]]
    datas <- cbind(as.character(dat[, 1]), as.character(dat[, 2]), datas)
    colnames(datas) <- colnames(dat)
    datas <- list(datas, l)
  })

  # Output original data and index of outliers, without order
  matchOutliers <- reactive({
    x <- shwOutliers()
    y <- DT_tab()
    output <- x[[1]]
    outliers <- rep("0", dim(output)[1])
    output <- cbind(output, outliers)
    vals <- x[[2]]
    options (warn = -1)
    findRows <- unique((which(is.na(apply(as.data.frame(output[, 3:dim(output)[2]]), 2, as.numeric)), arr.ind = T)[, 1]))
    options (warn = 1)
    output[, ncol(output)][findRows] <- "1"
    output <- as.data.frame(output)
    output <- cbind(y, output$outliers[match(y[, 1], output[, 1])])
    names(output) <- c(colnames(x[[1]]), "outliers")
    # output <- output[order(-(as.numeric(output$outliers))), ]
    transformeddata <- x[[1]]
    transformeddata <- as.data.frame(transformeddata)
    transformeddata <- apply(transformeddata, 2, as.character)
    # transformeddata <- as.matrix(transformeddata)
    vals <- (which(transformeddata %in% x[[2]]))
    output <- list(output, vals)
  })

  # DT with outliers found
  outliersTab1 <- reactive({
    vals <- NULL
    datas <- matchOutliers()
    vals <- c(datas[[2]])
    datas <- as.matrix(datas[[1]])
    datas[is.na(datas)] <- "MISSING"
    if (!is.null(matchOutliers()[[2]])) {
      vals <- c(matchOutliers()[[2]])
      datas[as.numeric(vals)] <- paste(datas[as.numeric(vals)], "*", sep = "")
    }
    vals <- c(datas[vals], "MISSING") # Added this (helped to remove dependence!)
    datas <- datas[order(-(as.numeric(datas[, ncol(datas)]))), ]
    datas <- list(datas, vals)
  })

  # DT table output for outliers found
  output$outliersTab1 <- renderDT({
    if (is.null(input$"outliersADtab2_rows_selected")) {
      datas <- outliersTab1()
      vals <- datas[[2]]
      datas <- (datas [[1]])
    } else {
      datas <- outliersTab1()
      vals <- datas[[2]]
      datas <- (datas [[1]])
      datas[, ncol(datas)] <- 0
      datas[which(datas[, 1] %in% outliersADtab3()[, 1]), ncol(datas)] <- 1
    }

    DT::datatable(datas) %>% formatStyle(
      columns = "outliers",
      target = "row",
      backgroundColor = styleEqual(1, "lightsalmon")
    ) %>% formatStyle(
      columns = colnames(datas),
      backgroundColor = styleEqual(vals, rep("cyan", length(vals)))
    )
  })

  # outliers advanced, tab 1
  outliersADtab1 <- reactive({
    datas <- shwOutliers()
    datas2 <- matchOutliers()
    if (length(matchOutliers()[[2]]) == 0) {
      outliersTab1()
    } else {
      vals <- c(datas[[2]])
      datas <- as.matrix(datas[[1]])
      datas[is.na(datas)] <- "MISSING"
      vals <- c(vals, "MISSING")
      datas2 <- as.matrix(datas2[[1]])
      datas <- cbind(datas, datas2[, ncol(datas2)])
      colnames(datas) <- colnames(datas2)
      datas <- datas[order(-(as.numeric(datas[, ncol(datas)]))), ]
      datas <- list(datas, vals)
      datas
    }
  })

  # DT table output for advanced table 1
  output$outliersADtab1 <- renderDT({
    datas <- outliersADtab1()
    vals <- datas[[2]]
    datas <- datas [[1]]
    DT::datatable(datas) %>% formatStyle(
      columns = colnames(datas),
      backgroundColor = styleEqual(vals, rep("darksalmon", length(vals)))
    )
  })

  # compute data table of standard deviates
  showstd <- reactive({
    dfoutput <- getSubsetstd(getspQQval())
    outliers <- rep(0, dim(dfoutput)[1])
    output <- cbind(dfoutput, outliers)
    if (!is.null(input$nStd)) {
      cut <- input$nStd
    } else {
      cut <- 2.576
    }
    val <- cut
    output[which(dfoutput > val, arr.ind = T)[, 1], ncol(output)] <- 1
    output[which(dfoutput < -val, arr.ind = T)[, 1], ncol(output)] <- 1
    output <- output[order(-(as.numeric(output[, ncol(output)]))), ]

    output
  })

  # output std data tab
  output$stdTab1 <- renderDT({
    if (!is.null(input$nStd)) {
      cut <- input$nStd
    } else {
      cut <- 2.576
    }
    DT::datatable(showstd()) %>% formatStyle(
      c(colnames(showstd())),
      backgroundColor = styleInterval(c(-(cut), cut), c("lightsalmon", "white", "lightsalmon")), fontWeight = "bold"
    )
  })

  # outliers advanced, tab 2 -> selection and diselection of data rows
  output$outliersADtab2 <- renderDT({
    datas <- outliersADtab1()
    vals <- datas[[2]]
    datas <- datas [[1]]
    DT::datatable(datas) %>% formatStyle(
      columns = "outliers",
      target = "row",
      backgroundColor = styleEqual(1, "lightsalmon")
    )
  })

  # Observe function that will run on NULL values (observing selection of rows for removal)
  an_observe_func <- observe(suspended = T, {
    input$"outliersADtab2_rows_selected"
    isolate({
      # do stuff here
      print(input$"outliersADtab2_rows_selected")
    })
  })

  # start the observer, without "suspended=T" the observer
  #  will start on init instead of when needed
  an_observe_func$resume()

  # outliers advanced, tab 3
  outliersADtab3 <- reactive({
    datas <- outliersADtab1()
    datas <- datas [[1]]

    if (is.null(input$"outliersADtab2_rows_selected")) {
      datas <- datas[which(datas[, ncol(datas)] == 1), ]
    } else {
      if (length(input$"outliersADtab2_rows_selected") == 1) {
        datas <- datas[input$"outliersADtab2_rows_selected", ]
        t(datas)
      }
      else {
        datas <- datas[input$"outliersADtab2_rows_selected", ]
      }
    }
    data.frame(datas)
  })

  # outliers advanced, tab 3 Output
  output$outliersADtab3 <- renderDT({
    outliersADtab3()
  })

  # outliers advanced, tab 4 Output (final result)
  outliersADtab4 <- reactive({
    datas <- DT_tab()
    datas[which(!datas[, 1] %in% outliersADtab3()[, 1]), ]
  })

  output$outliersADtab4 <- renderDT({
    outliersADtab4()
  })
  
  options (warn = 1)
  return (outliersADtab4)
}


