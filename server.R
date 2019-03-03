

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
source("modules_navtab2.R")
source("input_mod.R")
library("shiny")
library("DT")
library(parallel)
numCores <- detectCores() - 1


server <- function(input, output, session) {
  options(shiny.maxRequestSize = 70 * 1024^2) # Max csv data limit set to 60 mb

  source_data <- callModule(csvFile, "file1")

  source_columns <- callModule(columnChooser, "dat1", source_data)

  output$source_origin <- renderDT({
    source_columns()
  })

  # function for plots
  source_module <- callModule(inputMod, "dat1", source_columns)

  target_data <- callModule(csvFile, "file2")

  target_columns <- callModule(columnChooser, "dat2", target_data)

  output$trgs_origin <- renderDT({
    target_columns()
  })

  trgs_module <- callModule(inputMod, "dat2", target_columns)

  output$trgs_origin_mdl <- renderDT({
    target_function()
  })

  target_function <- reactive({
    trgs_module()[, which(names(trgs_module()) %in% names(src_ref_function()))]
  })

  # transformation output function
  src_origin_mdl <- callModule(checkD, "dat1", source_columns)

  src_ref_function <- reactive({
    src_origin_mdl()[, which(names(src_origin_mdl()) %in% names(trgs_module()))]
  })

  output$src_ref_output <- renderDT({
    src_ref_function()
  })

  output$ui_src_adjustfor <- renderUI({
    req(input$dat_text)
    selectInput("slcC", "Select columns with size, toc or together", choices = names(src_ref_function())[-c(1, 2)], multiple = TRUE, width = "100%")
  })

  # option to pick columns that will not be used for correction
  output$ui_src_remove <- renderUI({
    req(input$dat_text)
    selectInput("slcR", "Select columns to remove from correction", choices = names(src_ref_function())[-c(1, 2)], multiple = TRUE, width = "100%")
  })

  # select p value threshold of normality for residuals
  output$ui_src_shapiro <- renderUI({
    sliderInput("ui_src_shapiro", "Shapiro-Wilk Univariate Normality Test p-value:", value = 0.05, min = 0.001, max = 1, step = 0.01)
  })

  # sliderinput for corrplot R threshold
  output$ui_src_cor <- renderUI({
    sliderInput("ui_src_cor", "R", value = 0.6, min = 0, max = 0.99, step = 0.1, animate = F)
  })

  output$ui_src_split <- renderUI({
    sliderInput("ui_src_split", "Split proportion", value = 0.9, min = 0.1, max = 0.99, step = 0.05, animate = F)
  })

  output$ui_src_applycor <- renderUI({
    actionButton("ui_src_applycor", "Apply")
  })

  output$ui_src_applymix <- renderUI({
    actionButton("ui_src_applymix", "Use mixing")
  })

  # option to pick columns that will not be used for correction
  output$ui_dfa_remove <- renderUI({
    selectInput("dfaR", "Select columns to remove from DFA", choices = targets_brackets_output_func()[-c(1, 2)], multiple = TRUE, width = "100%")
  })



  src_corr_function <- eventReactive(input$ui_src_applycor, {
    if (input$corBut == "Cor") {
      datas <- src_ref_function()
      datas <- corect.func(as.data.frame(datas), input$ui_src_cor, input$ui_src_shapiro, input$slcC, input$slcR)
      datas <- datas[, -which(names(datas) %in% c("formula"))]
      datas <- datas[!duplicated(datas[, 5]), ]
      datas$id <- seq(1, nrow(datas))
      if (!is.null(datas)) {
        datas
      } else {

      }
    } else {
      NULL
    }
  })

  output$src_corr_function <- renderDT({
    if (!is.null(src_corr_function())) {
      datas <- src_corr_function()
      colnames(datas)[5] <- "Formula"
      datas[, 5] <- (gsub("I(datas$", replacement = "(", datas[, 5], fixed = T))
      datas[, 5] <- (gsub("datas$", replacement = "", datas[, 5], fixed = T))
      datas[, 5] <- (gsub("I", replacement = "", datas[, 5], fixed = T))
      datas[, 5] <- (gsub("logI", replacement = "log", datas[, 5], fixed = T))
      if (!is.null(datas)) {
        datas
      } else {

      }
    } else {
      NULL
    }
  })

  output$dat_text <- renderUI({
    textInput("dat_text", "", value = "Active")
  })

  corr_formulas_function <- reactive({
    if (!is.null(src_corr_function())) {
      dat <- src_corr_function()
      val_table <- dat
      varr <- val_table %>% group_by(source, element)
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
    } else {
      NULL
    }
  })

  corr_formulas_validtab_formula <- reactive({
    if (!is.null(src_corr_function())) {
      dat1 <- src_corr_function()
      dat2 <- correct_formulas_function()
      dat <- dat1[which(dat1$id %in% dat2$id), ]
      dat <- cbind(dat[, 1], 0, dat[, 2:ncol(dat)])
      dat
    }
    else {
      NULL
    }
  })

  corrected_function <- reactive({
    if (input$corBut == "Cor") {
      req(src_ref_function()) # source data
      req(corr_formulas_validtab_formula()) # formulas format
      req(target_function()) # target data
      req(input$ui_src_targets)

      datas <- src_ref_function()
      datas <- datas[, which(!colnames(datas) %in% input$slcR)]
      target <- as.data.frame(target_function())
      target <- target[, which(!colnames(target) %in% input$slcR)]

      y <- as.data.frame(corr_formulas_validtab_formula())
      x <- as.data.frame(datas)
      output <- list()
      slopes <- list()
      drops <- list()

      for (i in seq(1, dim(target)[1])) {
        dat <- correct(x, target[i, ], y, input$slcC)
        output[[i]] <- dat[[1]]
        #print ('output1 success')
        slopes[[i]] <- dat[[2]]
        #print ('output2 success')
        if (length(dat[[3]]) > 0) {
          drops [[i]] <- dat[[3]]
          #print ('output3 success')
        } else {
          drops[[i]] <- NA
        }
      }
      
      #print ('loop 1 success')

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
      
      #print ('slopes success')

      colnames(slopes.DT) <- c("var1", "var2", "var1*var2")
      slopes.DT_names <- rep(NA, nrow(slopes.DT))
      slopes.DT_names[seq(1, length(as.character(correct_formulas_function()[, 1])))] <- as.character(correct_formulas_function()[, 1])
      slopes.DT <- cbind(slopes.DT_names, slopes.DT)

      cleanT <- function(x) {
        x <- (gsub("I(datas$", replacement = "(", x, fixed = T))
        x <- (gsub("datas$", replacement = "", x, fixed = T))
        x <- (gsub("I", replacement = "", x, fixed = T))
        x <- (gsub("logI", replacement = "log", x, fixed = T))
      }
      drops <- lapply(drops, cleanT)

      drops <- plyr::ldply(drops, rbind)
      rownames(drops) <- target[, 1]

      list(output, slopes.DT, drops)
    } else {
      req(src_ref_function())
      req(target_function())
      output <- list()
      for (i in seq(1, nrow(target_function()))) {
        datas <- src_ref_function()
        output [[i]] <- datas
      }
      print("function outputcorrected")
      list(output, NULL, NULL)
    }
  })

  output$ui_src_targets <- renderUI({
    selectInput("ui_src_targets", "Target ID", choices = seq(1, nrow(target_function())), selected = 1)
  })

  correct_src_selected_function <- reactive({
    req(corrected_function())
    print (input$ui_src_targets)
    #try (print (corrected_function()[[1]]))
    corrected_function() [[1]][[as.numeric(input$ui_src_targets)]]
  })

  output$correct_src_selected_output <- renderDT({
    correct_src_selected_function()
  })

  output$correct_src_slopes_output <- renderDT({
    corrected_function() [[2]]
  })

  output$correct_src_drops_output <- renderDT({
    corrected_function() [[3]]
  })

  output$corr_formulas_output <- renderDT({
    if (!is.null(corr_formulas_function())) {
      dat <- corr_formulas_function()[, c(3, 2, 1, 6, 4, 5)]
      datatable(dat, filter = "top", options = list(
        pageLength = 5, autoWidth = F
      ))
    } else {
      NULL
    }
  })


  # radio button serve as indicator of user choice or default choice
  output$ui_formulas_selected <- renderUI({
    radioButtons(
      "ui_formulas_selected", "Select data :",
      c(
        "Top rank" = "def",
        "User choice" = "sel"
      )
    )
  })

  # get the top ranked formulas as function
  correct_formulasdef_function <- reactive({
    dat <- corr_formulas_function()
    datas <- dat[which(dat$rank == 1), ]
  })

  # if the radio button on default, in the function output top picks, else selected user choice
  correct_formulas_function <- reactive({
    req(input$ui_formulas_selected)
    if (input$ui_formulas_selected == "def") { # if it is default, have as an output top table
      datas <- correct_formulasdef_function()
    } else {
      datas <- corr_formulas_selected_function() # if it is user choice, have an output selected
    }
    datas
  })

  # table of top picks
  output$correct_formulasdef_output <- renderDT({
    datatable(correct_formulasdef_function(), filter = "top", options = list(
      pageLength = 5, autoWidth = F
    ))
  })

  # selected formulas after corretion #take the output of top choices, and apply selection
  corr_formulas_selected_function <- reactive({
    # req (input$corr_formulas_output_rows_selected)
    corr_formulas_function()[input$corr_formulas_output_rows_selected, ] # check
  })

  # table of final picks
  output$correct_formulas_output <- renderDT({
    datatable(correct_formulas_function(), filter = "top", options = list(
      pageLength = 5, autoWidth = F
    ))
  })

  output$ui_corvar_plot <- renderUI({
    if (!is.null(input$ui_formulas_selected)) {
      if (input$ui_formulas_selected == "def") {
        selectInput("ui_corvar_plot", label = "Select Element", choices = correct_formulas_function()[, 3])
      } else {
        selectInput("ui_corvar_plot", label = "Select Element", choices = correct_formulas_function()[, 3])
      }
    } else {}
  })

  output$ui_corvar_plot2 <- renderUI({
    if (!is.null(input$ui_formulas_selected)) {
      if (input$ui_formulas_selected == "def") {
        selectInput("ui_corvar_plot2", label = "Select Element", choices = correct_formulas_function()[, 3])
      } else {
        selectInput("ui_corvar_plot2", label = "Select Element", choices = correct_formulas_function()[, 3])
      }
    } else {}
  })


  # regressions plot function
  regression_plot_function <- reactive({
    req(correct_formulasdef_function())
    datas <- src_ref_function()
    if (!is.null(input$ui_formulas_selected)) {
      if (input$ui_formulas_selected == "def") {
        if (!is.null(input$ui_corvar_plot)) {
          sclass <- correct_formulasdef_function()[which(correct_formulasdef_function()[, 3] == input$ui_corvar_plot), 4]
          sclass <- as.character(sclass)
          datas <- datas[which(datas[, 2] == sclass), ]
        } else {}
      } else {
        if (!is.null(input$ui_corvar_plot)) {
          sclass <- corr_formulas_selected_function()[which(corr_formulas_selected_function()[, 3] == input$ui_corvar_plot), 4]
          sclass <- as.character(sclass)
          datas <- datas[which(datas[, 2] == sclass), ]
          datas
        } else {}
      }
    }
    datas
  })

  # regressions plot output
  output$regression_plot_output <- renderPlot({
    req(input$ui_formulas_selected)
    req(input$ui_corvar_plot)
    req(corr_formulas_selected_function())
    req(input$slcC)
    req(correct_formulas_function())
    req(regression_plot_function())
    dats <- NULL
    dats <- corr_formulas_selected_function()
    slcC <- input$slcC

    if (dim(dats)[1] < 1 & input$ui_formulas_selected == "sel") {} else {
      # print ('before error')
      datas <- regression_plot_function()

      if (length(slcC) > 1) {
        if (dim(datas)[1] == 0) {
          return(NULL)
        } # fixes error when we switch default to user
        var1 <- as.numeric(datas[, slcC[1]])
        var2 <- as.numeric(datas[, slcC[2]])
      } else {
        # print("else")
        # print(dim(datas))
        var1 <- as.numeric(datas[, slcC])
        var2 <- NULL
      }

      datas[, 3:dim(datas)[2]] <- apply(datas[, 3:dim(datas)[2]], 2, as.numeric)

      f <- input$ui_corvar_plot
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

  # Interaction and regression plots
  output$regression_plot_output2 <- renderPlotly({
    req(input$ui_formulas_selected)
    req(input$ui_corvar_plot2)
    req(corr_formulas_selected_function())
    req(input$slcC)
    req(correct_formulas_function())
    req(regression_plot_function())

    dats <- NULL
    dats <- corr_formulas_selected_function()
    slcC <- input$slcC

    if (dim(dats)[1] < 1 & input$ui_formulas_selected == "sel") {} else {
      # print ('before error')
      datas <- regression_plot_function()

      if (length(slcC) > 1) {
        if (dim(datas)[1] == 0) {
          return(NULL)
        } # fixes error when we switch default to user
        var1 <- as.numeric(datas[, slcC[1]])
        var2 <- as.numeric(datas[, slcC[2]])
      } else {
        #print("else")
        # print(dim(datas))
        var1 <- as.numeric(datas[, slcC])
        var2 <- NULL
      }

      datas[, 3:dim(datas)[2]] <- apply(datas[, 3:dim(datas)[2]], 2, as.numeric)

      # print (input$ui_corvar_plot2)
      f <- input$ui_corvar_plot2
      f <- sub("^.", "I(datas$", f)
      f <- sub("(var1", "I(var1", f, fixed = T)
      f <- sub("(var2", "I(var2", f, fixed = T)
      f <- sub("logI", "log", f, fixed = T)
      f <- sub("I(datas$og(", "log(datas$", f, fixed = T)
      fit <- lm(eval(parse(text = f)))

      dat_names <- names(datas)


      for (i in seq(1, length(dat_names))) {
        ind <- grep(dat_names[i], f)
        if (length(ind) > 0) {
          elementName <- dat_names[i]
        }
      }
      response <- datas[elementName]
      response <- response[, 1]
      intr <- grep("\\*", f)

      if (!is.null(var2)) {
        p <- ggplot(datas, aes(y = response, x = var1, color = var2)) + geom_point() + stat_smooth(method = "lm", se = FALSE)
        ggplotly(p, height = 600, width = 800)
      } else {
        p <- ggplot(datas, aes(y = response, x = var1)) + geom_point() + stat_smooth(method = "lm", se = FALSE)
        ggplotly(p, height = 600, width = 800)
      }
    }
  })

  output$xvar <- renderUI({
    selectInput("xvar", "Observed, Xvar", choices = names(correct_src_selected_function())[-c(1, 2)])
  })

  output$xvar2 <- renderUI({
    selectInput("xvar2", "Corrected, elements", choices = names(correct_src_selected_function())[-c(1, 2)])
  })

  output$yvar <- renderUI({
    selectInput("yvar", "Corrected, for", choices = c(input$slcC)) # 1234
  })

  output$scatterplot1 <- renderScatterD3({
    req(input$xvar)
    # req(input$yvar)
    mtdf <- correct_src_selected_function()
    dtorg <- src_ref_function()
    
    #print (mtdf)
    #print (dtorg)
    #print (input$xvar)
    
    
    x <- dtorg[[input$xvar]]
    Sclass <- dtorg [, 2]
    y <- mtdf[[input$xvar]]
    scatterD3(
      x = x, y = y,
      labels_size = 9, point_opacity = 1, col_var = Sclass,
      lines = data.frame(slope = 1, intercept = Inf),
      transitions = T
    )
  })

  output$scatterplot2 <- renderScatterD3({
    req(input$xvar2)
    req(input$yvar)
    mtdf <- correct_src_selected_function()
    dtorg <- src_ref_function()
    x <- dtorg[[input$xvar2]]
    Sclass <- dtorg [, 2]
    y <- mtdf[[input$yvar]]
    scatterD3(
      x = x, y = y,
      labels_size = 9, point_opacity = 1, col_var = Sclass,
      # lines = data.frame(slope = 1, intercept = Inf),
      transitions = T
    )
  })

  output$brack_range <- renderUI({
    sliderInput("brack_rng", "Bracket range parameter", value = 0.1, min = 0, max = 3, step = 0.05)
  })

  # output$targetDrop <- renderText ({
  #   print ('These targets will be removed, unless selected one will not be used for DFA')
  # })

  targets_brackets_function <- reactive({
    req(corrected_function())
    req(target_function())
    req(input$brack_rng)
    x <- corrected_function()[[1]] # list of sources corrected
    datas <- src_ref_function()
    datas <- datas[, which(!names(datas) %in% names(x[[1]]))]
    y <- trgs_module() # list of targets corrected
    y <- as.data.frame(y)

    # datas <- datas[,which(!names(datas) %in% colnames(y))]
    dat <- NULL
    trg <- NULL
    for (i in seq(1, nrow(y))) {
      x[[i]] <- cbind(x[[i]], datas)
      l <- bracketT(x[[i]], y[i, ], input$brack_rng, input$slcR)
      dat <- rbind(dat, l[[1]])
      trg <- c(trg, l[[2]])
    }

    dat[, 1] <- factor(as.character(dat[, 1]), levels = c(as.character(dat[, 1])))
    list(dat, trg)
  })

  output$targets_brackets_output <- renderDT({
    req(targets_brackets_function())
    dat <- targets_brackets_function()
    selection <- dat[[2]]

    # print (selection)
    # print (names (dat[[1]]))

    DT::datatable(dat[[1]]) %>% formatStyle(
      c(colnames(dat[[1]])),
      backgroundColor = styleEqual(selection, rep("lightsalmon", length(selection)))
    )
  })

  targets_brackets_output_func <- reactive({
    dat <- targets_brackets_function()
    names(dat[[1]])
  })



  target_droplist_function <- reactive({
    req(targets_brackets_function())

    dat <- targets_brackets_function()
    # print (dat)

    x <- dat[[1]]
    trg <- list()
    for (i in seq(1, nrow(x))) {
      k <- (colnames(x)[grepL(x[i, ])])

      if (length(k) == 0) {
        trg[[i]] <- "None"
      }
      else {
        trg[[i]] <- k
      }
    }
    names(trg) <- x$SampleName
    trg
  })

  output$target_droplist_output <- renderDT(
    target_removed_function(),
    selection = "multiple"
  )

  target_removed_function <- reactive({
    req(target_droplist_function())
    trg <- target_droplist_function()
    dat <- NULL
    for (i in seq(1, length(trg))) {
      tname <- names(trg)[i]
      for (j in seq(1, length(trg[i]))) {
        val <- trg[i][[j]]
        dat <- rbind(dat, cbind(tname, val))
      }
    }
    dat
  })

  dfa_apply_function <- reactive({
    req(corrected_function())
    sourceList <- corrected_function()[[1]]
    d <- data.frame(target_removed_function())
    # if (!is.null(input$target_droplist_output_rows_selected)) {
    #   d <- d[-input$target_droplist_output_rows_selected, ]
    # }
    d[, 1] <- as.character(d[, 1])
    d[, 2] <- as.character(d[, 2])
    x <- targets_brackets_function()[[1]]
    x[, 1] <- as.character(x[, 1])

    # probably target
    tdrops <- c(input$dfaR)
    x <- x[, !(names(x) %in% tdrops)]

    y <- sourceList

    for (i in seq(1, nrow(d))) {
      drops <- d[, 2][which(d[, 1] == d[, 1][i])]
      drops <- c(drops, tdrops)
      dat <- y[[which(x[, 1] == d[, 1][i])]]
      y[[which(x[, 1] == d[, 1][i])]] <- dat[, which(!names(dat) %in% c(drops))]
    }
    dat <- stepwiseDFA(y)
    dat
  })


  # Pick up here
  ####################################################################################################################################
  ####################################################################################################################################
  ####################################################################################################################################
  ####################################################################################################################################


  dfaList_x <- eventReactive(input$applyDFA, {

    # input$dfaR


    req(dfa_apply_function())
    req(target_function())
    datas <- dfa_apply_function()
    dfaList <- NULL
    sourceList <- datas
    targetList <- target_function()

    for (i in seq(1, nrow(targetList))) {
      dfa <- data.frame(matrix(data = 0, nrow = 1, ncol = (length(colnames(targetList)) - 2)))
      names(dfa) <- names(targetList)[-c(1, 2)]
      datas_i <- datas[[i]]
      datas_i <- t(datas_i)
      dfa[, match(datas_i[1, ], names(dfa))] <- as.numeric(datas_i[3, ])
      dfaList <- rbind(dfaList, dfa)
    }
    dfaList
  })

  output$dfaList <- renderDT(
    if (input$rbDFA == "def") {
      dfaList_x()
    } else {
      dat <- target_function()
      d <- dim(dat)
      dat_output <- data.frame(matrix(100, nrow = d[1], ncol = (d[2] - 2)))
      colnames(dat_output) <- names(dat)[-c(1, 2)]
      dat_output <- dat_output[, -c(which(colnames(dat_output) %in% target_removed_function()[, 2]))]
      rownames(dat_output) <- NULL
      dat_output
    }
  )

  mixingOutput <- eventReactive(input$ui_src_applymix, {
    req(dfaList_x())
    req(input$rbMix)
    l <- corrected_function()[[1]]
    DFA_l <- dfaList_x()
    targetD <- as.data.frame(target_function())
    finalDat <- NULL

    if (input$rbMix == "ALL") {

    } else {
      l <- l[which(target_function()[, 1] %in% input$selected_targets)]
      targetD <- targetD[which(targetD[, 1] %in% input$selected_targets), ]
      DFA_l <- DFA_l[which(targetD[, 1] %in% input$selected_targets), ]
    }

    for (i in seq(1, length(l))) {
      target <- targetD[i, -c(1, 2)]
      DFA <- DFA_l[i, ]
      x <- l[[i]]
      uniSource <- unique(x[, 2])
      uniSource <- as.character(uniSource)
      ui_src_split <- input$ui_src_split
      modelOutput <- NULL
      #print(paste0("Mixing source: ", i))

      for (j in seq(1, input$mcsimulations)) {
        inputTrain <- NULL
        inputValidate <- NULL
        for (i2 in seq(1, length(uniSource))) {
          dat <- x[which(x[, 2] == uniSource[i2]), ]
          train_index <- sample(1:nrow(dat), nrow(dat) * ui_src_split)
          training_dat <- dat[train_index, ]
          validate_dat <- dat[-train_index, ]
          inputTrain <- rbind(inputTrain, training_dat)
          inputValidate <- rbind(inputValidate, validate_dat)
        }
        datas <- getSubsetmean(inputTrain[, -1])

        DFA <- DFA[(which(colnames(DFA) %in% colnames(datas)))]
        DFA <- DFA[, colSums(DFA != 0) > 0]
        target <- target[, which(names(target) %in% colnames(DFA))]
        datas <- datas[, which(colnames(datas) %in% colnames(DFA))]


        dat <- inputValidate [, -c(1, 2)]
        dat <- dat[, which(names(dat) %in% colnames(DFA))]
        dat <- rbind(dat, target)

        if (any(dat == 0)) {
          dat[dat == 0] <- 0.001
        }

        rownames(dat) <- c(as.character(inputValidate[, 1]), as.character(targetD[i, 1]))
        # for (i3 in seq (1, nrow (dat))){
        # output <- UseUnMixing(dat[i3,], datas, DFA, method = "Nelder-Mead")
        # modelOutput <- rbind (modelOutput, output)
        # }

        cl <- makeCluster(numCores)

        optimMix <- function(x) {
          datas <- get("datas", envir = environment())
          DFA <- get("DFA", envir = environment())
          output <- UseUnMixing(x, datas, DFA, method = "Nelder-Mead")
        }
        clusterExport(cl, list("UseUnMixing", "datas", "DFA", "dat", "optimMix"), envir = environment())
        output <- t(parApply(cl, dat, 1, optimMix))
        output <- cbind(output, j)
        modelOutput <- rbind(modelOutput, output)
        stopCluster(cl)
      }
      modelOutput <- cbind(modelOutput, i)

      finalDat <- rbind(finalDat, modelOutput)
      finalDat <- round(finalDat, 3)
    }
    x <- as.data.frame(finalDat)
    x[, (ncol(x) - 1)] <- paste0("Monte Carlo: ", x[, (ncol(x) - 1)])
    x[, (ncol(x))] <- paste0("Target: ", x[, (ncol(x))])
    colnames(x) <- c(uniSource, "GOF", "Monte-Carlo", "Target")
    x
  })


  output$mixingOutput <- renderDT({
    x <- as.data.frame(mixingOutput())
    x
  })

  output$targetPlot <- renderUI({
    req(mixingOutput())
    dat <- as.data.frame(mixingOutput())
    dat <- dat[grep("target", rownames(dat)), ]

    un <- rownames(dat[which(rownames(dat) %in% as.character(target_function()[, 1])), ])

    if (!is.null(input$selected_targets)) {
      dat <- dat[which(rownames(dat) %in% input$selected_targets), ]
    }
    # selectInput('targetPlot', 'Select target', unique(dat[,ncol(dat)]), selected = NULL)
    selectInput("targetPlot", "Select target", un, selected = NULL)
  })

  output$srcPlot <- renderUI({


  })

  output$trg_mixing_plot <- renderPlot({
    req(mixingOutput())
    dat <- as.data.frame(mixingOutput())
    dat <- dat[grep("target", rownames(dat)), ]
    targNames <- rownames(dat)
    targetSelected <- dat$Target[which(rownames(dat) %in% input$targetPlot)]
    dat <- melt(dat)
    # datmeltGlobal <<- dat
    tryCatch({
      req(input$targetPlot)
      # dat <- dat[which (dat[,2] == input$targetPlot),]
      dat <- dat[which(dat[, 2] == targetSelected), ]
      ggplot(dat, aes(factor(variable), value, colour = variable)) +
        geom_violin(trim = FALSE) + geom_jitter(height = 0, width = 0.1) + geom_boxplot(width = 0.1, color = "black", alpha = 0.5) +
        facet_wrap(~variable, ncol = 2, scales = "free") #+ theme(legend.position="bottom")
    }, warning = function(cond) {}, error = function(cond) {})
  }, height = 400, width = 600)


  output$src_mixing_plot <- renderPlot({
    req(mixingOutput())
    dat <- as.data.frame(mixingOutput())
    targs <- dat[grep("target", rownames(dat)), ]
    dat <- dat[!which(rownames(dat) %in% rownames(targs))]
    dat <- melt(dat)
    # print(dat)
  })

  output$radioBut <- renderUI({
    radioButtons("rbMix", "Select mixing", choices = c("ALL", "SUBSET"), selected = "ALL")
  })

  output$selectTarget <- renderUI({
    req(input$rbMix)
    # print (target_function())
    if (input$rbMix == "ALL") { } else {
      checkboxGroupInput("selected_targets", "Targets",
        choices = unique(as.character(target_function()[, 1])), selected = unique(as.character(target_function()[, 1])), inline = FALSE,
        width = NULL
      )
    }
  })

  output$corBut <- renderUI({
    radioButtons(
      "corBut", "Correct for any of the elements?",
      c(
        "CORRECT" = "Cor",
        "NO CORRECTION" = "noCor"
      ),
      selected = "Cor"
    )
  })

  output$rbDFA <- renderUI({
    radioButtons(
      "rbDFA", "Apply DFA, default or uniform weights?",
      c(
        "DEFAULT" = "def",
        "UNIFORM" = "uni"
      ),
      selected = "def"
    )
  })

  output$verb1 <- renderText({
    "dat 1"
  })
  output$verb2 <- renderText({
    "dat 2"
  })
}


# User interface side of the user input
ui <- dashboardPage(
  skin = "blue",
  dashboardHeader(title = "SedSat_ShinyV2"),
  ## Sidebar content
  dashboardSidebar(
    sidebarMenu(
      menuItem("INPUT", tabName = "dataInput", icon = icon("upload")),
      # menuItem("Transformations & Outliers",
      #   icon = icon("upload"),
      #   menuSubItem("TRANSFORMATION", tabName = "Transformations"),
      #   menuSubItem("OUTLIERS", tabName = "Outliers")
      # ),
      menuItem("TRANSFORMATION", tabName = "Transformations", icon = icon("refresh")),
      menuItem("OUTLIERS", tabName = "Outliers", icon = icon("remove")),
      menuItem("SIZE & TOC CORRECTION", tabName = "regressions", icon = icon("random")),
      menuItem("BRACKET TEST", tabName = "Bracket", icon = icon("braille")),
      menuItem("DISCRIMINANT FUNCTION", tabName = "DFA", icon = icon("tasks")),
      menuItem("MIXING MODEL", tabName = "mixmod", icon = icon("rocket"))
      # ,
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
                  div(
                    style = "width: 50%;",
                    sidebarPanel(
                      fluidRow(
                        column(
                          width = 12,
                          csvFileInput("file1", "User data (.csv format)"),
                          columnChooserUI("dat1")
                        )
                      ),
                      width = 4
                    )
                  ),
                  mainPanel(
                    width = 10,
                    fluidRow(
                      tabsetPanel(
                        tabPanel(
                          "Data",
                          column(
                            width = 12,
                            DTOutput("source_origin"),
                            style = "height:'auto'; overflow-y: scroll;overflow-x: scroll;"
                          )
                        ),
                        navbarMenu(
                          "Plots",
                          tabPanel(
                            "Corplot",
                            fluidRow(
                              column(
                                width = 9,
                                br(),
                                srcCor("dat1")
                              ),
                              column(
                                width = 3,
                                mat_par("dat1")
                              )
                            )
                          ),
                          tabPanel(
                            "Distribution",
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
              )
            )
          ),
          tabPanel(
            "Target",
            fluidPage(
              fluidRow(
                sidebarLayout(
                  div(
                    style = "width: 50%;",
                    sidebarPanel(
                      fluidRow(
                        column(
                          width = 12,
                          csvFileInput("file2", "User data (.csv format)"),
                          columnChooserUI("dat2")
                        )
                      ),
                      width = 4
                    )
                  ),
                  mainPanel(
                    width = 10,
                    fluidRow(
                      tabsetPanel(
                        tabPanel(
                          "Data",
                          column(
                            width = 12,
                            DTOutput("trgs_origin"),
                            style = "height:'auto'; overflow-y: scroll;overflow-x: scroll;"
                          )
                        ),
                        navbarMenu(
                          "Plots",
                          tabPanel(
                            "Corplot",
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
                            )
                          ),
                          tabPanel(
                            "Distribution",
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
        )
      ),
      tabItem(
        tabName = "Transformations",
        box(
          title = ("Shapiro-Wilk test of normality, p values before and after transformations"), status = "success", height = "auto", width = 12, solidHeader = T,
          fluidPage(
            title = "Adjust p-value",
            sidebarLayout(
              div(
                style = "width: 50%;",
                sidebarPanel(
                  shapiroP("dat1"),
                  spPlotpick("dat1")
                ), width = 4
              ),
              mainPanel(
                width = 10,
                column(
                  width = 12,
                  tabsetPanel(
                    tabPanel(
                      "Normality test, Shapiro-Wilk p values",
                      box(
                        title = "Untransformed data p-values", status = "success", height = "630", width = 12, solidHeader = T,
                        srcSP("dat1"),
                        style = "height: 550px; overflow-y: scroll; overflow-x: scroll;"
                      )
                    ),
                    navbarMenu(
                      "Transformations",
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
                        "Achieved p-values",
                        box(
                          title = "Methods p-values ", status = "danger", height = "630", width = 12, solidHeader = T,
                          column(
                            width = 12,
                            getspPval("dat1"),
                            style = "height: 550px; overflow-y: scroll; overflow-x: scroll;"
                          )
                        )
                      )
                    ),
                    navbarMenu(
                      "QQ-Plots",
                      tabPanel(
                        "Transformed and untransformed plots",
                        box(
                          title = "QQ plot of Original Data", status = "success", height = "auto", width = 6, solidHeader = T,
                          getorigQQval("dat1"), style = "height:630px;overflow-y: scroll;overflow-x: scroll;"
                        ),
                        box(
                          title = "QQ plot of Transformed Data", status = "success", height = "auto", width = 6, solidHeader = T,
                          getspQQval("dat1"), style = "height:630px;overflow-y: scroll;overflow-x: scroll;" # 570
                        )
                      ),
                      tabPanel(
                        "Only untransformed",
                        box(
                          title = "QQ plot of Original Data", status = "success", height = "auto", width = 12, solidHeader = T,
                          getorigQQval2("dat1"), style = "height:630px;overflow-y: scroll;overflow-x: scroll;"
                        )
                      ),
                      tabPanel(
                        "Only transformed",
                        box(
                          title = "QQ plot of Transformed Data", status = "success", height = "auto", width = 12, solidHeader = T,
                          getspQQval2("dat1"), style = "height:630px;overflow-y: scroll;overflow-x: scroll;" # 570
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
      tabItem(
        tabName = "Outliers",
        box(
          title = ("Outliers within standard deviates"), status = "success", height = "auto", width = 12, solidHeader = T,
          fluidPage(
            title = "Adjust number of std",
            sidebarLayout(
              div(
                style = "width: 50%;",
                sidebarPanel(
                  nStd("dat1"),
                  rbSl("dat1"),
                  uiOutput("dat_text")
                ), width = 4
              ),
              mainPanel(
                width = 10,
                column(
                  width = 12,
                  tabsetPanel(
                    tabPanel(
                      "Data outliers",
                      box(
                        title = "Rows with outliers and selected outliers", status = "success", height =
                          "595", width = "12", solidHeader = T,
                        column(
                          width = 12,
                          tags$hr(),
                          outliersTab1("dat1"), style = "height:500px; overflow-y: scroll;overflow-x: scroll;"
                        )
                      )
                    ),
                    navbarMenu(
                      "Review",
                      tabPanel(
                        "Edit selection",
                        box(
                          # tags$style(HTML(".js-irs-1 .irs-single, .js-irs-1 .irs-bar-edge, .js-irs-1 .irs-bar {background: green;}")),
                          # tags$style(HTML(".js-irs-1 .irs-single, .js-irs-1 .irs-bar-edge, .js-irs-1 .irs-grid-text { font-size: 10pt;}")),
                          title = "Outliers ", status = "success", height =
                            "695", width = 12, solidHeader = T,
                          column(
                            width = 12,
                            outliersADtab2("dat1"), style = "height:500px; overflow-y: scroll;overflow-x: scroll;"
                          )
                        )
                      ),
                      tabPanel(
                        "Standard deviates",
                        box(
                          title = "Standard Normal Deviate", status = "success", height =
                            "695", width = 12, solidHeader = T,
                          column(
                            width = 12,
                            stdTab1("dat1"), style = "height:500px; overflow-y: scroll;overflow-x: scroll;"
                          )
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
            "Corrections",
            sidebarLayout(
              div(
                style = "width: 50%;",
                sidebarPanel(
                  uiOutput("corBut"),
                  uiOutput("ui_src_adjustfor"), # ,
                  uiOutput("ui_src_remove"),
                  uiOutput("ui_src_shapiro"),
                  uiOutput("ui_src_cor"),
                  uiOutput("ui_src_applycor")
                ), width = 4
              ),
              mainPanel(
                width = 10,
                box(
                  title = "Available options", status = "success", height =
                    "auto", solidHeader = T, width = "auto",
                  column(
                    width = 12,
                    withSpinner(DTOutput("src_corr_function")), style = "height:'auto'; overflow-y: scroll;overflow-x: scroll;"
                  )
                )
              )
            )
          ),
          tabPanel(
            "Target corrected",
            tabsetPanel(
              tabPanel(
                "Data Tables",
                fluidRow(
                  column(
                    width = 12,
                    box(
                      title = "Formulas", status = "success", height =
                        "auto", solidHeader = T,
                      DTOutput("corr_formulas_output"), style = "height:'auto'; overflow-y: scroll;overflow-x: scroll;"
                    ),
                    box(
                      title = "Selected", status = "primary", height =
                        "auto", solidHeader = T,
                      DTOutput("correct_formulas_output"), style = "height:'auto'; overflow-y: scroll;overflow-x: scroll;"
                    )
                  ),
                  column(
                    width = 12,
                    box(
                      title = "Top pick", status = "success", height =
                        "auto", solidHeader = T,
                      DTOutput("correct_formulasdef_output"), style = "height:'auto'; overflow-y: scroll;overflow-x: scroll;"
                    ),
                    box(
                      title = "Plots", status = "success", height =
                        "auto", solidHeader = T,
                      uiOutput("ui_formulas_selected"),
                      uiOutput("ui_corvar_plot"),
                      plotOutput("regression_plot_output")
                    )
                  )
                )
              ),
              tabPanel(
                "Plots",
                box(
                  title = "QQ plot of Original Data", status = "success", height = "600", width = 4, solidHeader = T,
                  uiOutput("ui_corvar_plot2")
                ),
                box(
                  title = "QQ plot of Original Data", status = "success", height = "670", width = 8, solidHeader = T,
                  column(
                    width = 12,
                    plotlyOutput("regression_plot_output2") # ,style = "height:'auto'; overflow-y: scroll;overflow-x: scroll;"
                  )
                )
              )
            )
          ),
          tabPanel(
            "Corrected Data",
            tabsetPanel(
              tabPanel(
                "Data",
                column(
                  width = 12,
                  uiOutput("ui_src_targets"),
                  DTOutput("correct_src_selected_output"), style = "height:'600'; overflow-y: scroll;overflow-x: scroll;"
                ),
                fluidRow(
                  box(
                    title = "Slopes", status = "primary", height = "800", solidHeader = T,
                    column(
                      width = 12,
                      DTOutput("correct_src_slopes_output"), style = "height:'600'; overflow-y: scroll;overflow-x: scroll;"
                    )
                  ),
                  box(
                    title = "Dropped", status = "warning", height = "800", solidHeader = T,
                    column(
                      width = 12,
                      DTOutput("correct_src_drops_output"), style = "height:'600'; overflow-y: scroll;overflow-x: scroll;"
                    )
                  )
                )
              ),
              navbarMenu(
                "Plots",
                tabPanel(
                  "Transformed~Original, 1-1 Plot",
                  uiOutput("xvar"),
                  scatterD3Output("scatterplot1")
                ),
                tabPanel(
                  "Adjusted Plot vs Correction",
                  uiOutput("xvar2"),
                  uiOutput("yvar"),
                  scatterD3Output("scatterplot2")
                )
              )
            )
          )
        )
      ),
      tabItem(
        "Bracket",
        # tabPanel(
        #   "Bracket test",
        fluidPage(
          uiOutput("brack_range"),
          sidebarLayout(
            div(
              style = "width: 50%;",
            sidebarPanel(
              # textOutput('targetDrop'),
              br(),
              br(),
              # uiOutput ('ui_dfa_remove'),
              DTOutput("target_droplist_output"),
              style = "height:'auto'; overflow-y: scroll;overflow-x: scroll;"
            ), width = 4),
            mainPanel(
              DTOutput("targets_brackets_output"),
              style = "height:'auto'; overflow-y: scroll;overflow-x: scroll;"
            )
          )
        )
        # )
      ),
      tabItem(
        tabName = "DFA",
        fluidPage(
          fluidRow(
            sidebarLayout(
              div(
                style = "width: 50%;",
              sidebarPanel(
                uiOutput("rbDFA"),
                uiOutput("ui_dfa_remove"),
                actionButton("applyDFA", "Apply")
              ), width = 4),
              mainPanel(
                column(
                  width = 12,
                  withSpinner(DTOutput("dfaList")), style = "height:auto; overflow-y: scroll;overflow-x: scroll;"
                )
              )
            )
          )
        )
      ),
      tabItem(
        tabName = "mixmod",
        fluidPage(
          fluidRow(
            sidebarLayout(
              div(
                style = "width: 50%;",
              sidebarPanel(
                uiOutput("ui_src_applymix"),
                uiOutput("ui_src_split"),
                uiOutput("radioBut"),
                uiOutput("selectTarget"),
                numericInput("mcsimulations", "Monte carlo simulations:", 2, min = 1, max = 1000)
              ), width = 4),
              mainPanel(
                tabsetPanel(
                  tabPanel(
                      'Table',
                      withSpinner(DTOutput("mixingOutput")), style = "height:auto; overflow-y: scroll;overflow-x: scroll;"
                    ),
                    tabPanel(
                      'Plot',
                      uiOutput("targetPlot"),
                      plotOutput("trg_mixing_plot", width = "100%")
                    )
                  )
                
                  # uiOutput ('srcPlot'),
                  # fluidRow (
                  # column(
                  #  width = 6,
                 # ,
                # column(
                # width =6,
                # plotOutput ('src_mixing_plot', width = '100%')
                # )
                # )
              )
            )
          )
        )
      )
    )
  )
)
# Run the app ----
shinyApp(ui, server)
