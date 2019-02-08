x <- trgDat[[1]] # target data table with marked values

y <- read.csv("SourceData.csv")
y <- na.omit(y)
y <- list(y, y, y, y, y, y, y, y, y, y)
# #1. Create a list for each target list column names where there is a symbol *
# #2. Show 1. target and all column-names selected.
# #3. Select all column names and targets with
# #4. List all targets and column names in the bracket test.
# #4.
# #2. Give an option to user choose column and update column selection
#
grepL <- function(x) {
  return(grep("*", x, fixed = T))
}
#
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

library(shinyWidgets)
library(shinyTree)


# Only run examples in interactive R sessions
if (interactive()) {
  ui <- fluidPage(
    sidebarLayout(
      sidebarPanel(
            shinyTree("tree", checkbox = TRUE, search = TRUE, theme = 'proton', themeIcons = F,
                      themeDots = F)
      ),
      mainPanel(
        DTOutput("dtf")
      )
    )
    # h4('Shiny hierarchical checkbox')

  )

  server <- function(input, output, session) {
    output$dtf <- renderDT({
      #print (input$tree)
      x2()
    })

    x2 <- reactive({
      x
    })
    output$clCh <- renderUI({
      req(x2())
      prettyCheckboxGroup("checkGroup",
        label = h3("Targets"),
        choices = as.character(x$SampleName),
        selected = as.character(x[which(unlist(lapply(trg, length)) > 0), 1]), shape = "curve"
      )
    })

    output$rwCh <- renderUI({
      req(x2())
      prettyCheckboxGroup("checkGroup",
        label = h3("Columns"),
        choices = names(x),
        selected = unique(na.omit(unlist(trg))), inline = T, width = "auto", shape = "curve"
      )
    })

    output$tree <- renderTree({
      atatrib <- function(x) {
        structure(x, stselected = TRUE)
      }

      #trg <- lapply(trg, atatrib)
      trg <- (list("bigtree" = (trg)))

      vL <- NULL
      k <- 1
      for (i in seq(1, length(trg$bigtree))) {
        v <- paste("trg$bigtree$target", paste(i, "$", sep = ""), sep = "")
        for (j in seq(1, length(trg$bigtree[[i]]))) {
          vName <- paste(v, trg$bigtree[[i]][j], sep = "")
          vName <- paste(vName, paste(" <- ", j, sep = ""), sep = "")
          vL[k] <- vName
          k <- k + 1
        }
      }


      for (i in seq(1, length(vL))) {
        options(warn = -1)
        eval(parse(text = vL[i]))
        options(warn = 1)
      }

      atatrib <- function(x) {
        x[which(names(x) == "")] <- NULL
        if (any(names(x) == "None")) {
          structure(x, stselected = FALSE)
        } else {
          structure(x, stselected = TRUE)
        }
      }

      trg <- lapply(trg$bigtree, atatrib)

      trg
    })
    
    observe (
      if (!is.null(input$tree)){
        l <- names(unlist(get_selected(input$tree, format = c("slices"))))
        d <- NULL
        
        if (!is.null(l)){
        for (i in seq (1, length (l))) {
          d <- rbind (d,unlist(strsplit(l[i], "[.]")) )
        }
        
          d<- d[-c(which (d[,1] == d[,2])),]
          print (data.frame(d))
        }
      }

    )
  }
  shinyApp(ui, server)
}
