#target bracket test function
tarBrackets <- reactive({
  # req(downloadFinal2())
  # req(targetData())
  # req(is.character(downloadFinal2()[, c(1:2)]) == TRUE)
  # x <- downloadFinal2()
  # y <- trg.Data()
  
  x <- read.csv ('SourceData.csv')
  x <- na.omit (x)
  x <- list(x,x,x,x,x,x,x,x,x,x)
  y <- read.csv ('targetData.csv')
  y[,3:ncol(y)] <- apply (y[,3:ncol(y)],2, as.numeric)
  y[,3:ncol(y)] <- y[,3:ncol(y)] *0.5

  bracketT<- function (x,y, r) {
    y <- y[,-c(1,2)]
    x <- x[, 3:ncol(x)]
    x <- apply(x, 2, as.numeric)
    x <- x[, -ncol(x)]
    colSrc <- match(colnames(x), colnames(y))
    yNum <- y[, na.omit(colSrc)]
    YNum_2 <- y[, na.omit(colSrc)]
    upperL <- x * (1 + r)
    lowerL <- x * (1 - r)
    l <- list()
    for (i in seq(1, ncol(upperL))) {
      yNum[(which(yNum [, i ] > max(upperL[, i]) | yNum[, i] < min(lowerL[, i]))), i] <- paste(yNum[(which(yNum [, i ] > max(upperL[, i]) | yNum[, i] < min(lowerL[, i]))), i], "*", sep = "")
      l[[i]] <- yNum[(which(YNum_2 [, i ] > max(upperL[, i]) | YNum_2[, i] < min(lowerL[, i]))), i]
    }
    dat <- cbind(y[, c(1, 2)], yNum) #take first and second column and combine them
    list(dat, c(na.omit(unlist(l))))
  }

})

#bracket test table output
output$tarBrackets <- renderDT({
  dat <- tarBrackets()
  selection <- dat[[2]]
  DT::datatable(dat[[1]]) %>% formatStyle(
    c(colnames(dat[[1]])),
    backgroundColor = styleEqual(selection, rep("lightsalmon", length(selection)))
  )
})





dat <- NULL
trg <- NULL
for (i in seq (1,nrow(y))){
  l <- bracketT (x[[i]],y[i,], 0.1)
  dat <- rbind (dat, l[[1]])
  trg <- c (trg, l[[2]])
}


