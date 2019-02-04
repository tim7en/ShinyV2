library(stringr)
library(plotly)
# By default, p value of equations is 0.05
corect.func <- function(x, R2, p_res, corFor, drops) {
  x <- x[, which(!colnames(x) %in% drops)]
  datas.main <- x

  # Create possible combination of formulas (7^3 possible combinations)
  eq1 <- c("(conc)", "log(conc)", "I(conc^2)", "I(conc^(1/2))", "I(conc^(1/3))", "I(conc^(-1))")
  eq2 <- c("var1", "log(var1)", "I(var1^2)", "I(var1^(1/2))", "I(var1^(1/3))", "I(var1^(-1))")
  eq3 <- c("var2", "log(var2)", "I(var2^2)", "I(var2^(1/2))", "I(var2^(1/3))", "I(var2^(-1))")
  eqComb <- expand.grid(eq1, eq2, eq3)
  indp <- paste(eqComb$Var2, eqComb$Var3, sep = "*")
  reg.eq <- paste(eqComb$Var1, indp, sep = "~")
  eqSize <- expand.grid(eq1, eq2)
  regSize.eq <- paste(eqSize$Var1, eqSize$Var2, sep = "~")

  # Calculate the number of cores
  no_cores <- detectCores() - 1

  # Initiate clustering
  cl <- makeCluster(no_cores)

  # Find unique sources
  un.source <- unique(datas.main[, 2])
  OutputTab <- NULL

  for (i in seq(1, length(un.source))) {

    # Extract only banks from the sources
    datas <- datas.main[which(datas.main[, 2] == un.source[i]), ]
    cooksDist <- function(x) {
      x <- lm(eval(parse(text = x)))
      cDist <- cooks.distance(x)
      cDist <- cDist / max(cDist)
      return(sum(cDist))
    }

    residSd <- function(x) {
      x <- lm(eval(parse(text = x)))
      fit <- lowess(fitted(x), residuals(x))
      return(sd(fit$y) / diff(quantile(fit$y, probs = c(.05, .95))))
    }

    lmp <- function(modelobject) { # From stackexchange
      if (class(modelobject) != "lm") stop("Not an object of class 'lm' ")
      f <- summary(modelobject)$fstatistic
      p <- pf(f[1], f[2], f[3], lower.tail = F)
      attributes(p) <- NULL
      return(p)
    }

    stepOver <- function(x) {
      out <- tryCatch({
        stepOut <- step(lm(eval(parse(text = x))), trace = F)
        if (length(stepOut$coefficients) > 1) {
          pval <- (lmp(stepOut))
          funcCall <- deparse(formula(stepOut))
        } else {
          pval <- 1
          funcCall <- deparse(formula(stepOut))
        }
        return(c(summary(stepOut)$adj, shapiro.test(stepOut$residuals)$p.value, pval, funcCall))
      },
      error = function(cond) {
        # Choose a return value in case of error
        return(c(NA, NA, NA, NA))
      }
      )
      return(out)
    }


    if (length(corFor) < 2) {
      var1 <- datas[, corFor]
      datas <- datas[!colnames(datas) %in% corFor]
      datas <- datas[, -c(1, 2)]

      # Create table that we need to fill in with the values
      table.out <- as.data.frame(matrix(nrow = (nrow(eqComb) * (dim(datas)[2])), ncol = 6))
      colnames(table.out) <- c("element", "formula", "R2", "p-resi", "p-eq", "stepV")
      table.out$element <- rep(colnames(datas), nrow(eqComb))
      table.out$element <- table.out$element[order(match(table.out$element, colnames(datas)))]
      table.out$element <- as.character(table.out$element)
      table.out$formula <- rep(regSize.eq, dim(datas)[2])

      # Parse each element (formula) of the list and apply evaluate equation
      table.out$formula <- str_replace(table.out$formula, "conc", paste("datas", table.out$element, sep = "$"))

      table.out <- table.out[!duplicated(table.out[, 2]), ]

      table.out[, 3:6] <- t(parSapply(cl, table.out$formula, stepOver))
      # table.out [, 3:6] <- t(sapply(table.out$formula, stepOver))
      table.out <- na.omit(table.out)

      table.out[, 3:5] <- apply(table.out[, 3:5], 2, as.numeric)
      table.out[, 3:5] <- round(table.out[, 3:5], 4)
      table.out$grade <- 0
      # Set threshold for p value of R2, p value of residuals and p of the model
      table.out$grade[which(table.out$R2 > R2 & table.out$`p-resi` > p_res & table.out$`p-eq` < 0.05)] <- 1
      table.out$source <- un.source[i]

      #
      el_left <- table.out[which(table.out$grade == 1), ]

      un_el <- unique(el_left$element)

      if (length(un_el) == 0) {

      } else {
        for (j in seq(1, length(un_el))) {
          p <- el_left[which(el_left$element == un_el[j]), ]
          p$Cooks <- sapply(p$stepV, cooksDist)
          p$Cooks <- round(p$Cooks, 3)
          p$residualsSD <- sapply(p$stepV, residSd)
          p$residualsSD <- round(p$residualsSD, 5)
          OutputTab <- rbind(OutputTab, p)
        }
      }
    } else if (length(corFor) == 2) {
      var1 <- datas [, corFor[1]]
      var2 <- datas [, corFor[2]]
      # Sort data out, extract concentration of elements that are not present in the list of size only
      datas <- datas[!colnames(datas) %in% corFor]
      datas <- datas[, -c(1, 2)]


      # Create table that we need to fill in with the values
      table.out <- as.data.frame(matrix(nrow = (nrow(eqComb) * (dim(datas)[2])), ncol = 6))
      colnames(table.out) <- c("element", "formula", "R2", "p-resi", "p-eq", "stepV")
      table.out$element <- rep(colnames(datas), nrow(eqComb))
      table.out$element <- table.out$element[order(match(table.out$element, colnames(datas)))]
      table.out$element <- as.character(table.out$element)
      table.out$formula <- rep(reg.eq, dim(datas)[2])

      # table.out$formula[which(table.out$element %in% size.ly)] <- paste(paste("datas", as.character(table.out$element[which(table.out$element %in% size.ly)]), sep = "$"), "size", sep = "~")

      # Parse each element (formula) of the list and apply evaluate equation
      table.out$formula <- str_replace(table.out$formula, "conc", paste("datas", table.out$element, sep = "$"))

      table.out[, 3:6] <- t(parSapply(cl, table.out$formula, stepOver)) # t(parSapply (cl, table.out$formula, stepOver))
      # table.out[, 3:6] <- t(sapply(table.out$formula, stepOver))
      table.out[, 3:5] <- apply(table.out[, 3:5], 2, as.numeric)
      table.out[, 3:5] <- round(table.out[, 3:5], 4)
      table.out$grade <- 0

      # Set threshold for p value of R2, p value of residuals and p of the model
      table.out$grade[which(table.out$R2 > R2 & table.out$`p-resi` > p_res & table.out$`p-eq` < 0.05)] <- 1
      table.out$source <- un.source[i]

      el_left <- table.out[which(table.out$grade == 1), ]
      un_el <- unique(el_left$element)
      if (length(un_el) == 0) {}
      else {
        for (j in seq(1, length(un_el))) {
          p <- el_left[which(el_left$element == un_el[j]), ]
          p$Cooks <- sapply(p$stepV, cooksDist)
          p$Cooks <- round(p$Cooks, 3)
          p$residualsSD <- sapply(p$stepV, residSd)
          p$residualsSD <- round(p$residualsSD, 5)
          OutputTab <- rbind(OutputTab, p)
        }
      }
    }
  }
  stopCluster(cl)
  closeAllConnections()
  if (!is.null(OutputTab)) {
    return(OutputTab)
  } else {return (NULL)}
}

# plot_ly(p, x=~p$`p-resi`, y=~p$Cooks, z=~p$R2, type="scatter3d", mode="markers", marker = list(size = p$Cooks*100)) %>%
#       add_markers(color=p$R2)
#
# fit<- lm (p$stepV[which(p$Cooks==max(p$Cooks))])
# cooksDist (fit)
# plot (residuals(fit))
# plot (fitted(fit),residuals(fit))
# acf(residuals(fit))
#
# #Proof
# set.seed(1);            #set the seed for reproducability
# N = 100;                #Sample size
# x = runif(N)            #Independant variable
# beta = 4;               #Regression coefficient
# epsilon = rnorm(N);     #Error with variance 1 and mean 0
# y = x * beta + epsilon  #Your generative model
# lin_mod <- lm(y ~x)  #Your linear model
#
# plot (residuals(fit), col = 'red', pch = 16, main = 'Our model') #Our model
# plot (residuals(lin_mod), col = 'blue', pch = 16, main = 'Random model') #Randmly generated model
# plot (fitted(fit),residuals(fit), col = 'red', pch = 16, main = 'Our model') #Our model
# plot(fitted(lin_mod), residuals(lin_mod), col = 'blue', pch = 16,  main = 'Random model') #Randmly generated model
# cooksDist (fit) #Our model
# cooksDist(lin_mod)#Randomly generated model
