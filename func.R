#extract shapiro wilk test p-values
extractpval <- function(y) {
  return(y$p.value)
}

#output data frame of shapiro wilk test p-values
rawShapiro <- function(x) {
  sourcefacts <- unique(x[, 2]) #column 2 should be factor (it used to find classes)
  l <- c(as.character(sourcefacts))
  dfShapiro <- function(y) {
    source <- x[which(x[, 2] == y), ]
    shapres <- apply(source[, 3:dim(source)[2]], 2, shapiro.test)
    shapiropvals <- unlist(lapply(shapres, extractpval))
  }
  lvals <- lapply(l, dfShapiro)
  names(lvals) <- l
  dfoutput <- t(round(as.data.frame(lvals), 3))
  return(dfoutput)
}

#transformation function for shapiro wilk test and best available normalization transformation
transform <- function(x, y, z) {
  shapiroP <- z
  rowcol <- which(y < shapiroP, arr.ind = T)
  pnames <- rownames(rowcol)
  rowcol <- as.data.frame(rowcol)
  rowcol$names <- pnames
  rownames(rowcol) <- NULL
  applytransform <- function(val) {
    xisolate <- x[which(x[, 2] == as.character(val[3])), as.numeric(val[2]) + 2]
    if (any(na.omit(xisolate) > 0) == FALSE) {
      xisolate <- xisolate * -1
      df <- data.frame(matrix(unlist(tfunc(xisolate)), nrow = length(xisolate), byrow = FALSE), stringsAsFactors = FALSE)
      names(df) <- paste(tnames, "*", sep = "")
    } else if (any(na.omit(xisolate) < 0) == FALSE) {
      df <- data.frame(matrix(unlist(tfunc(xisolate)), nrow = length(xisolate), byrow = FALSE), stringsAsFactors = FALSE)
      names(df) <- tnames
    } else if (any(na.omit(xisolate) == 0) == TRUE) {
      df <- data.frame(matrix(unlist(tfunczero(xisolate)), nrow = length(xisolate), byrow = FALSE), stringsAsFactors = FALSE)
      names(df) <- tnameszero
    } else {
      df <- data.frame(matrix(unlist(tfuncWneg(xisolate)), nrow = length(xisolate), byrow = FALSE), stringsAsFactors = FALSE)
      names(df) <- tnamesneg
    }
    shapresult <- apply(df, 2, shapiro.test)
    pshap <- unlist(lapply(shapresult, extractpval))
    xisolate <- df[, which(pshap == max(pshap))]
    pshap <- pshap[which(pshap == max(pshap))]
    output <- c(pshap, val, xisolate)
  }
  funcout <- apply(rowcol, 1, applytransform)
  extrp <- function(y) {
    return(y[1])
  }
  extrn <- function(y) {
    return(y[3])
  }
  extrval <- function(y) {
    return(y[-c(1, 2, 3, 4)])
  }
  
  pshapirotransformed <- unlist(lapply(funcout, extrp))
  pshapnames <- lapply(funcout, extrn)
  transformationtable <- y
  transformationtable [is.finite(transformationtable)] <- "None"
  transformationtable[which(y < shapiroP)] <- unlist(names(pshapirotransformed))
  names(pshapirotransformed) <- colnames(y)[as.numeric(unlist(pshapnames))]
  y[which(y < shapiroP)] <- as.numeric(unlist(as.numeric(pshapirotransformed)))
  vals <- lapply(funcout, extrval)
  return(list(round(y, 3), transformationtable, vals, rowcol))
}

#available functions for transformation applied for Shapiro-Wilk
tfunc <- function(x) {
  c(x, x^2, x^(1 / 2), x^(1 / 3), x^(-1), x^(-1 / 2), log10(x))
}
tnames <- c("x", "x^2", "x^(1/2)", "x^(1/3)", "x^(-1)", "x^(-1/2)", "log10(x)")
tfuncWneg <- function(x) {
  c(x, sign(x) * (abs(x))^(1 / 3), x^(-1))
}
tnamesneg <- c("x", "sign(x)*abs(x)^(1/3)", "x^(-1)")
tfunczero <- function(x) {
  c(x, sign(x) * (abs(x))^(1 / 3))
}
tnameszero <- c("x", "sign(x)*(abs(x)^(1/3)")

#needed for shapiro wilk test text return for QQ plots
replacevals <- function(x, y, z) {
  for (i in seq(1, length(z))) {
    x[which(x[, 2] == as.character(y[i, 3])), as.numeric(y[i, 2]) + 2] <- as.numeric(z[[i]])
  }
  return(x)
}

#tag asterix to outliers in the df
outlierRemoval <- function(y, N_outlier) {
  upper <- mean(y) + N_outlier * sd(y)
  lower <- mean(y) - N_outlier * sd(y)
  y[which(y > upper | y < lower)] <- paste(y[which(y > upper | y < lower)], "*", sep = "")
  return(y)
}

# Returns data frame of transformed data in original diminsion and location of outliers
detectoutliers <- function(x, y) {
  N_outlier <- y
  sourcefacts <- unique(x[, 2])
  l <- c(as.character(sourcefacts))
  x[, 3:dim(x)[2]] <- round(x[, 3:dim(x)[2]], 3)
  outliercheck <- function(y) {
    source <- x[which(x[, 2] == y), ]
    outliersFound <- apply(source[, 3:dim(source)[2]], 2, outlierRemoval, N_outlier = N_outlier)
  }
  lvals <- lapply(l, outliercheck)
  names(lvals) <- l
  dfoutput <- do.call("rbind", lvals)
  cond2 <- grep("\\*", dfoutput)
  # dfoutput <- apply(dfoutput,2,as.numeric)
  datas <- dfoutput
  datas <- list(datas, cond2)
}

findstd <- function(x) {
  x <- (x - mean(as.numeric(na.omit(x)))) / sd(as.numeric(na.omit(x)))
  return(x)
}

findMean <- function(x) {
  x <- mean(x)
  return(x)
}

getSubsetstd <- function(x) {
  sourcefacts <- unique(x[, 2])
  l <- c(as.character(sourcefacts))
  findstdval <- function(j) {
    source <- x[which(x[, 2] == j), ]
    stdfound <- apply(source[, 3:dim(source)[2]], 2, findstd)
  }
  lvals <- lapply(l, findstdval)
  names(lvals) <- l
  dfoutput <- do.call("rbind", lvals)
  dfoutput <- round(dfoutput, 3)
  dfoutput
}

getSubsetmean <- function(x) {
  sourcefacts <- unique(x[, 1])
  l <- c(as.character(sourcefacts))
  getMean <- function(j) {
    source <- x[which(x[, 1] == j), ]
    source[, 2:dim(source)[2]] <- apply(source[, 2:dim(source)[2]], 2, as.numeric)
    stdfound <- apply(source[, 2:dim(source)[2]], 2, findMean)
  }
  lvals <- lapply(l, getMean)
  names(lvals) <- l
  dfoutput <- do.call("rbind", lvals)
  dfoutput <- round(dfoutput, 3)
  dfoutput
}

untransform <- function(x, y, un) {
  loc1 <- which(y[(which(rownames(y) %in% un)), ] == "x^2")
  if (length(loc1) > 0) {
    x[, loc1] <- x[, loc1]^(1 / 2)
  }

  loc2 <- which(y[(which(rownames(y) == un)), ] == "x^(1/2)")
  if (length(loc2) > 0) {
    x[, loc2] <- x[, loc2]^(2)
  }

  loc3 <- which(y[(which(rownames(y) == un)), ] == "sign(x)*(abs(x))^(1/3)")
  if (length(loc3) > 0) {
    x[, loc3] <- x[, loc3]^3
  }

  loc4 <- which(y[(which(rownames(y) == un)), ] == "x^(-1)")
  if (length(loc4) > 0) {
    x[, loc4] <- x[, loc4]^(-1)
  }

  loc5 <- which(y[(which(rownames(y) == un)), ] == "x^(-1/2)")
  if (length(loc5) > 0) {
    x[, loc5] <- x[, loc5]^(-2)
  }

  loc6 <- which(y[(which(rownames(y) == un)), ] == "log10(x)")
  if (length(loc6) > 0) {
    x[, loc6] <- 10^(x[, loc6])
  }

  loc7 <- which(y[(which(rownames(y) == un)), ] == "x^(1/3)")
  if (length(loc7) > 0) {
    x[, loc7] <- (x[, loc7])^(3)
  }
  return(x)
}

CDF_mat <- function(Gen_mat, N, stat) { # Input, bootstrap numb, bin_size
  Gen_mat <- as.data.frame(Gen_mat)
  # print (dim(Gen_mat))
  x_mean <- array(0, dim(Gen_mat)[2])
  x_std <- array(0, dim(Gen_mat)[2])
  # boot_col_mean<-function (x) {boot(x,samplemean,100)[1][[1]]}
  # boot_col_std<-function (x) {boot(x,samplestd,100)[1][[1]]}
  # boot_col_median <- function (x) {boot (x, samplemedian,100) [1][[1]]}
  # boot_col_Robstd <- function (x) {boot (x, sampleRobustStd,100)[1][[1]]}

  x_mean <- colMeans(Gen_mat, na.rm = T)
  colSd <- function(x, na.rm = FALSE) apply(X = x, MARGIN = 2, FUN = sd, na.rm = na.rm)
  x_std <- colSd(Gen_mat, na.rm = T)
  colMedians <- function(x, na.rm = FALSE) apply(X = x, MARGIN = 2, FUN = median, na.rm = na.rm)
  x_median <- colMedians(Gen_mat, na.rm = TRUE) # boot_col_median)
  tryCatch(
    if (any(Gen_mat < 0) == F) {
      x_robstd <- apply(Gen_mat, 2, Qn)
    } else {
      stat <- "conv"
    },
    warning = function(cond) {
      stat <- "conv"
    }, error = function(cond) {
      stat <- "conv"
    }
  )


  if (stat == "lhs-conv") {
    Output <- randomLHS(N, dim(Gen_mat)[2])
    for (i in seq(1:dim(Gen_mat)[2])) {
      Output[, i] <- qnorm(Output[, i], x_mean[i], x_std[i])
    }
  } else if (stat == "lhs-rob") {
    Output <- randomLHS(N, dim(Gen_mat)[2])
    for (i in seq(1:dim(Gen_mat)[2])) {
      Output[, i] <- qnorm(Output[, i], x_median[i], x_robstd[i])
    }
  }

  if (stat == "conv") {
    Output <- matrix(NA, nrow = N, ncol = dim(Gen_mat)[2])
    for (i in seq(1:dim(Gen_mat)[2])) {
      Output[, i] <- rnorm(N, x_mean[i], x_std[i])
    }
  } else if (stat == "rob") {
    Output <- matrix(NA, nrow = N, ncol = dim(Gen_mat)[2])
    for (i in seq(1:dim(Gen_mat)[2])) {
      tryCatch(
        Output[, i] <- rnorm(N, x_median[i], x_robstd[i]),
        warning = function(cond) {
          Output[, i] <- qnorm(Output[, i], x_mean[i], x_std[i])
        }, error = function(cond) {
          Output[, i] <- qnorm(Output[, i], x_mean[i], x_std[i])
        }
      )
    }
  } else {
    print("Methods available (lhs-rob,lhs-conv,conv,rob)")
  }
  colnames(Output) <- colnames(Gen_mat)
  return(Output)
}


UseUnMixing <- function(samples, sources, weights, method = "Nelder-Mead") {
  clr <- function(x) {
    logx <- log(x)
    return(logx - mean(logx))
  }
  Iclr <- function(x) {
    expx <- exp(x)
    return(expx / sum(expx))
  }
  # Create the function to compute the "unmixing" model modified from Collins
  # Note that the parameters in par omit the first one, which is computed in this function
  Collins <- function(par, Conc, src, wt) {
    # Recover all source proportions
    par <- Iclr(c(0 - sum(par), par))
    SigmaPSi <- par %*% src
    return(sum(((Conc - SigmaPSi) / Conc)^2 * wt))
  }
  
  GOF <- function(par, Conc, src, wt) {
    # Recover all source proportions
    par <- Iclr(c(0 - sum(par), par))
    SigmaPSi <- par %*% src
    return(-(1 - (sqrt(sum(((Conc - SigmaPSi) / Conc)^2 * wt / 100))) / dim(src)[2]))
  }
  
  samples <- t(as.data.frame (samples))
  #print (samples)
  GOF_score <- double(nrow(samples))
  Pctgs <- matrix(0, nrow = nrow(samples), ncol = nrow(sources))
  # Step through each row in samples
  start <- rep(0, nrow(sources) - 1L)


  for (i in seq(nrow(samples))) {
    best <- 0
    WT <- weights
    for (j in seq(nrow(sources))) {
      src <- rbind(sources[j, , drop = FALSE], sources[-j, , drop = FALSE])
      ret <- optim(start, GOF, method = method, Conc = samples, src = src, wt = WT)
      ret$value <- sqrt(ret$value^2)
      if (ret$value > best) { # save the results
        best <- ret$value
        GOF_score[i] <- ret$value
        #GOF_score <- ret$value
        Pctgs[i, ] <- Iclr(append(ret$par, 0 - sum(ret$par), after = j - 1L))
        #Pctgs <- Iclr(append(ret$par, 0 - sum(ret$par), after = j - 1L))
      }
    }
  }
  
  colnames(Pctgs) <- rownames(sources)
  rownames(Pctgs) <- rownames (samples)
  return(cbind(Pctgs, GOF_score = GOF_score))
}





ranNumberGenUnmixing <- function(datas, way, Nsamples, stat, checkbox, origdf) {
  datas <- datas[, -ncol(datas)]
  l <- as.character(unique(datas[, 2]))
  way <- as.data.frame(way)
  generateRandoms <- function(j) {
    source <- datas[which(datas[, 2] == j), ]
    Gen_mat <- apply(source[, 3:dim(source)[2]], 2, as.numeric)
    transformedSource <- CDF_mat(Gen_mat, Nsamples, stat)
    untransformedDat <- untransform(transformedSource, way, j)
    untransformedDat <- round(untransformedDat, 3)
    untransformedDat <- untransformedDat [!is.infinite(rowSums(untransformedDat)), ]
    untransformedDat <- untransformedDat[complete.cases(untransformedDat), ]
    datas <- cbind(rep(j, Nsamples), untransformedDat)
  }
  lvals <- lapply(l, generateRandoms)
  names(lvals) <- l
  datasTemp <- do.call("rbind", lvals)
  colnames(datasTemp) <- c("SourceType", colnames(datasTemp)[-1])
  datas <- datasTemp
  negVal <- (unique(which(origdf < 0, arr.ind = T)[, 2])) - 1
  negatives <- (which(colnames(datas) %in% checkbox))
  negatives <- negatives[!negatives %in% negVal]
  datas[datas[, negatives] < 0] <- NA
  datas <- na.omit(datas)
}

powerTransform <- function(y, lambda1, lambda2 = NULL, method = "boxcox") {
  boxcoxTrans <- function(x, lam1, lam2 = NULL) {
    
    # if we set lambda2 to zero, it becomes a one parameter transformation (log(x))
    lam2 <- ifelse(is.null(lam2), 0, lam2)
    
    if (lam1 == 0L) {
      log(y + lam2)
    } else {
      (((y + lam2)^lam1) - 1) / lam1
    }
  }
  
  switch(method
         , boxcox = boxcoxTrans(y, lambda1, lambda2)
         , tukey = y^lambda1
  )
}

invBoxCox <- function(x, lambda) {
  if (lambda == 0) exp(x) else (lambda * x + 1)^(1 / lambda)
}


# lm_func <- function(y) {
#   lm(y ~ size * toc, data = input, na.action = "na.pass")
# }
# 
# 
# lm_size <- function(y) {
#   lm(y ~ size, na.action = "na.pass")
# }


step_lm <- function(y) {
  step(y, na.action = "na.pass")
}


extractcoeff <- function(y) {
  return(y$coefficients[-1])
}

extractcoefNames <- function(y) {
  return(names(y$coefficients)[-1])
}

extractrsquared <- function(y) {
  return(y$adj.r.squared)
}


extractpval <- function(y) {
  return(y$p.value)
}

# 
# lm_boxcox <- function(y) {
#   boxcox(y, lambda = seq(-4, 4, 0.005), plotit = F)
# }


getlambda <- function(y) {
  lambda <- y$x[which.max(y$y)]
}



getlmeq <- function(cc) {
  (eqn <- paste(paste(cc, names(cc), sep = " * (", collapse = ") + "),')',sep =""))
}


onlyslopes <- function(cc) {
  (eqn <- paste(paste(cc[-1], names(cc[-1]), sep = " * ", collapse = "+"), sep = " + "))
}

lm_func <- function(y) {
  lm(y ~ size * toc, data = input, na.action = "na.pass")
}

lm_size <- function(y) {
  lm(y ~ size, data = input, na.action = "na.pass")
}

step_lm <- function(y) {
  step(y, data = input, na.action = "na.pass")
}

lm_boxcox <- function(y) {
  boxcox(y, lambda = seq(-4, 4, 0.005), plotit = F)
}

lm_boxcoxNew <- function (y) {
  boxcox(y~size*toc, data= input, na.action = "na.pass")
}


#Working size & organic correction function
correct <- function(datas, trg, dff, corFor) {

  x <- as.data.frame(datas)
  x[,3:ncol(x)] <- apply (x[,3:ncol(x)],2, as.numeric)
  l <- NULL
  coefsOut <- list()
  k <- 1
  j <- 1

  for (i in seq(1:dim(dff)[1])) { #for each formula in the table
    f <- dff[i, ] #extract formula
    datas <- x #extract source
    datas <- datas[which(datas[, 2] == as.character(f[[8]])), ] #extract class
    
    datas[,3:ncol(datas)] <- apply (datas[,3:ncol(datas)],2, as.numeric) #make sure numeric
    vals <- strsplit(as.character(f[[6]]), "~") #split formula to get conversion of the initial concentration for src
    
    #print ('attempt to evaluate')
    element <- eval((parse(text = vals[[1]][1]))) #convert src element concentration
    
    if (length(corFor) > 1){
      var1 <- datas[,corFor[1]] # get source data (x)
      var2 <- datas[,corFor[2]] # get source data (y)

      #will be helpful in the future
      fit <- lm(eval((parse(text = as.character(f[[6]]))))) #evaluate equation
      #get conversions needed
      vals_n <- vals[[1]][2]
      vals_n <- strsplit(vals_n, " ")
      size_f <- vals_n[[1]][grep("var1", vals_n[[1]])]
      toc_f <- vals_n[[1]][grep("var2", vals_n[[1]])]
      
      #convert source data, size
      if (length(var1) == 0) {
      } else {
        sourceSize <- eval(parse(text = size_f))
      }
  
      #convert source data, toc
      if (length(var2) == 0) {} else {
        sourceToc <- eval(parse(text = toc_f))
      }
      
      var1 <- as.numeric(trg[corFor[1]]) #use values of only first target
      var2 <- as.numeric(trg[corFor[2]]) 
      
      targetSize <- eval(parse(text = size_f)) #convert them
      targetToc <- eval(parse(text = toc_f)) #convert values
  
      # Drop intercept if there any from coefficients
      coefs <- coefficients(fit)[names(coefficients(fit)) != "(Intercept)"]
      coefs <- as.numeric(coefs)
      Yi <- element
      Si <- sourceSize
      Ti <- sourceToc
      Sj <- targetSize
      Tj <- targetToc
    } else {
      var1 <- datas[,corFor] # get source data (x)
      
      #will be helpful in the future
      fit <- lm(eval((parse(text = as.character(f[[6]]))))) #evaluate equation
      #get conversions needed
      vals_n <- vals[[1]][2]
      vals_n <- strsplit(vals_n, " ")
      size_f <- vals_n[[1]][grep("var1", vals_n[[1]])]
      
      sourceSize <- eval(parse(text = size_f))
      var1 <- as.numeric(trg[corFor]) #use values of only first target
      targetSize <- eval(parse(text = size_f)) #convert them
      # Drop intercept if there any from coefficients
      coefs <- coefficients(fit)[names(coefficients(fit)) != "(Intercept)"]
      coefs <- list(as.numeric(coefs))
      
      Yi <- element
      Si <- sourceSize
      Ti <- NULL
      Sj <- targetSize
      Tj <- NULL
      
    }

    
    if (length(coefs) == 3) {
      Ss <- coefs[1] # Size slope
      Ts <- coefs[2] # Toc slope
      TjSj <- coefs[3] # Interaction slope
      Corrected <- (Yi - (((Si - Sj) * Ss) + ((Ti - Tj) * Ts) + (((Si * Ti) - (Sj * Tj)) * TjSj)))
    } else if (length(coefs) == 2) {
      if (length(size_f) == 0) { # toc and interaction
        Ts <- coefs[1]
        TjSj <- coefs[2]
        Corrected <- (Yi - (((Ti - Tj) * Ts)) + (((Si * Ti) - (Sj * Tj)) * TjSj))
      }
      else if (length(toc_f) == 0) { # size and interaction
        Ss <- coefs[1]
        TjSj <- coefs[2]
        Corrected <- (Yi - (((Si - Sj) * Ss)) + (((Si * Ti) - (Sj * Tj)) * TjSj))
      } else { # size and toc only
        Ss <- coefs[1]
        Ts <- coefs[2]
        Corrected <- (Yi - (((Si - Sj) * Ss) + ((Ti - Tj) * Ts)))
      }
    } else if (length(coefs) == 1) {
      if (length(size_f) == 0) {
        Ts <- coefs[[1]]
        Corrected <- (Yi - ((Ti - Tj) * Ts))
      } # toc
      else {
        Ss <- coefs[[1]]
        Corrected <- (Yi - ((Si - Sj) * Ss))
      }
    }
    
    dat <- as.character(f[6][1][[1]])
    to_listFormula <- as.list(scan(text = dat, what = "", sep = " "))
    
    varNames <- paste0(names(datas), collapse = "|")
    varElement <- gsub(varNames, "x", to_listFormula[[1]])
    
    initialConc <- datas[, which(names(datas) == as.character(f[1][[1]]))]
    
    options(warn=-1)
    
    if (grepl("x^2", varElement, fixed = T)) {
      Corrected <- sqrt(Corrected)
    } else if (grepl("x^(1/3)", varElement, fixed = T)) {
      Corrected <- Corrected^3
    } else if (grepl("x^(-1)", varElement, fixed = T)) {
      Corrected <- Corrected * initialConc
    } else if (grepl("x^(-1/2)", varElement, fixed = T)) {
      Corrected <- Corrected^2 * initialConc
    } else if (grepl("log10(x)", varElement, fixed = T)) {
      Corrected <- exp(Corrected)
    }
    options(warn=0)

    if (any(is.na(as.numeric(Corrected)))) {
      Corrected <- initialConc
      l[j] <- f[[6]]
      j<-j+1
    }
    
    if (any(Corrected < 0)) {
      Corrected <- initialConc
      l[j] <- f[[6]]
      j<-j+1
    }
    coefsOut[[k]] <- coefs
    k <- k +1
    x[which(x[, 2] == as.character(f[[8]])), which(names(x) == as.character(f[[1]])) ] <- Corrected
  }
  return(list(data.frame(x), coefsOut, l))
}


stepwiseDFA <- function(l) {
  DFA_l <- list()
  for (iter in seq(1, length(l))) {
    
    
    # read in data
    datas <- l[[iter]]
    rownames(datas) <- datas[, 1]
    datas <- datas[, -1]
    
    options(warn = -1)
    
    suppressMessages(attach(datas))
    
    options(warn = 0)
    
    
    # perform stepwise greedy wilks classification
    gw_obj <- greedy.wilks(SourceType ~ ., data = datas, grouping = SourceType, niveau = 0.1)
    
    gw_results <- as.data.frame(gw_obj$results)
    
    #print (formula(gw_obj))
    
    
    # equal priors#
    p <- 1 / length(unique(SourceType))
    p <- matrix(p, ncol = 1, nrow = length(unique(SourceType)))
    priors <- c(p)
    
    
    # LDA#
    testfit <- try({
      fit <- do.call(lda, list(
        formula = formula(gw_obj),
        data = quote(datas), prior = priors, CV = TRUE, tol = 0.001
      ))
    })
    
    
    # fit formulas
    fit <- do.call(lda, list(
      formula = formula(gw_obj, .env),
      data = quote(datas), prior = priors, CV = TRUE, tol = 0.001
    ))
    
    ## extract the formula in the environment of cv.step
    form <- as.formula(gw_obj$formula) 

    f <- (as.character(form))

    
    mod <- lda(formula(f), data = datas, prior = priors, tol = 0.001)
    pred <- as.data.frame(predict(mod)$x)
    
    
    # ggplot (biplot)
    lda.data <- cbind(datas, predict(mod)$x)
    ggplot(lda.data, aes(LD1, LD2, LD3)) +
      geom_point(aes(color = SourceType))
    
    
    
    # function for confusion matrix#
    errorRate <- function(object, ...) {
      if (!require(MASS)) stop("you need the MASS package installed")
      UseMethod("errorRate")
    }
    
    errorRate.lda <- function(object, data = eval.parent(object$call$data),
                              type = "plug-in") {
      pred <- predict(object, data, type = type)$class
      actu <- eval(formula(object)[[2]], data)
      conf <- table(pred, actu)
      1 - sum(diag(conf)) / sum(conf)
      #print(conf)
    }
    
    # confusion matrix and percentgood#
    confusionCVpr <- table(fit$class, SourceType)
    ConfusionMatrix <- as.data.frame(matrix(confusionCVpr, ncol = length(unique(SourceType))))
    
    
    # estimate percentage thats good
    percentgood <- as.data.frame(matrix(0, ncol = 1, nrow = 0))
    for (i in 1:length(ConfusionMatrix)) {
      good <- ConfusionMatrix[[i, i]]
      all <- as.data.frame(matrix(apply(ConfusionMatrix, 2, function(x) sum(x)), ncol = length(ConfusionMatrix)))
      percent <- good / all[[i]]
      percentgood <- rbind(percentgood, percent)
    }
    
    
    # percent good for all variables#
    percentgood <- sapply(percentgood, function(x) mean(x) * 100)
    percentgood <- as.data.frame(percentgood)
    names(percentgood) <- "PercentClassified"
    
    
    # results from greedy.wilks#
    results <- gw_obj$results
    
    # for Pis and Wi#
    Pitable <- as.data.frame(matrix(0, ncol = 2, nrow = 0))
    names(Pitable) <- c("Tracer", "Pi")
    
    # The big loop
    for (e in 1:(length(results[[1]]))) {
      indvformula <- paste("lda(SourceType~", as.character(results[[e, 1]]), sep = "")
      indvformula <- paste(indvformula, ",data=datas,prior=priors,CV=TRUE,tol=", 0.001, ")")
      #print(indvformula)
      # Sys.sleep(15)
      fitindv <- eval(parse(text = indvformula))
      confusionCVprINDV <- table(fitindv$class, SourceType)
      ConfusionMatrixINDV <- as.data.frame(matrix(confusionCVprINDV, ncol = length(unique(SourceType))))
      percentgoodINDV <- as.data.frame(matrix(0, ncol = 1, nrow = 0))
      for (j in 1:length(ConfusionMatrixINDV)) {
        good <- ConfusionMatrixINDV[[j, j]]
        all <- as.data.frame(matrix(apply(ConfusionMatrixINDV, 2, function(x) sum(x)), ncol = length(ConfusionMatrixINDV)))
        percent <- good / all[[j]]
        percentgoodINDV <- rbind(percentgoodINDV, percent)
      }
      # percent good for one variable#
      Pi <- sapply(percentgoodINDV, function(x) round(mean(x), digits = 3))
      Pi <- cbind.data.frame(results[[e, 1]], Pi)
      names(Pi)[1] <- "Tracer"
      Pitable <- rbind(Pitable, Pi)
    }
    
    # calculate weights#
    Popt <- min(Pitable[[2]])
    
    W <- as.data.frame(sapply(Pitable[c(2)], function(x) round(x / Popt, digits = 3)))
    
    
    names(W)[1] <- "Wi"
    Pitable <- cbind(Pitable, W)
    Pitable <- Pitable[order(-Pitable$Wi), ]
    Pitable$Wi <- Pitable$Wi * 100
    DFA_l[[iter]] <- Pitable
  }
  
  return(DFA_l)
}


bracketT<- function (x,y, r, nrm) {
  drops <- c(nrm)
  #DF[ , !(names(DF) %in% drops)]
  
  d <- y
  y <- y[,-c(1,2)]
  x <- x[, 3:ncol(x)]
  x <- apply(x, 2, as.numeric)

  #which(!names(datas) %in% names(x[[1]])
  x <- as.data.frame(x)
  x <- x[,!(names(x) %in% drops)]
  y <- y[,!(names(y) %in% drops)]
  
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
  dat <- cbind (d[,c(1,2)], yNum)
  dat <- cbind (dat, d[,!(names(d) %in% names(dat))])
  
  dat <- cbind(d[, c(1, 2)], yNum) #take first and second column and combine them
  list(dat, c(na.omit(unlist(l))))
}

grepL <- function(x) {
  return(grep("*", x, fixed = T))
}
