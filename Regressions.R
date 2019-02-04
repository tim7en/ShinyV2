# #* -------------------------------------------------------------------*/
# rm(list = ls())
# source("func.R")
#
# datainput <- read.csv("Source_Allen.csv")
# targetData <- read.csv("Targets_2.csv")
#
# sizeOnly <- NULL
# # Make sure data is data frame
# datainput <- as.data.frame(datainput) # Input source data
# target <- as.data.frame(targetData) # Input target data
#
# datainput[, 3:ncol(datainput)] <- apply(datainput[, 3:ncol(datainput)], 2, as.numeric)
#
# sizeonly <- colnames(datainput)[which(colnames(datainput)[-c(1, 2)] %in% sizeOnly) + 2] # Unselected tracers will be used only with size correction and not TOC
# R2 <- 0.5 # Slider input of R2 value used to assess correlation
# pVal <- 0.05 # Slider input of p value for normality test used to asssess normality of residuals
# constant <- 50

applycorrection <- function(datainput, target, sizeonly, R2, pVal, constant) {
  datainput <- as.data.frame(datainput) # Input source data
  datainput[, 3:ncol(datainput)] <- apply(datainput[, 3:ncol(datainput)], 2, as.numeric)
  target <- as.data.frame(target) # Input target data
  sizeonly <- colnames(datainput)[which(colnames(datainput)[-c(1, 2)] %in% sizeonly) + 2]

  Output_list <- list() # Final output list
  Onlyneg_list <- vector("list", length = dim(target)[1])
  Rsq_list <- list()
  Rsq_output <- list()
  pval_list <- vector("list", length = dim(target)[1])
  Pval_output <- list()
  finalvars_list <- list()
  modeleq_list <- list()
  modeleq_output <- list()


  for (k in seq(1, dim(target)[1])) {
    target_adj <- target[k, ] # Select element, target 1
    Result <- as.data.frame(matrix(nrow = 0, ncol = dim(datainput)[2])) # Create matrix to store results
    names(Result) <- names(datainput) # Give names of columns from input to newly created matrix

    # Make clusters
    cl <- makeCluster(getOption("cl.cores", detectCores() - 1))

    # convertNeg <- function(y) {
    #   if (any(y > 0) == TRUE) {
    #     y <- as.numeric(y) + constant
    #   } else {
    #     y <- y * -1
    #   }
    # }

    for (i in seq(1:length(unique(datainput[, 2])))) {
      input <- datainput[which(datainput[, 2] == unique(datainput[, 2])[i]), ] # Subset unique source
      if (any(input[, 3:dim(input)[2]] < 0) == TRUE) {
        negatives <- input[, unique(which(input[, 3:ncol(input)] < 0, arr.ind = T)[, 2]) + 2]
        negativecolumns <- list()

        j <- 1
        if (class(negatives) == "data.frame") {
          for (i_2 in seq(1, length(colnames(negatives)))) {
            if (any(negatives[, i_2] > 0) == TRUE) {
              negativecolumns[j] <- colnames(negatives)[i_2]
              j <- j + 1
            } else {

            }
          }
        } else if (class(negatives) == "numeric") {
          if (any(negatives > 0) == TRUE) {
            negativecolumns <- colnames(input[, unique(which(input[, 3:ncol(input)] < 0, arr.ind = T)[, 2]) + 2])
          } else {}
        }

        negativecolumns <- as.character(unlist(negativecolumns))
        if (length(negativecolumns) > 0) {
          Onlyneg_list[k] <- negativecolumns
          names(Onlyneg_list)[k] <- as.character(target_adj [, 1])
        }


        if (class(negatives) == "data.frame") {
          negatives <- apply(negatives, 2, convertNeg)
        } else {
          negatives <- convertNeg(negatives)
        }
        input[, unique(which(input[, 3:ncol(input)] < 0, arr.ind = T)[, 2]) + 2] <- negatives
      } else {}


      # Define dimensions of the unique source group data frame (nrows, ncols)
      dfdim <- dim(input)

      # Select column 3 and store it into variable size
      size <- as.numeric(input[, 3])

      # Select column 4 and store it into variable toc
      toc <- as.numeric(input[, 4])

      # Select tracers that are not in size only
      tracers <- which(!colnames(input) %in% sizeonly)

      # Select all columns that are above index 5, since we have 1-4 indexes in data frame for source id, source type, d50 , toc
      tracers <- tracers[tracers >= 5]

      # Load library MASS in clusters
      clusterEvalQ(cl, library(MASS))
      clusterExport(cl = cl, varlist = c("input", "size", "toc", "Result"), envir = environment())

      # Apply lm_func over selected columns of data frame (regs - regressions)
      regs <- parCapply(cl, input[, tracers], lm_func)

      # If there any size only selected element, apply size only function
      if (length(sizeonly) > 0) {
        if (length(sizeonly == 1)) {
          regsiso <- lm_size(input[, which(colnames(input) %in% sizeonly)])
          regsiso$call$formula <- as.formula(paste(sizeonly, "~size", sep = ""))
        } else {
          regsiso <- parCapply(cl, input[, which(colnames(input) %in% sizeonly)], lm_size)
          for (j in seq(1, length(sizeonly))) {
            regsiso[[j]]$call$formula <- as.formula(paste(names(regsiso)[j], "~size", sep = ""))
          }
        }
      } else {
        regsiso <- NULL
      }

      elementnames <- names(input [5:dfdim[2]])

      for (j in seq(1, length(regs))) {
        regs[[j]]$call$formula <- as.formula(paste(names(regs)[j], "~size*toc", sep = ""))
      }

      if (length(sizeonly) > 1) {
        regs <- c(regs, regsiso)
        names(regs)[length(regs)] <- sizeonly
      } else if (length(sizeonly) == 1) {
        library(rlist)
        regs <- list.append(regs, regsiso)
        names(regs)[length(regs)] <- sizeonly
      }


      # Do I need this step ? or I can just apply boxcox
      regsstep <- parLapply(cl, regs, step_lm)

      # BOXCOX transformation for linear models
      boxcoxtracers <- parLapply(cl, regsstep, lm_boxcox)
      # boxcoxtracers <- lapply (regs, boxcox)

      lambdatracers <- parLapply(cl, boxcoxtracers, getlambda)
      # lambdatracers <- lapply (boxcoxtracers, getlambda)

      if (length(sizeonly) > 0) {
        if (length(sizeonly) == 1) {
          boxcoxiso <- lm_boxcox(regsiso)
          lambdaiso <- getlambda(boxcoxiso)
        } else {
          boxcoxiso <- parLapply(cl, regsiso, lm_boxcox)
          lambdaiso <- parLapply(cl, boxcoxiso, getlambda)
        }
      } else {
        boxcoxiso <- NULL
        lambdaiso <- NULL
      }

      l <- c(names(lambdatracers), sizeonly)
      lambdares <- c(lambdatracers, lambdaiso)
      names(lambdares) <- l

      tracersloc <- which(!colnames(input) %in% sizeonly) # toc and size
      tracersloc <- tracersloc[which(tracersloc >= 5)]

      inputtransformed <- input
      dat <- as.matrix(inputtransformed[, tracersloc])
      dat <- apply(dat, 2, as.numeric)

      inputtransformed[, tracersloc] <- sapply(1:length(tracersloc), function(col) powerTransform(dat[, col], as.numeric(lambdatracers)[col]))

      clusterExport(cl = cl, varlist = c("input", "size", "toc", "Result"), envir = environment())

      # Apply mixing model to the tracers (not isotopes)
      regsnew <- parCapply(cl, inputtransformed[, tracersloc], lm_func)
      regsstep <- parLapply(cl, regsnew, step_lm)


      if (length(sizeonly) > 1) {
        regsstep <- c(regsstep, regsiso)
      } else if (length(sizeonly) == 1) {
        library(rlist)
        regsstep <- list.append(regsstep, regsiso)
        names(regsstep)[length(regsstep)] <- sizeonly
      }

      stepssummary <- parLapply(cl, regsstep, summary)
      # stepssummary <- lapply (regsstep, summary)
      coefs_r <- parLapply(cl, stepssummary, extractrsquared) #R2
      # coefs_r <- lapply (stepssummary, extractrsquared)
      coefs_r <- unlist(coefs_r)
      cors <- names(which(coefs_r > R2))
      residualsval <- NULL

      Rsq_list[[i]] <- coefs_r[which(coefs_r > R2)]
      names(Rsq_list)[i] <- as.character(unique(datainput[, 2])[i])


      # Check if the any variable above 0.5
      if (length(cors) == 0) {
        input[, which(colnames(input) %in% negativecolumns)] <- input[, which(colnames(input) %in% negativecolumns)] - constant
      } else {
        arr <- array(which(coefs_r > R2))
        for (j in seq(1, length(cors))) {
          length(regsstep)
          residualsval <- cbind(residualsval, regsstep[[arr[j]]]$residuals)
        }

        colnames(residualsval) <- cors

        shapiro_residuals <- tryCatch({
          parCapply(cl, residualsval, shapiro.test)
          # apply(residualsval, 2, shapiro.test)
        },
        error = function(cond) {
          shapiro_residuals <- apply(residualsval, 2, shapiro.test)
          return(shapiro_residuals)
        },
        finally = {
        }
        )

        p_vals <- parLapply(cl, shapiro_residuals, extractpval)
        # p_vals <- lapply (shapiro_residuals, extractpval)
        p_vals <- unlist(p_vals)
        # P value of residuals
        pval_list[[i]] <- p_vals
        names(pval_list)[i] <- as.character(unique(datainput[, 2])[i])

        drops <- which(p_vals < pVal)
        finalvars <- cors[which(!cors %in% names(drops))]
        finalvars_list[[i]] <- finalvars
        modslopes <- regsstep[which(names(regsstep) %in% finalvars)]

        modelslopes <- tryCatch({
          coefs <- parLapply(cl, modslopes, extractcoeff)
          # coefs <- lapply (modslopes, extractcoeff)
          modeleq <- parLapply(cl, coefs, getlmeq)
          # modeleq <- lapply (coefs, getlmeq)
          modelslopes <- parLapply(cl, coefs, onlyslopes)
          modeleq_list[[i]] <- unlist(modelslopes)
          names(modeleq_list)[i] <- as.character(unique(datainput[, 2])[i])
          # modelslopes <- lapply (coefs, onlyslopes)
        },
        error = function(cond) {
          coefs <- lapply(modslopes, extractcoeff)
          modeleq <- lapply(coefs, getlmeq)
          modelslopes <- lapply(coefs, onlyslopes)
          return(modelslopes)
        },
        finally = {
          inputmod <- inputtransformed[, which(names(inputtransformed) %in% finalvars)]
          lambdamod <- lambdares[which(names(lambdares) %in% finalvars)]
        }
        )

        if (any(dim(inputmod) == 0)) {
          outputcorrected <- as.data.frame(matrix(nrow = dim(inputmod)[1], ncol = dim(inputmod)[2]))
          names(outputcorrected) <- colnames(inputmod)
        } else if (is.null(dim(inputmod)) && length(inputmod) > 0) {
          outputcorrected <- as.data.frame(matrix(nrow = length(inputmod), ncol = 1))
          names(outputcorrected) <- finalvars
        } else {
          outputcorrected <- as.data.frame(matrix(nrow = dim(inputmod)[1], ncol = dim(inputmod)[2]))
          names(outputcorrected) <- finalvars
        }

        if (any(dim(inputmod) == 0)) {}
        else {
          # print (finalvars)
          for (j in seq(1, dim(outputcorrected)[2])) {
            varname <- names(outputcorrected)[j]
            len <- length(coefs[[which(names(coefs) %in% varname)]]) - 1
            if (is.null(dim(inputmod)) == FALSE) {
              Yi <- inputmod[, which(names(inputmod) %in% varname)] # Ini Concentration
            } else {
              Yi <- inputmod
            }
            Si <- size # Ini size
            Ti <- toc # Ini toc
            l <- lambdamod[which(names(lambdamod) %in% varname)] # lambda
            Sj <- target_adj[, 3] # Target size
            Tj <- target_adj [, 4] # Target toc

            if (len == 3) {
              Ss <- as.numeric(coefs[[which(names(coefs) %in% varname)]][-1][1]) # Size slope
              Ts <- as.numeric(coefs[[which(names(coefs) %in% varname)]][-1][2]) # Toc slope
              TjSj <- as.numeric(coefs[[which(names(coefs) %in% varname)]][-1][3]) # Interaction slope
              outputcorrected[, j] <- invBoxCox((Yi - (((Si - Sj) * Ss) + ((Ti - Tj) * Ts) + (((Si * Ti) - (Sj * Tj)) * TjSj))), as.numeric(l))
            } else if (len == 2) {
              comb1 <- c("size", "size:toc")
              comb2 <- c("toc", "size:toc")
              cnames <- names(coefs[[which(names(coefs) %in% varname)]])[-1]
              if (length(which(cnames %in% comb1) == TRUE) == 2) {
                Ss <- as.numeric(coefs[[which(names(coefs) %in% varname)]][-1][1])
                TjSj <- as.numeric(coefs[[which(names(coefs) %in% varname)]][-1][2])
                outputcorrected[, j] <- invBoxCox((Yi - (((Si - Sj) * Ss) + ((Si * Ti) - (Sj * Tj)) * TjSj)), as.numeric(l))
              } else if (length(which(cnames %in% comb2) == TRUE) == 2) {
                Ts <- as.numeric(coefs[[which(names(coefs) %in% varname)]][-1][1])
                TjSj <- as.numeric(coefs[[which(names(coefs) %in% varname)]][-1][2])
                outputcorrected[, j] <- invBoxCox((Yi - (((Ti - Tj) * Ts) + ((Si * Ti) - (Sj * Tj)) * TjSj)), as.numeric(l))
              }
              else {
                outputcorrected[, j] <- invBoxCox((Yi - (((Si - Sj) * Ss) + ((Ti - Tj) * Ts))), as.numeric(l))
              }
            } else {
              comb1 <- as.character("size")
              comb2 <- as.character("toc")
              cnames <- names(coefs[[which(names(coefs) %in% varname)]])[-1]
              if (length(which(cnames %in% comb1) == TRUE) == 1) {
                Ss <- as.numeric(coefs[[which(names(coefs) %in% varname)]][-1][1])
                outputcorrected[, j] <- invBoxCox((Yi - (((Si - Sj) * Ss))), as.numeric(l))
              }
              else if (length(which(cnames %in% comb2) == TRUE) == 1) {
                Ts <- as.numeric(coefs[[which(names(coefs) %in% varname)]][-1][1])
                outputcorrected[, j] <- invBoxCox((Yi - (((Ti - Tj) * Ts))), as.numeric(l))
              }
              else {
                TjSj <- as.numeric(coefs[[which(names(coefs) %in% varname)]][-1][1])
                outputcorrected[, j] <- invBoxCox((Yi - ((Si * Ti) - (Sj * Tj)) * TjSj), as.numeric(l))
              }
            }
          }
        }
        outputcorrected <- t(na.omit(t(outputcorrected)))
        outputcorrected <- as.data.frame(outputcorrected)
        # 3263
        input[, which(names(input) %in% colnames(outputcorrected))] <- outputcorrected
        input[, which(colnames(input) %in% negativecolumns)] <- input[, which(colnames(input) %in% negativecolumns)] - constant
      }
      Result <- rbind(Result, input)
    }

    Pval_output <- pval_list
    Rsq_output <- Rsq_list
    finalvars_output <- finalvars_list
    modeleq_output <- modeleq_list
    Output_list[[k]] <- Result
    stopCluster(cl)
    closeAllConnections()
  }
  return(list(Output_list, Rsq_output, Pval_output, finalvars_output, modeleq_output))
}
