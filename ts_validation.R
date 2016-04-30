TimeSeriesOutOfSampleValidation <- function(tsfrmla, tsdata, xregmat, methodname, k = 10, hrs.in.window = 168, logtransform = FALSE) {
  # Time series cross-validation loosely based on Rob Hyndman's example:
  # http://robjhyndman.com/hyndsight/tscvexample/
  #
  # Since November 1, 2011 is the TOU billing start in my data set, I want to be sure to 
  # include adequate amounts of training data on either side of that date. Because one of 
  # my models under consideration is a dynamic linear regression model with AR(1) error 
  # structure, I want to keep my out-of-sample timeseries forecasts in series with the 
  # rest of my data.
  #
  # This "rolling origin" time series validation function uses the tail observations as 
  # its out-of-sample candidates. A validation windows size is chosen (ie. number of observations 
  # used in a validation fold). By default the size of the window is 168 hours (ie. a week). 
  # The default number of validation folds is 10. These can both be changed.
  #
  # A series of observations floor(k * hrs.in.window) is reserved at the end of the sample 
  # for validation and held out of the fitting process. This validation function then uses 
  # the front portion of data floor(length(sample) - (k * hrs.in.window)) as the training set.
  # Forecasts of size floor(hrs.in.window) are then made for observations immediately following the 
  # training sample. Each iteration, the origin of the training sample rolls forward by floor(hrs.in.window)
  # and continues in this way until k validations have been run.
  #
  # Each iteration runs validation of a tslm model and a dynamic linear regression model with 
  # AR(1) errors.
  #
  # Args:
  #   tsfrmla: A formula(...) object to be used in the tslm function.
  #   tsdata: A ts(...) matrix suitable for use in the tslm function.
  #   xregmat: A matrix suitable for use in the Arima function's xreg parameter.
  #   methodname: A string used to identify the results, since calls to this function will happen 
  #               in batches. Used in figure labels and results directory structure.
  #   k: The number of validation window (default = 10)
  #   hrs.in.window: The number of observations in the validation window (default = 168, a week).
  #                  floor(...) within the body of the function handles this fraction well.
  #
  # Returns:
  #   A list with two named values (matrices) representing the linear regression validation results 
  #   and the dynamic regression validation results.  
  require(forecast)       # Package to model serially correlated error structure
  require(xlsx)           # Useful for exporting results to Excel
  
  trainsize <- floor(nrow(tsdata) - (k * hrs.in.window))
  trainsize.freq <- trainsize / 24
  windowsize.freq <- floor(hrs.in.window) / 24
  
  # Store accuracy results for each fold + average across all runs
  results.tslm <- results.ar1err <- matrix(NA, nrow = k+1, ncol = 10)
  colnames(results.tslm) <- colnames(results.ar1err) <- c("k", "trainsize", "trainstart", "trainend", 
                                                          "validsize", "validstart", "validend", 
                                                          "MAE", "MAPE", "ACF1")
  
  freq.step <- 1/24
  freqidx.start <- 1
  idx.start <- 1
  freqidx.stop <- trainsize.freq + 23/24
  idx.stop <- trainsize
  
  # Begin out-of-sample validation with rolling origin
  for(i in 1:k)
  {
    #####
    # Roll training/forecast window and fit models.
    #####
    
    # ts(...) version of data for tslm, with formula workaround so I don't have to type it all
    tsdata.train <- window(tsdata, start = freqidx.start, end = freqidx.stop)
    tsdata.valid <- window(tsdata, start = (freqidx.stop + freq.step), end = (freqidx.stop + windowsize.freq))
    
    # xreg data for Arima(...)
    xreg.train <- xregmat[c(idx.start:idx.stop),]
    xreg.valid <- xregmat[c((idx.stop+1):(idx.stop+floor(hrs.in.window))),]
    
    # Fit a standard multiple regression model using tslm
    tslm.trainfit <- tslm(formula = tsfrmla, data = tsdata.train)
    tslm.forecast <- forecast(object = tslm.trainfit, newdata = tsdata.valid)
    
    # Fit a dynamic regression model with AR(1) structured errors
    # frequency = 24 seems to fix, "Error in optim(init[mask], armaCSS, method = optim.method, hessian = FALSE,  : non-finite value supplied by optim "
    if (logtransform == TRUE) {
      ar1err.trainfit <- Arima(x = log(tsdata.train[,"kwh"]), xreg = xreg.train, order = c(1, 0, 0), include.drift = FALSE)
    } else {
      ar1err.trainfit <- Arima(x = tsdata.train[,"kwh"], xreg = xreg.train, order = c(1, 0, 0), include.drift = FALSE)
    }
    ar1err.forecast <- forecast(object = ar1err.trainfit, xreg = xreg.valid)
    
    
    #####
    # Whether fitting was done with logtransformed response or untransformed response,
    # I want to work with backtransformed fitted values and response.
    #####
    if (logtransform == TRUE) {
      tslm.forecast.ts <- exp(tslm.forecast$mean)
      tslm.residuals.ts <- tsdata.valid[,"kwh"] - tslm.forecast.ts
    } else {
      tslm.forecast.ts <- tslm.forecast$mean
      tslm.residuals.ts <- tsdata.valid[,"kwh"] - tslm.forecast.ts
    }
    
    if (logtransform == TRUE) {
      ar1err.forecast.ts <- exp(ar1err.forecast$mean)
      ar1err.residuals.ts <- tsdata.valid[,"kwh"] - ar1err.forecast.ts
    } else {
      ar1err.forecast.ts <- ar1err.forecast$mean
      ar1err.residuals.ts <- tsdata.valid[,"kwh"] - ar1err.forecast.ts
    }
    
    
    
    #####
    # Plot residuals over time and their density
    #####
    if (logtransform == TRUE) {
      tslm.resid.time.xlab <- "Backtransformed Residual (kWh)"
      tslm.resid.time.path <- paste0("./results/logtransform/", methodname, "/window/", floor(hrs.in.window), "/tslm_residuals", i, ".png")
    } else {
      tslm.resid.time.xlab <- "Residual (kWh)"
      tslm.resid.time.path <- paste0("./results/", methodname, "/window/", floor(hrs.in.window), "/tslm_residuals", i, ".png")
    }
    plot(tslm.residuals.ts,
         main = paste0("Validation: ", methodname, " Multiple Regression (k=", k, ", fold=", i, ", size=", floor(hrs.in.window), ")"),
         ylab = tslm.resid.time.xlab,
         xlab = "Validation Index (Day Number)",
         ylim = c(-2, 2),
         type = "p",
         pch = 16)
    dev.print(file = tslm.resid.time.path, device = png, height = 600, width = 800)
    
    if(logtransform == TRUE) {
      tslm.resid.density.title <- paste0("Validation Backtransformed Residuals Density: ", methodname, " Multiple Regression",
                                         "\n(k=", k, ", fold=", i, ")")
      tslm.resid.density.path <- paste0("./results/logtransform/", methodname, "/window/", floor(hrs.in.window), "/tslm_residual_density", i, ".png")
    } else {
      tslm.resid.density.title <- paste0("Validation Residuals Density: ", methodname, " Multiple Regression", 
                                         "\n(k=", k, ", fold=", i, ")")
      tslm.resid.density.path <- paste0("./results/", methodname, "/window/", floor(hrs.in.window), "/tslm_residual_density", i, ".png")
    }
    plot(density(as.vector(tslm.residuals.ts), bw = 0.06),
         main = tslm.resid.density.title,
         ylim = c(0,4),
         xlim = c(-2,2),
         lwd = 3)
    dev.print(file = tslm.resid.density.path, device = png, height = 600, width = 800)
    
    if (logtransform == TRUE) {
      ar1err.resid.time.xlab <- "Backtransformed Residual (kWh)"
      ar1err.resid.time.path <- paste0("./results/logtransform/", methodname, "/window/", floor(hrs.in.window), "/ar1err_residuals", i, ".png")
    } else {
      ar1err.resid.time.xlab <- "Residual (kWh)"
      ar1err.resid.time.path <- paste0("./results/", methodname, "/window/", floor(hrs.in.window), "/ar1err_residuals", i, ".png")
    }
    plot(ar1err.residuals.ts,
         main = paste0("Validation: ", methodname, " Multiple Regression with AR(1) Errors (k=", k, ", fold=", i, ", size=", floor(hrs.in.window), ")"),
         ylab = ar1err.resid.time.xlab,
         xlab = "Validation Day Index",
         ylim = c(-2, 2),
         type = "p",
         pch = 16)
    dev.print(file = ar1err.resid.time.path, device = png, height = 600, width = 800)
    
    if (logtransform == TRUE) {
      ar1err.resid.density.title <- paste0("Validation Backtransformed Residuals Density: ", methodname, " Mult. Reg. with AR(1) Errors", 
                                           "\n(k=", k, ", fold=", i, ")")
      ar1err.resid.density.path <- paste0("./results/logtransform/", methodname, "/window/", floor(hrs.in.window), "/ar1err_residual_density", i, ".png")
      
    } else {
      ar1err.resid.density.title <- paste0("Validation Residuals Density: ", methodname, " Mult. Reg. with AR(1) Errors", 
                                           "\n(k=", k, ", fold=", i, ")")
      ar1err.resid.density.path <- paste0("./results/", methodname, "/window/", floor(hrs.in.window), "/ar1err_residual_density", i, ".png")
    }
    plot(density(as.vector(ar1err.residuals.ts), bw = 0.06),
         main = ar1err.resid.density.title,
         ylim = c(0,4),
         xlim = c(-2,2),
         lwd = 3)
    dev.print(file = ar1err.resid.density.path, device = png, height = 600, width = 800)
    
    
    
    #####
    # Find the out-of-sample measures of accuracy
    #####
    tslm.accuracy <- accuracy(tslm.forecast.ts, x = tsdata.valid[,"kwh"])
    results.tslm[i,] <- c(i, (idx.stop - idx.start + 1), idx.start, idx.stop, 
                          floor(hrs.in.window), (idx.stop+1), (idx.stop+floor(hrs.in.window)), 
                          tslm.accuracy["Test set","MAE"], tslm.accuracy["Test set","MAPE"], tslm.accuracy["Test set","ACF1"])
    ar1err.accuracy <- accuracy(ar1err.forecast.ts, x = tsdata.valid[,"kwh"])
    results.ar1err[i,] <- c(i, (idx.stop - idx.start + 1), idx.start, idx.stop, 
                            floor(hrs.in.window), (idx.stop+1), (idx.stop+floor(hrs.in.window)), 
                            ar1err.accuracy["Test set","MAE"], ar1err.accuracy["Test set","MAPE"], ar1err.accuracy["Test set","ACF1"])
    
    # Sanity check, using array indeces to make sure I didn't screw of the time series frequency
    kwh.trimmed <- tsdata[,"kwh"]
    kwh.validation.vec <- kwh.trimmed[c((idx.stop+1):(idx.stop+hrs.in.window))]
    if (logtransform == TRUE) {
      tslm.forecast.vec <- exp(tslm.forecast$mean)
    } else {
      tslm.forecast.vec <- tslm.forecast$mean
    }
    tslm.accuracy.sanity <- accuracy(f = tslm.forecast.vec, x = kwh.validation.vec)
    print(paste0("TSLM Iteration #", i, " (Indexing/Frequency Check)"))
    print(tslm.accuracy.sanity["Test set",c(1:5)])
    print(tslm.accuracy["Test set",c(1:5)])
    print("...")
    if (logtransform == TRUE) {
      ar1err.forecast.vec <- exp(ar1err.forecast$mean)
    } else {
      ar1err.forecast.vec <- ar1err.forecast$mean
    }
    ar1err.accuracy.sanity <- accuracy(f = ar1err.forecast.vec, x = kwh.validation.vec)
    print(paste0("AR(1) Iteration #", i, " (Indexing/Frequency Check)"))
    print(ar1err.accuracy.sanity["Test set",c(1:5)])
    print(ar1err.accuracy["Test set",c(1:5)])
    print("======================================")
    
    # Move origin of training and validation time series
    idx.start <- idx.start + floor(hrs.in.window)
    idx.stop <- idx.stop + floor(hrs.in.window)
    freqidx.start <- freqidx.start + windowsize.freq
    freqidx.stop <- freqidx.stop + windowsize.freq
  }
  
  
  
  #####
  # Add row which averages accuracy measures over all folds
  #####
  results.tslm[k+1,] <- c("Average", NA, NA, NA, NA, NA, NA, 
                          mean(results.tslm[c(1:k),8]), mean(results.tslm[c(1:k),9]), mean(results.tslm[c(1:k),10]))
  results.ar1err[k+1,] <- c("Average", NA, NA, NA, NA, NA, NA, 
                            mean(results.ar1err[c(1:k),8]), mean(results.ar1err[c(1:k),9]), mean(results.ar1err[c(1:k),10]))
  
  
  #####
  # Write results to an Excel file
  #####
  if (logtransform == TRUE) {
    write.xlsx(results.tslm, paste0("./results/logtransform/", methodname, "/window/", floor(hrs.in.window), "/tslm_accuracy.xlsx"))
    write.xlsx(results.ar1err, paste0("./results/logtransform/", methodname, "/window/", floor(hrs.in.window), "/ar1err_accuracy.xlsx"))    
  } else {
    write.xlsx(results.tslm, paste0("./results/", methodname, "/window/", floor(hrs.in.window), "/tslm_accuracy.xlsx"))
    write.xlsx(results.ar1err, paste0("./results/", methodname, "/window/", floor(hrs.in.window), "/ar1err_accuracy.xlsx"))    
  }
  
  return(list("results.tslm" = results.tslm, "results.ar1err" = results.ar1err))
}