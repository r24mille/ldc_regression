TimeSeriesOutOfSampleValidation <- function(tsfrmla, tsdata, xregmat, respvec, k = 5, hrs.in.window = 730.484) {
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
  # used in a validation fold). By default the size of the window is 730.484 hours (ie. the number 
  # of hours in a month). The default number of validation folds is 5. These can both be changed.
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
  #   respvec: A response vector suitable for use in the Arima function's x parameter.
  #   k: The number of validation window (default = 5)
  #   hrs.in.window: The number of observations in the validation window (default = 730.484, a month).
  #                  floor(...) within the body of the function handles this fraction well.
  #
  # Returns:
  #   A list with two named values (matrices) representing the linear regression validation results 
  #   and the dynamic regression validation results.
  trainsize <- floor(length(respvec) - (k * hrs.in.window))
  
  # Store accuracy results for each fold
  results.tslm <- results.ar1err <- matrix(NA, nrow = k, ncol = 10)
  colnames(results.tslm) <- colnames(results.ar1err) <- c("k", "trainsize", "trainstart", "trainend", "validsize", "validstart", "validend", "MAE", "MAPE", "ACF1")
  
  idx.start <- 1
  idx.stop <- trainsize
  
  # Begin out-of-sample validation with rolling origin
  for(i in 1:k)
  {
    # ts(...) version of data for tslm, with formula workaround so I don't have to type it all
    tsdata.train <- window(tsdata, start = idx.start, end = idx.stop)
    tsdata.valid <- window(tsdata, start = (idx.stop+1), end = (idx.stop+floor(hrs.in.window)))
    
    # xreg data for Arima(...)
    xreg.train <- xregmat[c(idx.start:idx.stop),]
    kwh.train <- respvec[c(idx.start:idx.stop)]
    xreg.valid <- xregmat[c((idx.stop+1):(idx.stop+floor(hrs.in.window))),]
    kwh.valid <- respvec[c((idx.stop+1):(idx.stop+floor(hrs.in.window)))]
    
    # Fit a standard multiple regression model using tslm
    tslm.trainfit <- tslm(formula = tsfrmla, data = tsdata.train)
    tslm.forecast <- forecast(object = tslm.trainfit, newdata = tsdata.valid)
    
    # Fit a dynamic regression model with AR(1) structured errors
    ar1err.trainfit <- Arima(x = kwh.train, xreg = xreg.train, order = c(1, 0, 0), include.drift = FALSE)
    ar1err.forecast <- forecast(object = ar1err.trainfit, xreg = xreg.valid)
    
    # Find the out-of-sample measures of accuracy
    tslm.accuracy <- accuracy(tslm.forecast, x = tsdata.valid[,"kwh"])
    results.tslm[i,] <- c(i, (idx.stop - idx.start + 1), idx.start, idx.stop, floor(hrs.in.window), (idx.stop+1), (idx.stop+floor(hrs.in.window)), 
                     tslm.accuracy["Test set","MAE"], tslm.accuracy["Test set","MAPE"], tslm.accuracy["Test set","ACF1"])
    ar1err.accuracy <- accuracy(ar1err.forecast, x = ts(c(kwh.train, kwh.valid)))
    results.ar1err[i,] <- c(i, (idx.stop - idx.start + 1), idx.start, idx.stop, floor(hrs.in.window), (idx.stop+1), (idx.stop+floor(hrs.in.window)), 
                       ar1err.accuracy["Test set","MAE"], ar1err.accuracy["Test set","MAPE"], ar1err.accuracy["Test set","ACF1"])
    
    # Move origin of training and validation time series
    idx.start <- idx.start + floor(hrs.in.window)
    idx.stop <- idx.stop + floor(hrs.in.window)
  }
  
  return(list("results.tslm" = results.tslm, "results.ar1err" = results.ar1err))
}