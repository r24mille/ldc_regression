AdjustedRSquare <- function(r2, n, p) {
  # Because I want Adjusted R^2 from the timeseries model, I must define it 
  # myself. Also useful for getting an Adjusted R^2 from backtransformed 
  # residuals. I can't use values from the model itself.
  #
  # Args:
  #   r2: R^2 computed outside this function
  #   n: Number of observations
  #   p: Number of model parameters
  #
  # Return:
  #   Adjusted R^2
  adj.r2 <- 1 - ((1 - r2) * ((n - 1) / (n - p - 1)))
  return(adj.r2)
}

EvaluateTemperatureTransform <- function(datafr, temptransform, methodname, 
                                         trimsize = 0, logtransform = FALSE,
                                         nullmodel = FALSE, 
                                         includetemp = TRUE) {
  # Evaluate different transformations of the temperature explanatory variable. 
  # All other explanatory variables across methodologies remain fixed so that 
  # the effectiveness of temperature transformations can be evaluated.
  #
  # Args:
  #   datafr: The average household electricity demand dataframe, derived from 
  #           aggregate data.
  #   temptransform: The temperature transformation matrix.
  #   methodname: The temperature transformation methodology name, for use in 
  #               filesystem paths and figure titles.
  #   trimsize: Several transformations work with lagged temperature values, 
  #             requiring that the first few observations be trimmed so that 
  #             data used for fitting is full rank. This parameter controls how 
  #             many rows at the head of the dataframe will be trimmed before 
  #             model fitting.
  #   logtransform: Boolean indicating whether the model evaluated should 
  #                 logtransform the response (ie. log(kWh)). The default is 
  #                 FALSE.
  require(scales)         # For scaling axes (eg. datetime)
  require(forecast)       # Used for AR(1) structured errors
  require(ggplot2)        # Better plotting
  require(xlsx)           # Useful for exporting results to Excel
  require(lmtest)         # Durban-Watson test for serially correlated errors
  source("ts_validation.R")
  
  #####
  # Create matrices suitable for cbinding into a AR(1) xreg matrix
  #####
  hrstr_dummies <- model.matrix(~hrstr, data = datafr)[,-1]
  #price_dummies <- model.matrix(~price, data = datafr)[,-1]
  rateseason_dummies <- model.matrix(~rateseason, data = datafr)[,-1]
  
  # A little extra work to ensure the proper contrasts are used. All occurrances 
  # of hrstr0 and working_dayFALSE can be removed because those effects are 
  # captured by other contrasts.
  #
  # TODO(r24mille): Trimming out working_day:hrstr interactions that do not have p-value < 0.05
  #                 reported by Newey-West. However, rather than working with the fickle p-value,
  #                 I need to investigate what the effects are on Adj. R^2, BIC, and validation.
  #                 Was trimmed to -c(1:27,30:31,44:46)
  #
  # Also strip ":" out of column names so they can be parsed into a formula object 
  # without conflict.
  hrwrkdy_dummies <- model.matrix(~hrstr:working_day, data = datafr)[,-c(1:26)]
  colnames(hrwrkdy_dummies) <- gsub(":", "", colnames(hrwrkdy_dummies))
  # hrtou_dummies <- model.matrix(~hrstr:tou_active, data = datafr)[,-c(1:26)]
  # colnames(hrtou_dummies) <- gsub(":", "", colnames(hrtou_dummies))
  # hrseason_dummies <- model.matrix(~hrstr:rateseason, data = datafr)[,-c(1:26)]
  # colnames(hrseason_dummies) <- gsub(":", "", colnames(hrseason_dummies))
  # wrkdytou_dummies <- model.matrix(~working_day:tou_active, data = datafr)[,-c(1:4)]
  # wrkdyseason_dummies <- model.matrix(~working_day:rateseason, data = datafr)[,-c(1:4)]
  # touseason_dummies <- model.matrix(~tou_active:rateseason, data = datafr)[,-c(1:4)]
  
  
  
  #####
  # Create matrices suitable for ARIMA modeling
  #####
  ts_xvars <- cbind(datafr$working_day, 
                    datafr$tou_active, 
                    rateseason_dummies, 
                    # wrkdytou_dummies, 
                    # wrkdyseason_dummies, 
                    # touseason_dummies, 
                    hrstr_dummies, 
                    hrwrkdy_dummies
                    # hrtou_dummies,
                    # hrseason_dummies, 
  )
  # Leaving interactions out of chapter 4:
  # "working_dayTRUEtou_activeTRUE", "working_dayTRUErateseasonsummer", 
  # "tou_activeTRUErateseasonsummer"
  colnames(ts_xvars)[c(1:3)] <- c("working_day", "tou_active", 
                                  "rateseasonsummer")
  if (includetemp == TRUE) {
    if(ncol(temptransform) > 1) {
      ts_xvars <- cbind(ts_xvars, temptransform[,])
    } else {
      ts_xvars <- cbind(ts_xvars, temptransform)
    }
  }
  
  if (trimsize > 0) {
    ts_xvars.trimmed <- tail(ts_xvars, -trimsize)
  } else {
    ts_xvars.trimmed <- ts_xvars
  }
  
  
  
  #####
  # Create ts data suitable for use in tslm
  #####
  if (trimsize > 0) {
    tsdata.trimmed <- ts(data = cbind(tail(datafr$kwh, -trimsize), ts_xvars.trimmed), frequency = 24)
  } else {
    tsdata.trimmed <- ts(data = cbind(datafr$kwh, ts_xvars.trimmed), frequency = 24)
  }
  colnames(tsdata.trimmed)[1] <- c("kwh")
  
  
  #####
  # Total sum of squares for use later
  #####
  SStot <- sum((tsdata.trimmed[,"kwh"] - mean(tsdata.trimmed[,"kwh"]))^2)
  
  
  #####
  # First create a linear model with no error structure
  #####
  ts_xvar.names <- colnames(tsdata.trimmed)
  if(nullmodel == FALSE) {
    if (logtransform == TRUE) {
      frmla.workaround <- as.formula(paste("log(kwh) ~", 
                                           paste(ts_xvar.names[!ts_xvar.names %in% "kwh"], 
                                                 collapse = " + ")))
    } else {
      frmla.workaround <- as.formula(paste("kwh ~", 
                                           paste(ts_xvar.names[!ts_xvar.names %in% "kwh"], 
                                                 collapse = " + ")))
    }
  } else {
    if (logtransform == TRUE) {
      frmla.workaround <- as.formula("log(kwh) ~ 1")
    } else {
      frmla.workaround <- as.formula("kwh ~ 1")
    }
  }
  multreg <- tslm(formula = frmla.workaround, data = tsdata.trimmed)
  
  # Whether fitting was done with logtransformed response or untransformed response,
  # I want to work with backtransformed fitted values and response.
  if (logtransform == TRUE) {
    multreg.fitted.ts <- exp(multreg$fitted)
    multreg.resids.ts <- tsdata.trimmed[,"kwh"] - multreg.fitted.ts
  } else {
    multreg.fitted.ts <- multreg$fitted
    multreg.resids.ts <- multreg$residuals
  }
  
  # ACF of multiple regression residuals
  if (logtransform == TRUE) {
    multreg.acf.title <- paste0("Backtransformed Residuals ACF Plot: ", methodname, " Multiple Regression")
    multreg.acf.path <- paste0("./results/logtransform/", methodname, "/tslm_residuals_acf.png")
  } else {
    multreg.acf.title <- paste0("Residuals ACF Plot: ", methodname, " Multiple Regression")
    multreg.acf.path <- paste0("./results/", methodname, "/tslm_residuals_acf.png")
  }
  acf(multreg.resids.ts, main = multreg.acf.title)
  dev.print(file = multreg.acf.path, device = png, height = 600, width = 800)
  
  # Multiple regression residuals plotted over time
  if (logtransform == TRUE) {
    multreg.resid.time.title <- paste0("Backtransformed Residuals by Time: ", methodname, " Multiple Regression")
    multreg.resid.time.path <- paste0("./results/logtransform/", methodname, "/tslm_residuals_by_index.png")
  } else {
    multreg.resid.time.title <- paste0("Residuals by Time: ", methodname, " Multiple Regression")
    multreg.resid.time.path <- paste0("./results/", methodname, "/tslm_residuals_by_index.png")
  }
  if(trimsize > 0) {
    multreg.resid.time.datafr <- data.frame("posixtime" = tail(datafr, -trimsize)$timestamp_dst, 
                                            "resids" = as.vector(multreg.resids.ts))
  } else {
    multreg.resid.time.datafr <- data.frame("posixtime" = datafr$timestamp_dst, 
                                            "resids" = as.vector(multreg.resids.ts))
  }
  multreg.resid.idx.plot <- ResidualsOverTimePlot(thedata = multreg.resid.time.datafr, thetitle = multreg.resid.time.title)
  print(multreg.resid.idx.plot)
  dev.print(file = multreg.resid.time.path, device = png, height = 600, width = 800)
  
  # Multiple regression residuals plotted as a function of the fitted value
  if (logtransform == TRUE) {
    multreg.resid.resp.title <- paste0("Backtransformed Residuals: ", methodname, " Multiple Regression")
    multreg.resid.resp.path <- paste0("./results/logtransform/", methodname, "/tslm_residuals_by_fitted_response.png")
  } else {
    multreg.resid.resp.title <- paste0("Residuals: ", methodname, " Multiple Regression")
    multreg.resid.resp.path <- paste0("./results/", methodname, "/tslm_residuals_by_fitted_response.png")
  }
  multreg.resid.resp.plot <- ResidualsByResponseVal(thedata = data.frame("fittedvals" = multreg.fitted.ts, 
                                                                         "resids" = as.vector(multreg.resids.ts)), 
                                                    thetitle = multreg.resid.resp.title)
  print(multreg.resid.resp.plot)
  dev.print(file = multreg.resid.resp.path, device = png, height = 600, width = 800)
  
  # Multiple regression residuals plotted as a function of temperature
  if (logtransform == TRUE) {
    multreg.resid.temp.title <- paste0("Backtransformed Residuals as a Function of Observed Temperature:", 
                                       "\n", methodname, " Multiple Regression")
    multreg.resid.temp.path <- paste0("./results/logtransform/", methodname, "/tslm_residuals_by_temperature.png")
  } else {
    multreg.resid.temp.title <- paste0("Residuals as a Function of Observed Temperature:", 
                                       "\n", methodname, " Multiple Regression")
    multreg.resid.temp.path <- paste0("./results/", methodname, "/tslm_residuals_by_temperature.png")
  }
  if(trimsize > 0) {
    multreg.resid.temp.datafr <- data.frame("temperatures" = tail(datafr, -trimsize)$temperature, 
                                            "resids" = as.vector(multreg.resids.ts))
  } else {
    multreg.resid.temp.datafr <- data.frame("temperatures" = datafr$temperature, 
                                            "resids" = as.vector(multreg.resids.ts))
  }
  multreg.resid.temp.plot <- ResidualsByTemperature(thedata = multreg.resid.temp.datafr, thetitle = multreg.resid.temp.title)
  print(multreg.resid.temp.plot)
  dev.print(file = multreg.resid.temp.path, device = png, height = 600, width = 800)
  
  # Density plot of multiple regression residuals
  if (logtransform == TRUE) {
    multreg.density.title <- paste0("Backtransformed Residuals Density: ", methodname, " Multiple Regression")
    multreg.density.path <- paste0("./results/logtransform/", methodname, "/tslm_residual_density.png")
  } else {
    multreg.density.title <- paste0("Residuals Density: ", methodname, " Multiple Regression")
    multreg.density.path <- paste0("./results/", methodname, "/tslm_residual_density.png")
  }
  plot(density(as.vector(multreg.resids.ts), bw = 0.015),
       main = multreg.density.title,
       ylim = c(0, 12), 
       xlim = c(-2, 2),
       lwd = 3)
  dev.print(file = multreg.density.path, device = png, height = 600, width = 800)
  
  # Store multiple regression measures of explanatory power in an Excel file
  multreg.RSS <- sum(as.vector(multreg.resids.ts)^2)
  multreg.r2 <- 1 - (multreg.RSS / SStot)
  # TODO(r24mille): Why df[1] - 1 rather than just df[1]
  multreg.adjr2 <- AdjustedRSquare(r2 = multreg.r2, n = length(tsdata.trimmed[,"kwh"]), p = (summary(multreg)$df[1]-1))
  
  multreg.explpwr <- matrix(NA, nrow = 1, ncol = 4)
  colnames(multreg.explpwr) <- c("samplesize", "AdjR2", "BIC", "DW")
  if(logtransform == TRUE) {
    multreg.explpwr.path <- paste0("./results/logtransform/", methodname, "/tslm_explanatory_power.xlsx")
  } else {
    print("Do my Adj. R^2 and value fom lm() match?")
    print(paste0("Mine: ", multreg.adjr2))
    print(paste0("lm()'s: ", summary(multreg)$adj.r.squared))
    print("=================")
    
    multreg.explpwr.path <- paste0("./results/", methodname, "/tslm_explanatory_power.xlsx")
  }
  multreg.explpwr[1,] <- c(nrow(datafr), multreg.adjr2, BIC(multreg), dwtest(multreg)$statistic)
  write.xlsx(multreg.explpwr, multreg.explpwr.path)
  
  
  
  #####
  # AR(1) structured errors
  #####
  if (logtransform == TRUE) {
    ar1err <- Arima(x = log(tsdata.trimmed[,"kwh"]), 
                    xreg = ts_xvars.trimmed,
                    order = c(1, 0, 0), 
                    include.drift = FALSE)
  } else {
    ar1err <- Arima(x = tsdata.trimmed[,"kwh"], 
                    xreg = ts_xvars.trimmed,
                    order = c(1, 0, 0), 
                    include.drift = FALSE)
  }
  
  
  # Whether fitting was done with logtransformed response or untransformed response,
  # I want to work with backtransformed fitted values and response.
  if (logtransform == TRUE) {
    ar1err.fitted.ts <- exp(fitted(ar1err))
    ar1err.resids.ts <- tsdata.trimmed[,"kwh"] - ar1err.fitted.ts
  } else {
    ar1err.fitted.ts <- fitted(ar1err)
    ar1err.resids.ts <- ar1err$residuals
  }
  
  # ACF of multiple regression residuals remaining after fitting AR(1) error structure
  if(logtransform == TRUE) {
    ar1err.acf.title <- paste0("Backtransformed Residuals ACF Plot: ", methodname, " Multiple Regression with AR(1) Errors")
    ar1err.acf.path <- paste0("./results/logtransform/", methodname, "/ar1err_residuals_acf.png")
  } else {
    ar1err.acf.title <- paste0("Residuals ACF Plot: ", methodname, " Multiple Regression with AR(1) Errors")
    ar1err.acf.path <- paste0("./results/", methodname, "/ar1err_residuals_acf.png")
  }
  acf(ar1err.resids.ts, main = ar1err.acf.title)
  dev.print(file = ar1err.acf.path, device = png, height = 600, width = 800)
  
  # Multiple regression residuals remaining after fitting AR(1) error structure as a function of time
  if(logtransform == TRUE) {
    ar1err.resid.time.title <- paste0("Backtransformed Residuals by Time: ", methodname, " Multiple Regression with AR(1) Errors")
    ar1err.resid.time.path <- paste0("./results/logtransform/", methodname, "/ar1err_residuals_by_index.png")
  } else {
    ar1err.resid.time.title <- paste0("Residuals by Time: ", methodname, " Multiple Regression with AR(1) Errors")
    ar1err.resid.time.path <- paste0("./results/", methodname, "/ar1err_residuals_by_index.png")
  }
  if(trimsize > 0) {
    ar1err.resid.time.datafr <- data.frame("posixtime" = tail(datafr, -trimsize)$timestamp_dst, 
                                           "resids" = as.vector(ar1err.resids.ts))
  } else {
    ar1err.resid.time.datafr <- data.frame("posixtime" = datafr$timestamp_dst, 
                                           "resids" = as.vector(ar1err.resids.ts))
  }
  ar1err.resid.idx.plot <- ResidualsOverTimePlot(thedata = ar1err.resid.time.datafr, thetitle = ar1err.resid.time.title)
  print(ar1err.resid.idx.plot)
  dev.print(file = ar1err.resid.time.path, device = png, height = 600, width = 800)
  
  # Multiple regression residuals remaining after fitting AR(1) error structure as a function of the fitted response
  if(logtransform == TRUE) {
    ar1err.resid.resp.title <- paste0("Backtransformed Residuals: ", methodname, " Multiple Regression with AR(1) Errors")
    ar1err.resid.resp.path <- paste0("./results/logtransform/", methodname, "/ar1err_residuals_by_fitted_response.png")  
  } else {
    ar1err.resid.resp.title <- paste0("Residuals: ", methodname, " Multiple Regression with AR(1) Errors")
    ar1err.resid.resp.path <- paste0("./results/", methodname, "/ar1err_residuals_by_fitted_response.png")
  }
  ar1err.resid.resp.plot <- ResidualsByResponseVal(thedata = data.frame("fittedvals" = ar1err.fitted.ts, 
                                                                        "resids" = ar1err.resids.ts),
                                                   thetitle = ar1err.resid.resp.title)
  print(ar1err.resid.resp.plot)
  dev.print(file = ar1err.resid.resp.path, device = png, height = 600, width = 800)
  
  # Multiple regression residuals remaining after fitting AR(1) error structure plotted as a function of temperature
  if (logtransform == TRUE) {
    ar1err.resid.temp.title <- paste0("Backtransformed Residuals as a Function of Observed Temperature:", 
                                      "\n", methodname, " Multiple Regression with AR(1) Errors")
    ar1err.resid.temp.path <- paste0("./results/logtransform/", methodname, "/ar1err_residuals_by_temperature.png")
  } else {
    ar1err.resid.temp.title <- paste0("Residuals as a Function of Observed Temperature:", 
                                      "\n", methodname, " Multiple Regression with AR(1) Errors")
    ar1err.resid.temp.path <- paste0("./results/", methodname, "/ar1err_residuals_by_temperature.png")
  }
  if(trimsize > 0) {
    ar1err.resid.temp.datafr <- data.frame("temperatures" = tail(datafr, -trimsize)$temperature, 
                                           "resids" = as.vector(ar1err.resids.ts))
  } else {
    ar1err.resid.temp.datafr <- data.frame("temperatures" = datafr$temperature, 
                                           "resids" = as.vector(ar1err.resids.ts))
  }
  ar1err.resid.temp.plot <- ResidualsByTemperature(thedata = ar1err.resid.temp.datafr, thetitle = ar1err.resid.temp.title)
  print(ar1err.resid.temp.plot)
  dev.print(file = ar1err.resid.temp.path, device = png, height = 600, width = 800)
  
  # Density of multiple regression residuals remaining after fitting AR(1) error structure
  if(logtransform == TRUE) {
    ar1err.density.title <- paste0("Backtransformed Residuals Density: ", methodname, " Multiple Regression with AR(1) Errors")
    ar1err.density.path <- paste0("./results/logtransform/", methodname, "/ar1err_residual_density.png")
  } else {
    ar1err.density.title <- paste0("Residuals Density: ", methodname, " Multiple Regression with AR(1) Errors")
    ar1err.density.path <- paste0("./results/", methodname, "/ar1err_residual_density.png")
  }
  plot(density(as.vector(ar1err.resids.ts), bw = 0.015),
       main = ar1err.density.title,
       ylim = c(0, 12), 
       xlim = c(-2, 2),
       lwd = 3)
  dev.print(file = ar1err.density.path, device = png, height = 600, width = 800)
  
  # Save several measures of multiple regression with AR(1) errors explanatory power to an Excel spreadsheet
  if(logtransform == TRUE) {
    ar1err.explpwr.path <- paste0("./results/logtransform/", methodname, "/ar1err_explanatory_power.xlsx")
  } else {
    ar1err.explpwr.path <- paste0("./results/", methodname, "/ar1err_explanatory_power.xlsx")
  }
  ar1err.explpwr <- matrix(NA, nrow = 1, ncol = 3)
  colnames(ar1err.explpwr) <- c("samplesize", "AdjR2", "BIC")
  ar1err.RSS <- sum(as.vector(ar1err.resids.ts)^2)
  ar1err.r2 <- 1 - (ar1err.RSS / SStot)
  ar1err.adjr2 <- AdjustedRSquare(r2 = ar1err.r2, n = length(tsdata.trimmed[,"kwh"]), p = ncol(ts_xvars.trimmed))
  ar1err.explpwr[1,] <- c(nrow(datafr), ar1err.adjr2, ar1err$bic)
  write.xlsx(ar1err.explpwr, ar1err.explpwr.path)
  
  
  
  #####
  # Time series cross-validation based on Rob Hyndman's example:
  # http://robjhyndman.com/hyndsight/tscvexample/
  #####
  hrs.in.month <- 730.484
  hrs.in.week <- 168
  validation.results <- TimeSeriesOutOfSampleValidation(tsfrmla = frmla.workaround, 
                                                        tsdata = tsdata.trimmed, 
                                                        xregmat = ts_xvars.trimmed, 
                                                        methodname = methodname, 
                                                        k = 12, 
                                                        hrs.in.window = hrs.in.week,
                                                        logtransform = logtransform)
}

ResidualsOverTimePlot <- function(thedata, thetitle, 
                                  ylabel = "Residuals (kWh)"){
  require(ggplot2)
  
  p <- (ggplot() 
        + geom_point(mapping = aes(posixtime, resids),
                     data = thedata,
                     alpha = 0.5) 
        + labs(x = "Date (hourly)", 
               y = ylabel,
               title = thetitle)
        + ylim(-2, 2)
        + scale_x_datetime(labels = date_format(format = "%b %Y"), 
                           breaks = "1 month")
        + theme(axis.text.x = element_text(angle = 90, 
                                           vjust = 0.5, 
                                           size = 14),
                axis.text.y = element_text(size = 14),
                axis.title.x = element_text(size = 16), 
                axis.title.y = element_text(size = 16), 
                plot.title = element_text(size = 20)))
  return (p)
}

ResidualsByResponseVal <- function(thedata, thetitle,
                                   ylabel = "Residuals (kWh)") {
  require(ggplot2)
  
  p <- (ggplot() 
        + geom_point(mapping = aes(fittedvals, resids), 
                     data = thedata, 
                     alpha = 0.5) 
        + labs(x = "Estimated Response (kWh)", 
               y = ylabel,
               title = thetitle)
        + ylim(-2, 2)
        + theme(axis.text.x = element_text(size = 14),
                axis.text.y = element_text(size = 14),
                axis.title.x = element_text(size = 16), 
                axis.title.y = element_text(size = 16), 
                plot.title = element_text(size = 20)))
  return(p)
}

ResidualsByTemperature <- function(thedata, thetitle,
                                   xlabel = "Observed Temperature (Celsius)", 
                                   ylabel = "Residuals (kWh)") {
  require(ggplot2)
  
  p <- (ggplot() 
        + geom_point(mapping = aes(temperatures, resids), 
                     data = thedata, 
                     alpha = 0.5) 
        + labs(x = xlabel, 
               y = ylabel,
               title = thetitle)
        + ylim(-2, 2)
        + theme(axis.text.x = element_text(size = 14),
                axis.text.y = element_text(size = 14),
                axis.title.x = element_text(size = 16), 
                axis.title.y = element_text(size = 16), 
                plot.title = element_text(size = 20)))
  return(p)
}

SwitchingRegressionBasisMatrix <- function(explvar, breakpoint) {
  basismatrix <- matrix(NA, nrow = length(explvar), ncol = 2)
  
  basismatrix[,1] <- ifelse(test = explvar < breakpoint,
                            yes = breakpoint - explvar,
                            no = 0)
  basismatrix[,2] <- ifelse(test = explvar > breakpoint,
                            yes = explvar - breakpoint,
                            no = 0)
  
  colnames(basismatrix) <- c("below_brk", "above_brk")
  return(basismatrix)
}

DegreeHourBasisMatrix <- function(temperatures, balance, window) {
  require(zoo)
  
  sr.basismatrix <- SwitchingRegressionBasisMatrix(explvar = temperatures, breakpoint = balance)
  hdh <- rollsum(x = sr.basismatrix[,1], k = window)
  hdh <- c(rep(NA, window-1), hdh)
  cdh <- rollsum(x = sr.basismatrix[,2], k = window)
  cdh <- c(rep(NA, window-1), cdh)
  basismatrix <- cbind(hdh, cdh)
  
  colnames(basismatrix) <- c("HDH", "CDH")
  return(basismatrix)
}