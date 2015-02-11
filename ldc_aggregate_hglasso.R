library(ggplot2) # Better plotting!
library(scales) # For scaling axes (eg. datetime)
library(glinternet) # Hierarchical group-LASSO variable selection
library(car) # Companion to Applied Regression (better qqPlot)
library(reshape2) # For reshaping (ie. melting) data

# Source functions in other files
source("dataframe_processing.R")

results.idx <- 1
results.max <- 60
results.mat <- matrix(nrow = results.max,
                      ncol = 7)
colnames(results.mat) <- c("temp_lag", "poly_degree", "mape", "rmse", 
                           "mape_bktr", "rmse_bktr", "bad_resids")
for(temp.hrs in 0:5) {
  for(poly.deg in 1:10) {
    # Load SmartMeterReading data from CSV
    home <- Sys.getenv("HOME")
    fpath <- file.path(home, 
                       "../Dropbox/ISS4E/R", 
                       "aggregate_readings_01Mar2011_through_17Oct2012.csv")
    
    readings.aggregate <- InitReadingsDataFrame(fpath = fpath, 
                                                is.aggregate = TRUE)
    rm(fpath)
    
    # Load weather descriptions from CSV (for largest city in LDC)
    fpath2 <- file.path(home, 
                        "../Dropbox/ISS4E/R", 
                        "weather_desc_01Mar2011_through_17Oct2012.csv")
    weather <- read.csv(fpath2, 
                        na.strings = c("NULL", "NA", "NaN"), 
                        stringsAsFactors = FALSE)
    rm(fpath2)
    
    # Reduce weather descriptions to a simplified set of factors
    readings.aggregate$weather_reduced <- ReduceWeather(weather$weather_desc)
    rm(weather)
    
    # For clarity, reorder explanatory variables which are factors
    readings.aggregate <- OrderFactors(readings.aggregate)    
    
    # Past temperatures
    temps <- CreatePastTemperatureMatrix(nlags = temp.hrs, 
                                         df.readings = readings.aggregate)
    readings.aggregate$temperature <- temps[,(temp.hrs+1)]
    rm(temps)
    
    # Centered(?) polynomials
    poly.vars <- poly(readings.aggregate$temperature, degree = poly.deg)
    colnames(poly.vars) <- paste0("temp_poly", c(1:poly.deg))
    readings.aggregate <- cbind(readings.aggregate, poly.vars)
    rm(poly.vars)
    
    
    # Set the weather description to be the same as the past hour of temperature 
    # used (+1 due to weatherdesc_lag0)
    pastweather <- CreatePastWeatherDescDataFrame(nlags = temp.hrs, 
                                                  weather_reduced = readings.aggregate$weather_reduced)
    readings.aggregate$weather_reduced <- pastweather[,(temp.hrs+1)]
    rm(pastweather)
    
    # Trim up data frame
    readings.trimmed <- TrimExplanatoryVariables(readings.aggregate)
    
    # Don't include the first several rows due to missing past temperature info
    readings.trimmed <- tail(x= readings.trimmed, 
                             n = nrow(readings.trimmed) - temp.hrs)
    
    # Prep data into matrices for use in hierarchical group-lasso
    Y.mat <- as.matrix(log(readings.trimmed[,1]))
    X.mat <- NumericFactorCodedMatrix(readings.trimmed[,-1])
    numlvl.vec <- NumberFactorLevels(readings.trimmed[,-1])
    
    # From a similar package, glmnet, "The default depends on the sample size nobs
    # relative to the number of variables nvars. If nobs > nvars, the default is
    # 0.0001, close to zero. If nobs < nvars, the default is 0.01. A very small 
    # value of lambda.min.ratio will lead to a saturated fit in the 
    # nobs < nvars case.
    print(paste("Start of HG-LASSO:", Sys.time()))
    cv.hglasso <- glinternet.cv(X = X.mat, 
                                Y = Y.mat, 
                                numLevels = numlvl.vec)
    print(paste("End of HG-LASSO:", Sys.time()))    
    
    # Backtransform optimal, fitted values and compute more meaningful RMSE and MAPE
    Y.bktr <- exp(Y.mat)
    hglasso.fitted.bktr <- exp(cv.hglasso$fitted)
    hglasso.rmse <- sqrt(mean((cv.hglasso$fitted - Y.mat)^2, 
                              na.rm = TRUE))
    hglasso.rmse.bktr <- sqrt(mean((hglasso.fitted.bktr - Y.bktr)^2, 
                                   na.rm = TRUE))
    hglasso.mape <- (1/nrow(Y.mat)) * sum(abs((Y.mat - cv.hglasso$fitted)/Y.mat))
    hglasso.mape.bktr <- (1/nrow(Y.bktr)) * sum(abs((Y.bktr - hglasso.fitted.bktr)/Y.bktr))
    hglasso.residuals.bktr <- Y.bktr - hglasso.fitted.bktr
    hglasso.residuals <- Y.mat - cv.hglasso$fitted
    
    readings.aggregate.tail <- tail(x = readings.aggregate, 
                                    n = nrow(readings.trimmed))
    
    # Investigate the worst residuals
    bad.resids <- readings.aggregate.tail[abs(hglasso.residuals.bktr) > 0.6,]
    bad.resids$residual <- hglasso.residuals.bktr[abs(hglasso.residuals.bktr) > 0.6]
    bad.resids$timestamp_dst <- as.character(bad.resids$timestamp_dst)
    
    
    
    # Store results and increment index
    results.mat[results.idx, "temp_lag"] <- temp.hrs
    results.mat[results.idx, "poly_degree"] <- poly.deg
    results.mat[results.idx, "mape"] <- hglasso.mape
    results.mat[results.idx, "rmse"] <- hglasso.rmse
    results.mat[results.idx, "mape_bktr"] <- hglasso.mape.bktr
    results.mat[results.idx, "rmse_bktr"] <- hglasso.rmse.bktr
    results.mat[results.idx, "bad_resids"] <- nrow(bad.resids)
    results.idx <- results.idx + 1
    
    # Clean up before next loop iteration
    rm(home, X.mat, Y.mat, readings.aggregate, readings.trimmed, 
       readings.aggregate.tail, numlvl.vec, cv.hglasso, Y.bktr, 
       hglasso.rmse, hglasso.rmse.bktr, hglasso.mape, hglasso.mape.bktr, 
       hglasso.residuals, hglasso.residuals.bktr, hglasso.fitted.bktrbad.resids)
  }
}

# Plot out backtransformed MAPE from results
results.df <- as.data.frame(results.mat)
results.df$temp_lag <- as.factor(results.df$temp_lag)
mape.bktr.ggplot <- (ggplot(results.df, aes(y = mape_bktr, 
                                            x = poly_degree,
                                            group = temp_lag,
                                            colour = temp_lag))
                     + geom_line()
                     + scale_colour_brewer(palette="Dark2")
                     + scale_x_discrete()
                     + labs(x = "Polynomial Degrees", 
                            y = "Backtransformed Mean Absolute Percentage Error (MAPE)",
                            colour = "Temperature \nHours Ago",
                            title= "Change in Backtransformed MAPE as Temperature Polynomial Increases")
                     + theme(legend.title.align = 0.5,
                             legend.justification = "center")
)
mape.bktr.ggplot