library(dlnm)
library(ggplot2) # Better plotting!
library(scales) # For scaling axes (eg. datetime)
library(glinternet) # Hierarchical group-LASSO variable selection
library(car) # Companion to Applied Regression (better qqPlot)
library(reshape2) # For reshaping (ie. melting) data
library(lmtest) # Test for serially correlated errors
library(sandwich) # Newey West test for variable significance once serially correlated errors are accounted for

# Source functions in other files
source("dataframe_processing.R")

# Load SmartMeterReading data from CSV
home <- Sys.getenv("HOME")
fpath <- file.path(home, 
                   "./Dropbox/ISS4E/R", 
                   "aggregate_readings_01Mar2011_through_17Oct2012.csv")

readings.aggregate <- InitReadingsDataFrame(fpath = fpath, 
                                            is.aggregate = TRUE)
rm(fpath)

# Load weather descriptions from CSV (for largest city in LDC)
fpath2 <- file.path(home, 
                    "./Dropbox/ISS4E/R", 
                    "weather_desc_01Mar2011_through_17Oct2012.csv")
weather <- read.csv(fpath2, 
                    na.strings = c("NULL", "NA", "NaN"), 
                    stringsAsFactors = FALSE)
rm(fpath2, home)

# Reduce weather descriptions to a simplified set of factors
readings.aggregate$weather_reduced <- ReduceWeather(weather$weather_desc)

# For clarity, reorder explanatory variables which are factors
readings.aggregate <- OrderFactors(readings.aggregate)    

# Try replacing month with school year
readings.aggregate$schoolyr <- rep(TRUE, nrow(readings.aggregate))
readings.aggregate$schoolyr[readings.aggregate$month %in% c("m7", "m8")] <- FALSE

# Try replacing month with electricity rate seaons, even before TOU the flat rates were divided into two seasons
readings.aggregate$rateseason <- rep(NA, nrow(readings.aggregate))
readings.aggregate$rateseason[readings.aggregate$month %in% c("m5", "m6", "m7", "m8", "m9", "m10")] <- "summer"
readings.aggregate$rateseason[readings.aggregate$month %in% c("m11", "m12", "m1", "m2", "m3", "m4")] <- "winter"
readings.aggregate$rateseason <- factor(readings.aggregate$rateseason, c("winter", "summer"))



#####
# Start with a simple DLNM model in which the response variable is only predicted by temperature
#####
log_kwh <- log(readings.aggregate$kwh)
temps <- readings.aggregate$temperature

# Check a correlation matrix of past temperatures to see how many lags should be used, similar to Almon (1965)
temp.hrs <- 23
temp.lagmatrix <- CreatePastTemperatureMatrix(nlags = temp.hrs, df.readings = readings.aggregate)

# Trim out NAs from dependent and independent variables equally
log_kwh.trimmed <- log_kwh[(temp.hrs+1):length(log_kwh)]
temp.lagmatrix.trimmed <- tail(x = temp.lagmatrix, n = (nrow(temp.lagmatrix) - temp.hrs))

# Correlation matrix of the response variable and lags of temperature
cor(y = log_kwh.trimmed, x = temp.lagmatrix.trimmed)

# As per Almon (1965), look for lag that is as correlated with response as lag=0
# Viewing the above correlation matrix, this appears to be lag=5.
max.lag <- 5

cb <- crossbasis(x = temps, 
                 lag = max.lag, 
                 argvar = list(fun = "poly", degree = 3),
                 arglag = list(fun = "lin"))
#summary(cb)

# Prediction using the crossbasix matrix
model <- lm(log_kwh ~ cb, weights = readings.aggregate$agg_count)
summary(model)$r.squared

plot(x = temps[6:length(temps)], y = predict(model))
plot(model$residuals)

cb.predictions <- crosspred(cb, model)
plot(cb.predictions, ptype="3d", main="3D plot", xlab="Temperature", zlab="log(kWh)", theta = 200, ltheta = 180)
plot(cb.predictions, ptype="contour", key.title=title("log(kWh)"), plot.title=title("Contour plot", xlab="Temperature", ylab="Lag"))
par(mfrow = c(1, 1))

# A more complex linear model with other exogenous variables
model2 <- lm(log_kwh ~ cb + readings.aggregate$weather_reduced + readings.aggregate$hrstr + readings.aggregate$dayname + 
               readings.aggregate$holiday + readings.aggregate$price + readings.aggregate$schoolyr + readings.aggregate$rateseason)
summary(model2)
#plot(x = temps[6:length(temps)], y = predict(model2))
plot(model2$residuals)
#plot(model2)

# Compare R^2 to simple 3 degree polynomial
model3 <- lm(readings.aggregate$kwh ~ poly(readings.aggregate$temperature, degree = 3) + readings.aggregate$weather_reduced + 
               readings.aggregate$hrstr + readings.aggregate$dayname + readings.aggregate$holiday + readings.aggregate$price + 
               readings.aggregate$month)
summary(model3) # DNLM improves R^2 by ~.015


# Test for Variance inflation factor
vif(model2) # The VIF of lagged temperature crossbasis is high, related to hrstr, rateseason, and schoolyr (in that order)

# Durban Watson test for serially correlated errors
results.dwtest <- dwtest(formula = model2) 
results.dwtest # If residuals weren't obvious enough, the Durban Watson test shows they are serially correlated

# Newey-West standard errors
results.neweywest <- NeweyWest(x = model2,
                               lag = 5, 
                               diagnostics = TRUE, 
                               sandwich = TRUE,
                               verbose = TRUE)
results.neweywest
model2.nwadj <- coeftest(model2, vcov = NeweyWest(model2, lag = 5))
model2.nwadj # rateseason is not statistically signficant

# Remove statistically insignificant variable
model4 <- update(model2, "~ . - readings.aggregate$rateseason")
vif(model4) # VIF improved quite a bit, crossbasis is minimally worrysome
model4.nwadj <- coeftest(model4, vcov = NeweyWest(model4, lag = 5))
model4.nwadj # All remaining variables are statistically significant, even when accounting for serial correlation
summary(model4)
dwtest(formula = model4) # Still serially correlated though