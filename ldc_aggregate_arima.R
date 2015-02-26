library(forecast)
library(car)
library(lmtest)
library(sandwich)
library(ggplot2)

# Source functions in other files
source("dataframe_processing.R")

# Load SmartMeterReading data from CSV
readings_fpath <- file.path(Sys.getenv("HOME"), 
                   "./Dropbox/ISS4E/R", 
                   "aggregate_readings_01Mar2011_through_17Oct2012.csv")

readings.aggregate <- InitReadingsDataFrame(fpath = readings_fpath, 
                                            is.aggregate = TRUE)

# Load weather descriptions from CSV (for largest city in LDC)
weather_fpath <- file.path(Sys.getenv("HOME"), 
                    "./Dropbox/ISS4E/R", 
                    "weather_desc_01Mar2011_through_17Oct2012.csv")
weather <- read.csv(weather_fpath, 
                    na.strings = c("NULL", "NA", "NaN"), 
                    stringsAsFactors = FALSE)
rm(readings_fpath, weather_fpath)

# Reduce weather descriptions to a simplified set of factors
readings.aggregate$weather_reduced <- factor(ReduceWeather(weather$weather_desc))

# For clarity, reorder explanatory variables which are factors
readings.aggregate <- OrderFactors(readings.aggregate)    

# Past temperatures
temp.hrs <- 2
temps <- CreatePastTemperatureMatrix(nlags = temp.hrs, 
                                     df.readings = readings.aggregate)
temps <- tail(x = temps, n = -2) # Trim out first two rows of NAs

# Centered polynomials
poly.deg <- 3
poly.temperatures <- matrix(data = NA, nrow = nrow(temps), ncol = 0)
for(i in 0:temp.hrs) {
  poly.iter <- poly(temps[,i+1], degree = poly.deg)
  colnames(poly.iter) <- paste0("temp_lag", i, "_poly", c(1:poly.deg))
  poly.temperatures <- cbind(poly.temperatures, poly.iter)
  rm(poly.iter)
}

# Check the Variance Inflation Factors of polynomials and lags
# Obviously lags of poly1 vary together, poly2 vary together, and poly3 vary 
# together.
log_kwh <- log(tail(readings.aggregate$kwh, -2))
vif_check <- cbind(log_kwh, poly.temperatures)
vif_check.lm <- lm(log_kwh ~ .,
                     data = as.data.frame(vif_check))
plot(y = vif_check.lm$fitted, x = readings.trimmed$temperature)
vif(vif_check.lm)

# Try replacing month with school year
readings.aggregate$schoolyr <- rep(TRUE, nrow(readings.aggregate))
readings.aggregate$schoolyr[readings.aggregate$month %in% c("m7", "m8")] <- FALSE

# Try replacing month with electricity rate seaons, even before TOU the flat rates were divided into two seasons
readings.aggregate$rateseason <- rep(NA, nrow(readings.aggregate))
readings.aggregate$rateseason[readings.aggregate$month %in% c("m5", "m6", "m7", "m8", "m9", "m10")] <- "summer"
readings.aggregate$rateseason[readings.aggregate$month %in% c("m11", "m12", "m1", "m2", "m3", "m4")] <- "winter"
readings.aggregate$rateseason <- factor(readings.aggregate$rateseason, c("winter", "summer"))

# Create matrices suitable for ARIMA modeling
readings.trimmed <- tail(readings.aggregate, n = -2)
hrstr_dummies <- model.matrix(~hrstr, data = readings.trimmed)[,-1]
dayname_dummies <- model.matrix(~dayname, data = readings.trimmed)[,-1]
price_dummies <- model.matrix(~price, data = readings.trimmed)[,-1]
rateseason_summer <- model.matrix(~rateseason, data = readings.trimmed)[,-1]
weather_dummies <- model.matrix(~weather_reduced, data = readings.trimmed)[,-1]

arma_xvars <- cbind(readings.trimmed$holiday,
                    readings.trimmed$schoolyr,
                   hrstr_dummies, 
                   dayname_dummies, 
                   price_dummies,
                   rateseason_summer,
                   weather_dummies,
                   poly.temperatures[,1:9])
colnames(arma_xvars)[1:2] <- c("holiday", "schoolyr")
       
lm.arimaerr <- Arima(x = log(readings.trimmed$kwh), 
                     xreg = arma_xvars,
                     order = c(3, 0, 0), 
                     include.drift = TRUE)
summary(lm.arimaerr)
plot(lm.arimaerr$residuals, type = "p", pch = 20, col = rgb(red=0, green=0, blue=0, alpha=25, maxColorValue = 100),
     main = "Mutliple Regression Fit with AR(3) Error Structure",
     ylab = "Residuals",
     xlab = "Index")
tsdisplay(arima.errors(lm.arimaerr), main="ARIMA errors")
#autofit <- auto.arima(x = log(readings.aggregate$kwh),
#                      xreg = arma_xvars,
#                      seasonal = FALSE,
#                      stationary = TRUE,
#                      allowdrift = TRUE)
#tsdisplay(arima.errors(autofit), main="ARIMA errors")

# Play with dynlm a bit
#library(dynlm)
#readings.zoo <- zoo(x = data.matrix(readings.aggregate),
#                    order.by = readings.aggregate$timestamp_dst)

#dlm1 <- dynlm(log(kwh) ~ temperature,
#              data = readings.zoo)
#dlm2 <- dynlm(log(kwh) ~ temperature + lag(temperature, 2),
#              data = readings.zoo)

# Fit a simple multiple regression
readings.minimal <- cbind(subset(readings.trimmed, select = c("kwh", "holiday", "schoolyr", "hrstr", "dayname", "price", "rateseason", "weather_reduced")), poly.temperatures[,1:3])
lm1 <- lm(log(kwh) ~ .,
          data = readings.minimal)
summary(lm1)

# Durban Watson test for serially correlated errors
results.dwtest <- dwtest(formula = lm1) 
results.dwtest

# Newey-West standard errors
results.neweywest <- NeweyWest(x = lm1,
                               lag = 4, 
                               diagnostics = TRUE, 
                               sandwich = TRUE,
                               verbose = TRUE)
results.neweywest
lm1.adj <- coeftest(lm1, vcov = NeweyWest(lm1, lag = 4))
lm1.adj