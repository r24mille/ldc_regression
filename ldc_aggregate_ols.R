library(ggplot2) # Better plotting!
library(scales) # For scaling axes (eg. datetime)
library(car) # Companion to Applied Regression (better qqPlot)
library(reshape2) # For reshaping (ie. melting) data

# Source functions in other files
source("dataframe_processing.R")
source("lm_method_iteration.R")

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
readings.aggregate$weather_reduced <- factor(ReduceWeather(weather$weather_desc))
rm(weather)

# For clarity, reorder explanatory variables which are factors
readings.aggregate <- OrderFactors(readings.aggregate)    

# Past temperatures
temp.hrs <- 2
temps <- CreatePastTemperatureMatrix(nlags = temp.hrs, 
                                     df.readings = readings.aggregate)
readings.aggregate$temperature <- temps[,(temp.hrs+1)]
rm(temps)

# Centered polynomials
poly.deg <- 3
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


# Fit a saturated model
lm.saturated <- lm(formula = log(kwh) ~ (.)^2 - temp_poly1:temp_poly2 - temp_poly1:temp_poly3 - temp_poly2:temp_poly3, 
                     data = readings.trimmed)

