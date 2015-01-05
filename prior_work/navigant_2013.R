library(ggplot2) # Better plotting!
library(scales) # For scaling axes (eg. datetime)
library(car) # Companion to Applied Regression (better qqPlot)
library(reshape2) # For reshaping (ie. melting) data

# Source functions in other files
source("dataframe_processing.R")

# Load SmartMeterReading data from CSV
home <- Sys.getenv("HOME")
fpath <- file.path(home, 
                   "../Dropbox/ISS4E/R", 
                   "r24mille_lasso_25sample.csv")
readings.individual <- InitReadingsDataFrame(fpath = fpath, 
                                            is.aggregate = FALSE)


# A small minority ~0.013% of rows have NA values for dewpnt_temp. It is a 
# fairly simple variable to derive, so derive its value so that those rows 
# may be used.


# Reduce weather descriptions to a factor with 7 levels
readings.individual$weather_desc <- ReduceWeather(readings.individual$weather_desc)

# For clarity, reorder explanatory variables which are factors
readings.individual <- OrderFactors(readings.individual)

# Derive a few temperature humidity index (THI) variables according to Navigant 
# analysis white paper.
readings.individual$nvgnt_thi <- (17.5 
                                 + (0.55 * readings.individual$temperature) 
                                 + (0.2 * readings.individual$dewpnt_temp))
readings.individual$nvgnt_cool_thi <- sapply(readings.individual$nvgnt_thi, 
                                            function(x) max(x - 30, 0))
readings.individual$nvgnt_heat_thi <- sapply(readings.individual$nvgnt_thi, 
                                            function(x) max(25 - x, 0))

summary(readings.individual)
plot(density(log(readings.individual$kwh)))
plot(readings.individual$month)