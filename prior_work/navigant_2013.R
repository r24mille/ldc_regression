library(ggplot2) # Better plotting!
library(scales) # For scaling axes (eg. datetime)
library(car) # Companion to Applied Regression (better qqPlot)
library(reshape2) # For reshaping (ie. melting) data

# Source functions in other files
source("dataframe_processing.R")
source("individual_readings_quality_checks.R")

# Load SmartMeterReading data from CSV
home <- Sys.getenv("HOME")
fpath <- file.path(home, 
                   "../Dropbox/ISS4E/R", 
                   "r24mille_lasso_25sample.csv")
readings.individual <- InitReadingsDataFrame(fpath = fpath, 
                                            is.aggregate = FALSE)

# Reduce weather descriptions to a factor with 7 levels
readings.individual$weather_reduced <- ReduceWeather(readings.individual$weather_desc)

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

# Perform data quality checks
#expl_var_dist_checks(df = readings.individual)

