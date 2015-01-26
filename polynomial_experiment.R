library(segmented) # Find linear regression breakpoint
library(ggplot2) # Better plotting!
library(scales) # For scaling axes (eg. datetime)
library(glinternet) # Hierarchical group-LASSO variable selection
library(car) # Companion to Applied Regression (better qqPlot)
library(reshape2) # For reshaping (ie. melting) data

# Source functions in other files
source("dataframe_processing.R")

# Load SmartMeterReading data from CSV
home <- Sys.getenv("HOME")
fpath <- file.path(home, 
                   "../Dropbox/ISS4E/R", 
                   "aggregate_readings_01Mar2011_through_17Oct2012.csv")

readings.aggregate <- InitReadingsDataFrame(fpath = fpath, 
                                            is.aggregate = TRUE)

# Load weather descriptions from CSV (for largest city in LDC)
fpath2 <- file.path(home, 
                    "../Dropbox/ISS4E/R", 
                    "weather_desc_01Mar2011_through_17Oct2012.csv")
weather <- read.csv(fpath2, 
                    na.strings = c("NULL", "NA", "NaN"), 
                    stringsAsFactors = FALSE)

# Reduce weather descriptions to a simplified set of factors
readings.aggregate$weather_reduced <- ReduceWeather(weather$weather_desc)

# For clarity, reorder explanatory variables which are factors
readings.aggregate <- OrderFactors(readings.aggregate)



# Quadratic experiment
#poly.vars <- poly(readings.aggregate$temperature, degree = 2)
#colnames(poly.vars) <- c("temp_poly1", "temp_poly2")
#readings.aggregate <- cbind(readings.aggregate, poly.vars)
poly.lm1 <- lm(log(kwh) ~ month + hrstr + dayname + price + holiday + weather_reduced + poly(temperature, 1),
              data = readings.aggregate)

poly.lm2 <- lm(log(kwh) ~ month + hrstr + dayname + price + holiday + weather_reduced + poly(temperature, 2),
               data = readings.aggregate)

poly.lm3 <- lm(log(kwh) ~ month + hrstr + dayname + price + holiday + weather_reduced+ poly(temperature, 3),
                   data = readings.aggregate)

poly.lm4 <- lm(log(kwh) ~ month + hrstr + dayname + price + holiday + weather_reduced+ poly(temperature, 4),
               data = readings.aggregate)

poly.lm5 <- lm(log(kwh) ~ month + hrstr + dayname + price + holiday + weather_reduced+ poly(temperature, 5),
               data = readings.aggregate)

poly.lm6 <- lm(log(kwh) ~ month + hrstr + dayname + price + holiday + weather_reduced+ poly(temperature, 6),
               data = readings.aggregate)

poly.lm7 <- lm(log(kwh) ~ month + hrstr + dayname + price + holiday + weather_reduced+ poly(temperature, 7),
               data = readings.aggregate)

poly.lm8 <- lm(log(kwh) ~ month + hrstr + dayname + price + holiday + weather_reduced+ poly(temperature, 8),
               data = readings.aggregate)

anova(poly.lm1, poly.lm2, poly.lm3, poly.lm4, poly.lm5, poly.lm6, poly.lm7, poly.lm8)


fitted.bktr <- exp(poly.lm7$fitted.values)
qplot(x = readings.aggregate$temperature,
      y = readings.aggregate$kwh)
qplot(x = readings.aggregate$temperature, 
      y = fitted.bktr)
qplot(x = readings.aggregate$temperature,
      y = readings.aggregate$kwh - fitted.bktr)


# See what the polinomial line looks like without other vars
qplot(x = readings.aggregate$temperature, 
      y = log(readings.aggregate$kwh), 
      geom=c("point", "smooth"), 
      method="lm", 
      formula = y ~poly(x, 3)) 