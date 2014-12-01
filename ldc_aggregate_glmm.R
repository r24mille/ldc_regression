library(segmented) # Find linear regression breakpoint
library(glinternet) # Hierarchical group-LASSO variable selection

# Source functions in other files
source("dataframe_processing.R")

# Load SmartMeterReading data from CSV
home <- Sys.getenv("HOME")
fpath <- file.path(home, 
                   "../Dropbox/ISS4E/R/", 
                   "aggregate_readings_01Mar2011_through_17Oct2012.csv")
readings.aggregate <- read.csv(fpath)
readings.aggregate <- InitAggregateReadings(readings.aggregate)

##
# Use 'segmented' package to find the optimal temperature breakpoint.
# 
# TODO(r24mille: I attempted to use psi c(0,16) to find two temperature 
#                breakpoints (ie. Two-Threshold Regression from Moral-Carcedo 
#                et al.) but it does not converge. Looking at the plot, this 
#                seems sensible that only a Switching Regression applies.
model.readings.lm.presegment <- lm(log(kwh) ~ temperature, 
                                     data = readings.aggregate)
seg <- segmented(obj = model.readings.lm.presegment, 
                 seg.Z = ~ temperature,
                 psi = 16)

# plot(seg, 
#      res=TRUE, 
#      res.col = rgb(0, 0, 0, 10, maxColorValue=100),
#      pch = 16,
#      main = "Piecewise Linear Fit with Temperature Breakpoint",
#      xlab = "Summer Temperatures (Celsius)",
#      ylab = "Effect of Temperature",
#      col = "red",
#      lwd = 2)

temperature.break <- seg$psi[1,2]
readings.aggregate$temp_over_break <- ifelse(readings.aggregate$temperature > temperature.break, 
                                             readings.aggregate$temperature - temperature.break, 
                                             0)
readings.aggregate$temp_under_break <- ifelse(readings.aggregate$temperature < temperature.break, 
                                              readings.aggregate$temperature - temperature.break, 
                                              0)

# Clean up some unneeded variables from R environment
rm(model.readings.lm.presegment, seg)

# Create columns that contain the previous hours' temperature > breakpoint 
# and previous hours temperature < breakpoint.
temp.hrs <- 8
temps <- CreatePastTemperatureMatrix(nlags = temp.hrs, 
                                     df.readings = readings.aggregate)
readings.aggregate <- cbind(readings.aggregate, temps)
rm(temps)
readings.trimmed <- TrimExplanatoryVariables(readings.aggregate)

# Prep data into matrices for use in hierarchical group-lasso
Y.mat <- as.matrix(log(readings.trimmed[,1]))
X.mat <- NumericFactorCodedMatrix(readings.trimmed[,-1])
numlvl.vec <- NumberFactorLevels(readings.trimmed[,-1])
lasso <- glinternet(X = X.mat, 
                    Y = Y.mat, 
                    numLevels = numlvl.vec, 
                    family = "gaussian")