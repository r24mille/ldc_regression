library(segmented) # Find linear regression breakpoint
library(ggplot) # Better plots!
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
cv.hglasso <- glinternet.cv(X = X.mat, 
                            Y = Y.mat, 
                            numLevels = numlvl.vec, 
                            family = "gaussian")

# Backtransform optimal, fitted values and compute more meaningful RMSE and MAPE
Y.bktr <- exp(Y.mat)
hglasso.fitted.bktr <- exp(cv.hglasso$fitted)
hglasso.rmse.bktr <- sqrt(mean((hglasso.fitted.bktr - Y.bktr)^2, 
                                     na.rm = TRUE))
hglasso.mape.bktr <- (1/nrow(Y.bktr)) * sum(abs((Y.bktr - hglasso.fitted.bktr)/Y.bktr))
ghlasso.residuals.bktr <- Y.bktr - hglasso.fitted.bktr

# QQ Plot of fitted vs. observed values
library(car) # Companion to Applied Regression (better qqPlot)
qqPlot(x = ghlasso.residuals.bktr,
       distribution = "norm",
       envelope = .95,
       mean = 0,
       main = "Hierarchical Group-LASSO Residuals vs. Normal Distribution",
       ylab = "Backtransformed Residuals (kWh)",
       xlab = "Norm Quantiles")

# Plot the histogram of residuals
plot(density(ghlasso.residuals.bktr),
     main = "Histogram of Residuals",
     xlab = "Backtransformed Residuals (kWh)")
abline(v = mean(ghlasso.residuals.bktr),
       col = "darkgray",
       lty = 2)
axis(3, # top label
     padj = 1, # align to bottom of label bounding box
     at = round(mean(ghlasso.residuals.bktr), 3))

# Plot residuals to look for patterns
ghlasso.df.resids.bktr <- data.frame(ghlasso.residuals.bktr, 
                                     readings.aggregate$timestamp_dst)
names(ghlasso.df.resids.bktr) <- c("res", "dst")
ggresiduals <- (ggplot(ghlasso.df.resids.bktr, aes(y = res, x = dst))
                + geom_point(alpha = 1/3)
                + labs(x = "Date (hourly)", 
                       y = "Backtransformed Residuals (kWh)",
                       title= paste("Residuals of Optimal Hierarchical Group-LASSO Model",
                                    "\n(ie. Observed - Fitted)"))
                + scale_x_datetime(labels = date_format("%b %Y"), 
                                   breaks = "1 month")
                + theme(axis.text.x = element_text(angle=90, 
                                                   vjust=0.5, 
                                                   size=12),
                        axis.text.y = element_text(size=12))
)
ggresiduals