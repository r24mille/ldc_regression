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
                   "aggregate_readings_29Oct2010_through_17Oct2012.csv")

readings.aggregate <- InitReadingsDataFrame(fpath = fpath, 
                                            is.aggregate = TRUE)

# Load weather descriptions from CSV (for largest city in LDC)
fpath2 <- file.path(home, 
                   "../Dropbox/ISS4E/R", 
                   "weather_desc_29Oct2010_through_17Oct2012.csv")
weather <- read.csv(fpath2, 
                    na.strings = c("NULL", "NA", "NaN"), 
                    stringsAsFactors = FALSE)

# Reduce weather descriptions to a simplified set of factors
readings.aggregate$weather_reduced <- ReduceWeather(weather$weather_desc)

# For clarity, reorder explanatory variables which are factors
readings.aggregate <- OrderFactors(readings.aggregate)

# Use 'segmented' package to find the optimal temperature breakpoint(s)
model.readings.lm.presegment <- lm(log(kwh) ~ temperature + month + hrstr + price + holiday, 
                                     data = readings.aggregate)
seg <- segmented(obj = model.readings.lm.presegment, 
                 seg.Z = ~ temperature,
                 psi = c(17, 32))

plot(seg, 
     res = TRUE, 
     res.col = rgb(0, 0, 0, 10, maxColorValue=100),
     pch = 16,
     main = "Piecewise Linear Fit with Temperature Breakpoint",
     xlab = "Summer Temperatures (Celsius)",
     ylab = "Effect of Temperature",
     col = "red",
     lwd = 2)

temperature.break <- seg$psi[1,2]
temperature.exhaust <- seg$psi[2,2]

# Standard temperature over breakpoint metric that I've been using
readings.aggregate$temp_over_break <- ifelse(readings.aggregate$temperature > temperature.break, 
                                             readings.aggregate$temperature - temperature.break, 
                                             0)

# Humidex used as temperature over breakpoint (misnomer for now)
#readings.aggregate$temp_over_break <- readings.aggregate$humidex
#readings.aggregate$temp_over_break[is.na(readings.aggregate$temp_over_break)] <- 0
#readings.aggregate$temp_over_break <- readings.aggregate$humidex_diff

# Navigant cooling THI used as temperature over breakpoint (misnomer for now)
#readings.aggregate$temp_over_break <- readings.aggregate$nvgnt_cool_thi

# Limit temp_over_break to the exhaustion point ~32C
readings.aggregate$temp_over_break <- sapply(readings.aggregate$temp_over_break, 
                                             function(x) min(x, temperature.exhaust - temperature.break))


# Standard temperature under breakpoint metric that I've been using
readings.aggregate$temp_under_break <- ifelse(readings.aggregate$temperature < temperature.break, 
                                              readings.aggregate$temperature - temperature.break, 
                                              0)

# Wind chill used as temperature under breakpoint (misnomer for now)
#readings.aggregate$temp_under_break <- readings.aggregate$wind_chill
#readings.aggregate$temp_under_break[is.na(readings.aggregate$temp_under_break)] <- 0
#readings.aggregate$temp_under_break <- readings.aggregate$wind_chill_diff

# Navigant heating THI used as temperature under breakpoint (misnomer for now)
#readings.aggregate$temp_under_break <- readings.aggregate$nvgnt_heat_thi

# Clean up some unneeded variables from R environment
rm(model.readings.lm.presegment, seg)

# Create columns that contain the previous hours' temperature > breakpoint 
# and previous hours temperature < breakpoint.
temp.hrs <- 6
temps <- CreatePastTemperatureMatrix(nlags = temp.hrs, 
                                     df.readings = readings.aggregate)
readings.aggregate <- cbind(readings.aggregate, temps)
rm(temps)
readings.trimmed <- TrimExplanatoryVariables(readings.aggregate)
# Don't include the first several rows due to missing past temperature info
readings.trimmed <- tail(x= readings.trimmed, 
                         n = nrow(readings.trimmed) - temp.hrs) 

# Prep data into matrices for use in hierarchical group-lasso
Y.mat <- as.matrix(log(readings.trimmed[,1]))
X.mat <- NumericFactorCodedMatrix(readings.trimmed[,-1])
numlvl.vec <- NumberFactorLevels(readings.trimmed[,-1])
print(paste("Start of HG-LASSO:", Sys.time()))
cv.hglasso <- glinternet.cv(X = X.mat, 
                            Y = Y.mat, 
                            numLevels = numlvl.vec)
print(paste("End of HG-LASSO:", Sys.time()))

# Select the column names of main effects and interactions
cat.colnames <- names(readings.trimmed[c(FALSE, numlvl.vec > 1)])
cont.colnames <- names(readings.trimmed[c(FALSE, numlvl.vec == 1)])

# Get coefficients for all lambdas
coef.hglasso <- coef(cv.hglasso$glinternetFit)

# Column names for main effects at chosen lambda
#lmda.idx <- which(cv.hglasso$lambda == cv.hglasso$lambdaHat1Std)
for (lmda.idx in 1:which(cv.hglasso$lambda == cv.hglasso$lambdaHat1Std)) {
  me.colnames <- c(cat.colnames[coef.hglasso[[lmda.idx]]$mainEffects$cat],
                     cont.colnames[coef.hglasso[[lmda.idx]]$mainEffects$cont])
  # Column names for interactions at chosen lambda
  lmda.catcat <- coef.hglasso[[lmda.idx]]$interactions$catcat
  lmda.contcont <- coef.hglasso[[lmda.idx]]$interactions$contcont
  lmda.catcont <- coef.hglasso[[lmda.idx]]$interactions$catcont
  int.colnames <- c(paste(cat.colnames[lmda.catcat[,1]], cat.colnames[lmda.catcat[,2]], sep=":"),
                    paste(cont.colnames[lmda.contcont[,1]], cont.colnames[lmda.contcont[,2]], sep=":"),
                    paste(cat.colnames[lmda.catcont[,1]], cont.colnames[lmda.catcont[,2]], sep=":"))
  print(paste("Lmabda Index =", lmda.idx))
  print(c(me.colnames, int.colnames))
}


# Backtransform optimal, fitted values and compute more meaningful RMSE and MAPE
Y.bktr <- exp(Y.mat)
hglasso.fitted.bktr <- exp(cv.hglasso$fitted)
hglasso.rmse.bktr <- sqrt(mean((hglasso.fitted.bktr - Y.bktr)^2, 
                                     na.rm = TRUE))
hglasso.mape.bktr <- (1/nrow(Y.bktr)) * sum(abs((Y.bktr - hglasso.fitted.bktr)/Y.bktr))
hglasso.residuals.bktr <- Y.bktr - hglasso.fitted.bktr
hglasso.residuals <- Y.mat - cv.hglasso$fitted

# QQ Plot of fitted vs. observed values
qqPlot(x = hglasso.residuals,
       distribution = "norm",
       envelope = .95,
       mean = 0,
       main = "Hierarchical Group-LASSO Residuals vs. Normal Distribution",
       ylab = "Residuals of log(kWh)",
       xlab = "Norm Quantiles")

# Plot the histogram of residuals
plot(density(hglasso.residuals),
     main = "Histogram of Residuals",
     xlab = "Residuals of log(kWh)")
abline(v = mean(hglasso.residuals),
       col = "darkgray",
       lty = 2)
axis(3, # top label
     padj = 1, # align to bottom of label bounding box
     at = round(mean(hglasso.residuals), 3))

# Plot residuals to look for patterns
readings.aggregate.tail <- tail(x = readings.aggregate, 
                                n = nrow(readings.trimmed)) 
hglasso.df.resids <- data.frame(hglasso.residuals, 
                                readings.aggregate.tail$temperature,
                                readings.aggregate.tail$agg_count, 
                                X.mat[,"dayname"],
                                X.mat[,"holiday"], 
                                X.mat[,"price"],
                                X.mat[,"month"],
                                X.mat[,"hrstr"],
                                readings.aggregate.tail$timestamp_dst)
names(hglasso.df.resids) <- c("residuals", "temperature", "aggregate_count", 
                              "dayname", "holiday", "price", "month", "hour", 
                              "time")
ggresiduals <- (ggplot(hglasso.df.resids, aes(y = residuals, x = time))
                + geom_point(alpha = 1/3)
                + labs(x = "Date (hourly)", 
                       y = "Residuals of log(kWh)",
                       title= bquote(Residuals~of~Hierarchical~Group-LASSO~Model~at~hat(lambda)==.(cv.hglasso$lambdaHat, 2)))
                + scale_x_datetime(labels = date_format(format = "%b %Y"), 
                                   breaks = "1 month")
                + theme(axis.text.x = element_text(angle=90, 
                                                   vjust=0.5, 
                                                   size=12),
                        axis.text.y = element_text(size=12))
)
ggresiduals
hglasso.df.resids.long <- melt(hglasso.df.resids, id.vars="time")
ggresiduals.facet <- (ggplot() 
                      + geom_point(data = subset(hglasso.df.resids.long, 
                                                 variable == "residuals"), 
                                   aes(time, value)) 
                      + geom_line(data = subset(hglasso.df.resids.long, 
                                                variable == "temperature"),
                                  aes(time, value))
                      + geom_step(data = subset(hglasso.df.resids.long, 
                                                variable == "aggregate_count"), 
                                  aes(time, value)) 
                      + geom_step(data = subset(hglasso.df.resids.long, 
                                                variable == "holiday"), 
                                  aes(time, value)) 
                      + facet_grid(variable ~ ., scales = "free_y")
                      + labs(x = "Date (hourly)", 
                             y = "Facet Value (see right-hand label)",
                             title= "Residuals and Independent Variables of Interest")
                      + scale_x_datetime(labels = date_format(format = "%b %Y"), 
                                         breaks = "1 month")
                      + theme(axis.text.x = element_text(angle=90, 
                                                         vjust=0.5, 
                                                         size=12),
                              axis.text.y = element_text(size=12)))
ggresiduals.facet

# Plot pairs
hglasso.df.pairs <- data.frame(hglasso.residuals, 
                               hglasso.residuals.bktr, 
                               log(readings.aggregate.tail$kwh), 
                               readings.aggregate.tail$kwh,
                               readings.aggregate.tail$temperature,
                               readings.aggregate.tail$agg_count, 
                               readings.aggregate.tail$dayname,
                               readings.aggregate.tail$holiday, 
                               readings.aggregate.tail$price,
                               readings.aggregate.tail$month,
                               readings.aggregate.tail$hrstr,
                               readings.aggregate.tail$weather_reduced,
                               readings.aggregate.tail$timestamp_dst)
names(hglasso.df.pairs) <- c("residuals", "residuals_bktr", "log_kwh", "kwh", 
                             "temperature", "aggregate_count", "dayname", 
                             "holiday", "price", "month", "hour", "weather", 
                             "time")

ggresiduals.pairs.log_kwh <- (ggplot(hglasso.df.pairs, aes(y = residuals, 
                                                           x = log_kwh))
                              + geom_point(alpha = 1/2)
                              + labs(x = "Lognormal response, log(kWh)", 
                                     y = "Residuals of lognormal model",
                                     title= "Model Residuals as a Function of log(kWh)"))
ggresiduals.pairs.log_kwh

ggresiduals.pairs.kwh <- (ggplot(hglasso.df.pairs, aes(y = residuals_bktr, 
                                                        x = kwh))
                           + geom_point(alpha = 1/2)
                           + labs(x = "Observed Meter Reading (kWh)", 
                                  y = "Residuals",
                                  title= "Backtrasformed Residuals as a Function Observations"))
ggresiduals.pairs.kwh

ggresiduals.pairs.temp <- (ggplot(hglasso.df.pairs, aes(y = residuals, 
                                                        x = temperature))
                           + geom_point(alpha = 1/2)
                           + labs(x = "Outdoor Temperature (Celsius)", 
                                  y = "Residuals of log(kWh)",
                                  title= "Residuals as a Function of Temperature"))
ggresiduals.pairs.temp
ggresiduals.pairs.agg <- (ggplot(hglasso.df.pairs, aes(y = residuals, 
                                                       x = aggregate_count))
                          + geom_point(alpha = 1/2)
                          + labs(x = "Number of Meters Used to Determine Aggregate Average", 
                                 y = "Residuals of log(kWh)",
                                 title= "Residuals as a Function of the Number of Meters Used to Determine the Aggregate Average Reading"))
ggresiduals.pairs.agg
ggresiduals.pairs.holiday <- (ggplot(hglasso.df.pairs, aes(y = residuals, 
                                                           x = holiday))
                              + geom_boxplot()
                              + labs(x = "Holiday (as observed by Ontario Energy Board)", 
                                     y = "Residuals of log(kWh)",
                                     title= "Residuals as a Function of Holidays"))
ggresiduals.pairs.holiday
ggresiduals.pairs.hour <- (ggplot(hglasso.df.pairs, aes(y = residuals, 
                                                        x = hour))
                           + geom_boxplot()
                           + labs(x = "Hour of Day", 
                                  y = "Residuals of log(kWh)",
                                  title= "Residuals as a Function of Hour-of-Day"))
ggresiduals.pairs.hour
ggresiduals.pairs.dayname <- (ggplot(hglasso.df.pairs, aes(y = residuals, 
                                                           x = dayname))
                              + geom_boxplot()
                              + labs(x = "Day of Week", 
                                     y = "Residuals of log(kWh)",
                                     title= "Residuals as a Function of Day-of-Week"))
ggresiduals.pairs.dayname
ggresiduals.pairs.month <- (ggplot(hglasso.df.pairs, aes(y = residuals, 
                                                         x = month))
                            + geom_boxplot()
                            + labs(x = "Month", 
                                   y = "Residuals of log(kWh)",
                                   title= "Residuals as a Function of Month-of-Year"))
ggresiduals.pairs.month
ggresiduals.pairs.price <- (ggplot(hglasso.df.pairs, aes(y = residuals, 
                                                         x = price))
                            + geom_boxplot()
                            + labs(x = "Electricity Pricing Category", 
                                   y = "Residuals of log(kWh)",
                                   title= "Residuals as a Function of Pricing Sturcture"))
ggresiduals.pairs.price
ggresiduals.pairs.weather <- (ggplot(hglasso.df.pairs, aes(y = residuals, 
                                                         x = weather))
                            + geom_boxplot()
                            + labs(x = "Weather Description", 
                                   y = "Residuals of log(kWh)",
                                   title= "Residuals as a Function of Weather Description"))
ggresiduals.pairs.weather