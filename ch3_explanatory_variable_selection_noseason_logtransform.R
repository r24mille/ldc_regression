require(dlnm)           # Distributed Lag Non-linear Modeling for lagged effects of temperature
require(ggplot2)        # Better plotting!
require(scales)         # For scaling axes (eg. datetime)
require(car)            # Companion to Applied Regression (better qqPlot)
require(lmtest)         # Durban-Watson test for serially correlated errors
require(sandwich)       # Newey-West test for variable significance in the presence of serial correlation
require(cvTools)        # Used for reporting cross-validation results
require(segmented)      # For finding temperature breakpoint
require(xlsx)           # Useful for exporting results to analysis in Excel
require(reshape2)       # For reshaping (ie. melting) data
require(gridExtra)      # Placing multiple ggplot2 plots into a single viewport or image
require(splines)        # Natural cubic splines
require(visreg)         # Visualize multiple regression portions

# Source functions in other files
source("dataframe_processing.R")



#####
# Load aggregate smart meter readings and weighted average temperature data from a CSV.
# Similarly, load the Environment Canada Weather descriptions for the area from CSV.
#####
home <- Sys.getenv("HOME")
fpath <- file.path(home, 
                   "./Dropbox/ISS4E/R", 
                   "aggregate_readings_01Mar2011_through_17Oct2012.csv")

readings.aggregate <- InitReadingsDataFrame(fpath = fpath, 
                                            is.aggregate = TRUE)

fpath2 <- file.path(home, 
                    "./Dropbox/ISS4E/R", 
                    "weather_desc_01Mar2011_through_17Oct2012.csv")
weather <- read.csv(fpath2, 
                    na.strings = c("NULL", "NA", "NaN"), 
                    stringsAsFactors = FALSE)

# Add weather information to readings.aggregate
weather_desc.original <- unique(weather$weather_desc)
readings.aggregate$weather_reduced <- factor(ReduceWeatherCoarseTerms(weather$weather_desc))
readings.aggregate$severe_weather <- ReduceSevereWeather(weather$weather_desc)
summary(readings.aggregate$severe_weather) # 146 hours of severe conditions is still arguably low

# Add working_day boolean
readings.aggregate$working_day <- ((readings.aggregate$holiday == TRUE | 
                                      readings.aggregate$dayname == "Sat" | 
                                      readings.aggregate$dayname == "Sun") == FALSE)

# For clarity, reorder explanatory variables which are factors
readings.aggregate <- OrderFactors(readings.aggregate)
rm(fpath, fpath2, home)


#####
# Visualize the number of meters aggregated hourly
#####
agg_count.tsplot <- (ggplot() 
                     + geom_area(data = readings.aggregate, 
                                 aes(timestamp_dst, agg_count)) 
                     + labs(x = "Date (hourly)", 
                            y = paste0("Number of Observations, max=", max(readings.aggregate$agg_count)),
                            title= "Number of Households Reporting Each Hour")
                     + ylim(0, (max(readings.aggregate$agg_count) + 100))
                     + scale_x_datetime(labels = date_format(format = "%b %Y"), 
                                        breaks = "1 month")
                     + theme(axis.text.x = element_text(angle=90, 
                                                        vjust=0.5, 
                                                        size=12),
                             axis.text.y = element_text(size=12)))
agg_count.tsplot
dev.print(file="../../Figures/AggregateCountTimeseries.png", device=png, height = 600, width = 800)


#####
# Plot electricity consumption by hour of day and demonstrate that it needs to be 
# explicitly handed as the periodicity term. Hopefully it improves serially correlated 
# errors.
#####
hrstr.plot <- (ggplot(readings.aggregate, aes(hrstr, kwh)) + 
                 geom_boxplot() + 
                 theme(legend.position = "none") + 
                 labs(x = "Hour of Day", 
                      y = "Average Household Electricity Demand (kWh)",
                      title= "Distribution of Average Household Electricity Demand by Hour of Day"))
hrstr.plot
dev.print(file="../../Figures/HourOfDayDemand.png", device=png, height = 600, width = 800)



#####
# Plot electricity pricing method (ie. TOU vs. flat)
#####
price.palette <- c("#ffffff", "#c2da8b", "#fcdf6f", "#df9c7d") #white, green, yellow, red
price.plot <- (ggplot(readings.aggregate, aes(price, kwh, fill = price)) + 
                 geom_boxplot() + 
                 theme(legend.position = "none") + 
                 scale_x_discrete(labels = c("Flat (pre-TOU)","Off-Peak","Mid-Peak","On-Peak")) + 
                 scale_fill_manual(values = price.palette) + 
                 labs(x = "Price Period", 
                      y = "Average Household Electricity Demand (kWh)",
                      title= "Distribution of Average Household Electricity Demand by Price Period"))
price.plot
dev.print(file="../../Figures/PriceDemand.png", device=png, height = 600, width = 800)



#####
# Introduce working_day term, since distribution of weekdays seems similar.
#####
summary(subset(readings.aggregate, subset = holiday == TRUE)) # Holidays occur on Sun/Mon/Sat
dayname.plot <- (ggplot(readings.aggregate, aes(dayname, kwh)) + 
                   geom_boxplot(aes(fill = dayname)) + 
                   theme(legend.position = "none") + 
                   labs(x = "Day of the Week", 
                        y = "Average Household Electricity Demand (kWh)",
                        title= "Distribution of Response by Day of Week"))
holiday.plot <- (ggplot(readings.aggregate, aes(holiday, kwh)) + 
                   geom_boxplot() + 
                   labs(x = "Holiday", 
                        y = "Average Household Electricity Demand (kWh)",
                        title= "Distribution of Response by Holiday Indicator"))
grid.arrange(dayname.plot, holiday.plot, ncol = 2, nrow = 1)
dev.print(file="../../Figures/DaynameAndHoliday.png", device=png, height = 600, width = 1600)
nonholiday_dayname.plot <- (ggplot(subset(readings.aggregate, subset = holiday == FALSE), 
                                   aes(dayname, kwh)) + 
                              geom_boxplot(aes(fill = dayname)) + 
                              theme(legend.position = "none") + 
                              labs(x = "Day of Week", 
                                   y = "Average Household Electricity Demand (kWh)",
                                   title= "Distribution of Response by Day of Week (non-holidays)"))
working_day.plot <- (ggplot(readings.aggregate, aes(working_day, kwh)) + 
                       geom_boxplot() + 
                       labs(x = "Working Day", 
                            y = "Average Household Electricity Demand (kWh)",
                            title= "Distribution of Response by Working Day Indicator"))
grid.arrange(nonholiday_dayname.plot, working_day.plot, ncol = 2, nrow = 1)
dev.print(file="../../Figures/NonHolidayAndWorkingDay.png", device=png, height = 600, width = 1600)



#####
# Create an interrim model which reflects some temporal explanatory variables and the price 
# explanatory variable.
#####
readings.aggregate$working_day <- factor(readings.aggregate$working_day) # visreg likes "factor" not "logical"
readings.aggregate$working_day <- as.logical(readings.aggregate$working_day)
lm.interrim <- lm(kwh ~ hrstr + working_day + price + working_day:hrstr,
                  data = readings.aggregate)
summary(lm.interrim)
dwtest(lm.interrim)
interrim.visreg <- visreg(fit = lm.interrim, 
                          xvar = "hrstr", 
                          by = "working_day",
                          type = "conditional",
                          main = paste("Distribution of Average Household Electricity Demand (kWh) for Hour of Day,",
                                       "Split by \"Working Day\" Indicator\n",
                                       "Conditioned on Several Temporal and Price Explanatory Variables"),
                          ylab = "Electricity Demand (kWh)",
                          xlab = "Hour of Day")
dev.print(file="../../Figures/HourOfDayWorkingDayInteractionConditionalPlot.png", device=png, height = 600, width = 1600)



#####
# Simple plot of Temperature vs. kWh
#####
temp.vs.demand.plot <- (ggplot(readings.aggregate, aes(x = temperature, y = kwh)) + 
                          geom_point(alpha = 0.5) + 
                          labs(x = "Outdoor Temperature (Celsius)", 
                               y = "Average Household Demand (kWh)", 
                               title= "Electricity Demand as a Function of Temperature"))
temp.vs.demand.plot
dev.print(file="../../Figures/TemperatureVsDemand.png", device=png, height = 600, width = 800)



#####
# Demonstrate that the conditional effects of temperature are non-linear
#####
lm.lintemp <- lm(log(kwh) ~ temperature + hrstr + working_day + price + hrstr:working_day,
                     data = readings.aggregate,
                     weights = readings.aggregate$agg_count)
lintemp.adj.r2 <- summary(lm.lintemp)$adj.r.squared
lintemp.visreg <- visreg(fit = lm.lintemp, xvar = "temperature", type = "conditional",
                         alpha = 1, # no confidence interval
                         line.par = list(col = "red", lwd = 3),
                         points.par = list(col = scatter.black),
                         main = "Linear Fit of Temperature to Demand Conditioned on Other Explanatory Variables", 
                         ylab = "Conditional Electricity Demand, log(kWh)",
                         xlab = "Outdoor Temperature (Celsius)")
dev.print(file="../../Figures/LinearTemperatureConditionalPlot.png", device=png, height = 600, width = 800)



#####
# Temperature breakpoint plots
#####
par(mfrow = c(1, 2))
seg1 <- segmented(obj = lm.lintemp,
                 seg.Z = ~ temperature,
                 psi = 17)
temp.break <- seg1$psi[1,2]
seg1.adj.r2 <- summary(seg1)$adj.r.squared
plot(seg1, col = "red", lwd = 3, 
     res = TRUE, res.col = rgb(0, 0, 0, 25, maxColorValue = 100), pch = 20, 
     ylim = c(0.2, 3.45), 
     main = "Linear Breakpoint Fit of Temperature to Conditional Demand",
     xlab = "Outdoor Temperature (Celsius)",
     ylab = "Conditional Electricity Demand, log(kWh)",
     rug = FALSE)

seg2 <- segmented(obj = lm.lintemp,
                  seg.Z = ~ temperature,
                  psi = c(15, 22))
temp.switch <- seg2$psi[1,2]
temp.exhaust <- seg2$psi[2,2]
seg2.adj.r2 <- summary(seg2)$adj.r.squared
plot(seg2, col = "red", lwd = 3, 
     res = TRUE, res.col = rgb(0, 0, 0, 25, maxColorValue = 100), pch = 20, 
     ylim = c(0.2, 3.45), 
     main = "Three-Line Temperature Fit to Conditional Demand",
     xlab = "Outdoor Temperature (Celsius)",
     ylab = "Conditional Electricity Demand, log(kWh)",
     rug = FALSE)

dev.print(file="../../Figures/TemperatureBreakpointPlots.png", device=png, height = 600, width = 1600)
par(mfrow = c(1, 1))



#####
# Examine simple plot of polynomials to determine appropriate degree.
#####
scatter.black = rgb(0, 0, 0, 50, maxColorValue = 100)
par(mfrow = c(2, 2))
poly2 <- lm(log(kwh) ~ poly(temperature, degree = 2) + hrstr + working_day + price + hrstr:working_day, 
            data = readings.aggregate,
            weights = readings.aggregate$agg_count)
poly2.adj.r2 <- summary(poly2)$adj.r.squared
poly2.visreg <- visreg(fit = poly2, xvar = "temperature", type = "conditional",
                       alpha = 1, # no confidence interval
                       line.par = list(col = "red", lwd = 3),
                       points.par = list(col = scatter.black),
                       main = "Quadratic Polynomial Fit of Temperature to Conditional Demand", 
                       ylab = "Conditional Electricity Demand, log(kWh)",
                       xlab = "Outdoor Temperature (Celsius)")

poly3 <- lm(log(kwh) ~ poly(temperature, degree = 3) + hrstr + working_day + price + hrstr:working_day, 
            data = readings.aggregate,
            weights = readings.aggregate$agg_count)
poly3.adj.r2 <- summary(poly3)$adj.r.squared
poly3.visreg <- visreg(fit = poly3, xvar = "temperature", type = "conditional",
                       alpha = 1, # no confidence interval
                       line.par = list(col = "red", lwd = 3),
                       points.par = list(col = scatter.black),
                       main = "Cubic Polynomial Fit of Temperature to Conditional Demand", 
                       ylab = "Conditional Electricity Demand, log(kWh)",
                       xlab = "Outdoor Temperature (Celsius)")

poly4 <- lm(log(kwh) ~ poly(temperature, degree = 4) + hrstr + working_day + price + hrstr:working_day, 
            data = readings.aggregate,
            weights = readings.aggregate$agg_count)
poly4.adj.r2 <- summary(poly4)$adj.r.squared
poly4.visreg <- visreg(fit = poly4, xvar = "temperature", type = "conditional",
                       alpha = 1, # no confidence interval
                       line.par = list(col = "red", lwd = 3),
                       points.par = list(col = scatter.black),
                       main = "Degree 4 Polynomial Fit of Temperature to Conditional Demand", 
                       ylab = "Conditional Electricity Demand, log(kWh)",
                       xlab = "Outdoor Temperature (Celsius)")

poly5 <- lm(log(kwh) ~ poly(temperature, degree = 5) + hrstr + working_day + price + hrstr:working_day, 
            data = readings.aggregate,
            weights = readings.aggregate$agg_count)
poly5.adj.r2 <- summary(poly5)$adj.r.squared
poly5.visreg <- visreg(fit = poly5, xvar = "temperature", type = "conditional",
                       alpha = 1, # no confidence interval
                       line.par = list(col = "red", lwd = 3),
                       points.par = list(col = scatter.black),
                       main = "Degree 5 Polynomial Fit of Temperature to Conditional Demand", 
                       ylab = "Conditional Electricity Demand, log(kWh)",
                       xlab = "Outdoor Temperature (Celsius)")

dev.print(file="../../Figures/TemperaturePolynomialPlots.png", device=png, height = 1200, width = 1600)
par(mfrow = c(1, 1))
# Cubic polynomial (ie. degree = 3) is the highest, stable polynomial



#####
# Natural cubic splines
#####
lm.ns <- lm(log(kwh) ~ ns(temperature, knots = c(-1, 17, 32)) + hrstr + working_day + price + hrstr:working_day, 
            data = readings.aggregate,
            weights = readings.aggregate$agg_count)
ns.adj.r2 <- summary(lm.ns)$adj.r.squared
ns.adj.r2
ns.visreg <- visreg(fit = lm.ns, xvar = "temperature", type = "conditional",
                       alpha = 1, # no confidence interval
                       line.par = list(col = "red", lwd = 3),
                       points.par = list(col = scatter.black),
                       main = "Natural Cubic Splines Fit of Temperature to Conditional Demand", 
                       ylab = "Conditional Electricity Demand, log(kWh)",
                       xlab = "Outdoor Temperature (Celsius)")
ns.temp <- ns(readings.aggregate$temperature, knots = c(-1, 17, 32))
abline(v = attr(ns.temp, "knots"), col = "blue", lty = 2, lwd = 1)
legend("topleft", lty = c(1, 2), lwd = c(3, 1), col = c("red", "blue"), legend = c("Natural Cubic Spline Fit", "Knot Positions"))
dev.print(file="../../Figures/TemperatureNatrualSplinesPlots.png", device=png, height = 600, width = 800)



#####
# Examine correlation of past temperatures to determine number of lags for Distributed Lag Non-Linear Model (DLNM) 
# based on procedure described in Almon (1965).
#####
temps <- readings.aggregate$temperature
correlation.hrs <- 23
temp.lagmatrix <- CreatePastTemperatureMatrix(nlags = correlation.hrs, df.readings = readings.aggregate)

# Trim out NAs from dependent and independent variables equally
kwh.trimmed <- readings.aggregate$kwh[(correlation.hrs+1):nrow(readings.aggregate)]
temp.lagmatrix.trimmed <- tail(x = temp.lagmatrix, n = (nrow(temp.lagmatrix) - correlation.hrs))

# Correlation matrix of the response variable and lags of temperature
cor(y = kwh.trimmed, x = temp.lagmatrix.trimmed)

# As per Almon (1965), look for lag that is as correlated with response as lag=0
# Viewing the above correlation matrix, this appears to be lag=5.
max.lag <- 5
rm(kwh.trimmed, temp.lagmatrix.trimmed)



#####
# Distributed Lag Non-linear Model (DLNM) of lagged temperatures' effects 
# on electricity consumption.
#####
# lag = 12, 
# arglag = list(fun = "lin")
# adj.r2 = 0.9138482
#
# lag = 23, 
# arglag = list(fun = "ns", df = 3)
# adj.r2 = 0.9258107
#
# lag = 23, 
# arglag = list(fun = "ns", df = 4))
# adj.r2 = 0.9263082

# See if using breakpoints as threshold values improves external-response
cb <- crossbasis(x = readings.aggregate$temperature, 
                 argvar = list(fun = "ns", knots = c(-1, 17, 32)),
                 lag = 6, 
                 arglag = list(fun = "lin"))

# Prediction using the crossbasix matrix
lm.cb <- lm(log(kwh) ~ cb + hrstr + working_day + price + hrstr:working_day,
            data = readings.aggregate)
summary(lm.cb)$adj.r.squared
cb.pred <- crosspred(cb, lm.cb)
plot(cb.pred, 
     ptype = "3d", 
     main = "3D Plot of Exposure-Lag-Response Association", 
     xlab = "Temperature (Celsius)", 
     ylab = "Lag (hours)", 
     zlab = "log(kWh)", 
     theta = -40, 
     ltheta = 180)
#dev.print(file = "../../Figures/ExposureLagResponseAssociation3dPlot.png", device = png, height = 600, width = 800)
plot(cb.pred, ptype="contour", 
     key.title = title("log(kWh)"), 
     plot.title = title("Contour Plot of Exposure-Lag-Response Association", 
                      xlab = "Temperature (Celsius)", 
                      ylab = "Lag (hours)"))
#dev.print(file = "../../Figures/ExposureLagResponseAssociationContourPlot.png", device = png, height = 600, width = 800)



par(mfrow = c(1, 1))
plot(cb.pred, "overall", 
     main = "Overall Exposure-Lag-Response Association", 
     xlab = "Temperature", 
     ylab = "log(kWh)")

plot(cb.pred, "slices", 
     lag = 0, 
     ylim = c(-0.1, 0.6), 
     xlim = c(min(readings.aggregate$temperature), max(readings.aggregate$temperature)),
     main = "Effects of Temperature at Various Time Lags", 
     xlab = "Temperature (Celsius)", 
     ylab = "log(kWh)", 
     col = "blue", 
     ci = "n",
     ci.arg = list(lty=5))
lines(cb.pred, "slices",
      lag = 3, 
      type = "l",
      col = "green")
lines(cb.pred, "slices",
      lag = 5,
      type = "l",
      col = "red")

plot(cb.pred, "slices", 
     ylim = c(-0.1, 0.6),
     var = -10, 
     main = "Effects of Time Lag at Various Temperatures", 
     xlab = "Lag (hours)", 
     ylab = "log(kWh)", 
     type = "o", 
     pch = 19, 
     ci = "n",
     col = "blue")
lines(cb.pred, "slices",
      var = 13,
      type = "o",
      pch = 19,
      col = "green")
lines(cb.pred, "slices",
      var = 16,
      type = "o",
      pch = 19,
      col = "yellow")
lines(cb.pred, "slices",
      var = 19,
      type = "o",
      pch = 19,
      col = "orange")
lines(cb.pred, "slices",
      var = 35,
      type = "o",
      pch = 19,
      col = "red")



#####
# Run a Newey-West test to report heteroskedastically consistent standard errors and p-values
# Looking at the hrstr:working_day interaction, it appears hours 2, 3, 6-8, 10-17, and 21-23 
# are significant at p < 0.05.
#
# I may encorporate this into my model as 6-17 (work day) and 21-23 (bed time) interactions 
# rather than using them all. This may avoid some overfitting. Though it feels a bit 
# wrong ignoring the reported p-values for 2-3 in the morning. People seem to stay up later 
# on non-working days (party animals!).
#
# TODO(r24mille): Rather than just looking at p-values, examine changes to Adj. R^2, BIC, and 
#                 validation results when fewer hours of interaction are used.
#####
print(coeftest(lm.cb, vcov = NeweyWest(lm.cb, lag = 5)))



#####
# Plot residuals for analysis
#####
summary(lm.cb)$adj.r.squared
BIC(lm.cb)
dwtest(lm.cb)
residuals.bktr <- readings.aggregate[c(13:nrow(readings.aggregate)),]$kwh - exp(predict(lm.cb))
plot(density(residuals.bktr, bw = 0.014),
     main = "Distribution of Residuals \n(backtransformed, no seasonality term)",
     ylim = c(0,4.25),
     xlim = c(-1.6,0.5),
     lwd = 3)
dev.print(file = "../../Figures/ResidualComparison/NoSeasonBacktransformed/ResidualDistribution.png", device = png, height = 600, width = 800)

plot(x = exp(predict(lm.cb)), y = residuals.bktr,
     main = "Distribution of Residuals as a Function of Predicted Value \n(backtransformed, no seasonality term)",
     xlab = "Predicted Value (kWh)",
     ylab = "Residual",
     ylim = c(-1.8, 1.1))
dev.print(file = "../../Figures/ResidualComparison/NoSeasonBacktransformed/ResidualDiagnostics.png", device = png, height = 600, width = 800)

plot(residuals.bktr,
     main = "Serial Correlation of Residuals \n(backtransformed, no seasonality term)",
     ylab = "Residual",
     xlab = "Timeseries Index",
     ylim = c(-1.8, 1.1))
dev.print(file = "../../Figures/ResidualComparison/NoSeasonBacktransformed/ResidualSerialCorrelation.png", device = png, height = 600, width = 800)






#####
# Check variance inflation factor (VIF) of intuitive main effects. This test quantifies  
# multicolinearity in the data.
#####
# Test main effects model for multicolinearity with Variance Inflation Factor
lm.intuitive.me <- lm(log(kwh) ~ poly(temperature, degree = 3) + month + working_day + hrstr + price + severe_weather,
                      data = readings.aggregate)
vif(lm.intuitive.me) # Month has a high at VIF=9.18 and is colinear with temperature

# Introduce features implied by month that are not directly related to temperature. 
# Features are inspired by prior work. 
# TODO(r24mille): Find proper citations in notes.
readings.aggregate$schoolyr <- rep(TRUE, nrow(readings.aggregate))
readings.aggregate$schoolyr[readings.aggregate$month %in% c("m7", "m8")] <- FALSE

# Try replacing month with electricity rate seaons, even before TOU the flat rates were divided into two seasons
readings.aggregate$rateseason <- rep(NA, nrow(readings.aggregate))
readings.aggregate$rateseason[readings.aggregate$month %in% c("m5", "m6", "m7", "m8", "m9", "m10")] <- "summer"
readings.aggregate$rateseason[readings.aggregate$month %in% c("m11", "m12", "m1", "m2", "m3", "m4")] <- "winter"
readings.aggregate$rateseason <- factor(readings.aggregate$rateseason, c("winter", "summer"))

# Try replacing month with calendar seasons
readings.aggregate$season <- rep(NA, nrow(readings.aggregate))
readings.aggregate$season[readings.aggregate$month %in% c("m12", "m1", "m2")] <- "winter"
readings.aggregate$season[readings.aggregate$month %in% c("m3", "m4", "m5")] <- "spring"
readings.aggregate$season[readings.aggregate$month %in% c("m6", "m7", "m8")] <- "summer"
readings.aggregate$season[readings.aggregate$month %in% c("m9", "m10", "m11")] <- "autumn"
readings.aggregate$season <- factor(readings.aggregate$season, c("winter", "spring", "summer", "autumn"))