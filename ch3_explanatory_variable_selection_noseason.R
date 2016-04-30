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
require(pracma)         # De-trending experiments

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
readings.aggregate$working_day <- ((readings.aggregate$holiday == TRUE | 
                                      readings.aggregate$dayname == "Sat" | 
                                      readings.aggregate$dayname == "Sun") == FALSE)
summary(readings.aggregate$severe_weather) # 146 hours of severe conditions is still arguably low

# For clarity, reorder explanatory variables which are factors
readings.aggregate <- OrderFactors(readings.aggregate)
rm(fpath, fpath2, home)


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
     ylab = "Conditional Electricity Demand (kWh)",
     rug = FALSE)

dev.print(file="../../Figures/TemperatureBreakpointPlots.png", device=png, height = 600, width = 1600)
par(mfrow = c(1, 1))



#####
# Examine simple plot of polynomials to determine appropriate degree.
#####
scatter.black = rgb(0, 0, 0, 50, maxColorValue = 100)
par(mfrow = c(2, 2))
poly2 <- lm(kwh ~ poly(temperature, degree = 2) + hrstr + working_day + price + hrstr:working_day, 
            data = readings.aggregate,
            weights = readings.aggregate$agg_count)
poly2.adj.r2 <- summary(poly2)$adj.r.squared
poly2.visreg <- visreg(fit = poly2, xvar = "temperature", type = "conditional",
                       alpha = 1, # no confidence interval
                       line.par = list(col = "red", lwd = 3),
                       points.par = list(col = scatter.black),
                       main = "Quadratic Polynomial Fit of Temperature to Conditional Demand", 
                       ylab = "Conditional Electricity Demand (kWh)",
                       xlab = "Outdoor Temperature (Celsius)")

poly3 <- lm(kwh ~ poly(temperature, degree = 3) + hrstr + working_day + price + hrstr:working_day, 
            data = readings.aggregate,
            weights = readings.aggregate$agg_count)
poly3.adj.r2 <- summary(poly3)$adj.r.squared
poly3.visreg <- visreg(fit = poly3, xvar = "temperature", type = "conditional",
                       alpha = 1, # no confidence interval
                       line.par = list(col = "red", lwd = 3),
                       points.par = list(col = scatter.black),
                       main = "Cubic Polynomial Fit of Temperature to Conditional Demand", 
                       ylab = "Conditional Electricity Demand (kWh)",
                       xlab = "Outdoor Temperature (Celsius)")

poly4 <- lm(kwh ~ poly(temperature, degree = 4) + hrstr + working_day + price + hrstr:working_day, 
            data = readings.aggregate,
            weights = readings.aggregate$agg_count)
poly4.adj.r2 <- summary(poly4)$adj.r.squared
poly4.visreg <- visreg(fit = poly4, xvar = "temperature", type = "conditional",
                       alpha = 1, # no confidence interval
                       line.par = list(col = "red", lwd = 3),
                       points.par = list(col = scatter.black),
                       main = "Degree 4 Polynomial Fit of Temperature to Conditional Demand", 
                       ylab = "Conditional Electricity Demand (kWh)",
                       xlab = "Outdoor Temperature (Celsius)")

poly5 <- lm(kwh ~ poly(temperature, degree = 5) + hrstr + working_day + price + hrstr:working_day, 
            data = readings.aggregate,
            weights = readings.aggregate$agg_count)
poly5.adj.r2 <- summary(poly5)$adj.r.squared
poly5.visreg <- visreg(fit = poly5, xvar = "temperature", type = "conditional",
                       alpha = 1, # no confidence interval
                       line.par = list(col = "red", lwd = 3),
                       points.par = list(col = scatter.black),
                       main = "Degree 5 Polynomial Fit of Temperature to Conditional Demand", 
                       ylab = "Conditional Electricity Demand (kWh)",
                       xlab = "Outdoor Temperature (Celsius)")

dev.print(file="../../Figures/TemperaturePolynomialPlots.png", device=png, height = 1200, width = 1600)
par(mfrow = c(1, 1))
# Cubic polynomial (ie. degree = 3) is the highest, stable polynomial



#####
# Natural cubic splines
#####
lm.ns <- lm(kwh ~ ns(temperature, knots = c(3, 23, 30)) + hrstr + working_day + hrstr:working_day, 
            data = readings.aggregate)
ns.adj.r2 <- summary(lm.ns)$adj.r.squared
ns.adj.r2
ns.visreg <- visreg(fit = lm.ns, xvar = "temperature", type = "conditional",
                       alpha = 1, # no confidence interval
                       line.par = list(col = "red", lwd = 3),
                       points.par = list(col = ScatterBlack()),
                       main = "Natural Cubic Splines Fit of Temperature to Conditional Demand", 
                       ylab = "Conditional Electricity Demand (kWh)",
                       xlab = "Outdoor Temperature (Celsius)")
ns.temp <- ns(readings.aggregate$temperature, knots = c(2, 24, 30))
abline(v = attr(ns.temp, "knots"), col = "blue", lty = 2, lwd = 1)
legend("topleft", lty = c(1, 2), lwd = c(3, 1), col = c("red", "blue"), legend = c("Natural Cubic Spline Fit", "Knot Positions"))
dev.print(file="../../Figures/TemperatureTransformationNaturalSplines.png", device=png, height = 600, width = 800)


#####
# Natural Splines diagnostics to compare with Distributed Lag Nonliner Models
#####
summary(lm.ns)$adj.r.squared
BIC(lm.ns)
dwtest(lm.ns)
plot(density(lm.ns$residuals, bw = 0.014),
     main = "Distribution of Residuals \n(no DLNM, no kWh transform, no seasonality term)",
     ylim = c(0,4.75),
     xlim = c(-1.6,0.5),
     lwd = 3)
dev.print(file = "../../Figures/ResidualComparison/NoSeasonNaturalSplines/ResidualDistribution.png", device = png, height = 600, width = 800)

#par(mfrow = c(2, 2))
#plot(lm.ns)
#par(mfrow = c(1, 1))
plot(x = lm.ns$fitted, y = lm.ns$residuals,
     main = "Distribution of Residuals as a Function of Predicted Value \n(no DLNM, no kWh transform, no seasonality term)",
     xlab = "Predicted Value (kWh)",
     ylab = "Residual",
     ylim = c(-1.8, 1.1))
dev.print(file = "../../Figures/ResidualComparison/NoSeasonNaturalSplines/ResidualDiagnostics.png", device = png, height = 600, width = 800)

plot(lm.ns$residuals,
     main = "Serial Correlation of Residuals \n(no DLNM, no kWh transform, no seasonality term)",
     ylab = "Residual",
     xlab = "Timeseries Index",
     ylim = c(-1.8, 1.1))
dev.print(file = "../../Figures/ResidualComparison/NoSeasonNaturalSplines/ResidualSerialCorrelation.png", device = png, height = 600, width = 800)



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
                 argvar = list(fun = "ns", knots = c(2, 24, 30)),
                 lag = 10, 
                 arglag = list(fun = "poly", degree=3))

# Prediction using the crossbasix matrix
lm.cb <- lm(kwh ~ cb + hrstr + working_day + price + hrstr:working_day,
            data = readings.aggregate)
dwtest(lm.cb, alt="two.sided")
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
# Plot residuals for analysis
#####
summary(lm.cb)$adj.r.squared
BIC(lm.cb)
dwtest(lm.cb)
plot(density(lm.cb$residuals, bw = 0.014),
     main = "Distribution of Residuals \n(no kWh transform, no seasonality term)",
     ylim = c(0,4.75),
     xlim = c(-1.6,0.5),
     lwd = 3)
dev.print(file = "../../Figures/ResidualComparison/NoSeason/ResidualDistribution.png", device = png, height = 600, width = 800)

#par(mfrow = c(2, 2))
#plot(lm.cb)
#par(mfrow = c(1, 1))
plot(x = lm.cb$fitted, y = lm.cb$residuals,
     main = "Distribution of Residuals as a Function of Predicted Value \n(no kWh transform, no seasonality term)",
     xlab = "Predicted Value (kWh)",
     ylab = "Residual",
     ylim = c(-1.8, 1.1))
dev.print(file = "../../Figures/ResidualComparison/NoSeason/ResidualDiagnostics.png", device = png, height = 600, width = 800)

plot(lm.cb$residuals,
     main = "Serial Correlation of Residuals \n(no kWh transform, no seasonality term)",
     ylab = "Residual",
     xlab = "Timeseries Index",
     ylim = c(-1.8, 1.1))
dev.print(file = "../../Figures/ResidualComparison/NoSeason/ResidualSerialCorrelation.png", device = png, height = 600, width = 800)






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