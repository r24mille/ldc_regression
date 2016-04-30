require(ggplot2)        # Better plotting!
require(car)            # Companion to Applied Regression (better qqPlot)
require(segmented)      # For finding temperature breakpoint
require(gridExtra)      # Placing multiple ggplot2 plots into a single viewport or image

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

# Extract working_day factor
readings.aggregate$working_day <- ((readings.aggregate$holiday == TRUE | 
                                      readings.aggregate$dayname == "Sat" | 
                                      readings.aggregate$dayname == "Sun") == FALSE)

# For clarity, reorder explanatory variables which are factors
readings.aggregate <- OrderFactors(readings.aggregate)
rm(fpath, fpath2, home)



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
# Temperature breakpoint plots
#####
par(mfrow = c(1, 3))
lm.lintemp <- lm(kwh ~ temperature + hrstr + working_day + hrstr:working_day, 
                 data = readings.aggregate)
lintemp.adj.r2 <- summary(lm.lintemp)$adj.r.squared
seg1 <- segmented(obj = lm.lintemp,
                 seg.Z = ~ temperature,
                 psi = 17)
temp.break <- seg1$psi[1,2]
seg1.adj.r2 <- summary(seg1)$adj.r.squared
plot(seg1, col = "red", lwd = 3, 
     res = TRUE, res.col = rgb(0, 0, 0, 25, maxColorValue = 100), pch = 20, 
     ylim = c(0.2, 3.45), 
     main = "Switching Regression",
     xlab = "Outdoor Temperature (Celsius)",
     ylab = "Electricity Demand (kWh)",
     rug = FALSE)

plot(seg1, col = rgb(0, 0, 0, 0, maxColorValue = 100), lwd = 0, 
     res = TRUE, res.col = rgb(0, 0, 0, 25, maxColorValue = 100), pch = 20, 
     ylim = c(0.2, 3.45), 
     main = "Threshold Regression",
     xlab = "Outdoor Temperature (Celsius)",
     ylab = "Electricity Demand (kWh)",
     rug = FALSE)

plot(seg1, col = rgb(0, 0, 0, 0, maxColorValue = 100), lwd = 0, 
     res = TRUE, res.col = rgb(0, 0, 0, 25, maxColorValue = 100), pch = 20, 
     ylim = c(0.2, 3.45), 
     main = "Threshold Regression with Saturation at Extremes",
     xlab = "Outdoor Temperature (Celsius)",
     ylab = "Electricity Demand (kWh)",
     rug = FALSE)
dev.print(file="../../Figures/ConditionedPiecewiseLinearModels.png", device=png, height = 600, width = 2400)
par(mfrow = c(1, 1))