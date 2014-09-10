# Use ODBC for connection
library(RODBC)
channel <- odbcConnect("essex_annotated")

# Set up common color scheme
transRed <- rgb(100, 0, 0, 25, maxColorValue=100) # Red with 25% opacity
darkRed <- rgb(100, 0, 0, 100, maxColorValue=100) # Red with 100% opacity
transGreen <- rgb(0, 33, 0, 25, maxColorValue=100) # Green with 25% opacity
darkGreen <- rgb(0, 33, 0, 100, maxColorValue=100) # Green with 100% opacity

# Get weather data
summer2011 <- sqlQuery(channel, "select observation_datetime_dst, temp_metric 
                        from weathertables.wunderground_observation 
                        where location_id = 13 
                        and observation_datetime_dst >= '2011-05-01' 
                        and observation_datetime_dst <= '2011-10-31' 
                        order by observation_datetime_dst asc")

summer2012 <- sqlQuery(channel, "select observation_datetime_dst, temp_metric 
                        from weathertables.wunderground_observation 
                        where location_id = 13 
                        and observation_datetime_dst >= '2012-05-01' 
                        and observation_datetime_dst <= '2012-10-31' 
                        order by observation_datetime_dst asc")

# Plot raw/smooth timeseries data
# Before and after TOU years
plot(summer2011$temp_metric ~ summer2011$observation_datetime_dst, 
     xlab = "Date (Hourly)", 
     ylab = "Outdoor Temperature (Celsius)", 
     ylim = c(0, 40),    # Y axis limit, 0-40
     type = "l",         # Type, l is for line
     col  = transRed
     )
with(summer2011, 
     lines(loess.smooth(observation_datetime_dst, temp_metric, span=1/40), # Loess Smoothing
           col = darkRed)
     )

par(new = TRUE) # Don't redraw canvas with new plot

plot(summer2012$observation_datetime_dst, 
     summer2012$temp_metric, 
     xlab = "Date (Hourly)", 
     ylab = "Outdoor Temperature (Celsius)", 
     ylim = c(0, 40),    # Y axis limit, 0-40
     type = "l",     # Type, l is for line
     col  = transGreen
)
with(summer2012, 
     lines(loess.smooth(observation_datetime_dst, temp_metric, span=1/40), # Loess Smoothing
           col = darkGreen)
)

title(main = "Summer Temperature Comparison in LDC Region")

legend("topright", 
       legend = c("2011 Observed", 
                  "2011 Smoothed", 
                  "2012 Observed", 
                  "2012 Smoothed"), 
       lty  = 1,
       col  = c(transRed,
                darkRed,
                transGreen, 
                darkGreen)
)


# Histogram of data in two separate plots
par(mfrow=c(1,2))
hist(summer2011$temp_metric,
     xlab = "Hourly Outdoor Temperature (Celsius)",
     xlim = c(0, 40), 
     ylim = c(0, 550), 
     col = transRed,
     main = "Summer 2011 Temperature Histogram"
)
hist(summer2012$temp_metric,
     xlab = "Hourly Outdoor Temperature (Celsius)", 
     xlim = c(0, 40), 
     ylim = c(0, 550), 
     col = transGreen,
     main = "Summer 2012 Temperature Histogram"
)

# Combined boxplot of hourly temperature data
par(mfrow=c(1,1))
years <- numeric(8786)
for (i in 1:8786) {
  if (i <= 4393) {
    years[i] <- 2011
  } else {
    years[i] <- 2012
  } 
}

boxDf = data.frame(years, 
                   c(summer2011$temp_metric, summer2012$temp_metric)
)
colnames(boxDf) <- c("Year", "Temp")

boxplot(boxDf$Temp ~ boxDf$Year,
     ylim = c(0, 40),
     ylab = "Hourly Outdoor Temperature (Celsius)", 
     xlab = "Year",
     main = "Boxplot of Summer Temperatures",
     col  = c(transRed, transGreen)
)

# Written summary of temperatures
summary(summer2011$temp_metric)
summary(summer2012$temp_metric)