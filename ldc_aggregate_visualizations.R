home <- setwd(Sys.getenv("HOME"))
fpath <- file.path(home, 
                   "../Dropbox/ISS4E/R/", 
                   "full_aggregate_readings.csv")
readings.aggregate <- read.csv(fpath)

# Convert the observation_timestamp column to Date
readings.aggregate$timestamp_dst <- as.POSIXlt(readings.aggregate$timestamp_dst)

##
# Make a simple plot of number of households involved in each hourly measurement
plot(readings.aggregate$timestamp_dst, 
     readings.aggregate$agg_count,
     main = "Count of Observations Each Hour", 
     ylab = paste("Number of Households with a Reading (max=", 
                  max(readings.counts$number_observations), 
                  ")"),
     xlab = "Date (hourly)", 
     type = "p", # Points
     pch = 20, # Plot character
     col = rgb(0, 0, 0, 50, maxColorValue=255), 
     xaxt = "n" # x-axis type "n" turns of ticks
     )
tsrange <- as.POSIXct(round(range(readings.counts$observation_timestamp), 
                      "days"))
axis.POSIXct(1, 
             at = seq(tsrange[1], tsrange[2], by="month"), 
             format="%m/%Y")

##
# Make a simple plot of the timeseries data
plot(readings.aggregate$timestamp_dst, 
     readings.aggregate$kwh,
     main = "Timeseries of Aggregate Average Reading", 
     ylab = paste("Average Readings (kWh)"),
     xlab = "Date (hourly)", 
     type = "l", # Line
     col = rgb(0, 0, 0, 1, maxColorValue=1), 
     xaxt = "n" # x-axis type "n" turns of ticks
)
tsrange <- as.POSIXct(round(range(readings.counts$observation_timestamp), 
                            "days"))
axis.POSIXct(1, 
             at = seq(tsrange[1], tsrange[2], by="month"), 
             format="%m/%Y")