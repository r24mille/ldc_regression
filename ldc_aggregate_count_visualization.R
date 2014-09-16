# Load count of SmartMeterReadings data from CSV
home <- setwd(Sys.getenv("HOME"))
fpath <- file.path(home, 
                   "../Dropbox/ISS4E/R/", 
                   "observation_count_series.csv")
readings.counts <- read.csv(fpath)

# Convert the observation_timestamp column to Date
readings.counts$observation_timestamp <- as.POSIXlt(readings.counts$observation_timestamp)

# Make a simple plot of the data
plot(readings.counts$observation_timestamp, 
     readings.counts$number_observations,
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
