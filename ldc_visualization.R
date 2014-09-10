options(scipen=5)

# Load SmartMeterReading data from CSV
home <- setwd(Sys.getenv("HOME"))
fpath <- file.path(home, 
                   "../Dropbox/ISS4E/R/250_meterid_sample/", 
                   "250_random_sample.csv")
readings <- read.csv(fpath)
readings$tou_period <- factor(readings$tou_period, c("off_morning", "mid_morning", "on_peak", "mid_evening", "off_evening", "off_weekend"))

# Summarize Data
summary(readings$kwh)

# Boxplot separated by TOU Period
boxplot(readings$kwh ~ readings$tou_period,
        main = "Boxplot of Demand during TOU Period Across all Households",
        ylab = "Hourly Meter Reading (kWh)",
        xlab = "TOU Period",
        pch  = "-")

# Histogram of Hourly meter readings
hist(readings$kwh,
     main = "Histogram of Hourly Readings (2M observations, 150 bins)", 
     xlab = "Hourly Meter Reading (kWh)", 
     breaks = 150
)