library(lme4)
library(plyr) # Summarize dataframe
library(nlme) # groupedData
library(lattice) # Plotting
library(stargazer) # LaTeX tables

# Load SmartMeterReading data from CSV
home <- setwd(Sys.getenv("HOME"))
fpath <- file.path(home, 
                   "../Dropbox/ISS4E/R/25_meterid_sample/", 
                   "25_random_sample.csv")
readings <- read.csv(fpath)

# Re-orders TOU Period levels so that graphs are sorted accordingly
readings$tou_period <- factor(readings$tou_period, 
                              c("off_weekend", "off_morning", "mid_morning", 
                                "on_peak", "mid_evening", "off_evening"))

# Restructure meters into aggregate sample (ie. mean kWh for each hourindex).
readings.aggregate <- ddply(readings, 
                            ~hourindex, 
                            summarize, 
                            kwh = mean(kwh), 
                            daynum = unique(daynum), 
                            hour = unique(hour), 
                            temperature = unique(temperature), 
                            tou_period = unique(tou_period),
                            billing_active = unique(billing_active), 
                            agg_count = length(subject))
