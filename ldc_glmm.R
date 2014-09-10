library(lme4)


# Load SmartMeterReading data from CSV
home <- setwd(Sys.getenv("HOME"))
fpath <- file.path(home, 
                   "../Dropbox/ISS4E/R/25_meterid_sample/", 
                   "25_random_sample.csv")
readings <- read.csv(fpath)

# Re-orders TOU Period levels so that graphs are sorted accordingly
readings$tou_period <- factor(readings$tou_period, c("off_morning", "mid_morning", "on_peak", "mid_evening", "off_evening", "off_weekend"))
attach(readings)

# Playing with simple linear model
the_model <- lmer(kwh ~ temperature 
                  + (1|subject))
summary(the_model)
anova(the_model)
plot(the_model)


# Null hypothesis is that TOU billing has no effect
readings.null = lmer(kwh ~ temperature * tou_period 
                     + (1+billing_active|subject), 
                     data = readings, 
                     REML = FALSE)
summary(readings.null)

# Alternate hypothesis is that TOU billing has a significant effect
readings.model <- lmer(kwh ~ temperature * tou_period * billing_active 
                     + (1+billing_active|subject), 
                     data = readings, 
                     REML = FALSE)
summary(readings.model)

# Analysis of variance
anova(readings.model)