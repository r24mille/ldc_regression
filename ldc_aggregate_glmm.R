library(lme4)
library(plyr) # Summarize dataframe
library(nlme) # groupedData
library(lattice) # Plotting
library(stargazer) # LaTeX tables

# Load SmartMeterReading data from CSV
home <- setwd(Sys.getenv("HOME"))
fpath <- file.path(home, 
                   "../Dropbox/ISS4E/R/", 
                   "full_aggregate_readings.csv")
readings.aggregate <- read.csv(fpath)

# Re-orders TOU Period levels so that graphs are sorted accordingly
readings.aggregate$tou_period <- factor(readings.aggregate$tou_period, 
                              c("off_weekend", "off_morning", "mid_morning", 
                                "on_peak", "mid_evening", "off_evening"))

##
# A linear model will be build using lm(...) for each model variant so that 
# we can find an R-squared value. Also, a general linear model 
# will be fit through maximum likelihood estimation using glm(...). This 
# provides additional metrics such as Akaike Information Criterion (AIC) to 
# measure model quality.
##

# Temperature, TOU Period, and billing method modeled as fixed effects using 
# using coefficients for main effects only.
model.lm.fe.main.tmp_prd_bill <- lm(kwh ~ temperature + tou_period + billing_active, 
                                    data = readings.aggregate)
model.glm.fe.main.tmp_prd_bill <- glm(kwh ~ temperature + tou_period + billing_active, 
                                      data = readings.aggregate)
summary(model.lm.fe.main.tmp_prd_bill)
summary(model.glm.fe.main.tmp_prd_bill)

# Temperature, TOU Period, and billing method modeled as fixed effects using 
# using coefficients of main effects and all interactions.
model.lm.fe.inter.tmp_prd_bill <- lm(kwh ~ temperature * tou_period * billing_active, 
                                         data = readings.aggregate)
model.glm.fe.inter.tmp_prd_bill <- glm(kwh ~ temperature * tou_period * billing_active, 
                                           data = readings.aggregate)
summary(model.lm.fe.inter.tmp_prd_bill)
summary(model.glm.fe.inter.tmp_prd_bill)

# Temperature, hour of day, and billing method modeled as fixed effects using 
# using coefficients for main effects only.
model.lm.fe.main.tmp_hr_bill <- lm(kwh ~ temperature + hour + billing_active, 
                                    data = readings.aggregate)
model.glm.fe.main.tmp_hr_bill <- glm(kwh ~ temperature + hour + billing_active, 
                                      data = readings.aggregate)
summary(model.lm.fe.main.tmp_hr_bill)
summary(model.glm.fe.main.tmp_hr_bill)

# Temperature, hour of day, and billing method modeled as fixed effects using 
# using coefficients of main effects and all interactions.
model.lm.fe.inter.tmp_hr_bill <- lm(kwh ~ temperature * hour * billing_active, 
                                     data = readings.aggregate)
model.glm.fe.inter.tmp_hr_bill <- glm(kwh ~ temperature * hour * billing_active, 
                                       data = readings.aggregate)
summary(model.lm.fe.inter.tmp_hr_bill)
summary(model.glm.fe.inter.tmp_hr_bill)