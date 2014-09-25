library(lme4)
library(plyr) # Summarize dataframe
library(nlme) # groupedData
library(lattice) # Plotting
library(stargazer) # LaTeX tables
library(bestglm) # Compare models using leaps
library(glmulti)


# Source the function in another file
source('piecewise_cdh_setpoint.R')
source('reg_subset_visualization.R')

# Load SmartMeterReading data from CSV
home <- Sys.getenv("HOME")
fpath <- file.path(home, 
                   "../Dropbox/ISS4E/R/", 
                   "full_aggregate_readings.csv")
readings.aggregate <- read.csv(fpath)

# Re-orders TOU Period levels so that graphs are sorted accordingly
readings.aggregate$tou_period <- factor(readings.aggregate$tou_period, 
                              c("off_weekend", "off_morning", "mid_morning", 
                                "on_peak", "mid_evening", "off_evening"))

# Represent hours as levels rather than integers
readings.aggregate$hrstr <- paste0("h", readings.aggregate$hour)
readings.aggregate$hrstr <- factor(readings.aggregate$hrstr, 
                                   c("h0", "h1", "h2", "h3", "h4", "h5", "h6", 
                                     "h7", "h8", "h9", "h10", "h11", "h12", 
                                     "h13", "h14", "h15", "h16", "h17", "h18", 
                                     "h19", "h20", "h21", "h22", "h23"))

##
# Add yes/no weekend flag
readings.aggregate$weekend <- ifelse(readings.aggregate$tou_period %in% c("off_weekend"),
                                     "Yes",
                                     "No")

readings.aggregate$weekend <- factor(readings.aggregate$weekend, 
                                     unique(readings.aggregate$weekend))

##
# Add column which converts TOU Period to its price level
readings.aggregate$price <- readings.aggregate$tou_period
readings.aggregate$price <- gsub("off_weekend", "off_peak", readings.aggregate$price)
readings.aggregate$price <- gsub("off_morning", "off_peak", readings.aggregate$price)
readings.aggregate$price <- gsub("off_evening", "off_peak", readings.aggregate$price)
readings.aggregate$price <- gsub("mid_morning", "mid_peak", readings.aggregate$price)
readings.aggregate$price <- gsub("mid_evening", "mid_peak", readings.aggregate$price)
readings.aggregate$price <- factor(readings.aggregate$price, 
                                   c("off_peak", "mid_peak", "on_peak"))

##
# Add column which represents "month" as a categorical factor
readings.aggregate$timestamp_dst <- as.POSIXlt(readings.aggregate$timestamp_dst)
readings.aggregate$month <- paste0("m", (readings.aggregate$timestamp_dst$mon + 1))
readings.aggregate$month <- factor(readings.aggregate$month, 
                                   c("m5", "m6", "m7", "m8", "m9", "m10"))

## 
# Section finds the Cooling-Degree Hour (CDH) baseline for 10th, middle, and 
# 90th percentile observations for each temperature "bin" similar to the 
# three-line model.
##

# Determine kwh quantils for each temperature "bin"
quants <- ddply(readings.aggregate,
                   .(temperature), 
                   function(x) quantile(x$kwh, c(.1, .5, .9))
)
# Hacky due to temperatures matching quant indeces (by chance).
# TODO: Fix with a proper subset(...) I'm tired of fighting it for today.
readings.aggregate.tenth <- subset(readings.aggregate, 
                                   kwh <= quants[temperature, 2])
readings.aggregate.middle <- subset(readings.aggregate, 
                                    kwh > quants[temperature, 2] & kwh < quants[temperature, 4])
readings.aggregate.ninetieth <- subset(readings.aggregate, 
                                      kwh >= quants[temperature, 4])

# Determine the A/C setpoint threshold
piecewise.tenth <- findSingleSetpoint(readings.aggregate.tenth$temperature, 
                                      readings.aggregate.tenth$kwh)
piecewise.middle <- findSingleSetpoint(readings.aggregate.middle$temperature, 
                                       readings.aggregate.middle$kwh)
piecewise.ninetieth <- findSingleSetpoint(readings.aggregate.ninetieth$temperature, 
                                          readings.aggregate.ninetieth$kwh)
##
# Plot the temperature breakpoint data.
##

# plot(readings.aggregate$temperature, 
#      readings.aggregate$kwh,
#      pch = 4, 
#      main = "Cooling Degree Hour Breakpoints", 
#      ylab = "Average Meter Reading (kWh)", 
#      xlab = "Outdoor Temperature (Celsius)", 
#      col = rgb(0.5, 0.5, 0.5, 0.5, maxColorValue=1))
# abline(piecewise.tenth[1] + piecewise.tenth[4] * piecewise.tenth[2], 
#        -piecewise.tenth[2], 
#        lwd = 2, 
#        col = rgb(0.17, 0.61, 0.22, 1, maxColorValue=1)) #lhs  
# abline(piecewise.tenth[1] - piecewise.tenth[4] * piecewise.tenth[3], 
#        piecewise.tenth[3], 
#        lwd = 2, 
#        col = rgb(0.17, 0.61, 0.22, 1, maxColorValue=1))  #rs
# abline(piecewise.middle[1] + piecewise.middle[4] * piecewise.middle[2], 
#        -piecewise.middle[2], 
#        lwd = 2, 
#        col = rgb(0.59, 0.24, 0.17, 1, maxColorValue=1)) #lhs  
# abline(piecewise.middle[1] - piecewise.middle[4] * piecewise.middle[3], 
#        piecewise.middle[3], 
#        lwd = 2, 
#        col = rgb(0.59, 0.24, 0.17, 1, maxColorValue=1))  #rs
# abline(piecewise.ninetieth[1] + piecewise.ninetieth[4] * piecewise.ninetieth[2], 
#        -piecewise.ninetieth[2], 
#        lwd = 2, 
#        col = rgb(0.19, 0.22, 0.60, 1, maxColorValue=1)) #lhs  
# abline(piecewise.ninetieth[1] - piecewise.ninetieth[4] * piecewise.ninetieth[3], 
#        piecewise.ninetieth[3], 
#        lwd = 2, 
#        col = rgb(0.19, 0.22, 0.60, 1, maxColorValue=1))  #rs

##
# Prior work determined several Cooling-Degree Hour (CDH) breakpoints for three
# quantiles. Use the CDH breakpoint from the 90th percentile, meaning that 
# this is the lowest breakpoint that people may start reacting.
##
# cdhbreak <- floor(piecewise.middle[4]) # TODO: floor or round?
cdhbreak <- 18 # Fixing CDH breakpoint at 18C since it gave better results
readings.aggregate$cdh <- ifelse(readings.aggregate$temperature > cdhbreak, 
                                  readings.aggregate$temperature - cdhbreak, 
                                  0)

##
# CDH with 1-3 hour lags.
##
for(i in 1:nrow(readings.aggregate)) {
  if (i == 1) {
    readings.aggregate[i, "cdh_1lag"] <- readings.aggregate[i, "cdh"]
    readings.aggregate[i, "cdh_2lag"] <- readings.aggregate[i, "cdh"]
    readings.aggregate[i, "cdh_3lag"] <- readings.aggregate[i, "cdh"]
  } else if (i == 2) {
    readings.aggregate[i, "cdh_1lag"] <- readings.aggregate[i, "cdh"] + readings.aggregate[i-1, "cdh"]
    readings.aggregate[i, "cdh_2lag"] <- readings.aggregate[i, "cdh"] + readings.aggregate[i-1, "cdh"]
    readings.aggregate[i, "cdh_3lag"] <- readings.aggregate[i, "cdh"] + readings.aggregate[i-1, "cdh"]
  } else if (i == 3) {
    readings.aggregate[i, "cdh_1lag"] <- readings.aggregate[i, "cdh"] + readings.aggregate[i-1, "cdh"]
    readings.aggregate[i, "cdh_2lag"] <- readings.aggregate[i, "cdh"] + readings.aggregate[i-1, "cdh"] + readings.aggregate[i-2, "cdh"]
    readings.aggregate[i, "cdh_3lag"] <- readings.aggregate[i, "cdh"] + readings.aggregate[i-1, "cdh"] + readings.aggregate[i-2, "cdh"]
  }else {
    readings.aggregate[i, "cdh_1lag"] <- readings.aggregate[i, "cdh"] + readings.aggregate[i-1, "cdh"]
    readings.aggregate[i, "cdh_2lag"] <- readings.aggregate[i, "cdh"] + readings.aggregate[i-1, "cdh"] + readings.aggregate[i-2, "cdh"]
    readings.aggregate[i, "cdh_3lag"] <- readings.aggregate[i, "cdh"] + readings.aggregate[i-1, "cdh"] + readings.aggregate[i-2, "cdh"] + readings.aggregate[i-3, "cdh"]
  }
}

##
# A linear model will be build using lm(...) for each model variant.
##

# Temperature, TOU Period, and billing method modeled as fixed effects using 
# using coefficients for main effects only.
model.lm.fe.main.tmp_prd_bill <- lm(kwh ~ temperature + tou_period + billing_active, 
                                    data = readings.aggregate)
summary(model.lm.fe.main.tmp_prd_bill)
AIC(model.lm.fe.main.tmp_prd_bill)
anova(model.lm.fe.main.tmp_prd_bill)

# Temperature, TOU Period, and billing method modeled as fixed effects using 
# using coefficients of main effects and all interactions.
model.lm.fe.inter.tmp_prd_bill <- lm(kwh ~ temperature * tou_period * billing_active, 
                                         data = readings.aggregate)
summary(model.lm.fe.inter.tmp_prd_bill)
AIC(model.lm.fe.inter.tmp_prd_bill)
anova(model.lm.fe.inter.tmp_prd_bill)

# Temperature, hour of day, and billing method modeled as fixed effects using 
# using coefficients for main effects only.
model.lm.fe.main.tmp_hr_bill <- lm(kwh ~ temperature + hrstr + billing_active, 
                                    data = readings.aggregate)
summary(model.lm.fe.main.tmp_hr_bill)
AIC(model.lm.fe.main.tmp_hr_bill)
anova(model.lm.fe.main.tmp_hr_bill)

# Temperature, hour of day, and billing method modeled as fixed effects using 
# using coefficients of main effects and all interactions.
model.lm.fe.inter.tmp_hr_bill <- lm(kwh ~ temperature * hrstr * billing_active, 
                                     data = readings.aggregate)
summary(model.lm.fe.inter.tmp_hr_bill)
AIC(model.lm.fe.inter.tmp_hr_bill)
anova(model.lm.fe.inter.tmp_hr_bill)

# Cooling-Degree Hour, TOU Period, and billing method modeled as fixed effects
# using using coefficients for main effects only.
model.lm.fe.main.cdh_prd_bill <- lm(kwh ~ cdh + tou_period + billing_active, 
                                    data = readings.aggregate)
summary(model.lm.fe.main.cdh_prd_bill)
AIC(model.lm.fe.main.cdh_prd_bill)
anova(model.lm.fe.main.cdh_prd_bill)

# Cooling-Degree Hour, TOU Period, and billing method modeled as fixed effects
# using using coefficients of main effects and all interactions.
model.lm.fe.inter.cdh_prd_bill <- lm(kwh ~ cdh * tou_period * billing_active, 
                                     data = readings.aggregate)
summary(model.lm.fe.inter.cdh_prd_bill)
AIC(model.lm.fe.inter.cdh_prd_bill)
anova(model.lm.fe.inter.cdh_prd_bill)

# Cooling-Degree Hour, hour of day, and billing method modeled as fixed effects
# using using coefficients for main effects only.
model.lm.fe.main.cdh_hr_bill <- lm(kwh ~ cdh + hrstr + billing_active, 
                                   data = readings.aggregate)
summary(model.lm.fe.main.cdh_hr_bill)
AIC(model.lm.fe.main.cdh_hr_bill)
anova(model.lm.fe.main.cdh_hr_bill)

# Cooling-Degree Hour, hour of day, and billing method modeled as fixed effects
# using using coefficients of main effects and all interactions.
model.lm.fe.inter.cdh_hr_bill <- lm(kwh ~ cdh * hrstr * billing_active, 
                                    data = readings.aggregate)
summary(model.lm.fe.inter.cdh_hr_bill)
AIC(model.lm.fe.inter.cdh_hr_bill)
anova(model.lm.fe.inter.cdh_hr_bill)




# Experimentation area to try and fit as many fixed effects as possible for the
# highest descriptive power (regardless of AIC or meaning)
supamodel.lm <- lm(kwh ~ cdh + hrstr + weekend + tou_period + billing_active + cdh:hrstr + hrstr:weekend + cdh:billing_active + weekend:billing_active + tou_period:billing_active + cdh:weekend:billing_active,
                   data = readings.aggregate)
summary(supamodel.lm)
AIC(supamodel.lm)
anova(supamodel.lm)

# Experimentation area to try and fit as many fixed effects as possible for the
# highest descriptive power (regardless of AIC or meaning)
supamodel.lm.1lag <- lm(kwh ~ cdh_1lag + hrstr + weekend + tou_period + billing_active + cdh_1lag:hrstr + hrstr:weekend + cdh_1lag:billing_active + weekend:billing_active + tou_period:billing_active + cdh_1lag:weekend:billing_active,
                   data = readings.aggregate)
summary(supamodel.lm.1lag)
AIC(supamodel.lm.1lag)
anova(supamodel.lm.1lag)

# Experimentation area to try and fit as many fixed effects as possible for the
# highest descriptive power (regardless of AIC or meaning)
supamodel.lm.2lag <- lm(kwh ~ cdh_2lag + hrstr + weekend + tou_period + billing_active + cdh_2lag:hrstr + hrstr:weekend + cdh_2lag:billing_active + weekend:billing_active + tou_period:billing_active + cdh_2lag:weekend:billing_active,
                   data = readings.aggregate)
summary(supamodel.lm.2lag)
AIC(supamodel.lm.2lag)
anova(supamodel.lm.2lag)

# Experimentation area to try and fit as many fixed effects as possible for the
# highest descriptive power (regardless of AIC or meaning)
supamodel.lm.3lag <- lm(kwh ~ cdh_3lag*hrstr*weekend*tou_period*billing_active*price,
                        data = readings.aggregate)
summary(supamodel.lm.3lag)
AIC(supamodel.lm.3lag)
anova(supamodel.lm.3lag)



##
# Strip two data frames down to TOU Period (+others) and 
# "time components" (+others). Then use bestglm to determine the best 
# combination of fixed effects for each class of model.
readings.aggregate.timetou <- subset(readings.aggregate, 
                                    select = c(cdh_1lag, 
                                               month, 
                                               tou_period, 
                                               billing_active,
                                               kwh))

readings.aggregate.timecomponents <- subset(readings.aggregate, 
                                            select = c(cdh_1lag, 
                                                       month, 
                                                       hrstr, 
                                                       weekend, 
                                                       price, 
                                                       billing_active, 
                                                       kwh))
regss.fe.timetou <- regsubsets(kwh ~ cdh_1lag*month*tou_period*billing_active, 
                            data = readings.aggregate.timetou,
                            method = "backward", 
                            nvmax = 1000)

regss.fe.timetoue.summary <- summary(regss.fe.timetou)
regss.fe.timetoue.summary$bic
regss.fe.timetoue.summary$adjr2
plotRegSubSets(regss.fe.timetoue.summary$bic, 
               regss.fe.timetoue.summary$adjr2, 
               "Automated 'TOU Period' Model Identification")

regss.fe.timecomponents <- regsubsets(kwh ~ cdh_1lag*month*hrstr*weekend*price*billing_active, 
                               data = readings.aggregate.timecomponents,
                               method = "backward", 
                               nvmax = 1000)
regss.fe.timecomponents.summary <- summary(regss.fe.timecomponents)
plotRegSubSets(regss.fe.timecomponents.summary$bic, 
               regss.fe.timecomponents.summary$adjr2, 
               "Automated 'TOU/Time Components' Model Identification")
