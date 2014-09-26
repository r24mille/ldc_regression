library(lme4)
library(plyr) # Summarize dataframe
library(nlme) # groupedData
library(lattice) # Plotting
library(stargazer) # LaTeX tables
library(bestglm) # Compare models using leaps
library(glmulti)
library(segmented)
library(BMA)


# Source the function in another file
source('piecewise_cdh_setpoint.R')
source('reg_subset_visualization.R')
source('createCdhLag.R')

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
# Use 'segmented' package rather than the Stack Overflow post to determine 
# the optimal cooling degree hour breakpoint (though they give the same 
# result).
model.readings.lm.presegment <- lm(kwh ~ temperature*tou_period*billing_active,
                                   data = readings.aggregate)
seg <- segmented(obj = model.readings.lm.presegment, 
                 seg.Z = ~temperature,
                 psi = list(temperature = c(15)))
cdhbreak <- floor(seg$psi[1,2]) # TODO floor vs. round?
readings.aggregate$cdh <- ifelse(readings.aggregate$temperature > cdhbreak, 
                                 readings.aggregate$temperature - cdhbreak, 
                                 0)

# ## 
# # Section finds the Cooling-Degree Hour (CDH) baseline for 10th, middle, and 
# # 90th percentile observations for each temperature "bin" similar to the 
# # three-line model.
# ##
# 
# # Determine kwh quantils for each temperature "bin"
# quants <- ddply(readings.aggregate,
#                    .(temperature), 
#                    function(x) quantile(x$kwh, c(.1, .5, .9))
# )
# # Hacky due to temperatures matching quant indeces (by chance).
# # TODO: Fix with a proper subset(...) I'm tired of fighting it for today.
# readings.aggregate.tenth <- subset(readings.aggregate, 
#                                    kwh <= quants[temperature, 2])
# readings.aggregate.middle <- subset(readings.aggregate, 
#                                     kwh > quants[temperature, 2] & kwh < quants[temperature, 4])
# readings.aggregate.ninetieth <- subset(readings.aggregate, 
#                                       kwh >= quants[temperature, 4])
# 
# # Determine the A/C setpoint threshold
# piecewise.tenth <- findSingleSetpoint(readings.aggregate.tenth$temperature, 
#                                       readings.aggregate.tenth$kwh)
# piecewise.middle <- findSingleSetpoint(readings.aggregate.middle$temperature, 
#                                        readings.aggregate.middle$kwh)
# piecewise.ninetieth <- findSingleSetpoint(readings.aggregate.ninetieth$temperature, 
#                                           readings.aggregate.ninetieth$kwh)
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
# cdhbreak <- 18 # Fixing CDH breakpoint at 18C since it gave better results
# readings.aggregate$cdh <- ifelse(readings.aggregate$temperature > cdhbreak, 
#                                   readings.aggregate$temperature - cdhbreak, 
#                                   0)

##
# CDH with 1-3 hour lags.
##
nlags = 1
cdhlag <- createCdhLag(nlags,
                       readings.aggregate)
readings.aggregate$cdhlag <- cdhlag

##
# Strip two data frames down to TOU Period (+others) and 
# "time components" (+others). Then use bestglm to determine the best 
# combination of fixed effects for each class of model.
readings.aggregate.timetou <- subset(readings.aggregate, 
                                    select = c(cdhlag, 
                                               month, 
                                               tou_period, 
                                               billing_active,
                                               kwh))

readings.aggregate.timecomponents <- subset(readings.aggregate, 
                                            select = c(cdhlag, 
                                                       month, 
                                                       hrstr, 
                                                       weekend, 
                                                       price, 
                                                       billing_active, 
                                                       kwh))

regss.fe.timetou <- regsubsets(kwh ~ cdhlag*month*tou_period*billing_active, 
                            data = readings.aggregate.timetou,
                            method = "backward", 
                            nvmax = 500)

regss.fe.timetou.summary <- summary(regss.fe.timetou, 
                                    matrix = TRUE,
                                    matrix.logical = TRUE)
regss.fe.timetou.optidx <- which.min(regss.fe.timetou.summary$bic)
plotRegSubSets(regss.fe.timetou.summary$bic, 
               regss.fe.timetou.summary$adjr2, 
               "Automated 'TOU Period' Model Identification")

# Use the TRUE/FALSE summary from .which to get the optimal factor values as 
# a model string.
optwhich <- as.data.frame(regss.fe.timetou.summary$which[regss.fe.timetou.optidx,])
colnames(optwhich) <- c("included")
optwhich <- subset(optwhich,
                   included == TRUE)
optwhich <- tail(optwhich, regss.fe.timetou.optidx) # Trim (Intercept)
# model.fe.timetou.festring <- paste(rownames(optwhich), 
#                                    collapse = " + ")
# model.fe.timetou.optimal.formula <- as.formula(paste("kwh ~ ", 
#                                                      model.fe.timetou.festring))
model.fe.timetou.a <- lm(kwh ~ cdhlag + month + tou_period + cdhlag:month + cdhlag:tou_period + month:tou_period + month:billing_active + cdhlag:month:tou_period + cdhlag:month:billing_active + cdhlag:tou_period:billing_active,
                       data = readings.aggregate)
summary(model.fe.timetou.a)
AIC(model.fe.timetou.a)

model.fe.timetou.b <- lm(kwh ~ cdhlag:month:tou_period + cdhlag:month:billing_active + cdhlag:tou_period:billing_active,
                         data = readings.aggregate)
summary(model.fe.timetou.b)
AIC(model.fe.timetou.b)
plot(model.fe.timetou.b)






##
# Stepwise deletion
model.timetou.maximal <- lm(kwh ~ cdhlag * month * tou_period * billing_active, 
                            data = readings.aggregate)
summary(model.timetou.maximal)
anova(model.timetou.maximal)

# Step 1
model.timetou.step1 <- update(model.timetou.maximal, 
                              ~.- cdhlag:billing_active)
summary(model.timetou.step1)
anova(model.timetou.step1)
anova(model.timetou.maximal, model.timetou.step1)

# Step 2
model.timetou.step2 <- update(model.timetou.step1, 
                              ~.- cdhlag:tou_period:billing_active)
summary(model.timetou.step2)
anova(model.timetou.step2)
anova(model.timetou.step1, model.timetou.step2)

# Step 3
model.timetou.step3 <- update(model.timetou.step2, 
                              ~.- cdhlag:month:billing_active)
summary(model.timetou.step3)
anova(model.timetou.step3)
anova(model.timetou.step2, model.timetou.step3)

# Step 4
model.timetou.step4 <- update(model.timetou.step3, 
                              ~.- cdhlag:month:tou_period)
summary(model.timetou.step4)
anova(model.timetou.step4)
anova(model.timetou.step3, model.timetou.step4)

# Step 5
model.timetou.step5 <- update(model.timetou.step4, 
                              ~.- month:tou_period)
summary(model.timetou.step5)
anova(model.timetou.step5)
anova(model.timetou.step4, model.timetou.step5)

# Step 6
model.timetou.step6 <- update(model.timetou.step5, 
                              ~.- cdhlag:tou_period)
summary(model.timetou.step6)
anova(model.timetou.step6)
anova(model.timetou.step5, model.timetou.step6)

# Step 7
model.timetou.step7 <- update(model.timetou.step6, 
                              ~.- cdhlag:month)
summary(model.timetou.step7)
anova(model.timetou.step7)
anova(model.timetou.step6, model.timetou.step7)

# Step 8
model.timetou.step8 <- update(model.timetou.step7, 
                              ~.- tou_period:billing_active)
summary(model.timetou.step8)
anova(model.timetou.step8)
anova(model.timetou.step7, model.timetou.step8)

# Step 9
model.timetou.step9 <- update(model.timetou.step8, 
                              ~.- month:billing_active)
summary(model.timetou.step9)
anova(model.timetou.step9)
anova(model.timetou.step8, model.timetou.step9)

# Step 10
model.timetou.step10 <- update(model.timetou.step9, 
                              ~.- billing_active)
summary(model.timetou.step10)
anova(model.timetou.step10)
anova(model.timetou.step9, model.timetou.step10)

# Step 11
model.timetou.step11 <- update(model.timetou.step10, 
                               ~.- tou_period)
summary(model.timetou.step11)
anova(model.timetou.step11)
anova(model.timetou.step10, model.timetou.step11)

# Step 12
model.timetou.step12 <- update(model.timetou.step11, 
                               ~.- month)
summary(model.timetou.step12)
anova(model.timetou.step12)
anova(model.timetou.step11, model.timetou.step12)

# Step 13
model.timetou.step13 <- update(model.timetou.step12, 
                               ~.- cdhlag)
summary(model.timetou.step13)
anova(model.timetou.step13)
anova(model.timetou.step12, model.timetou.step13)


####
# HOLD UP...
####
glm.timetou.maximal <- glm(kwh ~ cdhlag * month * tou_period * billing_active, 
                           data = readings.aggregate)
summary(glm.timetou.maximal)
anova(glm.timetou.maximal)

# Step 1
glm.timetou.step1 <- update(glm.timetou.maximal,
                            ~.- cdhlag:month:tou_period:billing_active)
summary(glm.timetou.step1)
anova(glm.timetou.step1)
anova(glm.timetou.maximal, glm.timetou.step1, test="Chi")
qchisq(.95, 23) 
# 3.3 increase in deviance is significantly below 35.2 allowed by 
# qchisq(...). This simplification is justified.

# Step 2
glm.timetou.step2 <- update(glm.timetou.step1,
                            ~.- month:tou_period:billing_active)
summary(glm.timetou.step2)
anova(glm.timetou.step2)
anova(glm.timetou.step1, glm.timetou.step2, test="Chi")
qchisq(.95, 25)
# 3.2 increase in deviance is significantly below 37.7 allowed by 
# qchisq(...). This simplification is justified.

#Step 3
glm.timetou.step3 <- update(glm.timetou.step2,
                            ~.- cdhlag:tou_period:billing_active)
summary(glm.timetou.step3)
anova(glm.timetou.step3)
anova(glm.timetou.step2, glm.timetou.step3)
qchisq(.95, 5)
# 0.6 increase in deviance is significantly below 11.1 allowed by 
# qchisq(...). This simplification is justified.

#Step 4
glm.timetou.step4 <- update(glm.timetou.step3,
                            ~.- cdhlag:month:billing_active)
summary(glm.timetou.step4)
anova(glm.timetou.step4)
anova(glm.timetou.step3, glm.timetou.step4)
qchisq(.95, 5)
# 1.9 increase in deviance is significantly below 11.1 allowed by 
# qchisq(...). This simplification is justified.

#Step 5
glm.timetou.step5 <- update(glm.timetou.step4,
                            ~.- cdhlag:month:tou_period)
summary(glm.timetou.step5)
anova(glm.timetou.step5)
anova(glm.timetou.step4, glm.timetou.step5)
qchisq(.95, 24)
# 10.4 increase in deviance is significantly below 36.4 allowed by 
# qchisq(...). This simplification is justified.

#Step 6
glm.timetou.step6 <- update(glm.timetou.step5,
                            ~.- tou_period:billing_active)
summary(glm.timetou.step6)
anova(glm.timetou.step6)
anova(glm.timetou.step5, month:billing_active)
qchisq(.95, 5)
# 1.4 increase in deviance is significantly below 11.1 allowed by 
# qchisq(...). This simplification is justified.

#Step 7
glm.timetou.step7 <- update(glm.timetou.step6,
                            ~.- month:billing_active)
summary(glm.timetou.step7)
anova(glm.timetou.step7)
anova(glm.timetou.step6, glm.timetou.step7)
qchisq(.95, 5)
# 1.4 increase in deviance is significantly below 11.1 allowed by 
# qchisq(...). This simplification is justified.

#Step 8
glm.timetou.step8 <- update(glm.timetou.step7,
                            ~.- cdhlag:billing_active)
summary(glm.timetou.step8)
anova(glm.timetou.step8)
anova(glm.timetou.step7, glm.timetou.step8)
qchisq(.95, 1)
# 0.1 increase in deviance is significantly below 3.8 allowed by 
# qchisq(...). This simplification is justified.

#Step 9
glm.timetou.step9 <- update(glm.timetou.step8,
                            ~.- cdhlag:tou_period)
summary(glm.timetou.step9)
anova(glm.timetou.step9)
anova(glm.timetou.step8, glm.timetou.step9)
# 17.5 increase in deviance is below 37.7 allowed by qchisq(...). This 
# simplification is justified.

plot(c(BIC(glm.timetou.maximal), 
       BIC(glm.timetou.step1), 
       BIC(glm.timetou.step2), 
       BIC(glm.timetou.step3), 
       BIC(glm.timetou.step4), 
       BIC(glm.timetou.step5), 
       BIC(glm.timetou.step6), 
       BIC(glm.timetou.step7),
       BIC(glm.timetou.step8),
       BIC(glm.timetou.step9)),
     ylab = "BIC")

glm.timetou.opt.FT <- bic.glm(f = kwh ~ cdhlag*tou_period*billing_active + month + cdhlag:month + tou_period:month + billing_active:month, 
                             data = readings.aggregate.timetou, 
                             glm.family = "gaussian", # TODO: I think Gaussian is good
                             factor.type = TRUE)
bma.timetou.summary <- summary(glm.timetou.opt.FT)
imageplot.bma(glm.timetou.opt.FT)


# Plotzzzz
par(mfrow=c(2,2)) 
plot(kwh~cdhlag, col="red", data=readings.aggregate) 
plot(kwh~month, col="blue", data=readings.aggregate) 
plot(kwh~tou_period, col="green", data=readings.aggregate)
plot(kwh~billing_active, col="yellow", data=readings.aggregate)

regss.fe.timecomponents <- regsubsets(kwh ~ cdhlag*month*hrstr*weekend*price*billing_active, 
                               data = readings.aggregate.timecomponents,
                               method = "backward", 
                               nvmax = 1000)
regss.fe.timecomponents.summary <- summary(regss.fe.timecomponents)
plotRegSubSets(regss.fe.timecomponents.summary$bic, 
               regss.fe.timecomponents.summary$adjr2, 
               "Automated 'TOU/Time Components' Model Identification")
