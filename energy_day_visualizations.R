library(plyr) # Summarize dataframe
library(nlme)

readings.aggregate.2012 <- subset(readings.aggregate, billing_active == "Yes")

# Plot demand vs. hour-of-day
plot(x = readings.aggregate.2012$hour, 
     y = readings.aggregate.2012$kwh,
     main = "Average Hourly Household Demand by Time of Day, May-October 2012", 
     xlab = "Hour of Day",
     ylab = "Average Hourly Household Demand (kWh)", 
     pch = 16, # Plot character
     col = rgb(0, 0, 0, 25, maxColorValue=100))


# Plot demand vs. temperature
plot(x = readings.aggregate.2012$temperature, 
     y = readings.aggregate.2012$kwh,
     main = "Average Hourly Household Demand by Temperature, May-October 2012", 
     xlab = "Temperature (Celsius)",
     ylab = "Average Hourly Household Demand (kWh)", 
     pch = 16, # Plot character
     col = rgb(0, 0, 0, 25, maxColorValue=100))


# Working GLM for results
df.trimmed$hour <- readings.aggregate$hour
working.glm <- glm(data = df.trimmed, 
                   formula = "kwh ~ billing_active + month + hrstr + weekend + price + lag0 + lag1 + lag2 + lag3 + lag4 + lag5 + billing_active:month + billing_active:weekend + billing_active:price + billing_active:lag0 + billing_active:lag1 + billing_active:lag2 + billing_active:lag3 + billing_active:lag4 + billing_active:lag5 + month:weekend + month:price + month:lag0 + month:lag1 + month:lag5 + hrstr:lag0 + hrstr:lag5 + weekend:lag0 + weekend:lag1 + weekend:lag2 + weekend:lag3 + weekend:lag4 + weekend:lag5 + price:lag0 + price:lag1 + price:lag5 + lag0:lag1 + lag0:lag2 + lag0:lag3 + lag0:lag4 + lag0:lag5 + lag1:lag2 + lag1:lag3 + lag1:lag4 + lag1:lag5 + lag2:lag3 + lag2:lag4 + lag2:lag5 + lag3:lag4 + lag3:lag5 + lag4:lag5", 
                   family = Gamma(link = "log"))
df.trimmed.2012 <- subset(df.trimmed, billing_active == "Yes")
df.trimmed.2012.whatif <- df.trimmed.2012
df.trimmed.2012.whatif$billing_active <- "No"

##
# From https://www.mail-archive.com/r-help@r-project.org/msg39030.html

# compute the predicted values on the scale of the link
preds <- predict.glm(working.glm,
                     newdata = df.trimmed.2012.whatif,
                     se.fit = TRUE)

## get inverse of link as a function
ilog <- family(working.glm)$linkinv


## compute standard confidence intervals
crit.t <- qt(0.975, df = length(preds)-1)
ci.u <- with(preds, fit + (crit.t * se.fit))
ci.l <- with(preds, fit - (crit.t * se.fit))


## compute predictions on scale of response by applying ilogit
preds <- ilog(preds$fit)

## transform conf int on to scale of response
ci.u <- ilog(ci.u)
ci.l <- ilog(ci.l)



# Aggregate weekdays to get a summarized day
df.trimmed.2012.weekday <- subset(df.trimmed.2012, weekend == "No")
df.trimmed.2012.weekend <- subset(df.trimmed.2012, weekend == "Yes")
df.trimmed.2012.weekday.avg <- ddply(df.trimmed.2012.weekday,
                                     ~hrstr,
                                     summarize,
                                     kwh = mean(kwh),
                                     hour = unique(hour))
df.trimmed.2012.weekday.tou <- ddply(df.trimmed.2012.weekday,
                                     ~price,
                                     summarize,
                                     kwh = mean(kwh))
df.trimmed.2012.weekend.tou <- ddply(df.trimmed.2012.weekend,
                                     ~price,
                                     summarize,
                                     kwh = mean(kwh))

# Update the whatif dataframe with predictions
df.trimmed.2012.whatif$kwh <- preds
df.trimmed.2012.whatif$conf.upper <- ci.u
df.trimmed.2012.whatif$conf.lower <- ci.l

df.trimmed.2012.whatif.weekday <- subset(df.trimmed.2012.whatif, weekend == "No")
df.trimmed.2012.whatif.weekend <- subset(df.trimmed.2012.whatif, weekend == "Yes")
df.trimmed.2012.whatif.weekday.avg <- ddply(df.trimmed.2012.whatif.weekday,
                                     ~hrstr,
                                     summarize,
                                     conf.lower = mean(conf.lower),
                                     kwh = mean(kwh),
                                     conf.upper = mean(conf.upper),
                                     hour = unique(hour))
df.trimmed.2012.whatif.weekday.tou <- ddply(df.trimmed.2012.whatif.weekday,
                                     ~price,
                                     summarize,
                                     conf.lower = mean(conf.lower), 
                                     kwh = mean(kwh),
                                     conf.upper = mean(conf.upper))
df.trimmed.2012.whatif.weekend.tou <- ddply(df.trimmed.2012.whatif.weekend,
                                     ~price,
                                     summarize,
                                     conf.lower = mean(conf.lower), 
                                     kwh = mean(kwh),
                                     conf.upper = mean(conf.upper))

off_lwr_pcnt <- 100 - 100 * (df.trimmed.2012.whatif.weekday.tou[1, "conf.lower"] / df.trimmed.2012.weekday.tou[1, "kwh"]) # off-peak lower percent change
off_avg_pcnt <- 100 - 100 * (df.trimmed.2012.whatif.weekday.tou[1, "kwh"] / df.trimmed.2012.weekday.tou[1, "kwh"]) # off-peak average percent change
off_upp_pcnt <- 100 - 100 * (df.trimmed.2012.whatif.weekday.tou[1, "conf.upper"] / df.trimmed.2012.weekday.tou[1, "kwh"]) # off-peak upper percent change

mid_lwr_pcnt <- 100 - 100 * (df.trimmed.2012.whatif.weekday.tou[2, "conf.lower"] / df.trimmed.2012.weekday.tou[2, "kwh"]) # mid-peak lower percent change
mid_avg_pcnt <- 100 - 100 * (df.trimmed.2012.whatif.weekday.tou[2, "kwh"] / df.trimmed.2012.weekday.tou[2, "kwh"]) # mid-peak average percent change
mid_upp_pcnt <- 100 - 100 * (df.trimmed.2012.whatif.weekday.tou[2, "conf.upper"] / df.trimmed.2012.weekday.tou[2, "kwh"]) # mid-peak upper percent change

on_lwr_pcnt <- 100 - 100 * (df.trimmed.2012.whatif.weekday.tou[3, "conf.lower"] / df.trimmed.2012.weekday.tou[3, "kwh"]) # on-peak lower percent change
on_avg_pcnt <- 100 - 100 * (df.trimmed.2012.whatif.weekday.tou[3, "kwh"] / df.trimmed.2012.weekday.tou[3, "kwh"]) # on-peak average percent change
on_upp_pcnt <- 100 - 100 * (df.trimmed.2012.whatif.weekday.tou[3, "conf.upper"] / df.trimmed.2012.weekday.tou[3, "kwh"]) # on-peak upper percent change

wknd_lwr_pcnt <- 100 - 100 * (df.trimmed.2012.whatif.weekend.tou[1, "conf.lower"] / df.trimmed.2012.weekend.tou[1, "kwh"]) # weekend off-peak lower percent change
wknd_avg_pcnt <- 100 - 100 * (df.trimmed.2012.whatif.weekend.tou[1, "kwh"] / df.trimmed.2012.weekend.tou[1, "kwh"]) # weekend off-peak average percent change
wknd_upp_pcnt <- 100 - 100 * (df.trimmed.2012.whatif.weekend.tou[1, "conf.upper"] / df.trimmed.2012.weekend.tou[1, "kwh"]) # weekend off-peak upper percent change

overall_lwr_pcnt <- 100 - 100 * (mean(df.trimmed.2012.whatif$conf.lower) / mean(df.trimmed.2012$kwh))
overall_avg_pcnt <- 100 - 100 * (mean(df.trimmed.2012.whatif$kwh) / mean(df.trimmed.2012$kwh))
overall_upp_pcnt <- 100 - 100 * (mean(df.trimmed.2012.whatif$conf.upper) / mean(df.trimmed.2012$kwh))

# Plot the average summer weekday
yrng <- c((min(df.trimmed.2012.whatif.weekday.avg$conf.lower) - 0.25), 
          (max(df.trimmed.2012.whatif.weekday.avg$conf.upper) + 0.25))
obscol <- "black"
predcol <- "blue"
confcol <- "grey"
plot(x = df.trimmed.2012.whatif.weekday.avg$hour, 
     y = df.trimmed.2012.whatif.weekday.avg$kwh,
     ylim = yrng,
     ylab = "Average Household Demand (kWh)",
     xaxt = "n", 
     xlab = "Hour of Day",
     main = "Average Household Electricity Demand during a Weekday, Summer 2012",
     type = "l",
     lwd = 2, 
     lty = 2, 
     col = predcol)
)
axis(side = 1, # x-axis
     at = c(0:23), 
     labels = c(0, "", 2, "", 4, "", 6, "", 8, "", 10, "", 12, "", 14, "", 16, "", 18, "", 20, "", 22, ""))
lines(x = df.trimmed.2012.whatif.weekday.avg$hour,
      y = df.trimmed.2012.whatif.weekday.avg$conf.lower,
      type = "l",
      lwd = 2, 
      lty = 3, # dashed
      col = confcol)
lines(x = df.trimmed.2012.whatif.weekday.avg$hour,
      y = df.trimmed.2012.whatif.weekday.avg$conf.upper,
      type = "l",
      lwd = 2, 
      lty = 3, # dashed
      col = confcol)
lines(x = df.trimmed.2012.weekday.avg$hour, 
      y = df.trimmed.2012.weekday.avg$kwh,
      type = "l", 
      lwd = 1, 
      col = obscol)
legend("topleft", 
       col = c(obscol, predcol, confcol), 
       lty = c(1, 2, 3),
       lwd = c(1, 2, 2), # line width
       legend = c("Observed Demand with TOU Billing",
                  "Predicted Demand without TOU Billing", 
                  "95% Confidence Interval"))