require(lmtest)         # Durban-Watson test for serially correlated errors
require(splines)        # Natural cubic splines
require(zoo)            # For rollmean used in MA(q) transform of temperature
require(effects)        # For plotting effect sizes
require(sandwich)       # Heterskedastically and autocorrelation consistent SE
require(car)            # For Anova
require(xtable)         # For generatig LaTeX tables from data.frames
require(xlsx)           # Useful for exporting results to Excel

# Source functions in other files
source("dataframe_processing.R")
source("temperature_transformation.R")
source("thesis_standardization.R")



#####
# Load aggregate smart meter readings and weighted average temperature data from 
# a CSV. Similarly, load the Environment Canada Weather descriptions for the 
# area from CSV.
#####
home <- Sys.getenv("HOME")
fpath <- file.path(home, 
                   "../Dropbox/ISS4E/R", 
                   "aggregate_readings_01Mar2011_through_17Oct2012.csv")

readings.aggregate <- InitReadingsDataFrame(fpath = fpath, 
                                            is.aggregate = TRUE)

fpath2 <- file.path(home, 
                    "../Dropbox/ISS4E/R", 
                    "weather_desc_01Mar2011_through_17Oct2012.csv")
weather <- read.csv(fpath2, 
                    na.strings = c("NULL", "NA", "NaN"), 
                    stringsAsFactors = FALSE)

# Add weather information to readings.aggregate
weather_desc.original <- unique(weather$weather_desc)
readings.aggregate$weather_reduced <- factor(ReduceWeatherCoarseTerms(weather$weather_desc))
readings.aggregate$severe_weather <- ReduceSevereWeather(weather$weather_desc)
summary(readings.aggregate$severe_weather)

readings.aggregate$working_day <- ((readings.aggregate$holiday == TRUE | 
                                      readings.aggregate$dayname == "Sat" | 
                                      readings.aggregate$dayname == "Sun") == FALSE)

# Add column to identify utility rate seasons
readings.aggregate$rateseason <- rep(NA, nrow(readings.aggregate))
readings.aggregate$rateseason[readings.aggregate$month %in% c("m5", "m6", "m7", 
                                                              "m8", "m9", 
                                                              "m10")] <- "summer"
readings.aggregate$rateseason[readings.aggregate$month %in% c("m11", "m12", 
                                                              "m1", "m2", "m3", 
                                                              "m4")] <- "winter"
readings.aggregate$rateseason <- factor(readings.aggregate$rateseason, 
                                        c("winter", "summer"))

# Add a period of day to each observation for quantifying effects for each.
readings.aggregate$dayperiod <- rep(NA, nrow(readings.aggregate))
readings.aggregate$dayperiod[readings.aggregate$hrstr %in% c("h0", "h1", "h2", 
                                                             "h3", "h4", "h5", 
                                                             "h6", "h19", "h20", 
                                                             "h21", "h22", 
                                                             "h23")] <- "overnight"
readings.aggregate$dayperiod[readings.aggregate$hrstr %in% c("h7", "h8", "h9", 
                                                             "h10", "h17", 
                                                             "h18")] <- "transition"
readings.aggregate$dayperiod[readings.aggregate$hrstr %in% c("h11", "h12", 
                                                             "h13", "h14", 
                                                             "h15", "h16")] <- "afternoon"
readings.aggregate$dayperiod[readings.aggregate$working_day == FALSE] <- "nonwd"

# For clarity, reorder explanatory variables which are factors
readings.aggregate <- OrderFactors(readings.aggregate)
rm(fpath, fpath2, home)



#####
# Add "feels like" temperature to data.frame
#####
readings.aggregate$feelslike <- FeelsLike(readings.aggregate)



#####
# Transform temperature to a moving average of the past 4 hours of temperature.
# Optimal MA window is 6 hours, discovered empirically through an iteration of 
# 1-24 hour windows.
#####
ma.order <- 6
ma.trimsize <- (ma.order - 1)

# rollingmean(...) averages elements [i, ..., i+ma.order] and then trims 
# ma.order-1 from the end of the vector
rolling.avg <- rollmean(x = readings.aggregate$temperature, k = ma.order)

# To align the rollingmean(...) vector to be elements [i-1-ma.order, ..., i], 
# the first ma.order-1 elements should be trimmed from the head of the 
# observation data.frame
readings.trimmed <- tail(readings.aggregate, -ma.trimsize)

# Add the rolling mean as a column to the dataframe
readings.trimmed$rollingavg <- rolling.avg



#####
# Because Nov, Dec, Jan, and Feb (the most extreme parts of winter) occur only 
# once in the data set, it is skewing the effects of TOU pricing.
#
# Because chapter 5 focuses on the effects of TOU pricing, this unbalanced data 
# design is undesirable. Because I am not using autocorrelated errors, I can 
# trim the middle months out of the data set, model seasonality with "month" and 
# periodicity with "hrstr".
#
# Durbin-Watson test for serial correlation of residuals is still a bit high, 
# so long as I use heteroskedastically and autocorrelation consistent standard 
# errors, my reported errors will still have value.
#####
# Empirical, optimal natural spline knots were found for the temperature MA 
# vector in the presence of the "month" term and different subset of data.
rollingavg.ns <- ns(readings.trimmed$rollingavg, knots = c(6, 26, 28))
# Improve column names beyond V1, ..., V4
colnames(rollingavg.ns) <- c("ns_sec1", "ns_sec2", "ns_sec3", "ns_sec4")
readings.trimmed <- cbind(readings.trimmed, rollingavg.ns)

tou.lm <- lm(kwh ~ ns_sec1 + ns_sec2 + ns_sec3 + ns_sec4 + month + 
               (hrstr + working_day + rateseason + tou_active)^3,
             data = readings.trimmed)
summary(tou.lm)$adj.r.squared


# Heteroskedasticity and autocorrelation consistent covariance matrix
hac.vcov <- vcovHAC(tou.lm)

# 3-way interaction with the least statistical signficance is
# working_day:rateseason:tou_active; df=1, F=0.2144, pval=0.6433377.
Anova(mod = tou.lm, vcov = hac.vcov)

# t-test of coefficients
# working_dayTRUE:rateseasonsummer:tou_activeTRUE coef=0.011776
coeftest(tou.lm, vcov = hac.vcov)

# The effect size is fairly small, reasonably high pval, but only one degree of 
# freedom. Intuitively, does one expect the entire work day to be shifted 
# up/down uniquely based on its rateseason and pre-/post-TOU rates?
# 
# I don't expect this to be significantly different based on those criteria, 
# especially in the presence of other three-way interactions. I believe
# that removal of the term is justified.
tou.lm <- update(tou.lm, "~.-working_day:rateseason:tou_active")
summary(tou.lm)

# Update the covariance matrix.
hac.vcov <- vcovHAC(tou.lm)

# 3-way interaction with the least statistical significance is 
# hrstr:working_day:tou_active; df=23, F=1.0260, pval=0.4264493
Anova(mod = tou.lm, vcov = hac.vcov)

# t-test of coefficients
# All coefficients for hrstrh__:rateseasonsummer:tou_activeTRUE are all negative
# and abs(coef) < 0.037kWh, which is reasonably small. This suggests an 
# across-the-board reduction by summer TOU prices, which could be handled by 
# the two-way interaction of rateseasonsummer:tou_activeTRUE.
#
# I feel that the removal of this 3-way interaction is justified both 
# intuitively and statistically.
coeftest(tou.lm, vcov = hac.vcov)
tou.lm <- update(tou.lm, "~.-hrstr:working_day:tou_active")
summary(tou.lm)

# Update the covariance matrix.
hac.vcov <- vcovHAC(tou.lm)

# 2-way interaction with the least statistical significance is 
# working_day:tou_active; df=1, F=1.3445, pval=0.2462646
Anova(mod = tou.lm, vcov = hac.vcov)

# t-test of coefficient
# working_dayTRUE:tou_activeTRUE coef=-0.01658320
# I would err on the side of _not_ giving an effect to TOU pricing (ie. I'd 
# rather make a type I error).
coeftest(tou.lm, vcov = hac.vcov)
tou.lm <- update(tou.lm, "~.-working_day:tou_active")
summary(tou.lm)

# Update the covariance matrix.
hac.vcov <- vcovHAC(tou.lm)

# All remaining terms are either required by higher-order interactions (ie. 
# principle of marginality) or are statistically significant themselves.
Anova(mod = tou.lm, vcov = hac.vcov)
coeftest(tou.lm, vcov = hac.vcov)



#####
# Evaluate the effects of TOU pricing by price level periods. Fixing rateseason,
# working_day, and a temperature typical to those fixed values allows for 
# conditioning the model on those values and allowing tou_active and hrstr
# to vary.
#
# Building a lm(...) object by hand again so that I don't need to run the 
# hac.vcov each time I want to run this later section of code.
#####
lm.tou <- lm(kwh ~ ns_sec1 + ns_sec2 + ns_sec3 + ns_sec4 + 
               month + hrstr + working_day + rateseason + tou_active + 
               hrstr:working_day + hrstr:rateseason + hrstr:tou_active + 
               working_day:rateseason + rateseason:tou_active +
               hrstr:working_day:rateseason + hrstr:rateseason:tou_active, 
             data = readings.trimmed)

tou.start.idx <- min(which(readings.trimmed$tou_active == TRUE, 
                           arr.ind = TRUE))
counterfactual.trimmed <- readings.trimmed
# Apply inverse of tou_active to counterfactual data.frame
counterfactual.trimmed$tou_active = factor(readings.trimmed$tou_active == FALSE, 
                                           levels = c("FALSE", "TRUE"))
counterfactual.predict <- predict(object = lm.tou, 
                                  newdata = counterfactual.trimmed,
                                  se.fit = TRUE,
                                  interval="confidence")

# Find the on- to off-peak ratio for TOU case study
onoff.ratio.results <- OnToOffPeakRatio(observed.datafr = readings.trimmed,
                                        counterfactual.datafr = counterfactual.trimmed,
                                        prediction = counterfactual.predict,
                                        tou.idx = tou.start.idx)
write.xlsx(onoff.ratio.results, "./results/nosubset_toustudy_onoff_ratio.xlsx")

# Find the peak-to-average ratio for TOU case study
peakavg.ratio.results <- PeakToAverageRatio(observed.datafr = readings.trimmed,
                                            counterfactual.datafr = counterfactual.trimmed,
                                            prediction = counterfactual.predict,
                                            tou.idx = tou.start.idx)
View(peakavg.ratio.results)
write.xlsx(peakavg.ratio.results, "./results/nosubset_toustudy_peakavg_ratio.xlsx")

# Create tabular results of TOU effects
# TODO(r24mille): Is there some way to use heterscedastic and autocorrelation 
#                 consistent standard errors with the predict(...) function?
#                 vcov.=vcovHAC or something similar?
tou.effects <- EffectsTable(observed.datafr = readings.trimmed,
                            counterfactual.datafr = counterfactual.trimmed,
                            prediction = counterfactual.predict,
                            tou.idx = tou.start.idx)
print(xtable(tou.effects))
write.xlsx(tou.effects, "./results/nosubset_toustudy.xlsx")



# Create plots!
winter.wd.hourly <- GetConfident(observed.datafr = readings.trimmed,
                                 counterfactual.datafr = counterfactual.trimmed,
                                 prediction = counterfactual.predict,
                                 tou.idx = tou.start.idx,
                                 rate.season = "winter",
                                 working.day = TRUE)
winter.wd.plot <- PlotTouEffectResult(effect.datafr = winter.wd.hourly, 
                                      tou.datafr = WinterWorkingDayTouSequenceDataframe(),
                                      ggp.title = "Effects of TOU Pricing (Winter, Working Day)")
print(winter.wd.plot)
dev.print(file="../../Figures/TouEffectsWinterWorkingDay.png", device=png, height = 600, width = 900)

winter.nonwd.hourly <- GetConfident(observed.datafr = readings.trimmed,
                                    counterfactual.datafr = counterfactual.trimmed,
                                    prediction = counterfactual.predict,
                                    tou.idx = tou.start.idx,
                                    rate.season = "winter",
                                    working.day = FALSE)
winter.nonwd.plot <- PlotTouEffectResult(effect.datafr = winter.nonwd.hourly, 
                                         tou.datafr = NonWorkingDayTouSequenceDataframe(),
                                         ggp.title = "Effects of TOU Pricing (Winter, Non-Working Day)")
print(winter.nonwd.plot)
dev.print(file="../../Figures/TouEffectsWinterNonWorkingDay.png", device=png, height = 600, width = 900)

summer.wd.hourly <- GetConfident(observed.datafr = readings.trimmed,
                                 counterfactual.datafr = counterfactual.trimmed,
                                 prediction = counterfactual.predict,
                                 tou.idx = tou.start.idx,
                                 rate.season = "summer",
                                 working.day = TRUE)
summer.wd.plot <- PlotTouEffectResult(effect.datafr = summer.wd.hourly, 
                                      tou.datafr = SummerWorkingDayTouSequenceDataframe(),
                                      ggp.title = "Effects of TOU Pricing (Summer, Working Day)")
print(summer.wd.plot)
dev.print(file="../../Figures/TouEffectsSummerWorkingDay.png", device=png, height = 600, width = 900)

summer.nonwd.hourly <- GetConfident(observed.datafr = readings.trimmed,
                                    counterfactual.datafr = counterfactual.trimmed,
                                    prediction = counterfactual.predict,
                                    tou.idx = tou.start.idx,
                                    rate.season = "summer",
                                    working.day = FALSE)
summer.nonwd.plot <- PlotTouEffectResult(effect.datafr = summer.nonwd.hourly, 
                                         tou.datafr = NonWorkingDayTouSequenceDataframe(),
                                         ggp.title = "Effects of TOU Pricing (Summer, Non-Working Day)")
print(summer.nonwd.plot)
dev.print(file="../../Figures/TouEffectsSummerNonWorkingDay.png", device=png, height = 600, width = 900)