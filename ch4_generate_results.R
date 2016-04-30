require(segmented)  # For finding temperature breakpoint
require(splines)    # Natural cubic splines
require(dlnm)       # Package to model serially correlated error structure
require(zoo)        # Fore rollmean used in MA(q) transform of temperature
require(scales)         # For scaling axes (eg. datetime)
require(gridExtra)      # Placing multiple ggplot2 plots into a single viewport or image


# Source functions in other files
source("dataframe_processing.R")
source("temperature_transformation.R")



#####
# Load aggregate smart meter readings and weighted average temperature data from
# a CSV. Similarly, load the Environment Canada Weather descriptions for the
# area from CSV.
#####
home <- Sys.getenv("HOME")
fpath <- file.path(home,
                   "./Dropbox/ISS4E/R",
                   "aggregate_readings_01Mar2011_through_17Oct2012.csv")

readings.aggregate <- InitReadingsDataFrame(fpath = fpath,
                                            is.aggregate = TRUE)

fpath2 <- file.path(home,
                    "./Dropbox/ISS4E/R",
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

# For clarity, reorder explanatory variables which are factors
readings.aggregate <- OrderFactors(readings.aggregate)
rm(fpath, fpath2, home)



#####
# Create a vector where extreme temperatures have been replaced with wind chill
# and heat index accordingly.
#####
feelslike <- FeelsLike(readings.aggregate)



#####
# Moving average transformation
#####
ma.order <- 6
ma.trimsize <- (ma.order - 1)

# To align the rollingmean(...) vector to be elements [i-1-ma.order, ..., i],
# the first ma.order-1 elements should be trimmed from the head of the
# observation data.frame
readings.ma.pretrim <- tail(readings.aggregate, -ma.trimsize)

# rollingmean(...) averages elements [i, ..., i+ma.order] and then trims
# ma.order-1 from the end of the vector
ma.drybulb <- rollmean(x = readings.aggregate$temperature, k = ma.order)
ma.feelslike <- rollmean(x = feelslike, k = ma.order)

readings.ma.pretrim <- cbind(readings.ma.pretrim, ma.drybulb, ma.feelslike)
colnames(readings.ma.pretrim)[c(32, 33)] <- c("ma_db", "ma_fl")




#####
# Find the heating/cooling breakpoint for drybulb temperature and "feels like"
# temperature.
#####
lm.drybulb <- lm(kwh ~ temperature + hrstr + working_day + tou_active +
                   rateseason + hrstr:working_day,
                 data = readings.aggregate)
seg.drybulb <- segmented(obj = lm.drybulb,
                         seg.Z = ~ temperature,
                         psi = 17)
drybulb.break <- seg.drybulb$psi[1,2]

lm.logtrnfm.drybulb <- lm(log(kwh) ~ temperature + hrstr + working_day +
                            tou_active + rateseason + hrstr:working_day,
                          data = readings.aggregate)
seg.logtrnfm.drybulb <- segmented(obj = lm.logtrnfm.drybulb,
                                  seg.Z = ~ temperature,
                                  psi = 17)
drybulb.logtrnfm.break <- seg.logtrnfm.drybulb$psi[1,2]

lm.feelslike <- lm(kwh ~ feelslike + hrstr + working_day + tou_active +
                     rateseason + hrstr:working_day,
                   data = readings.aggregate)
seg.feelslike <- segmented(obj = lm.feelslike,
                           seg.Z = ~ feelslike,
                           psi = 17)
feelslike.break <- seg.feelslike$psi[1,2]

lm.logtrnfm.feelslike <- lm(log(kwh) ~ feelslike + hrstr + working_day +
                              tou_active + rateseason + hrstr:working_day,
                            data = readings.aggregate)
seg.logtrnfm.feelslike <- segmented(obj = lm.logtrnfm.feelslike,
                                    seg.Z = ~ feelslike,
                                    psi = 17)
feelslike.logtrnfm.break <- seg.logtrnfm.feelslike$psi[1,2]

lm.ma.drybulb <- lm(kwh ~ ma_db + hrstr + working_day + tou_active +
                      rateseason + hrstr:working_day,
                    data = readings.ma.pretrim)
seg.ma.drybulb <- segmented(obj = lm.ma.drybulb,
                            seg.Z = ~ ma_db,
                            psi = 18)
ma.drybulb.break <- seg.ma.drybulb$psi[1,2]

lm.ma.logtrnfm.drybulb <- lm(log(kwh) ~ ma_db + hrstr + working_day +
                               tou_active + rateseason + hrstr:working_day,
                             data = readings.ma.pretrim)
seg.ma.logtrnfm.drybulb <- segmented(obj = lm.ma.logtrnfm.drybulb,
                                     seg.Z = ~ ma_db,
                                     psi = 18)
ma.drybulb.logtrnfm.break <- seg.ma.logtrnfm.drybulb$psi[1,2]

lm.ma.feelslike <- lm(kwh ~ ma_fl + hrstr + working_day + tou_active +
                        rateseason + hrstr:working_day,
                      data = readings.ma.pretrim)
seg.ma.feelslike <- segmented(obj = lm.ma.feelslike,
                              seg.Z = ~ ma_fl,
                              psi = 18)
ma.feelslike.break <- seg.ma.feelslike$psi[1,2]

lm.ma.logtrnfm.feelslike <- lm(log(kwh) ~ ma_fl + hrstr + working_day +
                                 tou_active + rateseason + hrstr:working_day,
                               data = readings.ma.pretrim)
seg.ma.logtrnfm.feelslike <- segmented(obj = lm.ma.logtrnfm.feelslike,
                                       seg.Z = ~ ma_fl,
                                       psi = 18)
ma.feelslike.logtrnfm.break <- seg.ma.logtrnfm.feelslike$psi[1,2]



#####
# Create the intercept-only model (ie. null model)
#####
EvaluateTemperatureTransform(datafr = readings.aggregate,
                             temptransform = rep(NA, nrow(readings.aggregate)),
                             methodname = "NullModel",
                             trimsize = 0,
                             logtransform = FALSE,
                             nullmodel = TRUE,
                             includetemp = FALSE)

EvaluateTemperatureTransform(datafr = readings.aggregate,
                             temptransform = rep(NA, nrow(readings.aggregate)),
                             methodname = "NullModel",
                             trimsize = 0,
                             logtransform = TRUE,
                             nullmodel = TRUE,
                             includetemp = FALSE)



#####
# Create a model with no temperature effects
#####
EvaluateTemperatureTransform(datafr = readings.aggregate,
                             temptransform = rep(NA, nrow(readings.aggregate)),
                             methodname = "NoTemp",
                             trimsize = 0,
                             logtransform = FALSE,
                             includetemp = FALSE)

EvaluateTemperatureTransform(datafr = readings.aggregate,
                             temptransform = rep(NA, nrow(readings.aggregate)),
                             methodname = "NoTemp",
                             trimsize = 0,
                             logtransform = TRUE,
                             includetemp = FALSE)



#####
# Create a model with only linear, drybulb temperature effects
#####
temp.mat <- as.matrix(readings.aggregate$temperature)
colnames(temp.mat) <- c("temperature")
EvaluateTemperatureTransform(datafr = readings.aggregate,
                             temptransform = temp.mat,
                             methodname = "LinearDrybulb",
                             trimsize = 0,
                             logtransform = FALSE)

EvaluateTemperatureTransform(datafr = readings.aggregate,
                             temptransform = temp.mat,
                             methodname = "LinearDrybulb",
                             trimsize = 0,
                             logtransform = TRUE)



#####
# Create switching regression model
#####
sr.drybulb.basismatrix <- SwitchingRegressionBasisMatrix(explvar = readings.aggregate$temperature,
                                                         breakpoint = drybulb.break)
EvaluateTemperatureTransform(datafr = readings.aggregate,
                             temptransform = sr.drybulb.basismatrix,
                             methodname = "SR_Drybulb",
                             trimsize = 0,
                             logtransform = FALSE)

sr.logtrnfm.drybulb.basismatrix <- SwitchingRegressionBasisMatrix(explvar = readings.aggregate$temperature,
                                                                  breakpoint = drybulb.logtrnfm.break)
EvaluateTemperatureTransform(datafr = readings.aggregate,
                             temptransform = sr.logtrnfm.drybulb.basismatrix,
                             methodname = "SR_Drybulb",
                             trimsize = 0,
                             logtransform = TRUE)



#####
# Create switching regression model with wind chill + heat index replacements
# for temperature extremes
#####
sr.feelslike.basismatrix <- SwitchingRegressionBasisMatrix(explvar = feelslike,
                                                           breakpoint = feelslike.break)
EvaluateTemperatureTransform(datafr = readings.aggregate,
                             temptransform = sr.feelslike.basismatrix,
                             methodname = "SR_FeelsLike",
                             trimsize = 0,
                             logtransform = FALSE)

sr.logtrnfm.feelslike.basismatrix <- SwitchingRegressionBasisMatrix(explvar = feelslike,
                                                                    breakpoint = feelslike.logtrnfm.break)
EvaluateTemperatureTransform(datafr = readings.aggregate,
                             temptransform = sr.logtrnfm.feelslike.basismatrix,
                             methodname = "SR_FeelsLike",
                             trimsize = 0,
                             logtransform = TRUE)



#####
# Create multiple regression with natural cubic splines transformation of
# drybulb temperature. Knots chosen empirically by iteratively stepping through
# combinations of integers.
#####
drybulb.ns <- ns(readings.aggregate$temperature,
                 knots = c(3, 23, 30))
# Column names must be better than V1, ..., V4 so that my evaluation function doesn't break.
colnames(drybulb.ns) <- c("ns_sec1", "ns_sec2", "ns_sec3", "ns_sec4")
EvaluateTemperatureTransform(datafr = readings.aggregate,
                             temptransform = drybulb.ns,
                             methodname = "NS_Drybulb",
                             trimsize = 0,
                             logtransform = FALSE)

drybulb.logtrnfm.ns <- ns(readings.aggregate$temperature,
                          knots = c(3, 16, 31))
# Column names must be better than V1, ..., V4 so that my evaluation function doesn't break.
colnames(drybulb.logtrnfm.ns) <- c("ns_sec1", "ns_sec2", "ns_sec3", "ns_sec4")
EvaluateTemperatureTransform(datafr = readings.aggregate,
                             temptransform = drybulb.logtrnfm.ns,
                             methodname = "NS_Drybulb",
                             trimsize = 0,
                             logtransform = TRUE)



#####
# Create multiple regression with natural cubic splines transformation of
# "feels like" temperature.
#####
feelslike.ns <- ns(feelslike,
                   knots = c(4, 22, 27))
# Column names must be better than V1, ..., V4 so that my evaluation function doesn't break.
colnames(feelslike.ns) <- c("ns_sec1", "ns_sec2", "ns_sec3", "ns_sec4")
EvaluateTemperatureTransform(datafr = readings.aggregate,
                             temptransform = feelslike.ns,
                             methodname = "NS_FeelsLike",
                             trimsize = 0,
                             logtransform = FALSE)

feelslike.logtrnfm.ns <- ns(feelslike,
                            knots = c(4, 16, 27))
# Column names must be better than V1, ..., V4 so that my evaluation function doesn't break.
colnames(feelslike.logtrnfm.ns) <- c("ns_sec1", "ns_sec2", "ns_sec3", "ns_sec4")
EvaluateTemperatureTransform(datafr = readings.aggregate,
                             temptransform = feelslike.logtrnfm.ns,
                             methodname = "NS_FeelsLike",
                             trimsize = 0,
                             logtransform = TRUE)



#####
# Iterate through 1-6 lags of drybulb switching regression.
#####
for (lag in 1:6) {
  sr.drybulb.lag.basismatrix <- SwitchingRegressionBasisMatrix(explvar = head(x = c(rep(NA, lag), readings.aggregate$temperature),
                                                                              n = -lag),
                                                               breakpoint = drybulb.break)
  EvaluateTemperatureTransform(datafr = readings.aggregate,
                               temptransform = sr.drybulb.lag.basismatrix,
                               methodname = paste0("SR_Drybulb_Lag", lag),
                               trimsize = lag,
                               logtransform = FALSE)

  sr.logtrnfm.drybulb.lag.basismatrix <- SwitchingRegressionBasisMatrix(explvar = head(x = c(rep(NA, lag), readings.aggregate$temperature),
                                                                                       n = -lag),
                                                                        breakpoint = drybulb.logtrnfm.break)
  EvaluateTemperatureTransform(datafr = readings.aggregate,
                               temptransform = sr.logtrnfm.drybulb.lag.basismatrix,
                               methodname = paste0("SR_Drybulb_Lag", lag),
                               trimsize = lag,
                               logtransform = TRUE)
}



#####
# Iterate through 1-6 lags of drybulb natural splines.
#####
for (lag in 1:6) {
  drybulb.ns <- ns(head(c(rep(NA, lag), readings.aggregate$temperature), -lag),
                   knots = c(3, 23, 30))
  # Column names must be better than V1, ..., V4 so that my evaluation function doesn't break.
  colnames(drybulb.ns) <- c("ns_sec1", "ns_sec2", "ns_sec3", "ns_sec4")
  EvaluateTemperatureTransform(datafr = readings.aggregate,
                               temptransform = drybulb.ns,
                               methodname = paste0("NS_Drybulb_Lag", lag),
                               trimsize = lag,
                               logtransform = FALSE)

  drybulb.logtrnfm.ns <- ns(head(c(rep(NA, lag), readings.aggregate$temperature), -lag),
                            knots = c(3, 16, 31))
  # Column names must be better than V1, ..., V4 so that my evaluation function doesn't break.
  colnames(drybulb.logtrnfm.ns) <- c("ns_sec1", "ns_sec2", "ns_sec3", "ns_sec4")
  EvaluateTemperatureTransform(datafr = readings.aggregate,
                               temptransform = drybulb.logtrnfm.ns,
                               methodname = paste0("NS_Drybulb_Lag", lag),
                               trimsize = lag,
                               logtransform = TRUE)
}


#####
# Iterate through 1-6 lags of "feels like" switching regression.
#####
for (lag in 1:6) {
  sr.feelslike.lag.basismatrix <- SwitchingRegressionBasisMatrix(explvar = head(x = c(rep(NA, lag), feelslike),
                                                                                n = -lag),
                                                                 breakpoint = feelslike.break)
  EvaluateTemperatureTransform(datafr = readings.aggregate,
                               temptransform = sr.feelslike.lag.basismatrix,
                               methodname = paste0("SR_FeelsLike_Lag", lag),
                               trimsize = lag,
                               logtransform = FALSE)

  sr.logtrnfm.feelslike.lag.basismatrix <- SwitchingRegressionBasisMatrix(explvar = head(x = c(rep(NA, lag), feelslike),
                                                                                         n = -lag),
                                                                          breakpoint = feelslike.logtrnfm.break)
  EvaluateTemperatureTransform(datafr = readings.aggregate,
                               temptransform = sr.logtrnfm.feelslike.lag.basismatrix,
                               methodname = paste0("SR_FeelsLike_Lag", lag),
                               trimsize = lag,
                               logtransform = TRUE)
}



#####
# Iterate through 1-6 lags of "feels like" natural splines
#####
for (lag in 1:6) {
  feelslike.ns <- ns(head(c(rep(NA, lag), feelslike), -lag),
                     knots = c(4, 22, 27))
  # Column names must be better than V1, ..., V4 so that my evaluation function doesn't break.
  colnames(feelslike.ns) <- c("ns_sec1", "ns_sec2", "ns_sec3", "ns_sec4")
  EvaluateTemperatureTransform(datafr = readings.aggregate,
                               temptransform = feelslike.ns,
                               methodname = paste0("NS_FeelsLike_Lag", lag),
                               trimsize = lag,
                               logtransform = FALSE)

  feelslike.logtrnfm.ns <- ns(head(c(rep(NA, lag), feelslike), -lag),
                              knots = c(4, 16, 27))
  # Column names must be better than V1, ..., V4 so that my evaluation function doesn't break.
  colnames(feelslike.logtrnfm.ns) <- c("ns_sec1", "ns_sec2", "ns_sec3", "ns_sec4")
  EvaluateTemperatureTransform(datafr = readings.aggregate,
                               temptransform = feelslike.logtrnfm.ns,
                               methodname = paste0("NS_FeelsLike_Lag", lag),
                               trimsize = lag,
                               logtransform = TRUE)
}



#####
# Use a rolling sum of CDH and HDH modelled with a switching regression. The
# optimal degree-hour window is 6 hours, discovered empirically through an
# iteration of 1-24 hour windows.
#####
dh.window <- 6
dh.trimsize <- (dh.window - 1)

dh.drybulb.basismatrix <- DegreeHourBasisMatrix(temperatures = readings.aggregate$temperature,
                                                balance = drybulb.break,
                                                window = dh.window)
EvaluateTemperatureTransform(datafr = readings.aggregate,
                             temptransform = dh.drybulb.basismatrix,
                             methodname = "DH_Drybulb",
                             trimsize = dh.trimsize,
                             logtransform = FALSE)

dh.drybulb.logtrnfm.basismatrix <- DegreeHourBasisMatrix(temperatures = readings.aggregate$temperature,
                                                         balance = drybulb.logtrnfm.break,
                                                         window = dh.window)
EvaluateTemperatureTransform(datafr = readings.aggregate,
                             temptransform = dh.drybulb.logtrnfm.basismatrix,
                             methodname = "DH_Drybulb",
                             trimsize = dh.trimsize,
                             logtransform = TRUE)



#####
# Use a rolling sum of CDH and HDH using "feels like" temperatures, modelled
# with a switching regression. The optimal degree-hour window is 6 hours,
# discovered empirically through an iteration of 1-24 hour windows.
#####
dh.window <- 6
dh.trimsize <- (dh.window - 1)

dh.feelslike.basismatrix <- DegreeHourBasisMatrix(temperatures = feelslike,
                                                  balance = feelslike.break,
                                                  window = dh.window)
EvaluateTemperatureTransform(datafr = readings.aggregate,
                             temptransform = dh.feelslike.basismatrix,
                             methodname = "DH_FeelsLike",
                             trimsize = dh.trimsize,
                             logtransform = FALSE)

dh.feelslike.logtrnfm.basismatrix <- DegreeHourBasisMatrix(temperatures = feelslike,
                                                           balance = feelslike.logtrnfm.break,
                                                           window = dh.window)
EvaluateTemperatureTransform(datafr = readings.aggregate,
                             temptransform = dh.feelslike.logtrnfm.basismatrix,
                             methodname = "DH_FeelsLike",
                             trimsize = dh.trimsize,
                             logtransform = TRUE)



#####
# Transform temperature to a moving average of the past 4 hours of temperature.
# Optimal MA window is 6 hours, discovered empirically through an iteration of
# 1-24 hour windows.
#####
sr.ma.drybulb.basismatrix <- SwitchingRegressionBasisMatrix(explvar = ma.drybulb,
                                                            breakpoint = ma.drybulb.break)
EvaluateTemperatureTransform(datafr = readings.ma.pretrim,
                             temptransform = sr.ma.drybulb.basismatrix,
                             methodname = "MA_SR_Drybulb",
                             trimsize = 0,
                             logtransform = FALSE)

sr.ma.logtrnfm.drybulb.basismatrix <- SwitchingRegressionBasisMatrix(explvar = ma.drybulb,
                                                                     breakpoint = ma.drybulb.logtrnfm.break)
EvaluateTemperatureTransform(datafr = readings.ma.pretrim,
                             temptransform = sr.ma.logtrnfm.drybulb.basismatrix,
                             methodname = "MA_SR_Drybulb",
                             trimsize = 0,
                             logtransform = TRUE)



#####
# Transform temperature to a moving average of the past 4 hours of temperature.
# Optimal MA window is 6 hours, discovered empirically through an iteration of
# 1-24 hour windows.
#####
# Empirical, optimal natural spline knots were found for the temperature MA
# vector. Trimming has been done beforehand so that a full rank matrix may be
# used for the natural splines transform.
ma.drybulb.ns <- ns(ma.drybulb, knots = c(6, 24, 27))
# Column names must be better than V1, ..., V4 so that my evaluation function
# doesn't break.
colnames(ma.drybulb.ns) <- c("ns_sec1", "ns_sec2", "ns_sec3", "ns_sec4")
EvaluateTemperatureTransform(datafr = readings.ma.pretrim,
                             temptransform = ma.drybulb.ns,
                             methodname = "MA_NS_Drybulb",
                             trimsize = 0,
                             logtransform = FALSE)

ma.drybulb.logtrnfm.ns <- ns(ma.drybulb, knots = c(5, 18, 27))
# Column names must be better than V1, ..., V4 so that my evaluation function
# doesn't break.
colnames(ma.drybulb.logtrnfm.ns) <- c("ns_sec1", "ns_sec2", "ns_sec3", "ns_sec4")
EvaluateTemperatureTransform(datafr = readings.ma.pretrim,
                             temptransform = ma.drybulb.logtrnfm.ns,
                             methodname = "MA_NS_Drybulb",
                             trimsize = 0,
                             logtransform = TRUE)



#####
# Transform temperature to a moving average of the past 4 hours of temperature.
# Optimal MA window is 6 hours, discovered empirically through an iteration of
# 1-24 hour windows.
#####
sr.ma.feelslike.basismatrix <- SwitchingRegressionBasisMatrix(explvar = ma.feelslike,
                                                              breakpoint = ma.feelslike.break)
EvaluateTemperatureTransform(datafr = readings.ma.pretrim,
                             temptransform = sr.ma.feelslike.basismatrix,
                             methodname = "MA_SR_FeelsLike",
                             trimsize = 0,
                             logtransform = FALSE)

sr.ma.logtrnfm.feelslike.basismatrix <- SwitchingRegressionBasisMatrix(explvar = ma.feelslike,
                                                                       breakpoint = ma.feelslike.logtrnfm.break)
EvaluateTemperatureTransform(datafr = readings.ma.pretrim,
                             temptransform = sr.ma.logtrnfm.feelslike.basismatrix,
                             methodname = "MA_SR_FeelsLike",
                             trimsize = 0,
                             logtransform = TRUE)



#####
# Transform temperature to a moving average of the past 4 hours of temperature.
# Optimal MA window is 6 hours, discovered empirically through an iteration of
# 1-24 hour windows.
#####
# Empirical, optimal natural spline knots were found for the temperature MA
# vector. Trimming has been done beforehand so that a full rank matrix may be
# used for the natural splines transform.
ma.feelslike.ns <- ns(ma.feelslike, knots = c(7, 23, 26))
# Column names must be better than V1, ..., V4 so that my evaluation function
# doesn't break.
colnames(ma.feelslike.ns) <- c("ns_sec1", "ns_sec2", "ns_sec3", "ns_sec4")
EvaluateTemperatureTransform(datafr = readings.ma.pretrim,
                             temptransform = ma.feelslike.ns,
                             methodname = "MA_NS_FeelsLike",
                             trimsize = 0,
                             logtransform = FALSE)

ma.feelslike.logtrnfm.ns <- ns(ma.feelslike, knots = c(4, 18, 26))
# Column names must be better than V1, ..., V4 so that my evaluation function
# doesn't break.
colnames(ma.feelslike.logtrnfm.ns) <- c("ns_sec1", "ns_sec2", "ns_sec3", "ns_sec4")
EvaluateTemperatureTransform(datafr = readings.ma.pretrim,
                             temptransform = ma.feelslike.logtrnfm.ns,
                             methodname = "MA_NS_FeelsLike",
                             trimsize = 0,
                             logtransform = TRUE)



#####
# Create exposure-lag-response cross-basis matrix for "feels like" temperature
# and time.
#####
maxlag <- 6
cb.sr.drybulb <- crossbasis(x = readings.aggregate$temperature,
                            argvar = list(fun = "thr",
                                          thr.value = drybulb.break,
                                          side = "d"),
                            lag = maxlag,
                            arglag = list(fun = "poly", degree = 3))
EvaluateTemperatureTransform(datafr = readings.aggregate,
                             temptransform = cb.sr.drybulb,
                             methodname = "ELR_SR_Drybulb",
                             trimsize = maxlag,
                             logtransform = FALSE)

cb.sr.drybulb.logtrnfm <- crossbasis(x = readings.aggregate$temperature,
                                     argvar = list(fun = "thr",
                                                   thr.value = drybulb.break,
                                                   side = "d"),
                                     lag = maxlag,
                                     arglag = list(fun = "poly", degree = 3))
EvaluateTemperatureTransform(datafr = readings.aggregate,
                             temptransform = cb.sr.drybulb.logtrnfm,
                             methodname = "ELR_SR_Drybulb",
                             trimsize = maxlag,
                             logtransform = TRUE)



#####
# Create exposure-lag-response cross-basis matrix for temperature and time
#####
maxlag <- 6

cb.drybulb <- crossbasis(x = readings.aggregate$temperature,
                         argvar = list(fun = "ns", knots = c(2, 24, 30)),
                         lag = maxlag,
                         arglag = list(fun = "poly", degree = 3))
EvaluateTemperatureTransform(datafr = readings.aggregate,
                             temptransform = cb.drybulb,
                             methodname = "ELR_NS_Drybulb",
                             trimsize = maxlag,
                             logtransform = FALSE)

cb.drybulb.logtrnfm <- crossbasis(x = readings.aggregate$temperature,
                                  argvar = list(fun = "ns", knots = c(7, 17, 28)),
                                  lag = maxlag,
                                  arglag = list(fun = "poly", degree = 3))
EvaluateTemperatureTransform(datafr = readings.aggregate,
                             temptransform = cb.drybulb.logtrnfm,
                             methodname = "ELR_NS_Drybulb",
                             trimsize = maxlag,
                             logtransform = TRUE)



#####
# Create exposure-lag-response cross-basis matrix for "feels like" temperature
# and time.
#####
maxlag <- 6
cb.sr.feelslike <- crossbasis(x = feelslike,
                              argvar = list(fun = "thr",
                                            thr.value = feelslike.break,
                                            side = "d"),
                              lag = maxlag,
                              arglag = list(fun = "poly", degree = 3))
EvaluateTemperatureTransform(datafr = readings.aggregate,
                             temptransform = cb.sr.feelslike,
                             methodname = "ELR_SR_FeelsLike",
                             trimsize = maxlag,
                             logtransform = FALSE)

cb.sr.feelslike.logtrnfm <- crossbasis(x = feelslike,
                                       argvar = list(fun = "thr",
                                                     thr.value = feelslike.break,
                                                     side = "d"),
                                       lag = maxlag,
                                       arglag = list(fun = "poly", degree = 3))
EvaluateTemperatureTransform(datafr = readings.aggregate,
                             temptransform = cb.sr.feelslike.logtrnfm,
                             methodname = "ELR_SR_FeelsLike",
                             trimsize = maxlag,
                             logtransform = TRUE)



#####
# Create exposure-lag-response cross-basis matrix for "feels like" temperature
# and time.
#####
maxlag <- 6
cb.ns.feelslike <- crossbasis(x = feelslike,
                              argvar = list(fun = "ns", knots = c(9, 24, 25)),
                              lag = maxlag,
                              arglag = list(fun = "poly", degree = 3))
EvaluateTemperatureTransform(datafr = readings.aggregate,
                             temptransform = cb.ns.feelslike,
                             methodname = "ELR_NS_FeelsLike",
                             trimsize = maxlag,
                             logtransform = FALSE)

cb.ns.feelslike.logtrnfm <- crossbasis(x = feelslike,
                                       argvar = list(fun = "ns", knots = c(7, 17, 26)),
                                       lag = maxlag,
                                       arglag = list(fun = "poly", degree = 3))
EvaluateTemperatureTransform(datafr = readings.aggregate,
                             temptransform = cb.ns.feelslike.logtrnfm,
                             methodname = "ELR_NS_FeelsLike",
                             trimsize = maxlag,
                             logtransform = TRUE)



#####
# Create an illustrative countour plot for the 6-hour exposure-lag-response
# association.
#####
lm.cb.drybulb <- lm(kwh ~ cb.drybulb + hrstr + working_day + price + rateseason + tou_active +
                      hrstr:working_day,
                    data = readings.aggregate)
summary(lm.cb.drybulb)$adj.r.squared
cb.pred <- crosspred(cb.drybulb, lm.cb.drybulb)
plot(cb.pred, ptype="contour",
     key.title = title("kWh"),
     plot.title = title("Contour Plot of Exposure-Lag-Response Association",
                        xlab = "Outdoor, Dry-bulb Temperature (Celsius)",
                        ylab = "Lag (hours)"))
dev.print(file = "../../Figures/ExposureLagResponseAssociationContourPlot6Hours.png",
          device = png, height = 600, width = 800)



#####
# Create illustrative plots of rolling average residuals for untransformed and
# log transformed electricity demand.
#####
ma.order <- 6
ma.trimsize <- (ma.order - 1)
ma.drybulb <- rollmean(x = readings.aggregate$temperature, k = ma.order)

ma.drybulb.ns <- ns(ma.drybulb, knots = c(6, 24, 27))
colnames(ma.drybulb.ns) <- c("ns_sec1", "ns_sec2", "ns_sec3", "ns_sec4")
readings.ma.trim <- cbind(tail(readings.aggregate, -ma.trimsize), ma.drybulb.ns)
lm.ma <- lm(kwh ~ ns_sec1 + ns_sec2 + ns_sec3 + ns_sec4  +
              hrstr + working_day + price + rateseason + tou_active +
              hrstr:working_day,
            data = readings.ma.trim)


ma.drybulb.logtrnfm.ns <- ns(ma.drybulb, knots = c(5, 18, 27))
colnames(ma.drybulb.logtrnfm.ns) <- c("ns_sec1", "ns_sec2", "ns_sec3", "ns_sec4")
readings.ma.logtrnfm.trim <- cbind(tail(readings.aggregate, -ma.trimsize), ma.drybulb.logtrnfm.ns)
lm.ma.logtrnfm <- lm(log(kwh) ~ ns_sec1 + ns_sec2 + ns_sec3 + ns_sec4  +
                       hrstr + working_day + price + rateseason + tou_active +
                       hrstr:working_day,
                     data = readings.ma.logtrnfm.trim)
ma.bktrnfm.fitted <- exp(lm.ma.logtrnfm$fitted)
ma.bktrnfm.resids <- readings.ma.logtrnfm.trim$kwh - ma.bktrnfm.fitted


# Distribution of residuals
par(mfrow = c(1, 1))
resid.bw = 0.015
plot(density(lm.ma$residuals, bw = resid.bw),
     main = paste0("Regression of Aggregate Elec. Demand on Six-Hour Moving Average, Dry-bulb,",
                   "\nNatural Cubic Spline Temp. Transform: Residual Density (N=", 
                   nrow(readings.ma.trim), 
                   ", Bandwidth=", 
                   resid.bw, 
                   ")"),
     ylim = c(0, 4),
     xlim = c(-1.8, 1.8),
     xlab = "Residuals (kWh)",
     lwd = 3)
abline(v = 0, col = "blue", lty = 2, lwd = 1)
# plot(density(ma.bktrnfm.resids, bw = resid.bw),
#      main = paste0("Regression on log(kWh) with 6-Hour Rolling Average Dry-bulb Temperature Transform",
#                    "\nResidual Density (N=", nrow(readings.ma.logtrnfm.trim), ", Bandwidth=", resid.bw, ")"),
#      ylim = c(0, 4),
#      xlim = c(-1.8, 1.8),
#      xlab = "Backtransformed Residuals (kWh)",
#      lwd = 3)
# abline(v = 0, col = "blue", lty = 2, lwd = 1)
dev.print(file = "../../Figures/ResultsResidualDistributions.png", device = png, height = 600, width = 800)
par(mfrow = c(1, 1))

# Residuals as a function of time
lm.ma.resid.ts.datafr <- data.frame("posixtime" = readings.ma.trim$timestamp_dst,
                                    "resids" = lm.ma$residuals)
ma.resid.ts.plot <- ResidualsOverTimePlot(thedata = lm.ma.resid.ts.datafr,
                                          thetitle = paste0("Regression of Aggregate Elec. Demand on Six-Hour Moving Average, Dry-bulb,",
                                                            "\nNatural Cubic Spline Temp. Transform: Residuals Over Time"))
# lm.ma.bktrnfm.resid.ts.datafr <- data.frame("posixtime" = readings.ma.logtrnfm.trim$timestamp_dst,
#                                             "resids" = ma.bktrnfm.resids)
# ma.bktrnfm.resid.ts.plot <- ResidualsOverTimePlot(thedata = lm.ma.bktrnfm.resid.ts.datafr,
#                                                   thetitle = paste0("Regression on log(kWh) with 6-Hour Rolling Average Dry-bulb Temperature Transform",
#                                                                     "\nResiduals Plotted Over Time"),
#                                                   ylabel = "Backtransformed Residuals (kWh)")
# grid.arrange(ma.resid.ts.plot, ma.bktrnfm.resid.ts.plot, ncol = 2, nrow = 1)
print(ma.resid.ts.plot)
dev.print(file="../../Figures/ResultsResidualTimeseries.png", device=png, height = 600, width = 800)

# Residuals as a function of temperature
lm.ma.resid.temp.datafr <- data.frame("temperatures" = readings.ma.trim$temperature,
                                      "resids" = lm.ma$residuals)
ma.resid.temp.plot <- ResidualsByTemperature(thedata = lm.ma.resid.temp.datafr,
                                             thetitle = paste0("Regression of Aggregate Elec. Demand on Six-Hour Moving Average, Dry-bulb,",
                                                               "\nNatural Cubic Spline Temp. Transform: Residuals vs. Observed Dry-bulb Temperature"),
                                             xlabel = "Outdoor, Dry-bulb Temperature (Celsius)",
                                             ylabel = "Residuals (kWh)")
# lm.ma.bktrnfm.resid.temp.datafr <- data.frame("temperatures" = readings.ma.logtrnfm.trim$temperature,
#                                               "resids" = ma.bktrnfm.resids)
# ma.resid.bktrnfm.temp.plot <- ResidualsByTemperature(thedata = lm.ma.bktrnfm.resid.temp.datafr,
#                                                      thetitle = paste0("Regression on log(kWh) with 6-Hour Rolling Average Dry-bulb Temperature Transform",
#                                                                        "\nResiduals Plotted as a Function of Observed Dry-bulb Temperature"),
#                                                      xlabel = "Dry-bulb Temperature (Celsius)",
#                                                      ylabel = "Backtransformed Residuals (kWh)")
# grid.arrange(ma.resid.temp.plot, ma.resid.bktrnfm.temp.plot, ncol = 2, nrow = 1)
print(ma.resid.temp.plot)
dev.print(file="../../Figures/ResultsResidualByTemperature.png", device=png, height = 600, width = 800)

# Residuals as a function of estimated response
ma.resid.resp.plot <- ResidualsByResponseVal(thedata = data.frame("fittedvals" = lm.ma$fitted,
                                                                  "resids" = lm.ma$residuals),
                                             thetitle = paste0("Regression of Aggregate Elec. Demand on Six-Hour Moving Average, Dry-bulb,",
                                                               "\nNatural Cubic Spline Temp. Transform: Residuals vs. Estimated Response"))
# ma.bktrnfm.resid.resp.plot <- ResidualsByResponseVal(thedata = data.frame("fittedvals" = ma.bktrnfm.fitted,
#                                                                           "resids" = ma.bktrnfm.resids),
#                                                      thetitle = paste0("Regression on log(kWh) with 6-Hour Rolling Average Dry-bulb Temperature Transform",
#                                                                        "\nResiduals Plotted as a Function of Backtransformed Response Estimate"),
#                                                      ylabel = "Backtransformed Residuals (kWh)")
# grid.arrange(ma.resid.resp.plot, ma.bktrnfm.resid.resp.plot, ncol = 2, nrow = 1)
print(ma.resid.resp.plot)
dev.print(file="../../Figures/ResultsResidualByFittedResponse.png", device=png, height = 600, width = 800)