# See https://stats.stackexchange.com/questions/119835/correlation-between-a-nominal-iv-and-a-continuous-dv-variable/124618
# See https://stats.stackexchange.com/questions/108007/correlations-with-categorical-variables

# "MULTICOLINEARITY" is used to describe colinear relationship between 
# independent variables in a multiple regression.
#
# VARIANCE INFLATION FACTOR

# Source functions in other files
source("dataframe_processing.R")

# Load SmartMeterReading data from CSV
home <- Sys.getenv("HOME")
fpath <- file.path(home, 
                   "../Dropbox/ISS4E/R", 
                   "aggregate_readings_01Mar2011_through_17Oct2012.csv")

readings.aggregate <- InitReadingsDataFrame(fpath = fpath, 
                                            is.aggregate = TRUE)
rm(fpath)

# Load weather descriptions from CSV (for largest city in LDC)
fpath2 <- file.path(home, 
                    "../Dropbox/ISS4E/R", 
                    "weather_desc_01Mar2011_through_17Oct2012.csv")
weather <- read.csv(fpath2, 
                    na.strings = c("NULL", "NA", "NaN"), 
                    stringsAsFactors = FALSE)
rm(fpath2)

fpath3 <- file.path(home, 
                    "../Dropbox/ISS4E/R", 
                    "aggregate_exceptions_01Mar2011_through_17Oct2012.csv")
exceptions <- read.csv(fpath3, 
                       colClasses = c("integer", "POSIXct", "numeric"),
                       na.strings = c("NULL", "NA", "NaN"), 
                       stringsAsFactors = FALSE)
rm(fpath3, home)

# Summary of exceptions data
summary(exceptions)

# Reduce weather descriptions to a simplified set of factors
readings.aggregate$weather_reduced <- ReduceWeather(weather$weather_desc)
rm(weather)

# For clarity, reorder explanatory variables which are factors
readings.aggregate <- OrderFactors(readings.aggregate)    

# Past temperatures
temp.hrs <- 2
temps <- CreatePastTemperatureMatrix(nlags = temp.hrs, 
                                     df.readings = readings.aggregate)
readings.aggregate$temperature <- temps[,(temp.hrs+1)]
#rm(temps)

# Centered polynomials
poly.deg <- 5
poly.vars <- poly(readings.aggregate$temperature, degree = poly.deg)
colnames(poly.vars) <- paste0("temp_poly", c(1:poly.deg))
readings.aggregate <- cbind(readings.aggregate, poly.vars)
rm(poly.vars)


# Set the weather description to be the same as the past hour of temperature 
# used (+1 due to weatherdesc_lag0)
pastweather <- CreatePastWeatherDescDataFrame(nlags = temp.hrs, 
                                              weather_reduced = readings.aggregate$weather_reduced)
readings.aggregate$weather_reduced <- pastweather[,(temp.hrs+1)]
#rm(pastweather)

# Trim up data frame
readings.trimmed <- TrimExplanatoryVariables(readings.aggregate)

# Don't include the first several rows due to missing past temperature info
readings.trimmed <- tail(x= readings.trimmed, 
                         n = nrow(readings.trimmed) - temp.hrs)

# Contingency tables for all combinations of factors
contable.hrstr_month <- table(readings.trimmed$hrstr, readings.trimmed$month)
contable.hrstr_month
chisq.test(contable.hrstr_month) # Uncorrelated, as expected

contable.hrstr_dayname <- table(readings.trimmed$hrstr, readings.trimmed$dayname)
contable.hrstr_dayname
chisq.test(contable.hrstr_dayname) # Uncorrelated, as expected

contable.hrstr_holiday <- table(readings.trimmed$hrstr, readings.trimmed$holiday)
contable.hrstr_holiday
chisq.test(contable.hrstr_holiday) # Uncorrelated, as expected

contable.hrstr_price <- table(readings.trimmed$hrstr, readings.trimmed$price)
contable.hrstr_price
chisq.test(contable.hrstr_price) # Highly correlated, as expected

contable.hrstr_weather <- table(readings.trimmed$hrstr, readings.trimmed$weather_reduced)
contable.hrstr_weather
chisq.test(contable.hrstr_weather) # Highly correlated, reasonable

contable.month_dayname <- table(readings.trimmed$month, readings.trimmed$dayname)
contable.month_dayname
chisq.test(contable.month_dayname) # Correlated at p < 0.05, due to the low 
                                   # aggregate count pre-March 2011. Could be 
                                   # corrected with the expanded/weighted data.

contable.month_holiday <- table(readings.trimmed$month, readings.trimmed$holiday)
contable.month_holiday
chisq.test(contable.month_holiday) # Highly correlated, reasonable.
                                   # Apparently there are no vacations in April 
                                   # or November.

contable.month_price <- table(readings.trimmed$month, readings.trimmed$price)
contable.month_price
chisq.test(contable.month_price) # Highly correlated, as expected

contable.month_weather <- table(readings.trimmed$month, readings.trimmed$weather_reduced)
contable.month_weather
chisq.test(contable.month_weather) # Highly correlated, as expected

contable.dayname_holiday <- table(readings.trimmed$dayname, readings.trimmed$holiday)
contable.dayname_holiday
chisq.test(contable.dayname_holiday) # Highly correlated, reasonable.

contable.dayname_price <- table(readings.trimmed$dayname, readings.trimmed$price)
contable.dayname_price
chisq.test(contable.dayname_price) # Highly correlated, expected.

contable.dayname_weather <- table(readings.trimmed$dayname, readings.trimmed$weather_reduced)
contable.dayname_weather
chisq.test(contable.dayname_weather) # Highly correlated, unexpected. How can 
                                     # weather be correlated to the day of the 
                                     # week?

contable.holiday_price <- table(readings.trimmed$holiday, readings.trimmed$price)
contable.holiday_price
chisq.test(contable.holiday_price) # Highly correlated, expected.

contable.holiday_weather <- table(readings.trimmed$holiday, readings.trimmed$weather_reduced)
contable.holiday_weather
chisq.test(contable.holiday_weather) # Highly correlated, unexpected. This may 
                                     # just be spurious due to low amounts of 
                                     # data though.

contable.price_weather <- table(readings.trimmed$price, readings.trimmed$weather_reduced)
contable.price_weather
chisq.test(contable.price_weather) # Highly correlated, unexpected. Though, 
                                   # price reflects time-of-day, and weather 
                                   # patterns are correlated with hour-of-day.

# Correlation matrix for continuous variables
readings.continuous <- subset(x = readings.aggregate, 
                              select = c(kwh, temperature, agg_count))
# kwh and temperature correlation = 0.59 (fairly correlated)
# aggregate count and kwh correlation = 0.05 (uncorrelated)
# aggregate count and temperature correlation = 0.14 (fairly uncorrelated)
cor(readings.continuous)
pairs(readings.continuous)

# Check correlation matrix with temperature converted to polynomials
readings.continuous.temp_as_poly <- subset(x = readings.aggregate,
                                           select = c(kwh, agg_count, 
                                                      temp_poly1, temp_poly2,
                                                      temp_poly3, temp_poly4, 
                                                      temp_poly5))
# Correlation matrix seems to suggest that cubic polynomial is a good 
# stopping point.
cor(readings.continuous.temp_as_poly)
pairs(readings.continuous.temp_as_poly)

# Experiment with lognormal transformation of response variable and correlations
log_readings.continuous.temp_as_poly <- readings.continuous.temp_as_poly
log_readings.continuous.temp_as_poly$kwh <- log(readings.continuous.temp_as_poly$kwh)
cor(log_readings.continuous.temp_as_poly)
pairs(log_readings.continuous.temp_as_poly)


# Check correlation of exceptions with data
obs_exn.continuous <- data.frame(kwh = readings.aggregate$kwh,
                                agg_read_count = readings.aggregate$agg_count,
                                agg_exn_count = exceptions$agg_exn_count,
                                avg_exn_reading = exceptions$avg_exn_reading,
                                temperature = readings.aggregate$temperature)
# No correaltions between the exception data and the readings.
# However, there may be a small correlation between temperature and aggregate reading count.
cor(obs_exn.continuous)
pairs(obs_exn.continuous)

# Possibly lower aggregate counts in the temperature range -10C to 10C?
# I don't think it's that strong of a relationship
plot(x = readings.aggregate$temperature, y = readings.aggregate$agg_count)

# Plot reading count and exception count together
par(mfrow=c(2,1))
plot(readings.aggregate$agg_count)
plot(exceptions$agg_exn_count)


# Use ANOVA to determine whether categorical variables are correlated with 
# temperature
par(mfrow=c(1,1))
hist(readings.aggregate$temperature, breaks = 20)

# Temperature and hour-of-day are highly correlated, as expected
anova(lm(temperature ~ 1, data = readings.aggregate),
      lm(temperature ~ hrstr, data = readings.aggregate))

# Temperature and month are highly correlated, as expected
anova(lm(temperature ~ 1, data = readings.aggregate),
      lm(temperature ~ month, data = readings.aggregate))

# ANOVA indicates that temperature is correlated with dayname, at p<0.05.
# There should be no correlation
anova(lm(temperature ~ 1, data = readings.aggregate),
      lm(temperature ~ dayname, data = readings.aggregate))

# Temperature and holiday are not correlated, as expected
anova(lm(temperature ~ 1, data = readings.aggregate),
      lm(temperature ~ holiday, data = readings.aggregate))

# Temperature and price are correlated, which is reasonable. Prices are high 
# mid-day because electricity use is high. Electricity use is high because 
# temperature is high. Through that relationship, it is reasonable that 
# temperature and price are correlated.
anova(lm(temperature ~ 1, data = readings.aggregate),
      lm(temperature ~ price, data = readings.aggregate))