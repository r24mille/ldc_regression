library(ggplot2) # Better plotting!
library(scales) # For scaling axes (eg. datetime)
library(car) # Companion to Applied Regression (better qqPlot)
library(reshape2) # For reshaping (ie. melting) data
library(lme4) # Linear Mixed Effects Model

# Source functions in other files
source("dataframe_processing.R")
source("individual_readings_quality_checks.R")

# Load SmartMeterReading data from CSV
home <- Sys.getenv("HOME")
fpath <- file.path(home, 
                   "../Dropbox/ISS4E/R", 
                   "r24mille_lasso_25sample.csv")
readings.individual <- InitReadingsDataFrame(fpath = fpath, 
                                            is.aggregate = FALSE)

# Reduce weather descriptions to a factor with 7 levels
readings.individual$weather_reduced <- ReduceWeather(readings.individual$weather_desc)

# For clarity, reorder explanatory variables which are factors
readings.individual <- OrderFactors(readings.individual)

# Derive a few temperature humidity index (THI) variables according to Navigant 
# analysis white paper.
readings.individual$nvgnt_thi <- (17.5 
                                 + (0.55 * readings.individual$temperature) 
                                 + (0.2 * readings.individual$dewpnt_temp))
readings.individual$nvgnt_cool_thi <- sapply(readings.individual$nvgnt_thi, 
                                            function(x) max(x - 30, 0))
readings.individual$nvgnt_heat_thi <- sapply(readings.individual$nvgnt_thi, 
                                            function(x) max(25 - x, 0))

# With regard to the flags for the 1st and 2nd months of a season. If there are 
# two flags, then the 1st and 2nd months will be multiplied by some 
# coefficient. The 3rd month does not deviate from the mean. This is 
# functionally equivalent to using two boolean indicator terms. I'll take 
# "first" and "second" for the shoulder seasons to be the first and second 
# months as they are listed on page 4.
season_1stmnths <- c("m6", "m12", "m5", "m11")
readings.individual$season_1stmnth <- readings.individual$month %in% season_1stmnths
readings.individual$season_1stmnth <- as.factor(readings.individual$season_1stmnth)

season_2ndmnths <- c("m7", "m1", "m9", "m3")
readings.individual$season_2ndmnth <- readings.individual$month %in% season_2ndmnths
readings.individual$season_2ndmnth <- as.factor(readings.individual$season_2ndmnth)

# TOU flag indicates TRUE/FALSE whether TOU billing was active during 
# the observation.
readings.individual$tou_billing <- readings.individual$price != "flat"
readings.individual$tou_billing <- as.factor(readings.individual$tou_billing)

# Perform data quality checks
#expl_var_dist_checks(df = readings.individual)

# Break data into sets as described on pg. 5 of the report
# Rate Class:
#    All meters in my LDC dataset are residential, TOU-metered households. So 
#    there are no general service households to split out.

# Season:
#    A list (ie. R's map data type) using the season as a key, and value that 
#    is a vector of month strings corresponding to readings.individual$month.
seasons <- list(summer = c("m6", "m7", "m8"),
                winter = c("m12", "m1", "m2"),
                sumshoulder = c("m5", "m9", "m10"),
                winshoulder = c("m11", "m3", "m4"))

# Weekend:
#    A vector of values that weekend can take (ie. "Yes"/"No")
wknd_lvls <- attributes(readings.individual$weekend)$levels

# Hour of the Day:
#   A vector of values that hrstr can take (ie. "h0", "h1", ..., "h23")
hour_lvls <- attributes(readings.individual$hrstr)$levels

# Iterate through each combination of the three division terms
for(s in 1:length(seasons)) {
  for (w in 1:length(wknd_lvls)) {
    for (h in 1:length(hour_lvls)) {
      print(paste0(names(seasons)[s], "Wknd", wknd_lvls[w], hour_lvls[h]))
    }
  }
}

subdiv_lmms <- list()
i <- 1
for(s in 1:1) {
  for (w in 1:1) {
    for (h in 1:1) {
      uniqstr <- paste0(names(seasons)[s], 
                               "Wknd", 
                               wknd_lvls[w], 
                               hour_lvls[h])
      subdiv_lmms[[i]] <- list()
      names(subdiv_lmms)[i] <- uniqstr
      uniqdf <- readings.individual[readings.individual$month %in% seasons[[s]] & 
                                    readings.individual$weekend == wknd_lvls[w] & 
                                    readings.individual$hrstr == hour_lvls[h],]
      subdiv_lmms[[i]]$df <- uniqdf
      
      uniqlmm <- lmer(kwh ~ season_1stmnth + season_2ndmnth + nvgnt_cool_thi 
                      + nvgnt_heat_thi + tou_billing + tou_billing:nvgnt_cool_thi 
                      + tou_billing:nvgnt_heat_thi + (1|meterid),
                      data = uniqdf)
      subdiv_lmms[[i]]$lmm <- uniqlmm
      print(uniqstr)
      i <- i + 1
    }
  }
}
