CreatePastTemperatureMatrix <- function(nlags, datafr.readings) {
  # Creates a matrix of values such that each column looks an additional hour 
  # into the past at degrees over the temperature breakpoint.
  #
  # Args:
  #   nlags: The number of hours (ie. lags) to include in the matrix
  #   datafr.readings: A dataframe of smart meter readings.
  #
  # Returns:
  #   An [n x (nlags + 1)] matrix such that each hour into the past can have its
  #   own coefficient when modelled.
  lagnames <- c(paste0("temp_lag", c(0:nlags)))
  templags <- matrix(nrow = nrow(datafr.readings),
                     ncol = (nlags + 1))
  colnames(templags) <- lagnames
  
  # Unfortunately iterating over the dataframe seems to be the best method
  for(i in 1:nrow(datafr.readings)) {
    for(j in 0:nlags) {
      # Populate the past temperatures
      if (i-j > 0) {
        templags[i, (j+1)] <- datafr.readings[i-j, "temperature"]
      } else {
        templags[i, (j+1)] <- NA
      }
    }
  }
  
  return(templags)
}

CreatePastWeatherDescDataFrame <- function(nlags, weather_reduced) {
  # TODO(r24mille): Thiiiiis is not efficient.
  
  descnames <- c(paste0("weatherdesc_lag", c(0:nlags)))
  pastweather <- data.frame(weatherdesc_lag0 = weather_reduced,
                            stringsAsFactors = TRUE)
  
  # Iterate over the dataframe to create matrix of differences in temperature 
  # between hours.
  for (i in 1:length(weather_reduced)) {
    for (j in 0:nlags) {
      if (i-j > 0) {
        pastweather[i, j+1] <- as.character(weather_reduced[i-j])
      } else {
        pastweather[i, j+1] <- "unknown"
      }
    }
  }
  
  colnames(pastweather) <- descnames
  
  for (i in 1:(nlags+1)) {
    if (i == 1) {
      pastweather[,i] <- factor(pastweather[,i], 
                                c("clear", "cloudy_fog", "rain_tstorms", 
                                  "snow_ice"))
    } else {
      pastweather[,i] <- factor(pastweather[,i], 
                                c("clear", "cloudy_fog", "rain_tstorms", 
                                  "snow_ice", "unknown"))
    }
  }
  
  return(pastweather)
}


FeelsLike <- function(datafr) {
  # Converts dry-bulb temperatures to "feels like" temperatures that reflect 
  # what the apparent temperature would feel like to a person.
  #
  # Args:
  #   datafr: A data.frame with dry bulb temperature in Celsius, relative 
  #           humidity as a percentage, and wind speed in kph.
  #
  # Return:
  #   A vector of apparent temperature values.
  
  feels_like <- rep(NA, nrow(datafr))
  
  for (i in 1:nrow(datafr)){
    # Only calculate heat index and wind chill at valid values (as defined by 
    # NOAA or Environment Canada), otherwise the apparent temperature is the 
    # same as the dry bulb temperature.
    if (datafr$temperature[i] > 27 & datafr$rel_humidity_pct[i] > 40) {
      feels_like[i] <- HeatIndex(datafr$temperature[i], datafr$rel_humidity_pct[i])
    } else if (datafr$temperature[i] <= 10 & datafr$wind_speed_kph[i] > 4.8) {
      feels_like[i] <- WindChill(datafr$temperature[i], datafr$wind_speed_kph[i])
    } else {
      feels_like[i] <- datafr$temperature[i]
    }
  }
  
  return(feels_like)
}

HeatIndex <- function(temp, rel_humidity) {
  # Function computes the heat index as defined by NOAA. See 
  # http://www.srh.noaa.gov/images/ffc/pdf/ta_htindx.PDF
  #
  # Formula is valid for Fahrenheit values so they are converted in the body 
  # of this function but transformed back to celsius in the returns values.
  # 
  # Args:
  #   temp: Temperature value in Celsius
  #   rel_humidity: Relative humidity value as a percentage
  #
  # Return:
  #   The heat index for the provided values (ie. apparent temperature)
  
  # Relative humidity is only valid for temperatures > 27C AND relative 
  # humidity percentages > 40%.
  if (temp > 27 & rel_humidity > 40) {
    # Convert Celsius to Fahrenheit
    temp_f <- (temp * 9/5) + 32
    
    # Constants for heat index equation
    c1 <- -42.379
    c2 <- 2.04901523
    c3 <- 10.14333127
    c4 <- -0.22475541
    c5 <- -6.83783 * 10^-3
    c6 <- -5.481717 * 10^-2
    c7 <- 1.22874 * 10^-3
    c8 <- 8.5282 * 10^-4
    c9 <- -1.99 * 10^-6
    
    # Heat index equation is for fahrenheit values
    heat_index <- (c1 + 
                     (c2*temp_f) + 
                     (c3*rel_humidity) + 
                     (c4*temp_f*rel_humidity) + 
                     (c5*(temp_f^2)) + 
                     (c6*(rel_humidity^2)) + 
                     (c7*(temp_f^2)*rel_humidity) + 
                     (c8*temp_f*(rel_humidity^2)) + 
                     (c9*(temp_f^2)*(rel_humidity^2)))
    
    # Convert back to celsius
    heat_index_c <- (heat_index - 32) * 5/9
    
    return(heat_index_c)
  } else {
    return(temp)
  }
}

HumidexDiff <- function(datafr) {
  # Transform humidex to amount over a humidex threshold (ie. lowest recorded 
  # humidex is 25) using values reported by Environment Canada.
  # http://climate.weather.gc.ca/climate_normals/normals_documentation_e.html#humidex
  #
  # Args:
  #   datafr: A smart meter readings data.frame
  #
  # Return:
  #   A vector representing the amount that the humidex is over a minimum 
  #   threshold.
  humidex_threshold <- min(datafr$humidex, na.rm = TRUE) - 1
  humidex_diff <- datafr$humidex
  humidex_diff[is.na(humidex_diff)] <- humidex_threshold
  humidex_diff <- humidex_diff - humidex_threshold
  return(humidex_diff)
}

InitReadingsDataFrame <- function(fpath, is.aggregate = FALSE) {
  # Initializes a data.frame for the common readings and aggregate readings 
  # column structure. The data.frame is empty but columns have been typed 
  # appropriately with ordered factors.
  readings.colnames <- c("kwh", "sample_index", "daynum", "hrstr", "month", 
                         "year", "hourofday", "dayofweek", "dayofmonth", "dayofyear", 
                         "weekend", "timestamp_dst", "timestamp_std", "dayname", "holiday", 
                         "temperature", "dewpnt_temp", "rel_humidity_pct", 
                         "wind_speed_kph", "humidex", "wind_chill", "price", "tou_active")
  readings.colclasses <- c("numeric", "integer", "integer", "factor", "factor",
                           "integer", "integer", "integer", "integer", "integer", 
                           "factor", "POSIXct", "character", "factor", "factor", 
                           "numeric", "numeric", "numeric",  
                           "numeric", "numeric", "numeric", "factor", "factor")  
  
  if (is.aggregate == TRUE) {
    # Currently, the aggregate version of the CSV has column headers
    has.headers <- TRUE
    
    readings.colnames <- c(readings.colnames, c("agg_count"))
    readings.colclasses <- c(readings.colclasses, c("integer"))
  } else {
    # Currently, the meterid version of the CSV does not have column headers
    has.headers <- FALSE
    
    # The meterid version of the CSV also has two additional columns
    readings.colnames <- c(c("meterid"), readings.colnames, c("weather_desc"))
    readings.colclasses <- c(c("integer"), readings.colclasses, c("character"))
  }
  datafr <- read.csv(file = fpath, 
                     header = has.headers, 
                     col.names = readings.colnames,
                     na.strings = c("NULL", "NA", "NaN", "\\N"),
                     colClasses = readings.colclasses, 
                     stringsAsFactors = FALSE)
  datafr$how_frac <- (((datafr$dayofweek - 1) * 24) + (datafr$hourofday + 1)) / 168 # Hours in a week
  datafr$hom_frac <- (((datafr$dayofmonth - 1) * 24) + (datafr$hourofday + 1)) / 730 # Average number of hours in a month
  datafr$hoy_frac <- (((datafr$dayofyear - 1) * 24) + (datafr$hourofday + 1)) / (8760 + 6) # 2012 is a leap year, so add 6 hours to keep years consistent
  
  # Column type POSIXlt needs a manually-provided timezone for standard time
  datafr$timestamp_std <- as.POSIXlt(x = datafr$timestamp_std, 
                                     format = "%Y-%m-%d %H:%M:%S",
                                     tz = "EST")
  
  return(datafr)
}


NumberFactorLevels <- function(datafr) {
  # Args: 
  #   datafr: A smart meter reading data.frame that has been trimmed and is 
  #           ready for conversion to a matrix processed by ?glinternet.
  #
  # Return:
  #   A vector in which each value represents the number of factor levels for 
  #   each column of the data.frame passed in.
  numlevels <- rep(0, ncol(datafr))
  
  for (i in 1:ncol(datafr)) {
    if (is.factor(datafr[,i])) {
      numlevels[i] <- nlevels(datafr[,i])
    } else {
      numlevels[i] <- 1
    }
  }
  
  return(numlevels)
}

NumericFactorCodedMatrix <- function(datafr) {
  # Transforms an input data.frame to a matrix with factor levels represented 
  # as numeric. Returned value is suitable for ?glinternet.
  #
  # Args: 
  #   datafr: A smart meter reading data.frame that has been trimmed and is 
  #       ready for conversion to a matrix processed by ?glinternet.
  #
  # Return:
  #   The data.frame re-coded as a numeric matrix.
  numeric.m <- matrix(0, nrow = nrow(datafr), ncol = ncol(datafr))
  colnames(numeric.m) <- names(datafr)
  
  for (i in 1:ncol(datafr)) {
    if (is.factor(datafr[,i])) {
      numeric.m[,i] <- as.numeric(datafr[,i]) - 1
    } else {
      numeric.m[,i] <- datafr[,i]
    }
  }
  
  return(numeric.m)
}


OrderFactors <- function(datafr) {
  # Reorder factors so that model contrasts and graphs make a bit more sense
  # 
  # Args: 
  #   datafr: A smart meter reading data.frame that has been run through
  #       fixDataTypes(...) and addInferredInformation(...).
  # 
  # Return: 
  #   The dataframe with corrected/ordered factors.
  if("meterid" %in% colnames(datafr)) {
    datafr$meterid <- factor(datafr$meterid)
  }
  datafr$weekend <- factor(datafr$weekend, c("FALSE", "TRUE"))
  datafr$holiday <- factor(datafr$holiday, c("FALSE", "TRUE"))
  datafr$tou_active <- factor(datafr$tou_active, c("FALSE", "TRUE"))
  datafr$dayname <- factor(datafr$dayname, c("Sun", "Mon", "Tue", "Wed", "Thu", 
                                             "Fri", "Sat"))
  datafr$month <- factor(datafr$month, c("m1", "m2", "m3", "m4", "m5", "m6", 
                                         "m7", "m8", "m9", "m10", "m11", "m12"))
  datafr$hrstr <- factor(datafr$hrstr, 
                         c("h0", "h1", "h2", "h3", "h4", "h5", "h6", "h7", "h8", 
                           "h9", "h10", "h11", "h12", "h13", "h14", "h15", 
                           "h16", "h17", "h18", "h19", "h20", "h21", "h22", 
                           "h23"))
  
  datafr$price <- factor(datafr$price, c("flat", "off_peak", "mid_peak", 
                                         "on_peak"))
  datafr$weather_reduced <- factor(datafr$weather_reduced)
  
  if("severe_weather" %in% colnames(datafr)) {
    datafr$severe_weather <- factor(datafr$severe_weather, c("FALSE", "TRUE"))
  }
  
  if("working_day" %in% colnames(datafr)) {
    datafr$working_day <- factor(datafr$working_day, c("FALSE", "TRUE"))
  }
  
  return(datafr)
}

ReduceWeatherCoarseTerms <- function(vec){
  # Reduce the ~64 distinct combinations of weather terms to a set of 11 weather 
  # terms. The most severe weather description is used for reduction.
  #
  # Args:
  #   vec: A vector with character values (ie. weather_desc)
  #
  # Return:
  #   A vector of reduced weather descriptions corresponding to the indeces of 
  #   vec.
  weather_reduced <- rep(NA, length(vec))
  
  for (i in 1:length(vec)) {
    if (is.na(vec[i])) {
      weather_reduced[i] <- "unknown"
    } else if (grepl("Hail|Pellet", vec[i])) {
      # "Pellets" catches "Snow Pellets", "Ice Pellets" and "Ice Pellet Showers"
      weather_reduced[i] <- "hail"
    } else if (grepl("Freezing", vec[i])) {
      # "Freezing" to catch "Freezing Rain", "Freezing Drizzle", 
      weather_reduced[i] <- "freezing_rain"
    } else if (grepl("Heavy Snow|Moderate Snow", vec[i])) {
      # Heavy snow might keep an Ontarian in their house.
      weather_reduced[i] <- "heavy_snow"
    } else if (grepl("Heavy Rain", vec[i])) {
      # "Heavy Rain" to catch "Heavy Rain", "Heavy Rain Showers"
      weather_reduced[i] <- "heavy_rain"
    } else if (grepl("Moderate Rain|Rain", vec[i])) {
      # "Moderate Rain" to catch "Moderate Rain" and "Moderate Rain Showers"
      weather_reduced[i] <- "moderate_rain"
    } else if (grepl("Snow", vec[i])) {
      # "Snow" to catch "Snow Showers" and "Snow"
      weather_reduced[i] <- "snow"
    } else if (grepl("Thunderstorm", vec[i])) {
      # Thunderstorms
      weather_reduced[i] <- "tstorms"
    } else if (grepl("Overcast|Fog|Haze|Drizzle", vec[i])) {
      # "Rain|Drizzle" to catch "Rain", "Rain Showers", and "Drizzle"
      weather_reduced[i] <- "fog_drizzle"
    } else if (grepl("Cloudy", vec[i])) {
      # Fog|Haze would not change electricity use drastically compared to 
      # Cloudy, so I'm combinging them.
      weather_reduced[i] <- "cloudy"
    } else if (grepl("Clear", vec[i])) {
      weather_reduced[i] <- "clear"
    } else {
      weather_reduced[i] <- "unknown"
    }
  }
  
  return(weather_reduced)
}

ReduceSevereWeather <- function(vec){
  # Reduce the ~64 distinct combinations of weather terms to a boolean 
  # indicating whether there is severe weather in the area or not.
  #
  # Args:
  #   vec: A vector with character values (ie. weather_desc)
  #
  # Return:
  #   A vector of reduced weather descriptions corresponding to the indeces of 
  #   vec.
  severe_weather <- rep(NA, length(vec))
  
  for (i in 1:length(vec)) {
    if (is.na(vec[i])) {
      severe_weather[i] <- "unknown"
    } else if (grepl("Hail|Pellet|Freezing|Heavy Snow|Moderate Snow|Heavy Rain|Thunderstorms", vec[i])) {
      # "Hail" catches "Moderate Hail" and "Hail"
      # "Pellets" to catch "Snow Pellets", "Ice Pellets", and "Ice Pellet Showers"
      # "Freezing" to catch "Freezing Rain", "Freezing Drizzle"
      # "Heavy Snow" only occurs once in the dataset but seems severe
      # "Moderate Snow" might be considered severe enough to affect an Ontarian's behavior
      # "Heavy Rain" to catch "Heavy Rain", "Heavy Rain Showers"
      # "Thunderstorms" might be considered severe enough to affect an Ontarian's behavior
      severe_weather[i] <- TRUE
    } else {
      severe_weather[i] <- FALSE
    }
  }
  
  return(severe_weather)
}

ReduceWeatherFourTerms <- function(vec){
  # Reduce the ~64 distinct combinations of weather terms to a set of 4 weather 
  # terms. The most severe weather description is used for reduction.
  #
  # Args:
  #   vec: A vector with character values (ie. weather_desc)
  #
  # Return:
  #   A vector of reduced weather descriptions corresponding to the indeces of 
  #   vec.
  weather_reduced <- rep(NA, length(vec))
  
  for (i in 1:length(vec)) {
    if (is.na(vec[i])) {
      weather_reduced[i] <- "unknown"
    } else if (grepl("Snow|Ice|Freezing|Hail", vec[i])) {
      # Ice is under-represented in explanatory variables, even when it 
      # represents Ice|Freezing|Hail, so I'm combining it with Snow.
      weather_reduced[i] <- "snow_ice"
    } else if (grepl("Rain|Drizzle|Thunderstorm", vec[i])) {
      # Rain|Drizzle would not change electricity use drastically compared to 
      # Thunderstorms, so I'm combining them. Also, thunderstorms were arguably 
      # under-represented in explanatory variables as its own term.
      weather_reduced[i] <- "rain_tstorms"
    } else if (grepl("Cloudy|Overcast|Fog|Haze", vec[i])) {
      # Fog|Haze would not change electricity use drastically compared to 
      # Cloudy, so I'm combinging them.
      weather_reduced[i] <- "cloudy_fog"
    } else if (grepl("Clear", vec[i])) {
      weather_reduced[i] <- "clear"
    } else {
      weather_reduced[i] <- "unknown"
    }
  }
  
  return(weather_reduced)
}



WindChill <- function(temp, wind) {
  # Function computes the heat index as defined by Environment Canada. Code 
  # translated from the Javascript which backs their wind chill calculator:
  # http://www.ec.gc.ca/meteo-weather/default.js
  # 
  # Formulation is also easily seen at:
  # https://en.wikipedia.org/wiki/Wind_chill#North_American_and_United_Kingdom_wind_chill_index
  # 
  # Args:
  #   temp: Temperature value in Celsius
  #   wind: Wind speed in kph
  #
  # Return:
  #   The wind chill for provided values (ie. apparent temperature)
  
  # Windchill is defined only for temperatures at or below 10C and wind speeds 
  # above 4.8 kilometres per hour.
  if (temp <= 10 & wind > 4.8) {
    # Constants for wind chill equation
    c1 <- 13.12
    c2 <- 0.6215
    c3 <- 11.37
    c4 <- 0.3965
    
    # Wind chill equation for celsius values
    wind_chill <- c1 + (c2*temp) - (c3*(wind^0.16)) + (c4*temp*(wind^0.16));
    
    return(wind_chill)
  } else {
    return(temp)
  }
}



WindChillDiff <- function(datafr) {
  # Transform wind chill to amount under a windchill threshold (ie. highest 
  # recorded wind chill is -1) using values reported by Environment Canada in 
  # the original data.frame.
  #
  # Args:
  #   datafr: A smart meter readings data.frame
  #
  # Return:
  #   A vector representing the amount that the wind chill is under a maximum
  #   threshold.
  wind_chill_threshold <- max(datafr$wind_chill, na.rm = TRUE) + 1
  wind_chill_diff <- datafr$wind_chill
  wind_chill_diff[is.na(wind_chill_diff)] <- wind_chill_threshold
  wind_chill_diff <- wind_chill_diff - wind_chill_threshold
  return(wind_chill_diff)
}