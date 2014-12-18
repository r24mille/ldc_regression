AddInferredInformation <- function(df) {
  # Additional information can be inferred from existing dataframe column 
  # and placed into new columns that may be useful for modelling.
  #
  # Args:
  #   df: A smart meter reading data.frame
  #
  # Return:
  #   The dataframe with several new columns added, with information inferred 
  #   from columns that were originally in the CSV.
  # Add column which represents "month" as a categorical factor
  df$month <- paste0("m", (df$timestamp_dst$mon + 1))
  
  # Represent hours as levels rather than integers
  df$hrstr <- paste0("h", df$hour)
  
  # Add yes/no weekend flag
  df$weekend <- ifelse(df$timestamp_dst$wday == 0 | df$timestamp_dst$wday == 6,
                       "Yes",
                       "No")  
  
  return(df)
}

CreatePastTemperatureMatrix <- function(nlags, df.readings) {
  # Creates a matrix of values such that each column looks an additional hour 
  # into the past at degrees over the temperature breakpoint.
  #
  # Args:
  #   nlags: The number of hours (ie. lags) to include in the matrix
  #   df.readings: A dataframe of smart meter readings.
  #
  # Returns:
  #   An [n x (nlags + 1)] matrix such that each hour into the past can have its
  #   own coefficient when modelled.
  lagnames <- c(paste0("overbreak_lag", c(0:nlags)),
                paste0("underbreak_lag", c(0:nlags)))
  templags <- matrix(nrow = nrow(df.readings),
                     ncol = (2 * (nlags + 1)))
  colnames(templags) <- lagnames
  
  # Unfortunately iterating over the dataframe seems to be the best method
  for(i in 1:nrow(df.readings)) {
    for(j in 0:nlags) {
      if (i-j > 0) {
        templags[i, (j+1)] <- df.readings[i-j, "temp_over_break"]
        templags[i, (j+1+nlags+1)] <- df.readings[i-j, "temp_under_break"]
      } else {
        templags[i, (j+1)] <- 0
        templags[i, (j+1+nlags+1)] <- 0
      }
    }
  }
  
  return(templags)
}

FeelsLike <- function(df) {
  # Converts dry-bulb temperatures to "feels like" temperatures that reflect 
  # what the apparent temperature would feel like to a person.
  #
  # Args:
  #   df: A data.frame with dry bulb temperature in Celsius, relative humidity 
  #       as a percentage, and wind speed in kph.
  #
  # Return:
  #   A vector of apparent temperature values.
  
  feels_like <- rep(NA, nrow(df))
  
  for (i in 1:nrow(df)){
    # Only calculate heat index and wind chill at valid values (as defined by 
    # NOAA or Environment Canada), otherwise the apparent temperature is the 
    # same as the dry bulb temperature.
    if (df$temperature[i] > 27 & df$rel_humidity_pct[i] > 40) {
      feels_like[i] <- HeatIndex(df$temperature[i], df$rel_humidity_pct[i])
    } else if (df$temperature[i] <= 10 & df$wind_speed_kph[i] > 4.8) {
      feels_like[i] <- WindChill(df$temperature[i], df$wind_speed_kph[i])
    } else {
      feels_like[i] <- df$temperature[i]
    }
  }
  
  return(feels_like)
}

FixDataTypes <- function(df) {
  # Takes a data frame of smart meter reading data and fixes the column data types
  #
  # Args:
  #   df: A smart meter reading data.frame
  #
  # Return:
  #   The dataframe with corrected column types.
  
  # Coerses timestamp_dst to the POSIX datetime type
  df$timestamp_dst <- as.POSIXlt(df$timestamp_dst)
  
  # Enforce numeric column type on humidex and wind_chill
  df$humidex <- as.numeric(df$humidex)
  df$wind_chill <- as.numeric(df$wind_chill)
  
  return(df)
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

InitAggregateReadings <- function(df) {
  # Takes the data.frame of the parsed aggregate readings CSV and initializes 
  # datatypes, inferred information, factor levels, and tweaks a few terms 
  # to a feasible subset of interactions.
  #
  # Args:
  #   df: The initial aggregate smart meter readings from parsed CSV.
  #
  # Return:
  #   The dataframe initialized and reading for modelling.
  df <- FixDataTypes(df)
  df <- AddInferredInformation(df)
  df <- OrderFactors(df)
  
  return(df)
}

NumberFactorLevels <- function(df) {
  # Args: 
  #   df: A smart meter reading data.frame that has been trimmed and is ready 
  #       for conversion to a matrix processed by ?glinternet.
  #
  # Return:
  #   A vector in which each value represents the number of factor levels for 
  #   each column of the data.frame passed in.
  numlevels <- rep(0, ncol(df))
  
  for (i in 1:ncol(df)) {
    if (is.factor(df[,i])) {
      numlevels[i] <- nlevels(df[,i])
    } else {
      numlevels[i] <- 1
    }
  }
  
  return(numlevels)
}

NumericFactorCodedMatrix <- function(df) {
  # Transforms an input data.frame to a matrix with factor levels represented 
  # as numeric. Returned value is suitable for ?glinternet.
  #
  # Args: 
  #   df: A smart meter reading data.frame that has been trimmed and is ready 
  #       for conversion to a matrix processed by ?glinternet.
  #
  # Return:
  #   The data.frame re-coded as a numeric matrix.
  numeric.m <- matrix(0, nrow = nrow(df), ncol = ncol(df))
  colnames(numeric.m) <- names(df)
  
  for (i in 1:ncol(df)) {
    if (is.factor(df[,i])) {
      numeric.m[,i] <- as.numeric(df[,i]) - 1
    } else {
      numeric.m[,i] <- df[,i]
    }
  }
  
  return(numeric.m)
}

OrderFactors <- function(df) {
  # Reorder factors so that model contrasts and graphs make a bit more sense
  # 
  # Args: 
  #   df: A smart meter reading data.frame that has been run through
  #       fixDataTypes(...) and addInferredInformation(...).
  # 
  # Return: 
  #   The dataframe with corrected/ordered factors.
  #df$tou_period <- factor(df$tou_period, 
  #                        c("off_weekend", "off_morning", "mid_morning", 
  #                          "on_peak", "mid_evening", "off_evening"))
  df$dayname <- factor(df$dayname, c("Sun", "Mon", "Tue", "Wed", "Thu", "Fri", 
                                     "Sat"))
  df$month <- factor(df$month, c("m1", "m2", "m3", "m4", "m5", "m6", "m7", 
                                 "m8", "m9", "m10", "m11", "m12"))
  df$hrstr <- factor(df$hrstr, 
                     c("h0", "h1", "h2", "h3", "h4", "h5", "h6", "h7", "h8", 
                       "h9", "h10", "h11", "h12", "h13", "h14", "h15", "h16", 
                       "h17", "h18", "h19", "h20", "h21", "h22", "h23"))
  df$weekend <- factor(df$weekend, c("No", "Yes"))
  df$price <- factor(df$price, c("flat", "off_peak", "mid_peak", "on_peak"))
  df$weather_desc <- factor(df$weather_desc, 
                            c("clear", "fog", "cloudy", "rain", "snow", "ice", 
                              "thunderstorms"))
  
  return(df)
}

ReduceWeather <- function(df){
  # Reduce the ~64 distinct combinations of weather terms to a set of 7 weather 
  # terms. The most severe weather description is used for reduction.
  #
  # Args:
  #   df: A dataframe with ef$weather_desc column
  #
  # Return:
  #   A vector of reduced weather descriptions corresponding to the indeces of 
  #   df$weather_desc.
  weather_reduced <- rep(NA, nrow(df))
  
  for (i in 1:nrow(df)) {
    if (grepl("Thunderstorms", df$weather_desc[i])) {
      weather_reduced[i] <- "thunderstorms"
    } else if (grepl("Ice|Freezing|Hail", df$weather_desc[i])) {
      weather_reduced[i] <- "ice"
    } else if (grepl("Snow", df$weather_desc[i])) {
      weather_reduced[i] <- "snow"
    } else if (grepl("Rain|Drizzle", df$weather_desc[i])) {
      weather_reduced[i] <- "rain"
    } else if (grepl("Cloudy", df$weather_desc[i])) {
      weather_reduced[i] <- "cloudy"
    } else if (grepl("Fog|Haze", df$weather_desc[i])) {
      weather_reduced[i] <- "fog"
    } else {
      weather_reduced[i] <- "clear"
    } 
  }
  
  return(weather_reduced)
}

TrimExplanatoryVariables <- function(df) {
  # The dataframe formed from parsing the .csv file and manipulating past 
  # temperature breakpoint information has several extra columns. Trim the 
  # data.frame down to only those explanatory variables to be considered by the 
  # model.
  df.trimmed <- df[, ! colnames(df) %in% c("sample_index", 
                                           "daynum", 
                                           "hour", 
                                           "timestamp_dst", 
                                           "temperature", 
                                           "dewpoint_temp_c", 
                                           "rel_humidity_pct", 
                                           "wind_speed_kph", 
                                           "visibility_km", 
                                           "pressure_kpa",
                                           "wind_chill",
                                           "wind_chill_diff", 
                                           "agg_count", 
                                           "weekend", 
                                           "humidex",
                                           "humidex_diff",
                                           "temp_over_break",
                                           "temp_under_break",
                                           "nvgnt_thi",
                                           "nvgnt_cool_thi",
                                           "nvgnt_heat_thi")]
  return(df.trimmed)
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