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

HumidexDiff <- function(df) {
  # Transform humidex to amount over a humidex threshold (ie. lowest recorded 
  # humidex is 25) using values reported by Environment Canada.
  # http://climate.weather.gc.ca/climate_normals/normals_documentation_e.html#humidex
  #
  # Args:
  #   df: A smart meter readings data.frame
  #
  # Return:
  #   A vector representing the amount that the humidex is over a minimum 
  #   threshold.
  humidex_threshold <- min(df$humidex, na.rm = TRUE) - 1
  humidex_diff <- df$humidex
  humidex_diff[is.na(humidex_diff)] <- humidex_threshold
  humidex_diff <- humidex_diff - humidex_threshold
  return(humidex_diff)
}

InitReadingsDataFrame <- function(fpath, is.aggregate = FALSE) {
  # Initializes a data.frame for the common readings and aggregate readings 
  # column structure. The data.frame is empty but columns have been typed 
  # appropriately with ordered factors.
  readings.colnames <- c("kwh", "sample_index", "daynum", "hrstr", "month", 
                         "weekend", "timestamp_dst", "dayname", "holiday", 
                         "temperature", "dewpnt_temp", "rel_humidity_pct", 
                         "wind_speed_kph", "humidex", "wind_chill", "price")
  readings.colclasses <- c("numeric", "integer", "integer", "factor", "factor",
                           "factor", "POSIXct", "factor", "factor", 
                           "numeric", "numeric", "numeric",  
                           "numeric", "numeric", "numeric", "factor")  
  
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
  df <- read.csv(file = fpath, 
                 header = has.headers, 
                 col.names = readings.colnames,
                 na.strings = c("NULL", "NA", "NaN", "\\N"),
                 colClasses = readings.colclasses, 
                 stringsAsFactors = FALSE)
  
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
  if("meterid" %in% colnames(df)) {
    df$meterid <- factor(df$meterid)
  }
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
  df$weather_reduced <- factor(df$weather_reduced, 
                            c("clear", "cloudy_fog", "rain_tstorms", 
                              "snow_ice", "unknown"))
  
  return(df)
}

ReduceWeather <- function(vec){
  # Reduce the ~64 distinct combinations of weather terms to a set of 7 weather 
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
      # Rain|Drizzle would not change electricity use drrastically compared to 
      # Thunderstorms, so I'm combining them. Also, thunderstorms were arguably 
      # under-represented in explanatory variables as its own term.
      weather_reduced[i] <- "rain_tstorms"
    } else if (grepl("Cloudy|Overcast|Fog|Haze", vec[i])) {
      # Fog|Haze would not change electricity use drastically compared to 
      # Cloudy, so I'm combinging them.
      weather_reduced[i] <- "cloudy_fog"
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

WindChillDiff <- function(df) {
  # Transform wind chill to amount under a windchill threshold (ie. highest 
  # recorded wind chill is -1) using values reported by Environment Canada in 
  # the original data.frame.
  #
  # Args:
  #   df: A smart meter readings data.frame
  #
  # Return:
  #   A vector representing the amount that the wind chill is under a maximum
  #   threshold.
  wind_chill_threshold <- max(df$wind_chill, na.rm = TRUE) + 1
  wind_chill_diff <- df$wind_chill
  wind_chill_diff[is.na(wind_chill_diff)] <- wind_chill_threshold
  wind_chill_diff <- wind_chill_diff - wind_chill_threshold
  return(wind_chill_diff)
}