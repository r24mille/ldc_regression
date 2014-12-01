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
  
  return(df)
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
  
  return(df)
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
                                           "agg_count", 
                                           "temp_over_break",
                                           "temp_under_break")]
  return(df.trimmed)
}

NumericFactorCodedMatrix <- function(df) {
  # Transforms an input data.frame to a matrix with factor levels represented 
  # as numeric. Returned value is suitable for ?glinternet.
  #
  # Return:
  # The data.frame re-coded as a numeric matrix.
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

NumberFactorLevels <- function(df) {
  # Returns:
  # A vector in which each value represents the number of factor levels for 
  # each column of df
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