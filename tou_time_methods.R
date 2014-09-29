TrimColsTouPeriods <- function(readingsdf) {
  # Trim the columns down from full readings.aggregate to only those needed 
  # when working with "TOU Period" version of the model.
  #
  # Args:
  #   readingsdf: A dataframe object of aggregate readings
  #
  # Returns:
  #   A trimmed dataframe in which unneeded columns from readings.aggregate 
  #   dataframe have been trimmed out. Only column with data relevant to 
  #   representing TOU via its periods (ie. tou_period) are preserved. 
  #   Ontario's TOU periods are essentially a proxy for day, time, demand, and 
  #   price.  
  periodsdf <- readingsdf[, ! colnames(readingsdf) %in% c("daynum", 
                                                          "hour", 
                                                          "timestamp_dst", 
                                                          "temperature", 
                                                          "agg_count", 
                                                          "hrstr", 
                                                          "weekend", 
                                                          "price", 
                                                          "cdh")]
  return(periodsdf)
}


TrimColsTouTimeComponents <- function(readingsdf) {
  # Trim the columns down from full readings.aggregate to only those needed 
  # when working with "Components of TOU" version of the model.
  #
  # Args:
  #   readingdf: A dataframe object of aggregate readings
  #
  # Returns:
  #   A trimmed dataframe in which unneeded columns from readings.aggregate 
  #   dataframe have been trimmed out. Only columns with data relevant to 
  #   representing TOU via all of its components (ie. hour-of-day, 
  #   weekday/weekend, and price) are preserved.
  timecompdf <- readingsdf[, ! colnames(readingsdf) %in% c("daynum", 
                                                           "hour", 
                                                           "timestamp_dst", 
                                                           "temperature", 
                                                           "agg_count", 
                                                           "tou_period", 
                                                           "cdh")]
  return(timecompdf)
}