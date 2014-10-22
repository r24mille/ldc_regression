addInferredInformation <- function(df) {
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
  df$weekend <- ifelse(df$tou_period %in% c("off_weekend"), "Yes", "No")
  
  # Add column which converts TOU Period to its price level
  df$price <- readings.aggregate$tou_period
  df$price <- gsub("off_weekend", "off_peak", df$price)
  df$price <- gsub("off_morning", "off_peak", df$price)
  df$price <- gsub("off_evening", "off_peak", df$price)
  df$price <- gsub("mid_morning", "mid_peak", df$price)
  df$price <- gsub("mid_evening", "mid_peak", df$price)
  # If TOU billing isn't active, then price is flat rate
  df$price <- ifelse(df$billing_active %in% c("No"), "flat", df$price)
  
  return(df)
}

fixDataTypes <- function(df) {
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

initAggregateReadings <- function(df) {
  # Takes the data.frame of the parsed aggregate readings CSV and initializes 
  # datatypes, inferred information, factor levels, and tweaks a few terms 
  # to a feasible subset of interactions.
  #
  # Args:
  #   df: The initial aggregate smart meter readings from parsed CSV.
  #
  # Return:
  #   The dataframe initialized and reading for modelling.
  df <- fixDataTypes(df)
  df <- addInferredInformation(df)
  df <- orderFactors(df)
  df <- stripNonexistantInteractions(df)
  
  return(df)
}

orderFactors <- function(df) {
  # Reorder factors so that model contrasts and graphs make a bit more sense
  # 
  # Args: 
  #   df: A smart meter reading data.frame that has been run through
  #       fixDataTypes(...) and addInferredInformation(...).
  # 
  # Return: 
  #   The dataframe with corrected/ordered factors.
  df$tou_period <- factor(df$tou_period, 
                          c("off_weekend", "off_morning", "mid_morning", 
                            "on_peak", "mid_evening", "off_evening"))
  df$month <- factor(df$month, c("m5", "m6", "m7", "m8", "m9", "m10"))
  df$hrstr <- factor(df$hrstr, 
                     c("h0", "h1", "h2", "h3", "h4", "h5", "h6", "h7", "h8", 
                       "h9", "h10", "h11", "h12", "h13", "h14", "h15", "h16", 
                       "h17", "h18", "h19", "h20", "h21", "h22", "h23"))
  df$weekend <- factor(df$weekend, c("No", "Yes"))
  df$price <- factor(df$price, c("flat", "off_peak", "mid_peak", "on_peak"))
  
  return(df)
}

stripNonexistantInteractions <- function(df) {
  # Manually create a column representing the two-way interaction between 
  # several terms, with infeasible interaction between levels removed
  # (eg. hrstrh23:priceon_peak or weekendYes:pricemid_peak).
  #
  # Args:
  #   df: A smart meter reading data.frame
  #
  # Returns:
  #   The data.frame with several columns added which represent certain two-way 
  #   interactions. These special columns should not be considered as main 
  #   effects nor interacted with any other terms. They _are_ the interaction.
  
  # Some htrstr and price combinations are infeasible
  df$hrstr_price <- paste0(df$hrstr, df$price)
  df$hrstr_price <- as.factor(df$hrstr_price)
  
  # Some weekend and price combinations are infeasible
  df$wknd_price <- paste0(df$weekend, df$price)
  df$wknd_price <- as.factor(df$wknd_price)
  
  return (df)
}