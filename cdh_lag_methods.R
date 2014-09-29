CreateCdhLagSum <- function(nlags, readingsdf) {
  # For each observation, this function sums all temperature values above the
  # CDH breakpoint from the current time through a number of hours into the past
  # (nlags).
  # 
  # Args:
  #   nlags: The number of hours into the past (ie. lags) to sum for each 
  #          observation.
  #   readingsdf: A dataframe of smart meter readings.
  #
  # Returns:
  #   A vector of values based on the "conventional" definition of CDH.  
  cdhlag <- rep(0, nrow(readingsdf)) # Prealloate cdhlag vector
  
  # Unfortunately iterating over the dataframe seems to be the best method
  for(i in 1:nrow(readingsdf)) {
    cdhlag[i] <- readingsdf[i, "cdh"]
    tmplg <- 1
    
    # Sum what we can or all the values up to nlags
    while (i-tmplg > 0 & tmplg <= nlags) {
      cdhlag[i] <- cdhlag[i] + readingsdf[i-tmplg, "cdh"]
      tmplg <- tmplg + 1
    }
  }
  
  return(cdhlag)
}


CreateCdhLagMatrix <- function(nlags, readingsdf) {
  # Creates a matrix of values such that each column looks an additional hour 
  # into the past.
  #
  # Args:
  #   nlags: The number of hours (ie. lags) to include in the matrix
  #   readingsdf: A dataframe of smart meter readings.
  #
  # Returns:
  #   An [n x (nlags + 1)] matrix such that each hour into the past can have its
  #   own coefficient when modelled.
  lagnames <- paste0("lag", c(0:nlags))
  cdhlags <- matrix(nrow = nrow(readingsdf),
                    ncol = (nlags + 1))
  colnames(cdhlags) <- lagnames
  
  # Unfortunately iterating over the dataframe seems to be the best method
  for(i in 1:nrow(readingsdf)) {
    for(j in 0:nlags) {
      if (i-j > 0) {
        cdhlags[i, (j+1)] <- readingsdf[i-j, "cdh"]
      } else {
        cdhlags[i, (j+1)] <- 0
      }
    }
  }
  
  return(cdhlags)
}

CdhLagMaximalFormula <- function() {
  # Summed CDH lags can be incorporated into a GLM easily, testing the main 
  # effects of the column as well as two-way interractions with each other 
  # column. This function returns the forumla for summed
  #
  # Returns:
  #   A formula object, where each independent variable is modeled as a fixed 
  #   effect with main effects and two-way interactions.
  return(formula("kwh ~ (.)^2"))
}

CdhLagMaximalNestedFormula <- function(cdhlagmat.minimaldf, 
                                       cdhlagmat.colnames) {
  # Each CDH lag into the past may be provided as its own column in the 
  # dataframe and modelled as nested fixed effects and two-way interractions 
  # between all fixed effects.
  #
  # Args:
  #   cdhlagmat.minimaldf: A smart meter readings dataframe with CDH lags of 
  #                        previous hours are each a column of data (eg. "lag0",
  #                        "lag1", ..., "lagn").
  #   cdhlagmat.colnames: The ordered column names of previous hours' degrees 
  #                       over the CDH breakpoint.
  # 
  # Returns:
  #   A formula object, where the lagged obervation terms are nested, and 
  #   two-way interracted with other main effects.
  
  # TODO(r24mille): Optimize this so that an entire subset isn't created just 
  #                 to get me.colnames.
  me.colnames <- colnames(cdhlagmat.minimaldf[, ! colnames(cdhlagmat.minimaldf) %in% c(cdhlagmat.colnames, 
                                                                                       "kwh")]) 
  me.fmlastr <- paste(me.colnames,
                      collapse = " + ")
  nestedcdhlag.fmlastr <- paste(cdhlagmat.colnames,
                     collapse = "/")
  maximal.fmlastr <- paste0("kwh ~ (",
                            me.fmlastr,
                            " + ",
                            nestedcdhlag.fmlastr,  
                            ")^2")
  return(formula(maximal.fmlastr))
}