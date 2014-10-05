CdhLagMaximalFormula <- function(df.trimmed) {
  # This function returns the forumla for all main effects and two-way 
  # interactions of the columns in the trimmed dataframe provided.
  #
  # Args:
  #   trimmed.df: A trimmed version of the readings.aggregate dataframe.
  #
  # Returns:
  #   A formula object, where each independent variable is modeled as a fixed 
  #   effect with main effects and two-way interactions. This function 
  #   knows not to put the response variable or weights into the formula.
  explvars <- colnames(df.trimmed[, ! colnames(df.trimmed) %in% c("kwh",
                                                                  "weights",
                                                                  "wghts")]) 
  explvars.fmlastr <- paste(explvars, collapse = " + ")
  maximal.fmlastr <- paste0("kwh ~ (", explvars.fmlastr, ")^2")
  return(formula(maximal.fmlastr))
}

CdhLagMaximalNestedFormula <- function(df.trimmed, cdhlagmat.colnames) {
  # Each CDH lag into the past may be provided as its own column in the 
  # dataframe and modelled as nested fixed effects and two-way interractions 
  # between all fixed effects.
  #
  # Args:
  #   df.trimmed: A smart meter readings dataframe with CDH lags of previous
  #               hours are each a column of 
  #               data (eg. "lag0", "lag1", ..., "lagn").
  #   cdhlagmat.colnames: The ordered column names of previous hours' degrees 
  #                       over the CDH breakpoint.
  # 
  # Returns:
  #   A formula object, where the lagged obervation terms are nested, and 
  #   two-way interracted with other main effects. This function knows not to
  #   put the response variable or weights into the formula.
  
  # TODO(r24mille): Optimize this so that an entire subset isn't created just 
  #                 to get me.colnames.
  noncdh.colnames <- colnames(df.trimmed[, ! colnames(df.trimmed) %in% c(cdhlagmat.colnames, 
                                                                         "kwh",
                                                                         "weights",
                                                                         "wghts")]) 
  noncdh.fmlastr <- paste(noncdh.colnames,
                          collapse = " + ")
  nestedcdhlag.fmlastr <- paste(cdhlagmat.colnames,
                                collapse = "/")
  maximal.fmlastr <- paste0("kwh ~ (",
                            noncdh.fmlastr,
                            " + ",
                            nestedcdhlag.fmlastr,  
                            ")^2")
  return(formula(maximal.fmlastr))
}

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

GlmPowerResultsMatrix <- function(nlags) {
  # Function to standardize the creation of a matrix that stores model 
  # descriptive power results from multiple iterations of glm fitting.
  #
  # Args:
  #   nlags: The number of hours that will be lagged (iterated). Creates a 
  #          matrix big enough to store results from 0-nlags.
  #
  # Returns:
  #   A [3 x (nlags+1)] array that stores residual deviance, AICc, and BIC values. 
  #   Each column is a metric. Each row is an iteration.
  resmatrix <- matrix(nrow = (nlags + 1),
                      ncol = 3,
                      dimnames = list(c(0:nlags),
                                      c("ResidualDeviance", 
                                        "AICc",
                                        "BIC")))
  return(resmatrix)
}

IterativeGlmModel <- function(df.readings, wghts, nlags, 
                              is.touperiod, is.cdhlagsum, is.maxformula) {
  # Time-of-Use (TOU) can be represented as a categorical explanatory variable 
  # or its component parts can each be modelled as explanatory variables. 
  # Additionally, cooling degree-hours(CDH) can be modelled using the
  # conventional "sum" method or each hour can be modelled individually as
  # explanatory variables. Furthermore, because many hours of CDH history may
  # introduce a great deal of model complexity, there is one idea to represent
  # them as nested interactions.
  #
  # Args:
  #   df.readings: A dataframe of smart meter readings (likely aggregate)
  #   wghts: A vector of observation weights passed to the GLM model.
  #   nlags: The number of hours into the past that CDH will consider.
  #   is.touperiod: A boolean indicating whether the categorical tou_period
  #                 explanatory variable is used. If FALSE, then components of 
  #                 TOU will be used.
  #   is.cdhlagsum: A boolean indicating whether the CDH should be summed a 
  #                 number of hours into the past. If FALSE, then each 
  #                 historical cdh values will be represented as their own 
  #                 explanatory variables.
  #   is.maxformula: A boolean indicating whether a maximal two-way interaction 
  #                  formula should be used. If FALSE, then a nested formula 
  #                  (ie. the nested CDH lag experiment) will be created.
  #
  # Returns:
  #   A glm fitted according to the provided criteria.
  
  # CDH structure created (ie. vector or matrix)
  if(is.cdhlagsum == TRUE) {
    cdhlag <- CreateCdhLagSum(nlags, df.readings)
  } else {
    cdhlag <- CreateCdhLagMatrix(nlags, df.readings)
  }
  
  # Stitch the two data structures together into a dataframe. Also stitch 
  # weights (ie. wghts) as a column since glm(...) can't take the vector 
  # as a parameter. It can only read columns from its data=... dataframe. The 
  # later TrimColsTou functions know to ignore the weights/wghts column.
  readings.with.cdh <- cbind(df.readings, cdhlag, wghts)
  
  # Trim dataframe columns appropriate for desired representation of TOU
  if(is.touperiod == TRUE) {
    df.trimmed <- TrimColsTouPeriods(readings.with.cdh)
  } else {
    df.trimmed <- TrimColsTouTimeComponents(readings.with.cdh)
  }
  
  # Create appropriate GLM formula
  if (is.maxformula == TRUE) {
    fmla <- CdhLagMaximalFormula(df.trimmed)
  } else {
    fmla <- CdhLagMaximalNestedFormula(df.trimmed, colnames(cdhlag))
  }
  
  # Fit the GLM model
  iter.glm <- glm(formula = fmla, 
                  data = df.trimmed, 
                  weights = wghts, 
                  family = Gamma(link="log"))
  return(iter.glm)
}

IterativeGlmPower <- function(iter.glm) {
  # Creates a row to be placed into a matrix tracking the explanatory power over
  # model iterations of a CDH lag method.
  #
  # Args:
  #   iter.glm: A GLM object for one iteration.
  #
  # Returns:
  #   A 1x3 matrix representing a row of explanatory power information for the 
  #   provided model including residual deviance, AICc, and BIC.
  pwr.row <- matrix(nrow = 1, ncol = 3)
  pwr.row[1, 1] <- iter.glm$deviance
  pwr.row[1, 2] <- AICc(iter.glm)
  pwr.row[1, 3] <- BIC(iter.glm)
  return(pwr.row)
}

IterativeGlmPlot <- function(iter.glm.pwr, maintitle) {
  # Wraps the PlotGlmFitMeasures function since each iteration calls it in 
  # pretty much the same way. Encapsulating this call simplifies setting 
  # the many parameters required by the function.
  #
  # Args:
  #   iter.glm.pwr: An [n x 3] matrix created by the GlmPowerResultsMatrix 
  #                 function.
  #   maintitle: The "main" title passed along to the plot(...) function.
  PlotGlmFitMeasures(aiccs = iter.glm.pwr[, 2], 
                     bics = iter.glm.pwr[, 3], 
                     resdevs = iter.glm.pwr[, 1], 
                     xvals = c(0:(nrow(iter.glm.pwr) - 1)), 
                     xtitle = "Number of Past Hours Included", 
                     title = maintitle)
}

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