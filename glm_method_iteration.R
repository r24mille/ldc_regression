source('goodness_fit_visualization.R')
source('pseudo_rsquared.R')

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
  allcolnames <- colnames(df.trimmed)
  noncdh.colnames <- allcolnames[! allcolnames %in% c(cdhlagmat.colnames, 
                                               "kwh",
                                               "weights",
                                               "wghts")]
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

MaximumTractableInteractionsGlmModel <- function(df.readings, nlags, wghts) {
  # Certain two-way interactions can never interract all levels of each 
  # categorical explanatory variable (eg. weekend="Yes":price="on_peak"). So, 
  # the interacted term must be broken into its own column and then the 
  # formula updated to reflect the fact that
  #
  # Args:
  #   df.trimmed: A trimmed version of the readings.aggregate dataframe with 
  #               weights as a column.
  #
  # Returns:
  #   A formula object, where each independent variable is modeled as a fixed 
  #   effect with main effects and two-way interactions. This function 
  #   knows not to put the response variable or weights into the formula.
  cdhlag <- CreateCdhLagMatrix(nlags, df.readings)
  readings.with.cdh <- cbind(df.readings, cdhlag, wghts)
  df.trimmed <- TrimColsTouTimeComponents(readings.with.cdh)
  
  explvars <- colnames(df.trimmed[, ! colnames(df.trimmed) %in% c("kwh",
                                                                  "weights",
                                                                  "wghts")])
  explvars.fmlastr <- paste(explvars, collapse = " + ")
  maximal.fmlastr <- paste0("kwh ~ (", explvars.fmlastr, ")^2")
  tractable.fmlastr <- maximal.fmlastr
  
  # Create a column that represents only the feasible interactions of hrstr and
  # price columns. Then remove the hrstr:price interaction from the formula
  # string.
  # TODO(r24mille): The predict(...) and other functions all seem to baulk 
  #                 at the collinearity caused by the hrst:price term. Remove 
  #                 it from the model for now unless a solution can be found.
  #df.trimmed$hrstr_price <- paste0(df.trimmed$hrstr, df.trimmed$price)
  tractable.fmlastr <- paste(tractable.fmlastr, 
                            "- hrstr:price")
  
  # Create a column that represents only the feasible interaction of weekend and 
  # price columns. Then remove the weekend:price interaction from the formula 
  # string.
  # 
  # TODO(r24mille): Based on model simplification, weekend:price never becomes
  #                 linearly seprable from the main effects of price.
  #                   * Read on linear separability and what can lead to 
  #                     singularities in GLM model fitting.
  #df.trimmed$wknd_price <- paste0(df.trimmed$weekend, df.trimmed$price)
  tractable.fmlastr <- paste(tractable.fmlastr, 
                             "- weekend:price")

  
  # Fit the GLM model
  tractable.glm <- glm(formula = formula(tractable.fmlastr), 
                  data = df.trimmed, 
                  weights = wghts, 
                  family = Gamma(link="log"))
  return(tractable.glm)
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
  #   A [(nlags+1) x 4] array that stores residual deviance, AICc, BIC, and 
  #   pseudo R^2 values. Each column is a metric. Each row is an iteration.
  resmatrix <- matrix(nrow = (nlags + 1),
                      ncol = 4,
                      dimnames = list(c(0:nlags),
                                      c("ResidualDeviance", 
                                        "AICc",
                                        "BIC",
                                        "McFaddenPseudoR2")))
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
  #   A 1x4 matrix representing a row of explanatory power information for the 
  #   provided model including residual deviance, AICc, and BIC.
  pwr.row <- matrix(nrow = 1, ncol = 4)
  pwr.row[1, 1] <- iter.glm$deviance
  pwr.row[1, 2] <- AICc(iter.glm)
  pwr.row[1, 3] <- BIC(iter.glm)
  pwr.row[1, 4] <- McFaddenPseudoR2(iter.glm)
  return(pwr.row)
}

IterativeGlmPlots <- function(iter.glm.pwr, maintitle) {
  # Plot both the residual deviation and pseudo r-squared versions against two 
  # information criteria to convey the explanatory power of each GLM model 
  # iteration.
  #
  # Args:
  #   iter.glm.pwr: An [n x 4] matrix created by the GlmPowerResultsMatrix 
  #                 function.
  #   maintitle: The "main" title passed along to the plot(...) function.
  IterativeGlmPlotWithResidualDeviance(iter.glm.pwr, maintitle)
  IterativeGlmPlotWithPseudoRsquared(iter.glm.pwr, maintitle)
}

IterativeGlmPlotWithPseudoRsquared <- function(iter.glm.pwr, maintitle) {
  # Wraps the PlotGlmFitMeasures function since each iteration calls it in 
  # pretty much the same way. Encapsulating this call simplifies setting 
  # the many parameters required by the function.
  #
  # The y2-axis is McFadden pseudo r-squared values.
  #
  # Args:
  #   iter.glm.pwr: An [n x 4] matrix created by the GlmPowerResultsMatrix 
  #                 function.
  #   maintitle: The "main" title passed along to the plot(...) function.
  PlotGlmFitMeasures(aiccs = iter.glm.pwr[, 2], 
                     bics = iter.glm.pwr[, 3], 
                     y2vals = iter.glm.pwr[, 4], 
                     xvals = c(0:(nrow(iter.glm.pwr) - 1)), 
                     y2title = expression(paste("McFadden Pseudo ", R^2)), 
                     xtitle = "Number of Past Hours Included", 
                     title = maintitle)
}

IterativeGlmPlotWithResidualDeviance <- function(iter.glm.pwr, maintitle) {
  # Wraps the PlotGlmFitMeasures function since each iteration calls it in 
  # pretty much the same way. Encapsulating this call simplifies setting 
  # the many parameters required by the function.
  #
  # The y-axis is residual deviance values.
  #
  # Args:
  #   iter.glm.pwr: An [n x 4] matrix created by the GlmPowerResultsMatrix 
  #                 function.
  #   maintitle: The "main" title passed along to the plot(...) function.
  PlotGlmFitMeasures(aiccs = iter.glm.pwr[, 2], 
                     bics = iter.glm.pwr[, 3], 
                     y2vals = iter.glm.pwr[, 1], 
                     xvals = c(0:(nrow(iter.glm.pwr) - 1)), 
                     y2title = "Residual Deviance", 
                     xtitle = "Number of Past Hours Included", 
                     title = maintitle)
}

PerformTouCdhGlmIterations <- function(df.readings, weights, nhrs) {
  # Steps through 0-nhrs iterations of GLMs trying varioud combinations of 
  # structuring TOU and CDH explanatory variables. The explanatory power of each
  # model is plotted.
  #
  # Args:
  #   df.readings: A readings dataframe to be passed as the "data" parameter 
  #                to the GLMs.
  #   weights: Weights to be passed to the GLMs
  #   nhrs: Number of hours to include in the CDH component of the model
  cdhlagsum.touperiods.maxglm.pwr <- GlmPowerResultsMatrix(nhrs)
  cdhlagmat.touperiods.nestedglm.pwr <- GlmPowerResultsMatrix(nhrs)
  cdhlagmat.touperiods.maxglm.pwr <- GlmPowerResultsMatrix(nhrs)
  cdhlagsum.toucomps.maxglm.pwr <- GlmPowerResultsMatrix(nhrs)
  cdhlagmat.toucomps.nestedglm.pwr <- GlmPowerResultsMatrix(nhrs)
  cdhlagmat.toucomps.maxglm.pwr <- GlmPowerResultsMatrix(nhrs)
  
  for(i in 0:nhrs) {
    # 1. Summed CDH lags, TOU as periods
    cdhlagsum.touperiods.maxglm <- IterativeGlmModel(df.readings = df.readings, 
                                                     wghts = weights, 
                                                     nlags = i, 
                                                     is.touperiod = TRUE, 
                                                     is.cdhlagsum = TRUE, 
                                                     is.maxformula = TRUE)
    cdhlagsum.touperiods.maxglm.pwr[(i+1),] <- IterativeGlmPower(cdhlagsum.touperiods.maxglm)
    
    # 2. Matrix of CDH lags as nested interactions, TOU as periods
    cdhlagmat.touperiods.nestedglm <- IterativeGlmModel(df.readings = df.readings, 
                                                        wghts = weights, 
                                                        nlags = i, 
                                                        is.touperiod = TRUE, 
                                                        is.cdhlagsum = FALSE, 
                                                        is.maxformula = FALSE)
    cdhlagmat.touperiods.nestedglm.pwr[(i+1),] <- IterativeGlmPower(cdhlagmat.touperiods.nestedglm)
    
    # 3. Matrix of CDH lags as two-way interactions, TOU as periods
    cdhlagmat.touperiods.maxglm <- IterativeGlmModel(df.readings = df.readings, 
                                                     wghts = weights, 
                                                     nlags = i, 
                                                     is.touperiod = TRUE, 
                                                     is.cdhlagsum = FALSE, 
                                                     is.maxformula = TRUE)
    cdhlagmat.touperiods.maxglm.pwr[(i+1),] <- IterativeGlmPower(cdhlagmat.touperiods.maxglm)
    
    # 4. Summed CDH lags, TOU components of time
    cdhlagsum.toucomps.maxglm <- IterativeGlmModel(df.readings = df.readings, 
                                                   wghts = weights, 
                                                   nlags = i, 
                                                   is.touperiod = FALSE, 
                                                   is.cdhlagsum = TRUE, 
                                                   is.maxformula = TRUE)
    cdhlagsum.toucomps.maxglm.pwr[(i+1),] <- IterativeGlmPower(cdhlagsum.toucomps.maxglm)
    
    # 5. Matrix of CDH lags as nested interactions, TOU components of time
    cdhlagmat.toucomps.nestedglm <- IterativeGlmModel(df.readings = df.readings, 
                                                      wghts = weights, 
                                                      nlags = i, 
                                                      is.touperiod = FALSE, 
                                                      is.cdhlagsum = FALSE, 
                                                      is.maxformula = FALSE)
    cdhlagmat.toucomps.nestedglm.pwr[(i+1),] <- IterativeGlmPower(cdhlagmat.toucomps.nestedglm)
    
    # 6. Matrix of CDH lags as two-way interactions, TOU components of time
    cdhlagmat.toucomps.maxglm <- IterativeGlmModel(df.readings = df.readings, 
                                                   wghts = weights, 
                                                   nlags = i, 
                                                   is.touperiod = FALSE, 
                                                   is.cdhlagsum = FALSE, 
                                                   is.maxformula = TRUE)
    cdhlagmat.toucomps.maxglm.pwr[(i+1),] <- IterativeGlmPower(cdhlagmat.toucomps.maxglm)
  }
  
  IterativeGlmPlots(iter.glm.pwr = cdhlagsum.touperiods.maxglm.pwr, 
                    maintitle = expression(paste("TOU Period, Degrees>", 
                                                 CDH['break'], 
                                                 " Summed")))
  
  IterativeGlmPlots(iter.glm.pwr = cdhlagmat.touperiods.nestedglm.pwr, 
                    maintitle = expression(paste("TOU Period, Degrees>", 
                                                 CDH['break'], 
                                                 " as Coefficients w/ Nested Interaction")))
  
  IterativeGlmPlots(iter.glm.pwr = cdhlagmat.touperiods.maxglm.pwr, 
                    maintitle = expression(paste("TOU Period, Degrees>", 
                                                 CDH['break'], 
                                                 " as Coefficients w/ All 2-Way Interactions")))
  
  IterativeGlmPlots(iter.glm.pwr = cdhlagsum.toucomps.maxglm.pwr, 
                    maintitle = expression(paste("TOU Components, Degrees>", 
                                                 CDH['break'], 
                                                 " Summed")))
  
  IterativeGlmPlots(iter.glm.pwr = cdhlagmat.toucomps.nestedglm.pwr, 
                    maintitle = expression(paste("TOU Components, Degrees>", 
                                                 CDH['break'], 
                                                 " as Coefficients w/ Nested Interaction")))
  
  IterativeGlmPlots(iter.glm.pwr = cdhlagmat.toucomps.maxglm.pwr, 
                    maintitle = expression(paste("TOU Components, Degrees>", 
                                                 CDH['break'], 
                                                 " as Coefficients w/ All 2-Way Interactions")))
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