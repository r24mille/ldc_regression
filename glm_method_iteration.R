library(MASS) # For correcting Gamma glm AIC

source('dataframe_processing.R')
source('goodness_fit_visualization.R')
source('pseudo_rsquared.R')

BackwardStepwiseLinearRegression <- function(nlags, df.obs, formula.maximal,
                                             is.stopcriteria.bic = FALSE,
                                             pval.threshold = 0.05) {
  # Performs backward stepwise linear regression including the provided number 
  # of hours of temperature history in its search space.
  #
  # TODO(r24mille): Iteratively building a dataframe is frowned upon. Consider
  #                 some other, more optimal method of building results.
  #
  # Args:
  #   nlags: The number of hours into the past that temperature > breakpoint 
  #          information that can be incorporated into the model.
  #   df.obs: The data.frame of observations to be used.
  #   formula.maximal: The formula for the maximal model to be considered.
  #   is.stopcriteria.bic: Controls the stopping criteria for the forward 
  #                        stepwise explanatory variable exploration. If BIC 
  #                        stopping criteria is used, then search stops when 
  #                        BIC stops decreasing. Default value is FALSE and 
  #                        p-value stopping criteria is used.
  #   pval.threshold: Controls the p-value stopping criteria. A likelihood ratio
  #                   test (LRT) is performed and the term is kept so long as 
  #                   it's p-value is < 0.05 (default). Other thresholds 
  #                   (eg. 0.01) may be passed in via this parameter.
  #
  # Return:
  #   A data.frame of results fitting the df.stepresults convention established 
  #   elsewhere.
  df.stepresults <- InitStepresultsDataframe()
  glm.maximal <- glm(formula = formula.maximal, 
                     data = df.obs, 
                     family = gaussian)
  
  # I know... iteratively building a data.frame is bad
  glm.maximal.pwr <- IterativeGlmPower(glm.maximal)
  explvars <- attr(glm.maximal$terms, "term.labels")
  df.stepresults <- rbind(df.stepresults,
                          data.frame(num.lags = nlags,
                                     num.explvars = length(explvars),
                                     deviance.null = glm.maximal$null.deviance,
                                     deviance.residuals = glm.maximal$deviance,
                                     AICc = glm.maximal.pwr[1,2],
                                     BIC = glm.maximal.pwr[1,3],
                                     mcfadden.r2 = glm.maximal.pwr[1,4],
                                     formulastr = Reduce(paste, 
                                                         deparse(glm.maximal$formula,
                                                                 width.cutoff = 500)),
                                     variable.removed = "",
                                     probability.gt.chi = 0))
  
  # Print status update
  print(explvars)
  
  # Flag whether stepwise deletion of variables has completed
  is.stepwise.complete = FALSE
  while (is.stepwise.complete == FALSE) {
    drop1.results <- drop1(object = glm.maximal,
                           test = "LRT", 
                           k = log(nrow(df.obs)))
    
    if (is.stopcriteria.bic == TRUE) {
      # Sort according to p-value of LRT results
      explvar.todrop <- rownames(drop1.results[ order(drop1.results[,3], 
                                                      decreasing = TRUE), ][1,])
      pr.gt.chi <- drop1.results[ order(drop1.results[,3], 
                                        decreasing = TRUE), ][1,5]
      
      # If stepwise has reached the y-intercept, then no more model construction
      # can be done.
      if (explvar.todrop == "<none>") {
        is.stepwise.complete = TRUE
      }
    } else {
      # Sort according to p-value of LRT results
      explvar.todrop <- rownames(drop1.results[ order(drop1.results[,5], 
                                                      decreasing = TRUE), ][1,])
      pr.gt.chi <- drop1.results[ order(drop1.results[,5], 
                                        decreasing = TRUE), ][1,5]
      
      # If stepwise has reached the y-intercept, then no more model reduction can 
      # be done (linear seperability issue), stop.
      if (explvar.todrop == "<none>") {
        is.stepwise.complete = TRUE
      } else if (pr.gt.chi < pval.threshold) { # Dropped term is significant, stop.
        is.stepwise.complete = TRUE
      }
    }
    
    if (is.stepwise.complete == FALSE) {
      # Remove least significant variable and update the GLM
      glm.maximal <- update(glm.maximal, paste("~ . -", explvar.todrop))
      
      # Record results from simplified GLM
      stepwise.glm.pwr <- IterativeGlmPower(glm.maximal)
      explvars <- attr(glm.maximal$terms, "term.labels")
      df.stepresults <- rbind(df.stepresults,
                              data.frame(num.lags = i,
                                         num.explvars = length(explvars),
                                         deviance.null = glm.maximal$null.deviance,
                                         deviance.residuals = glm.maximal$deviance,
                                         AICc = stepwise.glm.pwr[1,2],
                                         BIC = stepwise.glm.pwr[1,3],
                                         mcfadden.r2 = stepwise.glm.pwr[1,4],
                                         formulastr = Reduce(paste, 
                                                             deparse(glm.maximal$formula,
                                                                     width.cutoff = 500)),
                                         variable.removed = explvar.todrop,
                                         probability.gt.chi = pr.gt.chi))
      
      # Print status update
      if (length(explvars) == 0) {
        is.stepwise.complete = TRUE
      } else {
        print(explvars)
      }
    }
  }
  
  return (df.stepresults)
}

CdhSumMaximalExplVars <- function(df.trimmed) {
  # This function returns the explanatory variable portion of a formula for all 
  # main effects and two-way interactions of the columns in a trimmed dataframe.
  # 
  # Args: 
  #   trimmed.df: A trimmed version of the meter readings dataframe.
  # 
  # Returns: 
  #   A string to be used in the explanatory variable portion of a formula which
  #   structures main effects and two-way interactions of every explanatory
  #   variable.
  explvars <- colnames(df.trimmed[, ! colnames(df.trimmed) %in% c("kwh")]) 
  explvars.collapsed <- paste(explvars, collapse = " + ")
  explvars.str <- paste0("(", explvars.collapsed, ")^2")
  return(explvars.str)
}

TempLagMaximalExplVars <- function(df.trimmed) {
  # This function returns the explanatory variable portion of a formula for all 
  # main effects and two-way interactions of the columns in a trimmed dataframe.
  # 
  # Args: 
  #   trimmed.df: A trimmed version of the meter readings dataframe.
  # 
  # Returns: 
  #   A string to be used in the explanatory variable portion of a formula which
  #   structures main effects and two-way interactions of every explanatory
  #   variable. Certain two-way interactions are removed and are replaced by 
  #   columns which represent the two-way interaction. (See function 
  #   StripNonexistantInteractions.)
  explvars <- colnames(df.trimmed[, ! colnames(df.trimmed) %in% c("kwh",
                                                                  "hrstr_price",
                                                                  "wknd_price")]) 
  explvars.collapsed <- paste(explvars, collapse = " + ")
  explvars.str <- paste0("(", explvars.collapsed, ")^2 - hrstr:price + hrstr_price - weekend:price + wknd_price")
  return(explvars.str)
}

TempLagNestedExplVars <- function(df.trimmed, templag.colnames) {
  # This function returns the explanatory variable portion of a formula for all
  # main effects of all dataframe columns, two-way interactions between 
  # categorical variables, and nested interactions of past hours' temperature >
  # breakpoint.
  #
  # Args:
  #   df.trimmed: A trimmed version of the meter readings dataframe.
  #   templag.colnames: The ordered column names of previous hours' degrees over
  #                     the temperature breakpoint.
  # 
  # Returns:
  #   A string to be used in the explanatory variable portion of a formula which
  #   structures main effects of all terms, two-way interactions of categorical 
  #   terms, and nested interactions of past hours' temperature > breakpoint.
  allcolnames <- colnames(df.trimmed)
  nontemp.colnames <- allcolnames[! allcolnames %in% c(templag.colnames, 
                                                       "kwh",
                                                       "hrstr_price",
                                                       "wknd_price")]
  main.str <- paste(c(nontemp.colnames, templag.colnames), collapse = " + ")
  nested.str <- paste(templag.colnames, collapse = "/")
  explvars.str <- paste0("(", main.str, ")^2 + ", nested.str)
  return(explvars.str)
}

CreateCdhSum <- function(nlags, df.readings) {
  # For each observation, this function sums all temperature values above the
  # CDH breakpoint from the current time through a number of hours into the 
  # past.
  # 
  # Args:
  #   nlags: The number of hours into the past (ie. lags) to sum for each 
  #          observation.
  #   df.readings: A dataframe of smart meter readings.
  #
  # Returns:
  #   A vector of values based on the conventional definition of CDH.  
  cdhvals <- rep(0, nrow(df.readings)) # Preallocate CDH vector
  
  # Unfortunately iterating over the dataframe seems to be the best method
  for(i in 1:nrow(df.readings)) {
    cdhvals[i] <- df.readings[i, "temp_over_break"]
    step <- 1
    
    # Sum what we can or all the values up to nlags
    while (i-step > 0 & step <= nlags) {
      cdhvals[i] <- cdhvals[i] + df.readings[i-step, "temp_over_break"]
      step <- step + 1
    }
  }
  
  return(cdhvals)
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
  lagnames <- paste0("templag", c(0:nlags))
  templags <- matrix(nrow = nrow(df.readings),
                     ncol = (nlags + 1))
  colnames(templags) <- lagnames
  
  # Unfortunately iterating over the dataframe seems to be the best method
  for(i in 1:nrow(df.readings)) {
    for(j in 0:nlags) {
      if (i-j > 0) {
        templags[i, (j+1)] <- df.readings[i-j, "temp_over_break"]
      } else {
        templags[i, (j+1)] <- 0
      }
    }
  }
  
  return(templags)
}

ForwardStepwiseLinearRegression <- function(nlags, df.obs, formula.maximal,
                                            is.stopcriteria.bic = FALSE,
                                            pval.threshold = 0.05) {
  # Performs forward stepwise linear regression including the provided number 
  # of hours of temperature history in its search space.
  #
  # TODO(r24mille): Iteratively building a dataframe is frowned upon. Consider
  #                 some other, more optimal method of building results.
  #
  # Args:
  #   nlags: The number of hours into the past that temperature > breakpoint 
  #          information that can be incorporated into the model.
  #   df.obs: The data.frame of observations to be used.
  #   formula.maximal: The formula for the maximal model to be considered.
  #   is.stopcriteria.bic: Controls the stopping criteria for the forward 
  #                        stepwise explanatory variable exploration. If BIC 
  #                        stopping criteria is used, then search stops when 
  #                        BIC stops decreasing. Default value is FALSE and 
  #                        p-value stopping criteria is used.
  #   pval.threshold: Controls the p-value stopping criteria. A likelihood ratio
  #                   test (LRT) is performed and the term is kept so long as 
  #                   it's p-value is < 0.05 (default). Other thresholds 
  #                   (eg. 0.01) may be passed in via this parameter.
  #
  # Return:
  #   A data.frame of results fitting the df.stepresults convention established 
  #   elsewhere.
  df.stepresults <- InitStepresultsDataframe()
  
  # Use the y-intercept GLM as a starting point
  stepwise.glm <- glm(formula = "(kwh) ~ 1", data = df.obs, family = gaussian)
  
  # I know... iteratively building a data.frame is bad  
  is.stepwise.complete = FALSE
  while (is.stepwise.complete == FALSE) {
    add1.results <- add1(object = stepwise.glm,
                         scope = formula.maximal,
                         test = "LRT", 
                         k = log(nrow(df.obs)))
    
    if (is.stopcriteria.bic == TRUE) {
      # Sort according to BIC
      explvar.toadd <- rownames(add1.results[ order(add1.results[,3],
                                                    decreasing = FALSE), ][1,])
      
      pr.gt.chi <- add1.results[ order(add1.results[,3], 
                                       decreasing = FALSE, 
                                       na.last = NA), ][1,5]
      
      # If stepwise has reached the y-intercept, then no more model construction
      # can be done.
      if (explvar.toadd == "<none>") {
        is.stepwise.complete = TRUE
      }
    } else {
      # Sort according to p-value of LRT results
      explvar.toadd <- rownames(add1.results[ order(add1.results[,5], 
                                                    decreasing = FALSE,
                                                    na.last = NA), ][1,])
      pr.gt.chi <- add1.results[ order(add1.results[,5], 
                                       decreasing = FALSE, 
                                       na.last = NA), ][1,5]
      
      # If stepwise has reached the y-intercept, then no more model construction
      # can be done.
      if (explvar.toadd == "<none>") {
        is.stepwise.complete = TRUE
      } else if (pr.gt.chi > pval.threshold) { # Term added isn't significant
        is.stepwise.complete = TRUE
      } 
    }
    
    if (is.stepwise.complete == FALSE) {
      # Remove least significant variable and update the GLM
      stepwise.glm <- update(stepwise.glm, 
                             paste("~ . + ", explvar.toadd))
      
      # Record results from simplified GLM
      stepwise.glm.pwr <- IterativeGlmPower(stepwise.glm)
      explvars <- attr(stepwise.glm$terms, "term.labels")
      df.stepresults <- rbind(df.stepresults,
                              data.frame(num.lags = nlags,
                                         num.explvars = length(explvars),
                                         deviance.null = stepwise.glm$null.deviance,
                                         deviance.residuals = stepwise.glm$deviance,
                                         AICc = stepwise.glm.pwr[1,2],
                                         BIC = stepwise.glm.pwr[1,3],
                                         mcfadden.r2 = stepwise.glm.pwr[1,4],
                                         formulastr = Reduce(paste, 
                                                             deparse(stepwise.glm$formula,
                                                                     width.cutoff = 500)),
                                         variable.added = explvar.toadd,
                                         probability.gt.chi = pr.gt.chi))
      
      print(explvars)
    }
  }
  
  return(df.stepresults)
}

GammaAIC <- function(fit) {
  # As per ?glm documentation, the default AIC value is wrong, "For gaussian,
  # Gamma and inverse gaussian families the dispersion is estimated from the
  # residual deviance, and the number of parameters is the number of
  # coefficients plus one. For a gaussian family the MLE of the dispersion is
  # used so this is a valid value of AIC, but for Gamma and inverse gaussian
  # families it is not."
  #
  # See supporting discussion at:
  # https://stackoverflow.com/questions/13405109/how-can-i-do-model-selection-by-aic-with-a-gamma-glm-in-r
  #
  # This function is inspired by that StackOverflow discussion. The conventional definition for AIC is:
  #
  # AIC = 2k − 2log(L) 
  # 
  # Where k is number of parameters and L is likelihood.
  #
  # Args:
  #   fit: A model fit by GLM, of family=Gamma()
  #
  # Return:
  #   A more accurate AIC value for the fitted model.
  
  # Dispersion is reciprocal of the estimate of the shape (eg. 1/shape) See
  # MASS:::gamma.dispersion
  # 
  # gamma.shape is estimated through maximum likelihood after fitting a Gamma
  # general linear model (ie. fit parameter). See ?MASS::gamma.shape
  # 
  # From documentation, "A glm fit for a Gamma family correctly calculates the
  # maximum likelihood estimate of the mean parameters but provides only a crude
  # estimate of the dispersion parameter. This function takes the results of the
  # glm fit and solves the maximum likelihood equation for the reciprocal of the
  # dispersion parameter, which is usually called the shape (or exponent) 
  # parameter."
  # 
  # So this method is a more accurate value for dispersion than the default
  # fit$deviance/length(fit$residuals).
  disp <- MASS::gamma.dispersion(fit) 
  # Rank is the number of parameters
  k <- fit$rank
  mu <- fit$fitted.values
  y <- fit$y
  return (2 * k - 2 * sum(dgamma(y, 1/disp, scale = mu * disp, log = TRUE)))
}

GammaAICc <- function(fit) {
  # Corrected AICc which uses the more accurate GammaAIC(...) function to 
  # computer the corrected AIC value. The equation for AICc is:
  #
  # AICc = AIC + (2k(k+1))/(n-k-1)
  #
  # Where k is the number of parameters and n denotes the sample size.
  #
  # See supporting discussion at:
  # https://stackoverflow.com/questions/13405109/how-can-i-do-model-selection-by-aic-with-a-gamma-glm-in-r
  # 
  # Args: 
  #   fit: A model fit by GLM, of family=Gamma()
  # 
  # Return: 
  #   A more accurate AICc value that incorporates the more accurate
  #   GammaAIC(...) function.
  val <- logLik(fit)
  k <- attributes(val)$df
  n <- attributes(val)$nobs
  return (GammaAIC(fit) + 2 * k * (k + 1) / (n - k - 1))
}

fixedGamma_extractAIC <- function(fit, scale=0, k=2, ...) {
  # A function that can fix the AIC measure for glm(family=Gamma()) fitted models
  # and does not break the AIC valu for other fitted models.
  #
  # Args:
  #   fit: A model fit by GLM, of family=Gamma()
  n <- length(fit$residuals)
  edf <- n - fit$df.residual  
  if (fit$family$family == "Gamma"){
    aic <- GammaAIC(fit)
  } else {
    aic <- fit$aic
  }
  return(c(edf, aic + (k - 2) * edf))
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

InitStepresultsDataframe <- function() {
  # Initializes an empty dataframe to store stepwise linear regression results 
  # in a common format.
  #
  # Return:
  #   An empty dataframe with columns structured appropriately for my iterative
  #   stepwise linear regression results.
  df.stepresults <- data.frame(num.lags = numeric(),
                               num.explvars = numeric(),
                               deviance.null = numeric(),
                               deviance.residuals = numeric(),
                               AICc = numeric(),
                               BIC = numeric(),
                               mcfadden.r2 = numeric(),
                               formulastr = character(),
                               variable.added = character(),
                               probability.gt.chi = numeric())
  
  return(df.stepresults)
}

IterativeGlmModel <- function(df.readings, nlags, 
                              is.touperiod, is.cdhsum, is.maximal) {
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
  #   nlags: The number of hours into the past that CDH will consider.
  #   is.touperiod: A boolean indicating whether the categorical tou_period
  #                 explanatory variable is used. If FALSE, then components of 
  #                 TOU will be used.
  #   is.cdhsum: A boolean indicating whether the CDH should be summed a number
  #              of hours into the past. If FALSE, then each historical 
  #              temperature > breakpoint value will be represented as its own 
  #              term.
  #   is.maximal: A boolean indicating whether a maximal two-way interaction 
  #               formula should be used. If FALSE, then a nested formula (ie. 
  #               past temperature nested interraction experiment) will be 
  #               created.
  #
  # Returns:
  #   A glm fitted according to the provided criteria.
  
  # CDH structure created (ie. vector or matrix)
  if(is.cdhsum == TRUE) {
    cdhsum <- CreateCdhSum(nlags, df.readings)
    readings.with.temp <- cbind(df.readings, cdhsum)
  } else {
    past_temps <- CreatePastTemperatureMatrix(nlags, df.readings)
    readings.with.temp <- cbind(df.readings, past_temps)
  }
  
  # Trim dataframe columns appropriate for desired representation of TOU
  if(is.touperiod == TRUE) {
    df.trimmed <- TrimColsTouPeriods(readings.with.temp)
  } else {
    df.trimmed <- TrimColsTouTimeComponents(readings.with.temp)
  }
  
  # Create appropriate GLM formula
  if (is.maximal == TRUE & is.touperiod == TRUE) {
    explvars.str <- CdhSumMaximalExplVars(df.trimmed)
  } else if (is.maximal == TRUE & is.touperiod == FALSE) {
    explvars.str <- TempLagMaximalExplVars(df.trimmed)
  } else {
    explvars.str <- TempLagNestedExplVars(df.trimmed, colnames(past_temps))
  }
  
  # Fit the GLM model
  fmla <- formula(paste0("log(kwh) ~ ", explvars.str))
  print(fmla)
  iter.glm <- glm(formula = fmla, 
                  data = df.trimmed, 
                  family = gaussian)
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

PerformTouCdhGlmIterations <- function(df.readings, nhrs) {
  # Steps through 0-nhrs iterations of GLMs trying varioud combinations of 
  # structuring TOU and CDH explanatory variables. The explanatory power of each
  # model is plotted.
  #
  # Args:
  #   df.readings: A readings dataframe to be passed as the "data" parameter 
  #                to the GLMs.
  #   nhrs: Number of hours to include in the CDH component of the model
  cdhsum.touperiods.maxglm.pwr <- GlmPowerResultsMatrix(nhrs)
  tempterms.touperiods.nestedglm.pwr <- GlmPowerResultsMatrix(nhrs)
  tempterms.touperiods.maxglm.pwr <- GlmPowerResultsMatrix(nhrs)
  cdhsum.timecomps.maxglm.pwr <- GlmPowerResultsMatrix(nhrs)
  tempterms.timecomps.nestedglm.pwr <- GlmPowerResultsMatrix(nhrs)
  tempterms.timecomps.maxglm.pwr <- GlmPowerResultsMatrix(nhrs)
  
  for(i in 0:nhrs) {
    # 1. Summed CDH lags, TOU as periods
    cdhsum.touperiods.maxglm <- IterativeGlmModel(df.readings = df.readings, 
                                                  nlags = i, 
                                                  is.touperiod = TRUE, 
                                                  is.cdhsum = TRUE, 
                                                  is.maximal = TRUE)
    cdhsum.touperiods.maxglm.pwr[(i+1),] <- IterativeGlmPower(cdhsum.touperiods.maxglm)
    
    # 2. Matrix of CDH lags as nested interactions, TOU as periods
    tempterms.touperiods.nestedglm <- IterativeGlmModel(df.readings = df.readings, 
                                                        nlags = i, 
                                                        is.touperiod = TRUE, 
                                                        is.cdhsum = FALSE, 
                                                        is.maximal = FALSE)
    tempterms.touperiods.nestedglm.pwr[(i+1),] <- IterativeGlmPower(tempterms.touperiods.nestedglm)
    
    # 3. Matrix of CDH lags as two-way interactions, TOU as periods
    tempterms.touperiods.maxglm <- IterativeGlmModel(df.readings = df.readings, 
                                                     nlags = i, 
                                                     is.touperiod = TRUE, 
                                                     is.cdhsum = FALSE, 
                                                     is.maximal = TRUE)
    tempterms.touperiods.maxglm.pwr[(i+1),] <- IterativeGlmPower(tempterms.touperiods.maxglm)
    
    # 4. Summed CDH lags, TOU components of time
    cdhsum.timecomps.maxglm <- IterativeGlmModel(df.readings = df.readings, 
                                                 nlags = i, 
                                                 is.touperiod = FALSE, 
                                                 is.cdhsum = TRUE, 
                                                 is.maximal = TRUE)
    cdhsum.timecomps.maxglm.pwr[(i+1),] <- IterativeGlmPower(cdhsum.timecomps.maxglm)
    
    # 5. Matrix of CDH lags as nested interactions, TOU components of time
    tempterms.timecomps.nestedglm <- IterativeGlmModel(df.readings = df.readings, 
                                                       nlags = i, 
                                                       is.touperiod = FALSE, 
                                                       is.cdhsum = FALSE, 
                                                       is.maximal = FALSE)
    tempterms.timecomps.nestedglm.pwr[(i+1),] <- IterativeGlmPower(tempterms.timecomps.nestedglm)
    
    # 6. Matrix of CDH lags as two-way interactions, TOU components of time
    tempterms.timecomps.maxglm <- IterativeGlmModel(df.readings = df.readings, 
                                                    nlags = i, 
                                                    is.touperiod = FALSE, 
                                                    is.cdhsum = FALSE, 
                                                    is.maximal = TRUE)
    tempterms.timecomps.maxglm.pwr[(i+1),] <- IterativeGlmPower(tempterms.timecomps.maxglm)
  }
  
  IterativeGlmPlots(iter.glm.pwr = cdhsum.touperiods.maxglm.pwr, 
                    maintitle = "TOU Period, Cooling Degree Hours")
  
  IterativeGlmPlots(iter.glm.pwr = tempterms.touperiods.nestedglm.pwr, 
                    maintitle = expression(paste("TOU Period, with Degrees>", 
                                                 Temp['break'], 
                                                 " Nested Interaction")))
  
  IterativeGlmPlots(iter.glm.pwr = tempterms.touperiods.maxglm.pwr, 
                    maintitle = expression(paste("TOU Period, with Degrees>", 
                                                 Temp['break'], 
                                                 " and All 2-Way Interactions")))
  
  IterativeGlmPlots(iter.glm.pwr = cdhsum.timecomps.maxglm.pwr, 
                    maintitle = "Time Components, Cooling Degree Hours")
  
  IterativeGlmPlots(iter.glm.pwr = tempterms.timecomps.nestedglm.pwr, 
                    maintitle = expression(paste("Time Components, with Degrees>", 
                                                 Temp['break'], 
                                                 " Nested Interaction")))
  
  IterativeGlmPlots(iter.glm.pwr = tempterms.timecomps.maxglm.pwr, 
                    maintitle = expression(paste("Time Components, with Degrees>", 
                                                 Temp['break'], 
                                                 " and All 2-Way Interactions")))
}

TrimColsTouPeriods <- function(df.readings) {
  # Trim the columns down from full readings.aggregate to only those needed 
  # when working with "TOU Period" version of the model.
  #
  # Args:
  #   df.readings: A dataframe object of aggregate readings
  #
  # Returns:
  #   A trimmed dataframe in which unneeded columns from readings.aggregate 
  #   dataframe have been trimmed out. Only column with data relevant to 
  #   representing TOU via its periods (ie. tou_period) are preserved. 
  #   Ontario's TOU periods are essentially a proxy for day, time, demand, and 
  #   price.  
  df.periods <- df.readings[, ! colnames(df.readings) %in% c("daynum", 
                                                             "hour", 
                                                             "timestamp_dst", 
                                                             "temperature", 
                                                             "agg_count", 
                                                             "hrstr", 
                                                             "weekend", 
                                                             "price", 
                                                             "temp_over_break",
                                                             "hrstr_price",
                                                             "wknd_price")]
  return(df.periods)
}


TrimColsTouTimeComponents <- function(df.readings) {
  # Trim the columns down from full readings.aggregate to only those needed 
  # when working with "Components of TOU" version of the model.
  #
  # Args:
  #   df.readings: A dataframe object of aggregate readings
  #
  # Returns:
  #   A trimmed dataframe in which unneeded columns from readings.aggregate 
  #   dataframe have been trimmed out. Only columns with data relevant to 
  #   representing TOU via all of its components (ie. hour-of-day, 
  #   weekday/weekend, and price) are preserved.
  df.timecomps <- df.readings[, ! colnames(df.readings) %in% c("daynum", 
                                                              "hour", 
                                                              "timestamp_dst", 
                                                              "temperature", 
                                                              "agg_count", 
                                                              "tou_period", 
                                                              "temp_over_break",
                                                              "billing_active")]
  return(df.timecomps)
}