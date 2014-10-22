library(MASS) # For correcting Gamma glm AIC

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
                          data.frame(num.cdhlags = nlags,
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
      # be done (linear seperability issue). Stop and iterate to the next cdhlag.
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
                              data.frame(num.cdhlags = i,
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

CdhLagMaximalFormula <- function(df.trimmed) {
  # This function returns the forumla for all main effects and two-way 
  # interactions of the columns in the trimmed dataframe provided.
  #
  # Args:
  #   trimmed.df: A trimmed version of the readings.aggregate dataframe.
  #
  # Returns:
  #   A formula object, where each independent variable is modeled as a fixed 
  #   effect with main effects and two-way interactions.
  explvars <- colnames(df.trimmed[, ! colnames(df.trimmed) %in% c("kwh")]) 
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
  #   two-way interracted with other main effects.
  allcolnames <- colnames(df.trimmed)
  noncdh.colnames <- allcolnames[! allcolnames %in% c(cdhlagmat.colnames, 
                                               "kwh")]
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
                              data.frame(num.cdhlags = nlags,
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
  df.stepresults <- data.frame(num.cdhlags = numeric(),
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
  
  # Stitch the two data structures together into a dataframe.
  readings.with.cdh <- cbind(df.readings, cdhlag)
  
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

PerformTouCdhGlmIterations <- function(df.readings, nhrs) {
  # Steps through 0-nhrs iterations of GLMs trying varioud combinations of 
  # structuring TOU and CDH explanatory variables. The explanatory power of each
  # model is plotted.
  #
  # Args:
  #   df.readings: A readings dataframe to be passed as the "data" parameter 
  #                to the GLMs.
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
                                                     nlags = i, 
                                                     is.touperiod = TRUE, 
                                                     is.cdhlagsum = TRUE, 
                                                     is.maxformula = TRUE)
    cdhlagsum.touperiods.maxglm.pwr[(i+1),] <- IterativeGlmPower(cdhlagsum.touperiods.maxglm)
    
    # 2. Matrix of CDH lags as nested interactions, TOU as periods
    cdhlagmat.touperiods.nestedglm <- IterativeGlmModel(df.readings = df.readings, 
                                                        nlags = i, 
                                                        is.touperiod = TRUE, 
                                                        is.cdhlagsum = FALSE, 
                                                        is.maxformula = FALSE)
    cdhlagmat.touperiods.nestedglm.pwr[(i+1),] <- IterativeGlmPower(cdhlagmat.touperiods.nestedglm)
    
    # 3. Matrix of CDH lags as two-way interactions, TOU as periods
    cdhlagmat.touperiods.maxglm <- IterativeGlmModel(df.readings = df.readings, 
                                                     nlags = i, 
                                                     is.touperiod = TRUE, 
                                                     is.cdhlagsum = FALSE, 
                                                     is.maxformula = TRUE)
    cdhlagmat.touperiods.maxglm.pwr[(i+1),] <- IterativeGlmPower(cdhlagmat.touperiods.maxglm)
    
    # 4. Summed CDH lags, TOU components of time
    cdhlagsum.toucomps.maxglm <- IterativeGlmModel(df.readings = df.readings, 
                                                   nlags = i, 
                                                   is.touperiod = FALSE, 
                                                   is.cdhlagsum = TRUE, 
                                                   is.maxformula = TRUE)
    cdhlagsum.toucomps.maxglm.pwr[(i+1),] <- IterativeGlmPower(cdhlagsum.toucomps.maxglm)
    
    # 5. Matrix of CDH lags as nested interactions, TOU components of time
    cdhlagmat.toucomps.nestedglm <- IterativeGlmModel(df.readings = df.readings, 
                                                      nlags = i, 
                                                      is.touperiod = FALSE, 
                                                      is.cdhlagsum = FALSE, 
                                                      is.maxformula = FALSE)
    cdhlagmat.toucomps.nestedglm.pwr[(i+1),] <- IterativeGlmPower(cdhlagmat.toucomps.nestedglm)
    
    # 6. Matrix of CDH lags as two-way interactions, TOU components of time
    cdhlagmat.toucomps.maxglm <- IterativeGlmModel(df.readings = df.readings, 
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