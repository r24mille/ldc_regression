BackwardStepwiseRemoval <- function(df.obs, 
                                    formula.saturated,
                                    pval.threshold = 0.05) {
  # Performs backward stepwise linear regression including the provided number 
  # of hours of temperature history in its search space.
  #
  # TODO(r24mille): Iteratively building a dataframe is frowned upon. Consider
  #                 some other, more optimal method of building results.
  #
  # Args:
  #   df.obs: The data.frame of observations to be used.
  #   formula.saturated: The formula for the maximal model to be considered.
  #   pval.threshold: Controls the p-value stopping criteria. A likelihood ratio
  #                   test (LRT) is performed and the term is kept so long as 
  #                   it's p-value is < 0.05 (default). Other thresholds 
  #                   (eg. 0.01) may be passed in via this parameter.
  #
  # Return:
  #   A data.frame of results fitting the df.stepresults convention established 
  #   elsewhere.  
  lm.fitted <- lm(formula = formula.saturated, 
                  data = df.obs)
  
  # I know... iteratively building a data.frame is bad
  explvars <- attr(lm.fitted$terms, "term.labels")
  lm.summary <- summary(lm.fitted)
  df.stepresults <- data.frame(num.explvars = length(explvars),
                               residual.se = lm.summary$sigma,
                               mult.r2 = lm.summary$r.squared,
                               adj.r2 = lm.summary$adj.r.squared,
                               fstat.val = lm.summary$fstatistic[["value"]],
                               fstat.numdf = lm.summary$fstatistic[["numdf"]],
                               fstat.dendf = lm.summary$fstatistic[["dendf"]],
                               formulastr = paste(lm.summary$call[[2]]),
                               variable.removed = "",
                               probability.gt.f = 0)
  
  # Print status update
  print(explvars)
  
  # Flag whether stepwise deletion of variables has completed
  is.stepwise.complete = FALSE
  while (is.stepwise.complete == FALSE) {
    drop1.results <- drop1(object = lm.fitted,
                           test = "F")
    
    # Sort according to p-value of F test
    explvar.todrop <- rownames(drop1.results[ order(drop1.results[,6], 
                                                    decreasing = TRUE), ][1,])
    pr.gt.f <- drop1.results[ order(drop1.results[,6], 
                                    decreasing = TRUE), ][1,6]
    
    # If stepwise has reached the y-intercept, then no more model reduction can 
    # be done (linear seperability issue), stop.
    if (explvar.todrop == "<none>") {
      is.stepwise.complete = TRUE
    } else if (pr.gt.f < pval.threshold) { # Dropped term is significant, stop.
      is.stepwise.complete = TRUE
    }
    
    if (is.stepwise.complete == FALSE) {
      # Remove least significant variable and update the GLM
      lm.fitted <- update(lm.fitted, paste("~ . -", explvar.todrop))
      
      # Record results from simplified GLM
      explvars <- attr(lm.fitted$terms, "term.labels")
      lm.summary <- summary(lm.fitted)
      df.stepresults <- rbind(df.stepresults,
                              data.frame(num.explvars = length(explvars),
                                         residual.se = lm.summary$sigma,
                                         mult.r2 = lm.summary$r.squared,
                                         adj.r2 = lm.summary$adj.r.squared,
                                         fstat.val = lm.summary$fstatistic[["value"]],
                                         fstat.numdf = lm.summary$fstatistic[["numdf"]],
                                         fstat.dendf = lm.summary$fstatistic[["dendf"]],
                                         formulastr = paste(lm.summary$call[[2]][3]),
                                         variable.removed = explvar.todrop,
                                         probability.gt.f = pr.gt.f))
      
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