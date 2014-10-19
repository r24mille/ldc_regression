library(stargazer) # LaTeX tables
library(segmented) # Find linear regression breakpoint
library(sme) # For AICc function
library(BMA) # Compare GLM models
library(reshape) # Reshape data.frame for the heatmap
library(rgl) # OpenGL functionality in R for better 3D plots

# Source the function in another file
source('glm_method_iteration.R')
# Fixes the AIC measure in the add1(...) and drop1(...) functions, also BIC
assignInNamespace(x="extractAIC.glm", value=fixedGamma_extractAIC, ns="stats")

# Load SmartMeterReading data from CSV
home <- Sys.getenv("HOME")
fpath <- file.path(home, 
                   "../Dropbox/ISS4E/R/", 
                   "full_aggregate_readings.csv")
readings.aggregate <- read.csv(fpath)

# Re-orders TOU Period levels so that graphs are sorted accordingly
readings.aggregate$tou_period <- factor(readings.aggregate$tou_period, 
                                        c("off_weekend", "off_morning", "mid_morning", 
                                          "on_peak", "mid_evening", "off_evening"))

##
# Add column which represents "month" as a categorical factor
readings.aggregate$timestamp_dst <- as.POSIXlt(readings.aggregate$timestamp_dst)
readings.aggregate$month <- paste0("m", (readings.aggregate$timestamp_dst$mon + 1))
readings.aggregate$month <- factor(readings.aggregate$month, 
                                   c("m5", "m6", "m7", "m8", "m9", "m10"))

# Represent hours as levels rather than integers
readings.aggregate$hrstr <- paste0("h", readings.aggregate$hour)
readings.aggregate$hrstr <- factor(readings.aggregate$hrstr, 
                                   c("h0", "h1", "h2", "h3", "h4", "h5", "h6", 
                                     "h7", "h8", "h9", "h10", "h11", "h12", 
                                     "h13", "h14", "h15", "h16", "h17", "h18", 
                                     "h19", "h20", "h21", "h22", "h23"))

##
# Add yes/no weekend flag
readings.aggregate$weekend <- ifelse(readings.aggregate$tou_period %in% c("off_weekend"),
                                     "Yes",
                                     "No")
readings.aggregate$weekend <- factor(readings.aggregate$weekend, 
                                     c("No", 
                                       "Yes"))

##
# Add column which converts TOU Period to its price level
readings.aggregate$price <- readings.aggregate$tou_period
readings.aggregate$price <- gsub("off_weekend", "off_peak", readings.aggregate$price)
readings.aggregate$price <- gsub("off_morning", "off_peak", readings.aggregate$price)
readings.aggregate$price <- gsub("off_evening", "off_peak", readings.aggregate$price)
readings.aggregate$price <- gsub("mid_morning", "mid_peak", readings.aggregate$price)
readings.aggregate$price <- gsub("mid_evening", "mid_peak", readings.aggregate$price)
readings.aggregate$price <- factor(readings.aggregate$price, 
                                   c("off_peak", "mid_peak", "on_peak"))

##
# Use 'segmented' package rather than my prior home-grown method of finding 
# the optimal cooling degree hour breakpoint (though they give the same 
# result).
#
# TODO(r24mille): segmented chooses very different cdhbreak numbers based on 
#                 the model chosen:
#                      lm(...) is 18
#                      glm(...) is 18
#                      glm(..., family=Gamma) is 13
#                      glm(...,family=Gamma(link="log")) is 16
#                 Maybe just making an empirical choice for CDH break within 
#                 the context of the final model is the best choice.
model.readings.glm.presegment <- glm(kwh ~ temperature*tou_period*billing_active,
                                     data = readings.aggregate,
                                     family = Gamma(link="log"))
seg <- segmented(obj = model.readings.glm.presegment, 
                 seg.Z = ~temperature,
                 psi = list(temperature = c(18)))
cdhbreak <- floor(seg$psi[1,2]) # TODO floor vs. round vs. real temps
readings.aggregate$cdh <- ifelse(readings.aggregate$temperature > cdhbreak, 
                                 readings.aggregate$temperature - cdhbreak, 
                                 0)

##
# Find the optimal number of hours lag, and the best method for incorporating 
# CDH. Is it best to sum up current and past hours, similar to traditional CDH 
# or should previous hours be nested under the current hour?
#
# TODO(r24mille): Resid. Dev is lower with link="log". Justify this choice.
#                   * This is the same as lm(log(kwh) ~ .) but allows for 
#                     easier interpretation of the expected value since it does 
#                     not need to be back transformed.
#                   * Also allows for different GLM models to be compared on the
#                     same terms using ANOVA (ie. analysis of deviance)

# Commenting out the iterative comparison, now that I have the results from it.
# PerformTouCdhGlmIterations(df.readings = readings.aggregate,
#                          nhrs = 5)

# TOU components and each CDH lag as its own coefficient turns out to have the 
# highest predictive power with 8 hours of history.
# cdhlagmat.toucomps.maxglm <- IterativeGlmModel(df.readings = readings.aggregate, 
#                                                nlags = 8, 
#                                                is.touperiod = FALSE, 
#                                                is.cdhlagsum = FALSE, 
#                                                is.maxformula = TRUE)

##
# STEP Work-in-Progress
nlags = 24

# Set up 3D results dataframe
df.stepresults <- data.frame(num.cdhlags = numeric(),
                             num.explvars = numeric(),
                             deviance.null = numeric(),
                             deviance.residuals = numeric(),
                             AICc = numeric(),
                             BIC = numeric(),
                             mcfadden.r2 = numeric(),
                             formulastr = character(),
                             variable.removed = character(),
                             probability.gt.chi = numeric())

# Iterate through all steps of the current nlag
for(i in 0:nlags) {
  maxtractable.glm <- MaximumTractableInteractionsGlmModel(df.readings = readings.aggregate, 
                                                           nlags = i)
  
  # Copy trimmed dataframe created within function for use later
  df.trimmed <- maxtractable.glm$model
  num.observations <- nrow(df.trimmed)
  
  # I know... iteratively building a data.frame is bad
  maxtractable.glm.pwr <- IterativeGlmPower(maxtractable.glm)
  explvars <- attr(maxtractable.glm$terms, "term.labels")
  df.stepresults <- rbind(df.stepresults,
                          data.frame(num.cdhlags = i,
                                     num.explvars = length(explvars),
                                     deviance.null = maxtractable.glm$null.deviance,
                                     deviance.residuals = maxtractable.glm$deviance,
                                     AICc = maxtractable.glm.pwr[1,2],
                                     BIC = maxtractable.glm.pwr[1,3],
                                     mcfadden.r2 = maxtractable.glm.pwr[1,4],
                                     formulastr = Reduce(paste, 
                                                         deparse(maxtractable.glm$formula,
                                                                 width.cutoff = 500)),
                                     variable.removed = "",
                                     probability.gt.chi = 0))
  
  # Print status update
  print(explvars)
  
  # Flag whether stepwise deletion of variables has completed
  is.stepwise.nlags.complete = FALSE
  while (is.stepwise.nlags.complete == FALSE) {
    drop1.results <- drop1(object = maxtractable.glm,
                           test = "LRT", 
                           k = log(num.observations))
    # Sort drop1.results according to BIC descending, remove the term resulting 
    # in the highest BIC.
    explvar.leastsig <- rownames(drop1.results[ order(drop1.results[,5], 
                                                      decreasing = TRUE), ][1,])
    pr.gt.chi <- drop1.results[ order(drop1.results[,5], 
                                      decreasing = TRUE), ][1,5]
    
    # If stepwise has reached the y-intercept, then no more model reduction can 
    # be done (linear seperability issue). Stop and iterate to the next cdhlag.
    #
    # TODO(r24mille): Can this be corrected in the MaximumTractableInteractionsGlmModel(...)
    #                 function?
    if (explvar.leastsig == "<none>") {
      is.stepwise.nlags.complete = TRUE
    } else if (pr.gt.chi < 0.05) { # Term is statistically significant, stop.
      is.stepwise.nlags.complete = TRUE
    } else {
      # Remove least significant variable and update the GLM
      maxtractable.glm <- update(maxtractable.glm, 
                                 paste("~ . -", explvar.leastsig))
      
      # Record results from simplified GLM
      maxtractable.glm.pwr <- IterativeGlmPower(maxtractable.glm)
      explvars <- attr(maxtractable.glm$terms, "term.labels")
      df.stepresults <- rbind(df.stepresults,
                              data.frame(num.cdhlags = i,
                                         num.explvars = length(explvars),
                                         deviance.null = maxtractable.glm$null.deviance,
                                         deviance.residuals = maxtractable.glm$deviance,
                                         AICc = maxtractable.glm.pwr[1,2],
                                         BIC = maxtractable.glm.pwr[1,3],
                                         mcfadden.r2 = maxtractable.glm.pwr[1,4],
                                         formulastr = Reduce(paste, 
                                                             deparse(maxtractable.glm$formula,
                                                                     width.cutoff = 500)),
                                         variable.removed = explvar.leastsig,
                                         probability.gt.chi = pr.gt.chi))
      
      # Print status update
      if (length(explvars) == 0) {
        is.stepwise.nlags.complete = TRUE
      } else {
        print(explvars)
      }
    }
  }
}



##
# 3D Plot goodness-of-fit

# The y-intercept is not interesting, the same for every model, and makes 
# the plot difficult to look at.
df.stepresults.trimmed <- subset(df.stepresults, num.explvars > 0)
wireframe(BIC ~ num.explvars * num.cdhlags, 
          data = df.stepresults.trimmed,
          xlab = list("Number of Explanatory Variables", rot=-13), 
          ylab = list(expression(paste("Hours of Temp.>", CDH['break'], " Included")), rot=65),
          zlab = list("Bayesian Information Criterion (BIC)", rot=92), 
          main = "Change in BIC During Stepwise Removal of Terms from Maximal Model",
          drape = TRUE,
          colorkey = TRUE,
          scales = list(arrows = FALSE), # Switches unlabelled arrows to ticks
          screen = list(z = -20, x = -60))
# Decent setting for BIC
# screen = list(z = -20, x = -60)

# Commenting out cloud points for now
# minbic <- min(df.stepresults.trimmed$BIC)
# df.testcloud <- data.frame(thex = c(40, 43, 45), they = c(3, 4, 5), thez = c(-5000, -1000, -2000))
# cloud(thez ~ thex * they, data = df.testcloud)

# That same information plotted with RGL
df.heat <- df.stepresults.trimmed
df.heat.reshape <- melt(data = df.heat,
                        id.vars = c("num.cdhlags", "num.explvars"),
                        measure.vars = c("BIC"))
cast.heat <- cast(data = df.heat.reshape, 
                  formula = num.cdhlags ~ num.explvars)
stepresults.matrix <- as.matrix(cast.heat)
bic.lim <- range(df.stepresults.trimmed$BIC)
bic.len <- bic.lim[2] - bic.lim[1] + 1
hgt.lookup <- terrain.colors(bic.len) # height color lookup table
bic.colors <- hgt.lookup[stepresults.matrix - bic.lim[1] + 1] # assign colors to heights
persp3d(x = seq(0, (nrow(stepresults.matrix) - 1), len = nrow(stepresults.matrix)),
        y = seq(1, ncol(stepresults.matrix), len = ncol(stepresults.matrix)),
        z = stepresults.matrix, 
        xlab = "Number of Past Hours Included",
        ylab = "Number of Explanatory Variables",
        zlab = "Bayesian Information Criterion (BIC)",
        main = "Change in BIC During Stepwise Removal of Terms from Maximal Model",
        col = bic.colors)

# Heatmap
#maxbic <- max(df.stepresults.trimmed$BIC)
#stepresults.matrix[which(is.na(stepresults.matrix))] = maxbic;
heatmap(stepresults.matrix, 
        Rowv=NA, 
        Colv=NA, 
        col = heat.colors(round(bic.len) * 2), 
        margins=c(5,10),
        na.rm = TRUE)

# ADD1 STEP work-in-process
# Iterate steps forward and store in a similar dataframe
nlags = 24

df.results.stepforward <- data.frame(num.cdhlags = numeric(),
                                     num.explvars = numeric(),
                                     deviance.null = numeric(),
                                     deviance.residuals = numeric(),
                                     AICc = numeric(),
                                     BIC = numeric(),
                                     mcfadden.r2 = numeric(),
                                     formulastr = character(),
                                     variable.added = character(),
                                     probability.gt.chi = numeric())

# Iterate through all steps of the current nlag
for(i in 0:nlags) {
  # Iterate through all steps of the current nlag
  maxtractable.glm <- IterativeGlmModel(df.readings = readings.aggregate, 
                                        nlags = i, 
                                        is.touperiod = FALSE, 
                                        is.cdhlagsum = FALSE, 
                                        is.maxformula = TRUE)
  maxformulastr <- Reduce(paste, 
                          deparse(maxtractable.glm$formula,
                                  width.cutoff = 500))
  maxformulastr <- paste(maxformulastr, " - hrstr:price - weekend:price")
  formula.maxtractable <- formula(maxformulastr)
  
  # Copy trimmed dataframe created within function for use later
  df.trimmed <- maxtractable.glm$model
  num.observations <- nrow(df.trimmed)
  
  # Use the y-intercept GLM as a starting point
  stepwise.glm <- glm(formula = "kwh ~ 1", 
                      data = df.trimmed, 
                      family = Gamma(link="log"))
  
  # I know... iteratively building a data.frame is bad  
  is.stepwise.nlags.complete = FALSE
  while (is.stepwise.nlags.complete == FALSE) {
    add1.results <- add1(object = stepwise.glm,
                         scope = formula.maxtractable,
                         test = "LRT", 
                         k = log(num.observations))
    # Sort add1.results according to information criteria and assign to a 
    # variable.
    explvar.lowestPVal <- rownames(add1.results[ order(add1.results[,5], 
                                                      decreasing = FALSE,
                                                      na.last = NA), ][1,])
    pr.gt.chi <- add1.results[ order(add1.results[,5], 
                                     decreasing = FALSE, 
                                     na.last = NA), ][1,5]
    
    # If stepwise has reached the y-intercept, then no more model construction
    # can be done.
    if (explvar.lowestPVal == "<none>") {
      is.stepwise.nlags.complete = TRUE
    } else if (pr.gt.chi > 0.05) { # Addition of term is not significant at p=0.05
      is.stepwise.nlags.complete = TRUE
    }else {
      # Remove least significant variable and update the GLM
      stepwise.glm <- update(stepwise.glm, 
                             paste("~ . + ", explvar.lowestPVal))
      
      # Record results from simplified GLM
      stepwise.glm.pwr <- IterativeGlmPower(stepwise.glm)
      explvars <- attr(stepwise.glm$terms, "term.labels")
      df.results.stepforward <- rbind(df.results.stepforward,
                                      data.frame(num.cdhlags = i,
                                                 num.explvars = length(explvars),
                                                 deviance.null = stepwise.glm$null.deviance,
                                                 deviance.residuals = stepwise.glm$deviance,
                                                 AICc = stepwise.glm.pwr[1,2],
                                                 BIC = stepwise.glm.pwr[1,3],
                                                 mcfadden.r2 = stepwise.glm.pwr[1,4],
                                                 formulastr = Reduce(paste, 
                                                                     deparse(stepwise.glm$formula,
                                                                             width.cutoff = 500)),
                                                 variable.added = explvar.lowestPVal,
                                                 probability.gt.chi = pr.gt.chi))
      
      print(explvars)
    }
  }
}

df.results.stepforward.trimmed <- subset(df.results.stepforward, num.explvars > 0)
wireframe(BIC ~ num.explvars * num.cdhlags, 
          data = df.results.stepforward.trimmed,
          xlab = list("Number of Explanatory Variables", 
                      rot=-13), 
          ylab = list("Number of Past Hours of Temp. > CDH_break Included", 
                      rot=65),
          zlab = list("Bayesian Information Criterion (BIC)", 
                      rot=90), 
          main = "Change in BIC During Stepwise Addition of Terms to GLM",
          drape = TRUE,
          colorkey = TRUE,
          scales = list(arrows = FALSE), # Switches unlabelled arrows to ticks
          screen = list(z = -20, x = -60))


##
# BIC.GLM leap variant work-in-progress
# bic.glm(...) isn't working due to singularity issues, I believe.
# averagingtractable.glm <- MaximumTractableInteractionsGlmModel(df.readings = readings.aggregate,
#                                                          nlags = 8)
# df.averaging.trimmed <- averagingtractable.glm$model
# formula.string.averaging = Reduce(paste, deparse(averagingtractable.glm$formula, 
#                                                  width.cutoff = 500))
# formula.averaging <- formula(formula.string.averaging)
# bma.results <- bic.glm(f = formula.averaging,
#                        data = df.averaging.trimmed,
#                        glm.family = Gamma(link="log"),
#                        maxCol = 50,
#                        nbest = 5)


Plot2DFitByExplVarCountWithMultiplePastHrsTemp(df.steps = df.stepresults,
                                               is.bic = FALSE,
                                               title = "BIC Change as Terms are Removed from Maximal Models",
                                               subtitle = "(p-value stopping criterion)")