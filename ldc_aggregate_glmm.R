library(segmented) # Find linear regression breakpoint
library(sme) # For AICc function
library(BMA) # Bayesian Model Averaging
library(bestglm) # For bestglm selection using leaps
library(lars) # A separate lasso package
library(glmnet) # Explanatory variable selection by lasso
library(boot) # Cross-validation

# Experiment with devtools and fastVAR
#library(devtools)
#install_github("fastVAR", "jeffwong", "master")
#instalibrary(fastVAR)

# Source functions in other files
source('dataframe_processing.R')
source('glm_method_iteration.R')
# Fixes the AIC measure in the add1(...) and drop1(...) functions, also BIC
#assignInNamespace(x="extractAIC.glm", value=fixedGamma_extractAIC, ns="stats")

# Load SmartMeterReading data from CSV
home <- Sys.getenv("HOME")
fpath <- file.path(home, 
                   "../Dropbox/ISS4E/R/", 
                   "full_aggregate_readings.csv")
readings.aggregate <- read.csv(fpath)
readings.aggregate <- InitAggregateReadings(readings.aggregate)

##
# Use 'segmented' package rather than my prior home-grown method of finding 
# the optimal cooling degree hour breakpoint (though they give the same 
# result).
#
# TODO(r24mille): segmented chooses very different temperature.break numbers 
#                 based on the model chosen. Maybe just making an empirical 
#                 choice for CDH break within the context of the final model is 
#                 the best choice.
model.readings.glm.presegment <- glm(log(kwh)~ month + hrstr + weekend + price + temperature,
                                     data = readings.aggregate,
                                     family = gaussian)
seg <- segmented(obj = model.readings.glm.presegment, 
                 seg.Z = ~ temperature,
                 psi = list(temperature = c(18)))

# plot(seg, 
#      res=TRUE, 
#      res.col = rgb(0, 0, 0, 10, maxColorValue=100),
#      pch = 16,
#      main = "Piecewise Linear Fit with Temperature Breakpoint",
#      xlab = "Summer Temperatures (Celsius)",
#      ylab = "Effect of Temperature",
#      col = "red",
#      lwd = 2)

temperature.break <- floor(seg$psi[1,2]) # TODO floor vs. round vs. real temps
readings.aggregate$temp_over_break <- ifelse(readings.aggregate$temperature > temperature.break, 
                                             readings.aggregate$temperature - temperature.break, 
                                             0)

plot(y = readings.aggregate$kwh,
     x = readings.aggregate$temp_over_break,
     col = rgb(0, 0, 0, 20, maxColorValue=100),
     pch = 16,
     main = "Average Household Demand as a Function of Current Degrees>Breakpoint",
     ylab = "Average Household Demand (kWh)",
     xlab = "Current Degrees>Breakpoint (Celsius)")

# Geometric mean of past temperatures
for(n in 1:5) {
  temps.geometricmean.order.natural <- rep(1, nrow(readings.aggregate))
  temps.geometricmean.order.correlation <- rep(1, nrow(readings.aggregate))
  temp.correlation.order <- c(2, 3, 1, 4, 0)
  for(j in 1:n) {
    for(i in 1:nrow(readings.aggregate)) {
      # If index of past reading is greater than 0, use it... otherwise insert 0
      ifelse(i-j+1 > 0,
             temps.geometricmean.order.natural[i] <- temps.geometricmean.order.natural[i] * readings.aggregate[i-j+1,]$temp_over_break,
             temps.geometricmean.order.natural[i] <- 0)
      
      # If index of past reading is greater than 0, use it... otherwise insert 0
      ifelse(i-temp.correlation.order[j] > 0, 
             temps.geometricmean.order.correlation[i] <- temps.geometricmean.order.correlation[i] * readings.aggregate[i - temp.correlation.order[j],]$temp_over_break,
             temps.geometricmean.order.correlation[i] <- 0)
    }
  }
  
  temps.geometricmean.order.natural <- temps.geometricmean.order.natural ^ (1/n)
  temps.geometricmean.order.correlation <- temps.geometricmean.order.correlation ^ (1/n)
  
  plot(y = readings.aggregate$kwh,
       x = temps.geometricmean.order.natural,
       col = rgb(0, 0, 0, 20, maxColorValue=100),
       pch = 16,
       main = paste0("Geometric Mean of Past Degrees>Breakpoint (n=", n, ")"),
       ylab = "Average Household Demand (kWh)",
       xlab = "Geometric Mean of Past Degrees>Breakpoint")
  
  plot(y = readings.aggregate$kwh,
       x = temps.geometricmean.order.correlation,
       col = rgb(0, 0, 0, 20, maxColorValue=100),
       pch = 16,
       main = paste0("Geometric Mean of Past Degrees>Breakpoint, Ordered by Correlation (n=", n, ")"),
       ylab = "Average Household Demand (kWh)",
       xlab = "Geometric Mean of Past Degrees>Breakpoint, Ordered According to Correlation with Residuals")
}

# Aritmetic mean of past temperatures
for(n in 1:5) {
  temps.arithmean.order.natural <- rep(0, nrow(readings.aggregate))
  temps.arithmean.order.correlation <- rep(0, nrow(readings.aggregate))
  temp.correlation.order <- c(2, 3, 1, 4, 0)
  for(j in 1:n) {
    for(i in 1:nrow(readings.aggregate)) {
      # If index of past reading is greater than 0, add it
      if (i-j+1 > 0) {
             temps.arithmean.order.natural[i] <- temps.arithmean.order.natural[i] + readings.aggregate[i-j+1,]$temperature
      }
      
      # If index of past reading is greater than 0, add it
      if (i-temp.correlation.order[j] > 0) {
             temps.arithmean.order.correlation[i] <- temps.arithmean.order.correlation[i] + readings.aggregate[i - temp.correlation.order[j],]$temperature
      }
    }
  }
  
  temps.arithmean.order.natural <- temps.arithmean.order.natural / n
  temps.arithmean.order.correlation <- temps.arithmean.order.correlation / n
  
  plot(y = readings.aggregate$kwh,
       x = temps.arithmean.order.natural,
       col = rgb(0, 0, 0, 20, maxColorValue=100),
       pch = 16,
       main = paste0("Arithmetic Mean of Past Temperature (n=", n, ")"),
       ylab = "Average Household Demand (kWh)",
       xlab = "Arithmetic Mean of Past Temperature")
  
  plot(y = readings.aggregate$kwh,
       x = temps.arithmean.order.correlation,
       col = rgb(0, 0, 0, 20, maxColorValue=100),
       pch = 16,
       main = paste0("Arithmetic Mean of Past Temperature, Ordered by Correlation (n=", n, ")"),
       ylab = "Average Household Demand (kWh)",
       xlab = "Arithmetic Mean of Past Temperature, Ordered According to Correlation with Residuals")
}



##
# Find the optimal number of hours lag, and the best method for incorporating 
# CDH. Is it best to sum up current and past hours, similar to traditional CDH 
# or should previous hours be nested under the current hour?
# (Commenting out the iterative comparison for now)
#PerformTouCdhGlmIterations(df.readings = readings.aggregate, nhrs = 5)


# Iterate backward stepwise linear regression through all combinations of 
# past hours' temperature > breakpoint.
#df.stepresults <- InitStepresultsDataframe()
# Iterate through all steps of the current nlag
#for(i in 0:24) {
#  df.bkwdseries <- BackwardStepwiseLinearRegression(nlags = i,
#                                                    df.obs = readings.trimmed,
#                                                    formula.maximal = formula("(kwh) ~ (.)^2"))
#  df.stepresults <- rbind(df.stepresults, df.bkwdseries)
#}
#Plot3DStepwiseResults(df = df.stepresults,
#                      title = "Change in BIC During Stepwise Removal of Terms from Maximal Model")

# Iterate forward stepwise linear regression through all combinations of 
# past hours' temperature > breakpoint.
#df.stepresults <- InitStepresultsDataframe()
#for (i in 0:24) {
#  df.fwdseries <- ForewardStepwiseLinearRegression(nlags = i, 
#                                                   df.obs = readings.trimmed, 
#                                                   formula.maximal = formula("(kwh) ~ (.)^2"))
#  df.stepresults <- rbind(df.stepresults, df.fwdseries)
#}
#Plot3DStepwiseResults(df = df.stepresults,
#                      title = "Change in BIC During Stepwise Addition of Terms to GLM")


# Create 2D version of stepwise linear regression results.
# (Commenting out for now)
#df.results.stepforward.trimmed <- subset(df.stepresults, num.cdhlags %% 2 == 0)
#
#Plot2DFitByExplVarCountWithMultiplePastHrsTemp(df.steps = df.stepresults,
#                                               is.bic = FALSE,
#                                               title = expression(paste("Pseudo-", 
#                                                                        R^2, 
#                                                                        " Change as Terms are Removed from Maximal Models")),
#                                               subtitle = "(p-value<0.05 stopping criterion)")
#
#Plot2DFitByExplVarCountWithMultiplePastHrsTemp(df.steps = df.stepresults,
#                                               is.bic = TRUE,
#                                               title = "BIC Change as Terms are Removed from Maximal Models",
#                                               subtitle = "(p-value<0.05 stopping criterion)")

# Some data manipulation before selecting a model with bestglm
temp.hrs = 8
readings.aggregate <- AddAllCategoricalInteractions(df = readings.aggregate)
temps <- CreatePastTemperatureMatrix(nlags = temp.hrs, 
                                     df.readings = readings.aggregate)

# Initial trim of readings
readings.trimmed <- TrimColsTouTimeComponents(readings.aggregate)
FormulaStringNoInterTemperatureInteractions <- function(df.readings, 
                                                        matrix.temps,
                                                        incl.y = TRUE) {
  # Two vectors of variable names
  catvars <- c("month", "hrstr", "weekend", "price")
  tempvars <- colnames(matrix.temps)
  
  # Create main effect portion of string
  frmla.str <- paste(cbind(paste(catvars, collapse = " + "),
                           paste(tempvars, collapse = " + ")),
                     collapse = " + ")
  
  # Create interactions between categorical variables
  frmla.str <- paste(c(frmla.str, "mnth_hrstr", "mnth_wknd", "mnth_price", 
                       "hrstr_wknd", "hrstr_price", "wknd_price"),
                     collapse = " + ")
  
  # Interact each categorical term with temperature, but do not interact 
  # temperature terms with one another.
  for (i in 1:length(catvars)) {
    frmla.str <- paste(cbind(c(frmla.str),
                             paste(catvars[i], tempvars, collapse = " + ", sep = ":")),
                       collapse = " + ")
  }
  
  # logtransform the response variable
  if (incl.y == TRUE) {
    frmla.str <- paste(c("log(kwh)", frmla.str),
                       collapse = " ~ ")
  }
  
  return(frmla.str)
}
frmla.str.no.intertemp.interac <- FormulaStringNoInterTemperatureInteractions(df.readings = readings.trimmed,
                                                                              matrix.temps = temps)

FormulaStringOnlyCategoricalInteractions<- function(df.readings, 
                                                    matrix.temps,
                                                    incl.y = TRUE) {
  # Two vectors of variable names
  catvars <- colnames(df.readings)
  catvars <- catvars[-1] # trim kwh
  tempvars <- colnames(matrix.temps)
  
  # Create main effect portion of string
  frmla.str <- paste(cbind(paste(catvars, collapse = " + "),
                           paste(tempvars, collapse = " + ")),
                     collapse = " + ")
  
  # logtransform the response variable
  if (incl.y == TRUE) {
    frmla.str <- paste(c("log(kwh)", frmla.str),
                       collapse = " ~ ")
  }
  return(frmla.str)
}
frmla.str.only.category.interac <- FormulaStringOnlyCategoricalInteractions(df.readings = readings.trimmed,
                                                                            matrix.temps = temps)

FormulaStringOnlyMainEffects <- function(df.readings, 
                                         matrix.temps,
                                         incl.y = TRUE) {
  # Two vectors of variable names
  readings.colnames <- colnames(df.readings)
  catvars <- subset(readings.colnames, 
                    readings.colnames %in% c("month", "hrstr", "weekend", "price"))
  tempvars <- colnames(matrix.temps)
  
  # Create main effect portion of string
  frmla.str <- paste(cbind(paste(catvars, collapse = " + "),
                           paste(tempvars, collapse = " + ")),
                     collapse = " + ")
  
  # logtransform the response variable
  if (incl.y == TRUE) {
    frmla.str <- paste(c("log(kwh)", frmla.str),
                       collapse = " ~ ")
  }
  
  return(frmla.str)
}
frmla.str.only.main.effects <- FormulaStringOnlyMainEffects(df.readings = readings.trimmed,
                                                            matrix.temps = temps)


#readings.aggregate <- cbind(readings.aggregate, temps)
#readings.trimmed <- TrimColsTouTimeComponents(readings.aggregate)
readings.trimmed <- cbind(readings.trimmed, temps)
#explvar.str <- TempLagMaximalExplVars(readings.trimmed)
#formula.maximal <- formula(paste0("log(kwh) ~ ", explvar.str))
frmla.no.intertemp.interac <- formula(frmla.str.no.intertemp.interac)
frmla.only.category.interac <- formula(frmla.str.only.category.interac)
frmla.only.main.effects <- formula(frmla.str.only.main.effects)

# Create a design matrix for reuse in several functions
kwh.logtransformed <- log(readings.trimmed$kwh)
only.category.interac.designmatrix <- model.matrix(object = frmla.only.category.interac,
                                                   data = readings.trimmed)

# Take advantage of LM and use regsubset
# readings.regsubsets <- regsubsets(x = frmla.only.category.interac,
#                                   data = readings.trimmed,
#                                   nbest = 1, # 1 best model for each number of predictors
#                                   method = "exhaustive",
#                                   really.big = TRUE)
# regsubsets.summary <- summary(readings.regsubsets)
# regsubsets.which.adjr2.max <- which.max(regsubsets.summary$adjr2)
# regsubsets.summary$which[regsubsets.which.adjr2.max,]
# 
# readings.lm <- lm(log(kwh) ~ month + hrstr + weekend + price + templag0 + templag1 + templag2 + templag3 + templag4 + templag5 + templag6 + templag8 + templag12,
#                   data = readings.trimmed)
# summary(readings.lm)

# LEAPs errors out?
# readings.X <- as.matrix(readings.trimmed[,-1])
# y <- log(readings.trimmed$kwh)
# readings.Xy <- as.data.frame(cbind(readings.X, y))
# readings.bestglm <- bestglm(Xy = readings.Xy,
#                             family = gaussian,
#                             IC = "BIC")

# BMA
# frmla.bma <- update(frmla.only.main.effects,
#                     . ~ . + month:hrstr + month:weekend + month:price)
# readings.bma <- bic.glm(f = frmla.bma,
#                         data = readings.trimmed,
#                         glm.family = gaussian,
#                         factor.type = TRUE)
# summary(readings.bma)
# imageplot.bma(readings.bma)

# Example LEAPs following the LISA short course
temps.scaled <- scale(temps) # Use scale(...) to z-score

# Visualize after z-score
#pairs(cbind(kwh.logtransformed, temps.scaled[,1:5]))
lasso1 <- lars(x = temps.scaled, 
               y = kwh.logtransformed,
               type = 'lasso',
               trace = TRUE,
               normalize = TRUE,
               intercept = TRUE)
plot(lasso1)
lasso1
round(coef(lasso1), 4)

# Order that terms were added
for (i in 1:length(lasso1$actions)) {
  print(names(lasso1$actions[[i]]))
}

predict.lars(object = lasso1, 
             s = .375,
             mode = "fraction",
             type = "coefficients")
predict.lars(object = lasso1,
             newx = temps.scaled, 
             s = .375,
             mode = "fraction",
             type = "fit")

PlotLasso(lars.obj = lasso1, y = kwh.logtransformed, design.mat = temps.scaled,
          xvar = "degf")

PlotLassoCrossValidation(design.mat = temps.scaled, y.vec = kwh.logtransformed, 
                         k = 10, backtransform.mse = "log", xvar = "step")


# 2nd iteration of LASSO method which includes contrasts of categorical 
# variables.
readings.categorical <- readings.trimmed[, colnames(readings.trimmed) %in% c("month",
                                                                             "hrstr",
                                                                             "weekend",
                                                                             "price",
                                                                             "hrstr_price",
                                                                             "wknd_price",
                                                                             "mnth_hrstr",
                                                                             "mnth_wknd",
                                                                             "mnth_price",
                                                                             "hrstr_wknd")]
# Design matrix representing categorical factors and categorical interactions
categorical.design <- model.matrix(~ ., 
                                   data=readings.categorical, 
                                   contrasts.arg = lapply(readings.categorical[,sapply(readings.categorical, is.factor)], 
                                                          contrasts, 
                                                          contrasts=FALSE))
# Trim incercept off (LASSO computes intercept itself)
categorical.design <- categorical.design[,-1]

# Stitch together categorical and continuous (scaled) variables for LASSO
readings.lasso2 <- cbind(categorical.design, temps.scaled)
lasso2 <- lars(x = readings.lasso2, 
               y = kwh.logtransformed,
               type = 'lasso',
               trace = FALSE,
               normalize = TRUE,
               intercept = TRUE)

# 3rd iteraction of LASSO including categorical variable interactions with temp
trimcat.scaletemp <- cbind(readings.categorical, as.data.frame(temps.scaled))
twoway.nointertemp.str <- FormulaStringNoInterTemperatureInteractions(readings.categorical, 
                                                                      temps.scaled,
                                                                      incl.y = FALSE)

# Design matrix representing main effects, categorical interactions, and temp
# interaction with categorical factors
twoway.nointertemp.str <- paste("~", twoway.nointertemp.str, collapse = " ")
twoway.nointertemp.design <- model.matrix(formula(twoway.nointertemp.str),
                                          data = trimcat.scaletemp,
                                          contrasts.arg = lapply(trimcat.scaletemp[,sapply(trimcat.scaletemp, is.factor)], 
                                                                 contrasts, 
                                                                 contrasts=FALSE))
twoway.nointertemp.design <- twoway.nointertemp.design[,-1] # Trim off intercept

# Stitch together categorical and continuous (scaled) variables for LASSO
lasso3 <- lars(x = twoway.nointertemp.design, 
               y = kwh.logtransformed,
               type = 'lasso',
               trace = FALSE,
               normalize = TRUE,
               intercept = TRUE)


# 4th example of LASSO, run with only pre-TOU data and no price factor
pretou.trimmed <- subset(readings.trimmed, price == "flat")
pretou.categorical.main <- pretou.trimmed[, colnames(pretou.trimmed) %in% c("month",
                                                                            "hrstr",
                                                                            "weekend")]
temps.pretou <- pretou.trimmed[, ! colnames(pretou.trimmed) %in% c("kwh", 
                                                                   "month",
                                                                   "hrstr",
                                                                   "weekend",
                                                                   "price",
                                                                   "hrstr_price",
                                                                   "wknd_price",
                                                                   "mnth_hrstr",
                                                                   "mnth_wknd",
                                                                   "mnth_price",
                                                                   "hrstr_wknd")]
temps.pretou.scaled <- scale(temps.pretou)
kwh.pretou.logtransformed <- log(pretou.trimmed$kwh)

pretou.main.str <- FormulaStringOnlyMainEffects(df.readings = pretou.categorical.main,
                                                matrix.temps = temps.pretou,
                                                incl.y = FALSE)
pretou.main.str <- paste("~", pretou.main.str, collapse = " ")

pretou.main.scaletemp <- cbind(pretou.categorical.main, as.data.frame(temps.pretou.scaled))
pretou.main.design <- model.matrix(formula(pretou.main.str),
                                   data = pretou.main.scaletemp,
                                   contrasts.arg = lapply(pretou.main.scaletemp[,sapply(pretou.main.scaletemp, is.factor)], 
                                                          contrasts, 
                                                          contrasts=FALSE))
pretou.main.design <- pretou.main.design[,-1] # Trim off intercept

# Stitch together categorical and continuous (scaled) variables for LASSO
step.1se <- PlotLassoCrossValidation(design.mat = pretou.main.design, 
                                     y.vec = kwh.pretou.logtransformed, 
                                     k = 10, 
                                     backtransform.mse = "log", 
                                     rmse = TRUE,
                                     xvar = "step")

lasso4 <- lars(x = pretou.main.design, 
               y = kwh.pretou.logtransformed,
               type = 'lasso',
               trace = TRUE,
               normalize = TRUE,
               intercept = TRUE)

# degf are 1 higher than LASSO step due to intercept
degf.1se <- unlist(step.1se, use.names = FALSE) + 1 
PlotLasso(lars.obj = lasso4, y = kwh.pretou.logtransformed, 
          design.mat = pretou.main.design, xvar = "degf", 
          term.labels = "catleg", line.marker = "1se", xvar.parsimonious = degf.1se)

# lasso4.cv <- cv.glmnet(x = pretou.main.design, 
#                        y = kwh.pretou.logtransformed,
#                        family = "gaussian")
# plot(x = lasso4.cv)
# 
# 
# lasso4.glmnet <- glmnet(x = pretou.main.design, 
#                         y = kwh.pretou.logtransformed,
#                         family = "gaussian")
# plot(x = lasso4.glmnet,
#      xvar = "dev",
#      label = TRUE)


lasso5 <- lars(x = pretou.main.design, 
               y = kwh.pretou.logtransformed,
               type = 'lasso',
               trace = TRUE,
               normalize = TRUE,
               intercept = TRUE,
               max.steps = unlist(step.1se, use.names = FALSE))

PlotLasso(lars.obj = lasso5, y = kwh.pretou.logtransformed, 
          design.mat = pretou.main.design, xvar = "degf", 
          term.labels = "catleg", line.marker = "1se", xvar.parsimonious = degf.1se)


# CV over more data
lasso6.step.1se <- PlotLassoCrossValidation(design.mat = twoway.nointertemp.design, 
                                            y.vec = kwh.logtransformed,
                                            k = 10, 
                                            backtransform.mse = "log", 
                                            rmse = TRUE,
                                            xvar = "step")

lasso6 <- lars(x = twoway.nointertemp.design, 
               y = kwh.logtransformed,
               type = 'lasso',
               trace = TRUE,
               normalize = TRUE,
               intercept = TRUE,
               max.steps = unlist(lasso6.step.1se, use.names = FALSE))

for (i in 1:length(lasso6$actions)) {
  print(names(lasso6$actions[[i]]))
}

PlotLasso(lars.obj = lasso6, y = kwh.logtransformed, 
          design.mat = twoway.nointertemp.design, xvar = "degf", 
          term.labels = "catleg", line.marker = "1se", xvar.parsimonious = lasso6.step.1se)