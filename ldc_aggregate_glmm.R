library(segmented) # Find linear regression breakpoint
library(sme) # For AICc function
library(BMA) # Bayesian Model Averaging
library(bestglm) # For bestglm selection using leaps
library(lars) # A separate lasso package
library(glmnet) # Explanatory variable selection by lasso
library(boot) # Cross-validation
library(lars) # Different LASSO library

# Experiment with devtools and fastVAR
#library(devtools)
#install_github("fastVAR", "jeffwong", "master")
library(fastVAR)

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
temperature.break <- floor(seg$psi[1,2]) # TODO floor vs. round vs. real temps
readings.aggregate$temp_over_break <- ifelse(readings.aggregate$temperature > temperature.break, 
                                             readings.aggregate$temperature - temperature.break, 
                                             0)

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
temp.hrs = 6
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
  catvars <- c("month", "hrstr", "weekend", "price")
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
readings.lasso1 <- as.data.frame(cbind(kwh.logtransformed, temps))
#pairs(readings.lasso1)

# Use scale(...) to z-score
readings.lasso1.scaled <- scale(readings.lasso1[,2:(temp.hrs+1)])

# Visualize after z-score
#pairs(readings.lasso1.scaled)
lasso1 <- lars(x = readings.lasso1.scaled, 
               y = kwh.logtransformed,
               type = 'lasso',
               trace = FALSE,
               normalize = TRUE,
               intercept = TRUE)
plot(lasso1)
lasso1
round(coef(lasso1), 4)

predict.lars(object = lasso1, 
             s = .375,
             mode = "fraction",
             type = "coefficients")
predict.lars(object = lasso1,
             newx = readings.lasso1.scaled, 
             s = .375,
             mode = "fraction",
             type = "fit")

# Examine the ordinary least squares fit of the full model using multiple regression(?)
OLS <- summary(lm(kwh.logtransformed~., data = as.data.frame(readings.lasso1.scaled)))
OLS
absum <- sum(abs(OLS$coeff[-1,1]))

# Manually create the LASSO step plot from the ground up
t <- apply(abs(coef(lasso1)), 1, sum)
s <- t/absum

coef.max = max(coef(lasso1))
coef.min = min(coef(lasso1))
coef.range = abs(coef.max - coef.min)
plot(x = s,
     y = coef(lasso1)[, 1], 
     ylim = c(coef.min, coef.max),
     type = "l",
     lwd = 2,
     xlab = "Shrinkage factor s",
     main = "LASSO path - coefficients as a function of shrinkage factor s",
     xlim = c(0,1.2),
     axes = FALSE,
     ylab = "Coefficient",
     cex.lab = 1.5,
     cex.axis = 1.4) # Plot one line to initiate things
axis(1, at = seq(0, 1, .2), 
     cex.axis = 1.1) # Control x-axis
axis(2, 
     at = seq(coef.min, coef.max, round((coef.range/10), 2)),
     cex.axis = 1.1, 
     las = 1) # Control y-axis

lines(s, coef(lasso1)[,2], lwd = 2)
lines(s, coef(lasso1)[,3], lwd = 2)
lines(s, coef(lasso1)[,4], lwd = 2)
lines(s, coef(lasso1)[,5], lwd = 2)
lines(s, coef(lasso1)[,6], lwd = 2)

text(1.07, 
     coef(lasso1)[nrow(coef(lasso1)), 1],
     colnames(readings.lasso1.scaled)[1])
text(1.07, 
     coef(lasso1)[nrow(coef(lasso1)), 2],
     colnames(readings.lasso1.scaled)[2])
text(1.07, 
     coef(lasso1)[nrow(coef(lasso1)), 3],
     colnames(readings.lasso1.scaled)[3])
text(1.07, 
     coef(lasso1)[nrow(coef(lasso1)), 4],
     colnames(readings.lasso1.scaled)[4])
text(1.07, 
     coef(lasso1)[nrow(coef(lasso1)), 5],
     colnames(readings.lasso1.scaled)[5])
text(1.07, 
     coef(lasso1)[nrow(coef(lasso1)), 6],
     colnames(readings.lasso1.scaled)[6])

abline(v=s,
       col = "darkgray",
       lty = 3)


# Use LASSO for better selection of explanatory variables
# transform dataframe to matrices as required by glmnet
#explvar.matrix <- model.matrix(formula.maximal, readings.trimmed)
frmla.str.without.y <- FormulaStringOnlyMainEffects(df.readings = readings.trimmed,
                                                    matrix.temps = temps,
                                                    incl.y = FALSE)
frmla.str.without.y <- paste(c("~", frmla.str.without.y), collapse = " ")
explvar.sparse.no.intertemp.interac <- model.matrix(object = frmla.str.without.y, 
                                                    data = readings.trimmed)

# Remove intercept, glmnet fits an intercept
explvar.sparse.nointercept <- explvar.sparse.no.intertemp.interac[,-1] 

# Use a penalty factor vector to ensure that the main effects are always kept 
# in the model. Because glmnet requires data in sparse matrix form, it does 
# not know that the earlier variables are main effects, and that main effects 
# are required if an interaction is to be used.
# 
# Vector of 1s by default.
explvar.sparse.pf <- rep(1, ncol(explvar.sparse.nointercept)) 
# Index of templag0 marks start of interactions temperatures, then interactions
idx.templag.start <- match("templag0", 
                                 colnames(explvar.sparse.nointercept))
# Every variable before hrstr_priceh0off_peak is a main effect and should have 
# penalty.factor=0 to ensure that no shrinkage is applied and variable is always
# included in the model.
explvar.sparse.pf[1:(idx.templag.start-1)] <- 0

# Cross-validate to select a value for lambda
readings.glmnet.cv <- cv.glmnet(x = explvar.sparse.nointercept, 
                                y = depvar.matrix,
                                family = "gaussian",
                                penalty.factor = explvar.sparse.pf)

print(coef(readings.glmnet.cv, s="lambda.1se"))
plot(readings.glmnet.cv)
lambda.1se.idx <- match(readings.glmnet.cv$lambda.1se,
                        readings.glmnet.cv$lambda)
