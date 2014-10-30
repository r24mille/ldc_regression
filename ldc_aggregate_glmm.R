library(segmented) # Find linear regression breakpoint
library(sme) # For AICc function
library(bestglm) # For bestglm selection using leaps
library(lars) # A separate lasso package
library(glmnet) # Explanatory variable selection by lasso

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
model.readings.glm.presegment <- glm(log(kwh)~(month+hrstr+weekend+price+temperature)^2 - hrstr:price - weekend:price,
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
temps <- CreatePastTemperatureMatrix(nlags = 6, 
                                     df.readings = readings.aggregate)

# Initial trim of readings
readings.trimmed <- TrimColsTouTimeComponents(readings.aggregate)
FormulaStringNoInterTemperatureInteractions <- function(df.readings, 
                                                        matrix.temps) {
  # Two vectors of variable names
  catvars <- colnames(df.readings[, ! colnames(df.readings) %in% c("kwh",
                                                                   "hrstr_price",
                                                                   "wknd_price")])
  tempvars <- colnames(matrix.temps)
  
  # Create main effect portion of string
  frmla.str <- paste(cbind(paste(catvars, collapse = " + "),
                           paste(tempvars, collapse = " + ")),
                     collapse = " + ")
  
  # Create interactions between categorical variables
  frmla.str <- paste(c(frmla.str, "month:hrstr", "month:weekend", "month:price", 
                       "hrstr:weekend", "hrstr_price", "wknd_price"),
                     collapse = " + ")
  
  # Interact each categorical term with temperature, but do not interact 
  # temperature terms with one another.
  for (i in 1:length(catvars)) {
    frmla.str <- paste(cbind(c(frmla.str),
                             paste(catvars[i], tempvars, collapse = " + ", sep = ":")),
                       collapse = " + ")
  }
  
  # logtransform the response variable
  frmla.str <- paste(c("log(kwh)", frmla.str),
                     collapse = " ~ ")
  
  return(frmla.str)
}
frmla.notempintr.str <- FormulaStringNoInterTemperatureInteractions(df.readings = readings.trimmed,
                                                                    matrix.temps = temps)

library(BMA)
readings.bma <- bic.glm(data = readings.trimmed,
                        f = formula(frmla.notempintr.str),
                        glm.family = gaussian,
                        nbest = 50)

# LEAPs
#xmatrix <- as.matrix(readings.trimmed[,-1])
#ymatrix <- as.matrix(readings.trimmed[,1])
#designmatrix <- as.data.frame(cbind(explvar.matrix.nointercept, 
#                                    depvar.matrix))
#readings.bestglm <- bestglm(Xy = designmatrix)

# Use LASSO for better selection of explanatory variables
# transform dataframe to matrices as required by glmnet
#readings.aggregate <- cbind(readings.aggregate, temps)
#readings.trimmed <- TrimColsTouTimeComponents(readings.aggregate)
readings.trimmed <- cbind(readings.trimmed, temps)
#explvar.str <- TempLagMaximalExplVars(readings.trimmed)
#formula.maximal <- formula(paste0("log(kwh) ~ ", explvar.str))
formula.notempintr <- formula(frmla.notempintr.str)
#explvar.matrix <- model.matrix(formula.maximal, readings.trimmed)
explvar.matrix <- model.matrix(formula.notempintr, readings.trimmed)
depvar.matrix <- as.matrix(log(readings.trimmed$kwh), ncol=1)

# Remove intercept, glmnet fits an intercept
explvar.matrix.nointercept <- explvar.matrix[,-1] 

# Use a penalty factor vector to ensure that the main effects are always kept 
# in the model. Because glmnet requires data in sparse matrix form, it does 
# not know that the earlier variables are main effects, and that main effects 
# are required if an interaction is to be used.
# 
# Vector of 1s by default.
explvar.matrix.nointercept.pf <- rep(1, ncol(explvar.matrix.nointercept)) 
# Index of templag0 marks start of interactions temperatures, then interactions
idx.ineteractions.start <- match("templag0", 
                                 colnames(explvar.matrix.nointercept))
# Every variable before hrstr_priceh0off_peak is a main effect and should have 
# penalty.factor=0 to ensure that no shrinkage is applied and variable is always
# included in the model.
explvar.matrix.nointercept.pf[1:(idx.ineteractions.start-1)] <- 0

# Cross-validate to select a value for lambda
readings.glmnet.cv <- cv.glmnet(x = explvar.matrix.nointercept, 
                                y = depvar.matrix,
                                family = "gaussian",
                                penalty.factor = explvar.matrix.nointercept.pf)

print(coef(readings.glmnet.cv, s="lambda.1se"))
plot(readings.glmnet.cv)