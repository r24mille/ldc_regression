library(segmented) # Find linear regression breakpoint
library(sme) # For AICc function
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
# TODO(r24mille): segmented chooses very different cdhbreak numbers based on 
#                 the model chosen. Maybe just making an empirical choice for 
#                 CDH break within the context of the final model is the best 
#                 choice.
model.readings.glm.presegment <- glm(log(kwh)~(month+hrstr+weekend+price+temperature)^2 - hrstr:price - weekend:price,
                                     data = readings.aggregate,
                                     family = gaussian)
seg <- segmented(obj = model.readings.glm.presegment, 
                 seg.Z = ~ temperature,
                 psi = list(temperature = c(18)))
cdhbreak <- floor(seg$psi[1,2]) # TODO floor vs. round vs. real temps
readings.aggregate$cdh <- ifelse(readings.aggregate$temperature > cdhbreak, 
                                 readings.aggregate$temperature - cdhbreak, 
                                 0)

##
# Find the optimal number of hours lag, and the best method for incorporating 
# CDH. Is it best to sum up current and past hours, similar to traditional CDH 
# or should previous hours be nested under the current hour?
# (Commenting out the iterative comparison for now)
#PerformTouCdhGlmIterations(df.readings = readings.aggregate,
#                           nhrs = 5)


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