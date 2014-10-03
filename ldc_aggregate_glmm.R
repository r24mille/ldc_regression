library(stargazer) # LaTeX tables
library(segmented) # Find linear regression breakpoint
library(BMA) # Compare GLM models
library(sme) # For AICc function

# Source the function in another file
source('cdh_lag_methods.R')
source('tou_time_methods.R')
source('goodness_fit_visualization.R')

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
# Weight terms based on the number of readings involved in the aggregate average
#
# TODO(r24mille): Resid. Dev is lower with weights. Justify this choice.
#                   * Weights really have more to do with dispersion and non-
#                     homogeneity of variance. 
#                   * From ?glm it states, "Non-NULL weights can be used to
#                     indicate that different observations have different
#                     dispersions (with the values in weights being inversely
#                     proportional to the dispersions); or equivalently, when the
#                     elements of weights are positive integers w_i, that each
#                     response y_i is the mean of w_i unit-weight observations. For
#                     a binomial GLM prior weights are used to give the number of
#                     trials when the response is the proportion of successes: they
#                     would rarely be used for a Poisson GLM."
#                   * The 'number of trials' doesn't quite fit with agg_count the 
#                     way it is written below.
maxagg <- max(readings.aggregate$agg_count)
wghts <- readings.aggregate$agg_count/maxagg

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
#                     not need to be back transformed
# TODO(r24mille): This whole section could be two nicely built functions
nlags <- 12

# 1. Summed CDH lags, TOU as periods
cdhlagsum.touperiods.maxglm.pwr <- matrix(nrow = (nlags + 1),
                                          ncol = 3,
                                          dimnames = list(c(0:nlags),
                                                          c("ResidualDeviance", 
                                                            "AICc",
                                                            "BIC")))
for(i in 0:nlags) {
  cdhlagsum <- CreateCdhLagSum(i, readings.aggregate)
  readings.aggregate.cdhlagsum <- cbind(readings.aggregate, 
                                        cdhlagsum)
  cdhlagsum.touperiods <- TrimColsTouPeriods(readings.aggregate.cdhlagsum)
  cdhlagsum.touperiods.maxfmla <- CdhLagMaximalFormula()
  cdhlagsum.touperiods.maxglm <- glm(formula = cdhlagsum.touperiods.maxfmla, 
                                     data = cdhlagsum.touperiods,
                                     weights = wghts, 
                                     family = Gamma(link="log")) 
  
  cdhlagsum.touperiods.maxglm.pwr[(i+1), 1] <- cdhlagsum.touperiods.maxglm$deviance
  cdhlagsum.touperiods.maxglm.pwr[(i+1), 2] <- AICc(cdhlagsum.touperiods.maxglm)
  cdhlagsum.touperiods.maxglm.pwr[(i+1), 3] <- BIC(cdhlagsum.touperiods.maxglm)
}
PlotGlmFitMeasures(aiccs = cdhlagsum.touperiods.maxglm.pwr[, 2], 
                   bics = cdhlagsum.touperiods.maxglm.pwr[, 3], 
                   resdevs = cdhlagsum.touperiods.maxglm.pwr[, 1], 
                   xvals = c(0:nlags), 
                   xtitle = "Number of Past Hours Included", 
                   title = expression(paste("TOU Period, Degrees>", 
                                            CDH['break'], 
                                            " Summed")))

# 2. Matrix of CDH lags as nested interactions, TOU as periods
cdhlagmat.touperiods.nestedglm.pwr <- matrix(nrow = (nlags + 1),
                                          ncol = 3,
                                          dimnames = list(c(0:nlags),
                                                          c("ResidualDeviance", 
                                                            "AICc",
                                                            "BIC")))
for(i in 0:nlags) {
  cdhlagmat <- CreateCdhLagMatrix(i, readings.aggregate)
  readings.aggregate.cdhlagmat <- cbind(readings.aggregate, 
                                        cdhlagmat)
  cdhlagmat.touperiods <- TrimColsTouPeriods(readings.aggregate.cdhlagmat)
  cdhlagmat.touperiods.nestedfmla <- CdhLagMaximalNestedFormula(cdhlagmat.touperiods,
                                                                colnames(cdhlagmat))
  cdhlagmat.touperiods.nestedglm <- glm(formula = cdhlagmat.touperiods.nestedfmla, 
                                     data = cdhlagmat.touperiods,
                                     weights = wghts, 
                                     family = Gamma(link="log")) 
  
  cdhlagmat.touperiods.nestedglm.pwr[(i+1), 1] <- cdhlagmat.touperiods.nestedglm$deviance
  cdhlagmat.touperiods.nestedglm.pwr[(i+1), 2] <- AICc(cdhlagmat.touperiods.nestedglm)
  cdhlagmat.touperiods.nestedglm.pwr[(i+1), 3] <- BIC(cdhlagmat.touperiods.nestedglm)
}
PlotGlmFitMeasures(aiccs = cdhlagmat.touperiods.nestedglm.pwr[, 2], 
                   bics = cdhlagmat.touperiods.nestedglm.pwr[, 3], 
                   resdevs = cdhlagmat.touperiods.nestedglm.pwr[, 1], 
                   xvals = c(0:nlags), 
                   xtitle = "Number of Past Hours Included", 
                   title = expression(paste("TOU Period, Degrees>", 
                                            CDH['break'], 
                                            " as Coefficients w/ Nested Interaction")))

# 3. Matrix of CDH lags as two-way interactions, TOU as periods
cdhlagmat.touperiods.maxglm.pwr <- matrix(nrow = (nlags + 1),
                                          ncol = 3,
                                          dimnames = list(c(0:nlags),
                                                          c("ResidualDeviance", 
                                                            "AICc",
                                                            "BIC")))
for(i in 0:nlags) {
  cdhlagmat <- CreateCdhLagMatrix(i, readings.aggregate)
  readings.aggregate.cdhlagmat <- cbind(readings.aggregate, 
                                        cdhlagmat)
  cdhlagmat.touperiods <- TrimColsTouPeriods(readings.aggregate.cdhlagmat)
  cdhlagmat.touperiods.maxfmla <- CdhLagMaximalFormula()
  cdhlagmat.touperiods.maxglm <- glm(formula = cdhlagmat.touperiods.maxfmla, 
                                          data = cdhlagmat.touperiods,
                                          weights = wghts, 
                                          family = Gamma(link="log")) 
  
  cdhlagmat.touperiods.maxglm.pwr[(i+1), 1] <- cdhlagmat.touperiods.maxglm$deviance
  cdhlagmat.touperiods.maxglm.pwr[(i+1), 2] <- AICc(cdhlagmat.touperiods.maxglm)
  cdhlagmat.touperiods.maxglm.pwr[(i+1), 3] <- BIC(cdhlagmat.touperiods.maxglm)
}
PlotGlmFitMeasures(aiccs = cdhlagmat.touperiods.maxglm.pwr[, 2], 
                   bics = cdhlagmat.touperiods.maxglm.pwr[, 3], 
                   resdevs = cdhlagmat.touperiods.maxglm.pwr[, 1], 
                   xvals = c(0:nlags), 
                   xtitle = "Number of Past Hours Included", 
                   title = expression(paste("TOU Period, Degrees>", 
                                            CDH['break'], 
                                            " as Coefficients w/ All 2-Way Interactions")))

# 4. Summed CDH lags, TOU components of time
cdhlagsum.toucomps.maxglm.pwr <- matrix(nrow = (nlags + 1),
                                          ncol = 3,
                                          dimnames = list(c(0:nlags),
                                                          c("ResidualDeviance", 
                                                            "AICc",
                                                            "BIC")))
for(i in 0:nlags) {
  cdhlagsum <- CreateCdhLagSum(i, readings.aggregate)
  readings.aggregate.cdhlagsum <- cbind(readings.aggregate, 
                                        cdhlagsum)
  cdhlagsum.toucomps <- TrimColsTouTimeComponents(readings.aggregate.cdhlagsum)
  cdhlagsum.toucomps.maxfmla <- CdhLagMaximalFormula()
  cdhlagsum.toucomps.maxglm <- glm(formula = cdhlagsum.toucomps.maxfmla, 
                                     data = cdhlagsum.toucomps,
                                     weights = wghts, 
                                     family = Gamma(link="log")) 
  
  cdhlagsum.toucomps.maxglm.pwr[(i+1), 1] <- cdhlagsum.toucomps.maxglm$deviance
  cdhlagsum.toucomps.maxglm.pwr[(i+1), 2] <- AICc(cdhlagsum.toucomps.maxglm)
  cdhlagsum.toucomps.maxglm.pwr[(i+1), 3] <- BIC(cdhlagsum.toucomps.maxglm)
}
PlotGlmFitMeasures(aiccs = cdhlagsum.toucomps.maxglm.pwr[, 2], 
                   bics = cdhlagsum.toucomps.maxglm.pwr[, 3], 
                   resdevs = cdhlagsum.toucomps.maxglm.pwr[, 1], 
                   xvals = c(0:nlags), 
                   xtitle = "Number of Past Hours Included", 
                   title = expression(paste("TOU Components, Degrees>", 
                                            CDH['break'], 
                                            " Summed")))

# 5. Matrix of CDH lags as nested interactions, TOU components of time
cdhlagmat.toucomps.nestedglm.pwr <- matrix(nrow = (nlags + 1),
                                          ncol = 3,
                                          dimnames = list(c(0:nlags),
                                                          c("ResidualDeviance", 
                                                            "AICc",
                                                            "BIC")))
for(i in 0:nlags) {
  cdhlagmat <- CreateCdhLagMatrix(i, readings.aggregate)
  readings.aggregate.cdhlagmat <- cbind(readings.aggregate, 
                                        cdhlagmat)
  cdhlagmat.toucomps <- TrimColsTouPeriods(readings.aggregate.cdhlagmat)
  cdhlagmat.toucomps.nestedfmla <- CdhLagMaximalNestedFormula(cdhlagmat.toucomps,
                                                           colnames(cdhlagmat))
  cdhlagmat.toucomps.nestedglm <- glm(formula = cdhlagmat.toucomps.nestedfmla, 
                                     data = cdhlagmat.toucomps,
                                     weights = wghts, 
                                     family = Gamma(link="log")) 
  
  cdhlagmat.toucomps.nestedglm.pwr[(i+1), 1] <- cdhlagmat.toucomps.nestedglm$deviance
  cdhlagmat.toucomps.nestedglm.pwr[(i+1), 2] <- AICc(cdhlagmat.toucomps.nestedglm)
  cdhlagmat.toucomps.nestedglm.pwr[(i+1), 3] <- BIC(cdhlagmat.toucomps.nestedglm)
}
PlotGlmFitMeasures(aiccs = cdhlagmat.toucomps.nestedglm.pwr[, 2], 
                   bics = cdhlagmat.toucomps.nestedglm.pwr[, 3], 
                   resdevs = cdhlagmat.toucomps.nestedglm.pwr[, 1], 
                   xvals = c(0:nlags), 
                   xtitle = "Number of Past Hours Included", 
                   title = expression(paste("TOU Components, Degrees>", 
                                            CDH['break'], 
                                            " as Coefficients w/ Nested Interaction")))

# 6. Matrix of CDH lags as two-way interactions, TOU components of time
cdhlagmat.toucomps.maxglm.pwr <- matrix(nrow = (nlags + 1),
                                          ncol = 3,
                                          dimnames = list(c(0:nlags),
                                                          c("ResidualDeviance", 
                                                            "AICc",
                                                            "BIC")))
for(i in 0:nlags) {
  cdhlagmat <- CreateCdhLagMatrix(i, readings.aggregate)
  readings.aggregate.cdhlagmat <- cbind(readings.aggregate, 
                                        cdhlagmat)
  cdhlagmat.toucomps <- TrimColsTouPeriods(readings.aggregate.cdhlagmat)
  cdhlagmat.toucomps.maxfmla <- CdhLagMaximalFormula()
  cdhlagmat.toucomps.maxglm <- glm(formula = cdhlagmat.toucomps.maxfmla, 
                                     data = cdhlagmat.toucomps,
                                     weights = wghts, 
                                     family = Gamma(link="log")) 
  
  cdhlagmat.toucomps.maxglm.pwr[(i+1), 1] <- cdhlagmat.toucomps.maxglm$deviance
  cdhlagmat.toucomps.maxglm.pwr[(i+1), 2] <- AICc(cdhlagmat.toucomps.maxglm)
  cdhlagmat.toucomps.maxglm.pwr[(i+1), 3] <- BIC(cdhlagmat.toucomps.maxglm)
}
PlotGlmFitMeasures(aiccs = cdhlagmat.toucomps.maxglm.pwr[, 2], 
                   bics = cdhlagmat.toucomps.maxglm.pwr[, 3], 
                   resdevs = cdhlagmat.toucomps.maxglm.pwr[, 1], 
                   xvals = c(0:nlags), 
                   xtitle = "Number of Past Hours Included", 
                   title = expression(paste("TOU Components, Degrees>", 
                                            CDH['break'], 
                                            " as Coefficients w/ All 2-Way Interactions")))
