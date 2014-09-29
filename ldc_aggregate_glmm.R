library(stargazer) # LaTeX tables
library(segmented) # Find linear regression breakpoint
library(BMA) # Compare GLM models

# Source the function in another file
source('reg_subset_visualization.R')
source('cdh_lag_methods.R')

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
                                     unique(readings.aggregate$weekend))

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
model.readings.lm.presegment <- glm(kwh ~ temperature*tou_period*billing_active,
                                   data = readings.aggregate)
seg <- segmented(obj = model.readings.lm.presegment, 
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

# Trim the columns down from full readings.aggregate to only those needed 
# when working with "TOU Period" version of the model.
trimColsToTimeTou <- function(readingdf) {
  timetoudf <- readingdf[, ! colnames(readingdf) %in% c("daynum", 
                                                        "hour", 
                                                        "timestamp_dst", 
                                                        "temperature", 
                                                        "agg_count", 
                                                        "hrstr", 
                                                        "weekend", 
                                                        "price", 
                                                        "cdh")]
  return(timetoudf)
}

# Trim the columns down from full readings.aggregate to only those needed 
# when working with "Components of TOU" version of the model.
trimColsToTimeComponents <- function(readingdf) {
  timecompdf <- readingdf[, ! colnames(readingdf) %in% c("daynum", 
                                                         "hour", 
                                                         "timestamp_dst", 
                                                         "temperature", 
                                                         "agg_count", 
                                                         "tou_period", 
                                                         "cdh")]
  return(timecompdf)
}

# Summed CDH lags can be incorporated into a GLM easily, testing the main 
# effects of the column as well as two-way interractions with each other column.
findBmaCdhLagSum <- function(df) {
  bma.res <- bic.glm(f = kwh ~ .,
                    data = df,
                    glm.family = "gaussian", 
                    factor.type = TRUE)
  return(bma.res)
}

# The matrix of CDH lags (ie. hours in the past) should be modeled as a nested 
# relationship and then incorporated into the main effects of the nested 
# terms and their two-way interractions with each other fixed effect.
#
# This gets a bit tricky to create in a reusable fashion for a changing number 
# of columns (ie. lag0, lag1, ..., lagn). So a string for the nested term 
# is constructed with the help of cdhlagmatcolnames. This string is then 
# stitched as a formula string and then passed to Bayesian Model 
# Averaging (BMA).
findBmaCdhLagMat <- function(df, cdhlagmatcolnames) {
  fecolnames <- colnames(df[, ! colnames(df) %in% c(cdhlagmatcolnames, "kwh")]) # Optimize
  mainstr <- paste(fecolnames,
                   collapse = " + ")
  nestedstr <- paste(cdhlagmatcolnames,
                     collapse = "/")
  nested.formulastr <- paste0("kwh ~ (",
                         mainstr,
                         ")^2 + ",
                         nestedstr)
  
  bma.res <- bic.glm(f = formula(nested.formulastr),
                    data = df,
                    glm.family = "gaussian", 
                    factor.type = TRUE)
  return(bma.res)
}

nlags <- 4
cdhsumbics <- rep(0, (nlags + 1))
cdhmatbics <- rep(0, (nlags + 1))
cdhmatdeviances <- rep(0, (nlags + 1))
for(i in 0:nlags) {
  # Summed CDH for nlag
  cdhlagsums <- createCdhLagSum(i, readings.aggregate)
  readings.aggregate.cdhsum <- cbind(readings.aggregate, 
                                     cdhlagsums)
  readings.aggregate.cdhsum.timetou <- trimColsToTimeTou(readings.aggregate.cdhsum)
  bma.cdhsum.timetou <- findBmaCdhLagSum(readings.aggregate.cdhsum.timetou)
  cdhsumbics[i + 1] <- mean(bma.cdhsum.timetou$bic)
  
  # Matrix of CDH lags
  cdhlagmat <- createCdhLagMatrix(i, readings.aggregate)
  readings.aggregate.cdhmat <- cbind(readings.aggregate, 
                                     cdhlagmat)
  readings.aggregate.cdhmat.timetou <- trimColsToTimeTou(readings.aggregate.cdhmat)
  bma.cdhmat.timetou <- findBmaCdhLagMat(readings.aggregate.cdhmat.timetou, 
                                         colnames(cdhlagmat))
  cdhmatbics[i + 1] <- mean(bma.cdhmat.timetou$bic)
  cdhmatdeviances[i + 1] <- mean(bma.cdhmat.timetou$deviance)
}

# Plot the evolution of BIC as different CDH lags are used.
plot(cdhmatbics,
     xlab = "CDH (hours lag)", 
     ylab = "Average BIC of Best GLMs",
     main = "BIC of Models as the Hours in CDH Increase",
     type = "l")

# Plot the evolution of BIC as different CDH lags are used.
plot(cdhmatdeviances,
     xlab = "CDH (hours lag)", 
     ylab = "Average Deviance of Best GLMs",
     main = "Deviances of Best GLM Models as the Hours in CDH Increase",
     type = "l")
