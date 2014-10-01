library(stargazer) # LaTeX tables
library(segmented) # Find linear regression breakpoint
library(BMA) # Compare GLM models

# Source the function in another file
source('reg_subset_visualization.R')
source('cdh_lag_methods.R')
source('tou_time_methods.R')

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
model.readings.glm.presegment <- glm(kwh ~ temperature*tou_period*billing_active,
                                   data = readings.aggregate,
                                   family = Gamma)
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
nlags <- 4
cdhlagmat.touperiods.maxglm.pwr <- matrix(nrow = (nlags + 1),
                                          ncol = 2,
                                          dimnames = list(c(0:nlags),
                                                          c("ResidualDeviance", 
                                                            "AIC")))
for(i in 0:nlags) {
  # Matrix of CDH lags, TOU as periods
  cdhlagmat <- CreateCdhLagMatrix(i, readings.aggregate)
  readings.aggregate.cdhlagmat <- cbind(readings.aggregate, 
                                        cdhlagmat)
  cdhlagmat.touperiods <- TrimColsTouPeriods(readings.aggregate.cdhlagmat)
  cdhlagmat.touperiods.maxfmla <- CdhLagMaximalFormula()
  # maxagg <- max(readings.aggregate.cdhlagmat$agg_count)
  # wghts <- readings.aggregate.cdhlagmat$agg_count/maxagg
  cdhlagmat.touperiods.maxglm <- glm(formula = cdhlagmat.touperiods.maxfmla, 
                                          data = cdhlagmat.touperiods,
                                          # weights = wghts, 
                                          family = Gamma)
  
  cdhlagmat.touperiods.maxglm.pwr[(i+1), 1] <- cdhlagmat.touperiods.maxglm$deviance
  cdhlagmat.touperiods.maxglm.pwr[(i+1), 2] <- cdhlagmat.touperiods.maxglm$aic
}