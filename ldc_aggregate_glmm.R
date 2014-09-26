library(stargazer) # LaTeX tables
library(segmented)
library(BMA)


# Source the function in another file
source('reg_subset_visualization.R')
source('createCdhLag.R')

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
# Add column which represents "month" as a categorical factor
readings.aggregate$timestamp_dst <- as.POSIXlt(readings.aggregate$timestamp_dst)
readings.aggregate$month <- paste0("m", (readings.aggregate$timestamp_dst$mon + 1))
readings.aggregate$month <- factor(readings.aggregate$month, 
                                   c("m5", "m6", "m7", "m8", "m9", "m10"))

##
# Use 'segmented' package rather than my prior home-grown method of finding 
# the optimal cooling degree hour breakpoint (though they give the same 
# result).
model.readings.lm.presegment <- lm(kwh ~ temperature*tou_period*billing_active,
                                   data = readings.aggregate)
seg <- segmented(obj = model.readings.lm.presegment, 
                 seg.Z = ~temperature,
                 psi = list(temperature = c(15)))
cdhbreak <- floor(seg$psi[1,2]) # TODO floor vs. round?
readings.aggregate$cdh <- ifelse(readings.aggregate$temperature > cdhbreak, 
                                 readings.aggregate$temperature - cdhbreak, 
                                 0)


##
# CDH with 1-24 hour lags. 
# TODO: 0-24 lags
# Convert to a coefficient for each hour into the past
##
cdhbics <- rep(0, 24)
for(nlags in 1:24) {
  cdhlag <- createCdhLag(nlags,
                         readings.aggregate)
  readings.aggregate$cdhlag <- cdhlag
  
  glm.timetou.opt.FT <- bic.glm(f = kwh ~ (cdhlag + month + tou_period + billing_active)^2, 
                                data = readings.aggregate, 
                                glm.family = "gaussian", # TODO: I think Gaussian is good
                                factor.type = TRUE)
  cdhbics[nlags + 1] <- glm.timetou.opt.FT$bic
  imageplot.bma(glm.timetou.opt.FT)
}

# Plot the evolution of BIC as different CDH lags are used.
plot(cdhbics,
     xlab = "CDH (hours lag)", 
     ylab = "BIC of Optimal GLM",
     main = "BIC of Models Using Varioud CDH Lags",
     type = "l")
