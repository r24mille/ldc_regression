library(lme4)

# Load SmartMeterReading data from CSV
home <- setwd(Sys.getenv("HOME"))
fpath <- file.path(home, 
                   "../Dropbox/ISS4E/R/25_meterid_sample/", 
                   "25_random_sample.csv")
readings <- read.csv(fpath)

# Re-orders TOU Period levels so that graphs are sorted accordingly
readings$tou_period <- factor(readings$tou_period, 
                              c("off_weekend", "off_morning", "mid_morning", 
                                "on_peak", "mid_evening", "off_evening"))

# Renumbers subjects 1 through n for simplicity.
meterids <- unique(readings$subject)
readings$subject_reindexed <- match(readings$subject, meterids)

# Attach to readings data frame to simplify following code.
attach(readings)

# Null model is just one parameter, the overall mean 
# (pg. 333 of "The R Book"). 1 is the intercept in 
# regression models, but here it is the overall mean 
# (ie. grand mean).
null_model <- lm(kwh ~ 1, readings)

# Suppose that there are no fixed effects, so that all of the 
# categorical variables are random effects. Then the fixed 
# effect simply estimates the intercept (parameter 1):


# Linear Mixed-Effects Regression
# Subject is a random effect because its factor levels have no meaning 
# (ie. convey no information). However, subjects are sampled from 
# a larger population and differ in many ways but we do not know exactly how.
# The subject being sampled affects the variance of the result, not the 
# mean of the result.
#
# Further refinement is made for temporal pseudoreplication as described 
# on pg. 642 of "The R Book". This considers daynum (a continuous random 
# effect) as pseudoreplication within each subject.
#
# lmeControl added as described in 
# https://stat.ethz.ch/pipermail/r-help/2008-June/164806.html to work around 
# an error.
pseudorep_model <- lme(kwh~temperature*tou_period, 
                       random=~hourindex|subject_reindexed, 
                       control=lmeControl(opt="optim")
)
