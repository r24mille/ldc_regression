library(lme4)

# Load SmartMeterReading data from CSV
home <- setwd(Sys.getenv("HOME"))
fpath <- file.path(home, 
                   "../Dropbox/ISS4E/R/25_meterid_sample/", 
                   "25_random_sample.csv")
readings <- read.csv(fpath)

# Re-orders TOU Period levels so that graphs are sorted accordingly
readings$tou_period <- factor(readings$tou_period, c("off_weekend", "off_morning", "mid_morning", "on_peak", "mid_evening", "off_evening"))
attach(readings)

# Null model is just one parameter, the overall mean 
# (pg. 333 of "The R Book"). 1 is the intercept in 
# regression models, but here it is the overall mean 
# (ie. grand mean).
null_model <- lm(kwh ~ 1, readings)

# The simplest useful equation is:
#
# y = a + bx
#
# A two-parameter model with one parameter for the 
# intercept, a, and another for the slope, b, of the graph 
# of the continuous response variable y against a continuous # explanatory 
# variable x.
#
# Modeled as:
simple_regression <- glm(kwh ~ temperature)

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
pseudorep_model1 <- lme(kwh~temperature*tou_period, 
                       random=~hourindex|subject, 
                       control=lmeControl(opt="optim")
)
pseudorep_model2 <- lme(kwh~temperature*tou_period*billing_active, 
                        random=~hourindex|subject, 
                        control=lmeControl(opt="optim")
)
anova(pseudorep_model1, pseudorep_model2)
summary(pseudorep_model2)
plot(pseudorep_model2,subject~fitted(.))
