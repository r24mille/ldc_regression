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
# fixed = kwh~1

# The fixed effect is just the overall mean value of the response variable. The 
# random effects show the indentities of the random variables and their relative 
# locations in the hierarchy. The random effects are specified like:
# random = ~ 1 | subject_reindexed

# An important detail to notice is that the name of the response variable (kwh) 
# is not repeated in the random-effects forumula: there is a blank space to the 
# left of the tilde. In most mixed-effects models we assume that the random 
# effects have a mean of zero and that we are interested in quantifying 
# variation in the intercept (ie. parameter 1) caused by differences between 
# the factor levels of the random effects.
#
# After the intercept comes the vertical bar | which is read as "given the 
# following spatial arrangement of the random variables". In this case, there 
# is ony random effect, subject_reindexed.
#
# TODO: random <- ~ 1 | subject_reindexed/daynum/hour
# TODO: random <- ~ 1 | subject_reindexed/hourindex
no_fixed_sub_only <- lme(fixed = kwh~1,
                 random = ~ 1 | subject_reindexed)
no_fixed_sub_nested_dayhour <- lme(fixed = kwh~1,
                         random = ~ 1 | subject_reindexed/daynum/hour,
                         control=lmeControl(opt="optim"))
no_fixed_sub_hourindex <- lme(fixed = kwh~1,
                         random = ~ 1 | subject_reindexed/hourindex,
                         control=lmeControl(opt="optim"))
no_fixed_sub_daynum <- lme(fixed = kwh~1,
                               random = ~ 1 | subject_reindexed/daynum,
                               control=lmeControl(opt="optim"))
no_fixed_day_nested_sub <- lme(fixed = kwh~1,
                           random = ~ 1 | daynum/subject_reindexed,
                           control=lmeControl(opt="optim"))

summary(no_fixed_sub_only)
summary(no_fixed_sub_nested_dayhour)
summary(no_fixed_sub_hourindex)
summary(no_fixed_sub_daynum)

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
                       control=lmeControl(opt="optim"))
