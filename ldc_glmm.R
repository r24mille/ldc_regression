library(lme4)
library(lattice) # plotting

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
readings$subject <- match(readings$subject, meterids)

# Null model is just one parameter, the overall mean 
# (pg. 333 of "The R Book"). 1 is the intercept in 
# regression models, but here it is the overall mean 
# (ie. grand mean).
model.null <- lm(kwh ~ 1, readings)

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
# is ony random effect, subject
model.nofixed.subonly <- lme(fixed = kwh~1, 
                             random = ~ 1 | subject,
                             data = readings)
model.nofixed.sub.nested.daynum <- lme(fixed = kwh~1,
                                       random = ~ 1 | subject/daynum,
                                       control=lmeControl(opt="optim"),
                                       data = readings)
#model.nofixed.sub.nested.daynum.hour <- lme(fixed = kwh~1, 
#                                            random = ~ 1 | subject/daynum/hour, 
#                                            control=lmeControl(opt="optim"),
#                                            data = readings)
#model.nofixed.sub.nested.hourindex <- lme(fixed = kwh~1, 
#                                          random = ~ 1 | subject/hourindex, 
#                                          control=lmeControl(opt="optim"), 
#                                          data = readings)
#model.nofixed.daynum.nested.sub <- lme(fixed = kwh~1, 
#                                       random = ~ 1 | daynum/subject, 
#                                       control=lmeControl(opt="optim"), 
#                                       data = readings)

summary(model.nofixed.subonly)
summary(model.nofixed.sub.nested.daynum) # Lowest AIC
#summary(model.nofixed.sub.nested.daynum.hour)
#summary(model.nofixed.sub.nested.hourindex) 
#summary(model.nofixed.daynum.nested.sub)

# Attempt to rewrite the successful model from above using lmer
model.lmer.nofixed.sub.nested.daynum <- lmer(kwh ~ 1 + (1 | subject/daynum),
                                             data = readings)
summary(model.lmer.nofixed.sub.nested.daynum)

# Linear Mixed-Effects Regression
# Subject is a random effect because its factor levels have no meaning 
# (ie. convey no information). However, subjects are sampled from 
# a larger population and differ in many ways but we do not know exactly how.
# The subject being sampled affects the variance of the result, not the 
# mean of the result.
#
# A common cause of temporal pseudoreplication in growth experiments with 
# fixed effects is when each individual is measured several times as it grows 
# during the course of an experiment.
#
# To use trellis plotting, begin by turning the readings dataframe into a  
# groupedData object. Specify the nesting structure of the random effects, and  
# indicate the fixed effect by defining temperature*tou_period*billing_active 
# as "outer" to the nesting (ie. a fixed effect):
#readings <- groupedData(kwh ~ daynum|subject,
#                        outer = ~ billing_active, 
#                        readings)
#plot(readings, outer=T)

# Further refinement is made for temporal pseudoreplication as described 
# on pg. 965 of "The R Book, 2nd Ed". This considers daynum (a continuous 
# random effect) as pseudoreplication within each subject.
# 
# lmeControl added as described in 
# https://stat.ethz.ch/pipermail/r-help/2008-June/164806.html to work around 
# an error.
#model.threefixed.psuedorep.hourindex <- lme(fixed = kwh~temperature*tou_period*billing_active, 
#                                  random = ~hourindex|subject, 
#                                  control = lmeControl(opt="optim"), 
#                                  data = readings)
#model.threefixed.psuedorep.daynum <- lme(fixed = kwh~temperature*tou_period*billing_active, 
#                                         random = ~daynum|subject, 
#                                         control = lmeControl(opt="optim"), 
#                                         data = readings)
model.threefixed.psuedorep.hourindex.nested.daynum <- lme(fixed = kwh~temperature*tou_period*billing_active, 
                                                          random = ~hourindex|subject/daynum, 
                                                          control = lmeControl(opt="optim"), 
                                                          data = readings)


#summary(model.threefixed.psuedorep.hourindex)
#summary(model.threefixed.psuedorep.daynum)
summary(model.threefixed.psuedorep.hourindex.nested.daynum) # Lowest AIC

# Time series analysis in mixed-effects models (pg. 699)
# It 
