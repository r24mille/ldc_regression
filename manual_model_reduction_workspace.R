##
# Manually construct model informed by iterative add1(...) and drop1(...) 
# results.
cdhlag <- CreateCdhLagMatrix(nlags = 11, readingsdf = readings.aggregate)
readings.with.cdh <- cbind(readings.aggregate, cdhlag)
readings.with.cdh.trimmed <- TrimColsTouTimeComponents(readings.with.cdh)

# With cdhlag=8, the iterative model proposes the following model: 
#        kwh ~ billing_active + month + hrstr + weekend + price + lag0 + lag1 + 
#        lag2 + lag3 +lag6 + lag8 + hrstr_price + billing_active:lag8 + 
#        month:hrstr + month:price + month:lag2 + month:lag8 + hrstr:lag0 + 
#        hrstr:lag3 + hrstr:lag6 + hrstr:lag8 + price:lag1 + lag0:lag3
#
# The past temperatures that were preserved via model reduction are difficult 
# to explain. It seems easier to explain interactions from lag0 through the 
# largest lag interacted with a categorical variable. So I have added all 
# interactions in between. From there I will proceed with manual model reduction 
# stopping according to t-test at p=0.05.
r24mille.alt1.glm <- glm(formula = "kwh ~ billing_active + month + hrstr + weekend + price + lag0 + lag1 + lag2 + lag3 + lag4 + lag5 + lag6 + lag7 + lag8 + hrstr:price + billing_active:lag0 + billing_active:lag1 + billing_active:lag2 + billing_active:lag3 + billing_active:lag4 + billing_active:lag5 + billing_active:lag6 + billing_active:lag7 + billing_active:lag8 + month:hrstr + month:price + month:lag0 + month:lag1 + month:lag2 + month:lag3 + month:lag4 + month:lag5 + month:lag6 + month:lag7 + month:lag8 + hrstr:lag0  + hrstr:lag1  + hrstr:lag2 + hrstr:lag3 + hrstr:lag4 + hrstr:lag5 + hrstr:lag6 + hrstr:lag7 + hrstr:lag8 + price:lag0 + price:lag1 + lag0:lag1 + lag0:lag2 + lag0:lag3 + lag1:lag2 + lag1:lag3 + lag2:lag3", 
                         data = readings.with.cdh.trimmed, 
                         weights = wghts, 
                         family = Gamma(link="log"))


# Start by trying to reduce temperature interacting with itself.
r24mille.alt2.glm <- update(r24mille.alt1.glm, ~ . - lag2:lag3)
anova(r24mille.alt1.glm, r24mille.alt2.glm, test = "LRT")
# lag2:lag3 is significant. I will not reduce temperature interactions further.

# Attempt to reduce price interactions with temperature
r24mille.alt3.glm <- update(r24mille.alt1.glm, ~ . - price:lag1)
anova(r24mille.alt1.glm, r24mille.alt3.glm, test = "LRT") # Accept reduction

r24mille.alt4.glm <- update(r24mille.alt3.glm, ~ . - price:lag0)
anova(r24mille.alt3.glm, r24mille.alt4.glm, test = "LRT")
# price:lag0 is significant. I will not reduce price and temp interactions further.

# Attempt to reduce hour-of-day interacted with temperature
r24mille.alt5.glm <- update(r24mille.alt3.glm, ~ . - hrstr:lag8)
anova(r24mille.alt3.glm, r24mille.alt5.glm, test = "LRT")
# hrstr:lag8 is significant. I will not reduce hour-of-day and temp interactions further.

# Attempt to reduce month interacted with temperature
r24mille.alt6.glm <- update(r24mille.alt3.glm, ~ . - month:lag8)
anova(r24mille.alt3.glm, r24mille.alt6.glm, test = "LRT") # Accept reduction

r24mille.alt7.glm <- update(r24mille.alt6.glm, ~ . - month:lag7)
anova(r24mille.alt6.glm, r24mille.alt7.glm, test = "LRT") # Accept reduction

r24mille.alt8.glm <- update(r24mille.alt7.glm, ~ . - month:lag6)
anova(r24mille.alt7.glm, r24mille.alt8.glm, test = "LRT") # Accept reduction

r24mille.alt9.glm <- update(r24mille.alt8.glm, ~ . - month:lag5)
anova(r24mille.alt8.glm, r24mille.alt9.glm, test = "LRT") 
# month:lag5 is significant. I will not reduce month and temp interactions further.

# Attempt to reduce billing interacted with temperature
r24mille.alt10.glm <- update(r24mille.alt8.glm, ~ . - billing_active:lag8)
anova(r24mille.alt8.glm, r24mille.alt10.glm, test = "LRT")
# billing_active:lag8 is significant. I will not reduce billing and temp interactions further.

BIC(r24mille.alt8.glm) # -10445.12


# Try fewer hours into the past according to the same criteria, that if a 
# past temperature is included all hours in between must be included as well 
# (unless I am confident I can explain/justify a gap).
#
# Model proposed by iterative reduction at cdhlag=7:
#     kwh ~ billing_active + month + hrstr + weekend + price + lag0 + lag2 + 
#     lag3 + lag5 + lag7 + hrstr_price + billing_active:lag7 + month:hrstr + 
#     month:price + month:lag2 + month:lag7 + hrstr:lag0 + hrstr:lag3 + 
#     hrstr:lag5 + hrstr:lag7 + price:lag0 + lag0:lag3
r24mille.alt11.glm <- glm(formula = "kwh ~ billing_active + month + hrstr + weekend + price + lag0 + lag1 + lag2 + lag3 + lag4 + lag5 + lag6 + lag7 + hrstr:price + billing_active:lag0 + billing_active:lag1 + billing_active:lag2 + billing_active:lag3 + billing_active:lag4 + billing_active:lag5 + billing_active:lag6 + billing_active:lag7 + month:hrstr + month:price + month:lag0 + month:lag1 + month:lag2 + month:lag3 + month:lag4 + month:lag5 + month:lag6 + month:lag7 + hrstr:lag0 + hrstr:lag1 + hrstr:lag2 + hrstr:lag3 + hrstr:lag4 + hrstr:lag5 + hrstr:lag6 + hrstr:lag7 + price:lag0 + lag0:lag1 + lag0:lag2 + lag0:lag3 + lag1:lag2 + lag1:lag3 + lag2:lag3", 
                         data = readings.with.cdh.trimmed, 
                         weights = wghts, 
                         family = Gamma(link="log"))

# Start by trying to reduce temperature interacting with itself.
r24mille.alt12.glm <- update(r24mille.alt11.glm, ~ . - lag2:lag3)
anova(r24mille.alt11.glm, r24mille.alt12.glm, test = "LRT") # Reject reduction

# Hour-of-day interacting with temperature
r24mille.alt13.glm <- update(r24mille.alt11.glm, ~ . - hrstr:lag7)
anova(r24mille.alt11.glm, r24mille.alt13.glm, test = "LRT") # Reject reduction

# Month interacting with temperature
r24mille.alt14.glm <- update(r24mille.alt11.glm, ~ . - month:lag7)
anova(r24mille.alt11.glm, r24mille.alt14.glm, test = "LRT") # Accept reduction

r24mille.alt15.glm <- update(r24mille.alt14.glm, ~ . - month:lag6)
anova(r24mille.alt14.glm, r24mille.alt15.glm, test = "LRT") # Reject reduction

# Attempt to reduce billing interacted with temperature
r24mille.alt16.glm <- update(r24mille.alt14.glm, ~ . - billing_active:lag7)
anova(r24mille.alt14.glm, r24mille.alt16.glm, test = "LRT") # Reject reduction

BIC(r24mille.alt14.glm) # -10275.78


# Try fewer hours into the past according to the same criteria, that if a 
# past temperature is included all hours in between must be included as well 
# (unless I am confident I can explain/justify a gap).
#
# Model proposed by iterative reduction at cdhlag=6:
#     kwh ~ billing_active + month + hrstr + weekend + price + lag0 + lag1 + 
#     lag2 + lag3 + lag6 + hrstr_price + billing_active:lag6 + month:hrstr + 
#     month:price + month:lag2 + month:lag6 + hrstr:lag0 + hrstr:lag3 + 
#     hrstr:lag6 + weekend:lag1 + lag0:lag3
r24mille.alt17.glm <- glm(formula = "kwh ~ billing_active + month + hrstr + weekend + price + lag0 + lag1 + lag2 + lag3 + lag4 + lag5 + lag6 + hrstr:price + billing_active:lag0 + billing_active:lag1 + billing_active:lag2 + billing_active:lag3 + billing_active:lag4 + billing_active:lag5 + billing_active:lag6 + month:hrstr + month:price + month:lag0 + month:lag1 + month:lag2 + month:lag3 + month:lag4 + month:lag5 + month:lag6 + hrstr:lag0 + hrstr:lag1 + hrstr:lag2 + hrstr:lag3 + hrstr:lag4 + hrstr:lag5 + hrstr:lag6 + weekend:lag0 + weekend:lag1 + lag0:lag3 + lag0:lag1 + lag0:lag2 + lag0:lag3 + lag1:lag2 + lag1:lag3 + lag2:lag3", 
                          data = readings.with.cdh.trimmed, 
                          weights = wghts, 
                          family = Gamma(link="log"))

# Start by trying to reduce temperature interacting with itself.
r24mille.alt18.glm <- update(r24mille.alt17.glm, ~ . - lag2:lag3)
anova(r24mille.alt17.glm, r24mille.alt18.glm, test = "LRT") # Reject reduction

# Hour-of-day interacting with temperature
r24mille.alt19.glm <- update(r24mille.alt17.glm, ~ . - hrstr:lag6)
anova(r24mille.alt17.glm, r24mille.alt19.glm, test = "LRT") # Reject reduction

# Month interacting with temperature
r24mille.alt20.glm <- update(r24mille.alt17.glm, ~ . - month:lag6)
anova(r24mille.alt17.glm, r24mille.alt20.glm, test = "LRT") # Reject reduction

# Attempt to reduce billing interacted with temperature
r24mille.alt21.glm <- update(r24mille.alt17.glm, ~ . - billing_active:lag6)
anova(r24mille.alt17.glm, r24mille.alt21.glm, test = "LRT") # Reject reduction

BIC(r24mille.alt17.glm) # -10088.3


# Try fewer hours into the past according to the same criteria, that if a 
# past temperature is included all hours in between must be included as well 
# (unless I am confident I can explain/justify a gap).
#
# Model proposed by iterative reduction at cdhlag=5:
#     kwh ~ billing_active + month + hrstr + weekend + price + lag0 + lag1 + 
#     lag3 + lag5 + hrstr_price + billing_active:price + billing_active:lag5 + 
#     month:hrstr + month:price + month:lag0 + month:lag5 + hrstr:lag0 + 
#     hrstr:lag3 + hrstr:lag5 + price:lag1 + lag0:lag3
r24mille.alt22.glm <- glm(formula = "kwh ~ billing_active + month + hrstr + weekend + price + lag0 + lag1 + + lag2 + lag3 + lag4 + lag5 + hrstr:price + billing_active:price + billing_active:lag0 + billing_active:lag1 + billing_active:lag2 + billing_active:lag3 + billing_active:lag4 + billing_active:lag5 + month:hrstr + month:price + month:lag0 + month:lag1 + month:lag2 + month:lag3 + month:lag4 + month:lag5 + hrstr:lag0 + hrstr:lag1 + hrstr:lag2 + hrstr:lag3 + hrstr:lag4 + hrstr:lag5 + price:lag0 + price:lag1 + lag0:lag3 + lag0:lag1 + lag0:lag2 + lag0:lag3 + lag1:lag2 + lag1:lag3 + lag2:lag3", 
                          data = readings.with.cdh.trimmed, 
                          weights = wghts, 
                          family = Gamma(link="log"))

# Start by trying to reduce temperature interacting with itself.
r24mille.alt23.glm <- update(r24mille.alt22.glm, ~ . - lag2:lag3)
anova(r24mille.alt22.glm, r24mille.alt23.glm, test = "LRT") # Reject reduction

# Reduce price interacting with temperature
r24mille.alt24.glm <- update(r24mille.alt22.glm, ~ . - price:lag1)
anova(r24mille.alt22.glm, r24mille.alt24.glm, test = "LRT") # Accept reduction

r24mille.alt25.glm <- update(r24mille.alt24.glm, ~ . - price:lag0)
anova(r24mille.alt24.glm, r24mille.alt25.glm, test = "LRT") # Reject reduction

# Hour-of-day interacting with temperature
r24mille.alt26.glm <- update(r24mille.alt24.glm, ~ . - hrstr:lag5)
anova(r24mille.alt24.glm, r24mille.alt26.glm, test = "LRT") # Reject reduction

# Month interacting with temperature
r24mille.alt27.glm <- update(r24mille.alt24.glm, ~ . - month:lag5)
anova(r24mille.alt24.glm, r24mille.alt27.glm, test = "LRT") # Reject reduction

# Attempt to reduce billing interacted with temperature
r24mille.alt28.glm <- update(r24mille.alt24.glm, ~ . - billing_active:lag5)
anova(r24mille.alt24.glm, r24mille.alt28.glm, test = "LRT") # Reject reduction

BIC(r24mille.alt24.glm) # -9925.804


# Try fewer hours into the past according to the same criteria, that if a 
# past temperature is included all hours in between must be included as well 
# (unless I am confident I can explain/justify a gap).
#
# Model proposed by iterative reduction at cdhlag=4:
#     kwh ~ billing_active + month + hrstr + weekend + price + lag0 + 
#     lag1 + lag3 + lag4 + hrstr_price + billing_active:month + 
#     billing_active:price + month:hrstr + month:price + month:lag0 + 
#     month:lag4 + hrstr:lag0 + hrstr:lag3 + hrstr:lag4 + weekend:lag1 + 
#     lag0:lag3
r24mille.alt29.glm <- glm(formula = "kwh ~ billing_active + month + hrstr + weekend + price + lag0 + lag1 + lag2 + lag3 + lag4 + hrstr:price + billing_active:month + billing_active:price + month:hrstr + month:price +  + month:lag0 + month:lag1 + month:lag2 + month:lag3 + month:lag4 + hrstr:lag0 + hrstr:lag1 + hrstr:lag2 + hrstr:lag3 + hrstr:lag4 + weekend:lag0 + weekend:lag1 + lag0:lag3 + lag0:lag1 + lag0:lag2 + lag0:lag3 + lag1:lag2 + lag1:lag3 + lag2:lag3", 
                          data = readings.with.cdh.trimmed, 
                          weights = wghts, 
                          family = Gamma(link="log"))

# Start by trying to reduce temperature interacting with itself.
r24mille.alt30.glm <- update(r24mille.alt29.glm, ~ . - lag2:lag3)
anova(r24mille.alt29.glm, r24mille.alt30.glm, test = "LRT") # Accept reduction

r24mille.alt31.glm <- update(r24mille.alt30.glm, ~ . - lag1:lag3)
anova(r24mille.alt30.glm, r24mille.alt31.glm, test = "LRT") # Accept reduction

r24mille.alt32.glm <- update(r24mille.alt31.glm, ~ . - lag0:lag3)
anova(r24mille.alt31.glm, r24mille.alt32.glm, test = "LRT") # Accept reduction

r24mille.alt33.glm <- update(r24mille.alt32.glm, ~ . - lag1:lag2)
anova(r24mille.alt32.glm, r24mille.alt33.glm, test = "LRT") # Accept reduction

r24mille.alt34.glm <- update(r24mille.alt33.glm, ~ . - lag0:lag2)
anova(r24mille.alt33.glm, r24mille.alt34.glm, test = "LRT") # Reject reduction

# Reduce weekend interacting with temperature
r24mille.alt35.glm <- update(r24mille.alt33.glm, ~ . - weekend:lag1)
anova(r24mille.alt33.glm, r24mille.alt35.glm, test = "LRT") # Accept reduction

r24mille.alt36.glm <- update(r24mille.alt35.glm, ~ . - weekend:lag0)
anova(r24mille.alt35.glm, r24mille.alt36.glm, test = "LRT") # Reject reduction

# Hour-of-day interacting with temperature
r24mille.alt37.glm <- update(r24mille.alt35.glm, ~ . - hrstr:lag4)
anova(r24mille.alt35.glm, r24mille.alt37.glm, test = "LRT") # Reject reduction

# Month interacting with temperature
r24mille.alt38.glm <- update(r24mille.alt35.glm, ~ . - month:lag4)
anova(r24mille.alt35.glm, r24mille.alt38.glm, test = "LRT") # Reject reduction

BIC(r24mille.alt35.glm) # -9730.674


# Try fewer hours into the past according to the same criteria, that if a 
# past temperature is included all hours in between must be included as well 
# (unless I am confident I can explain/justify a gap).
#
# Model proposed by iterative reduction at cdhlag=3:
#     kwh ~ billing_active + month + hrstr + weekend + price + lag0 + lag2 + 
#     lag3 + hrstr_price + billing_active:month + billing_active:price + 
#     billing_active:lag3 + month:hrstr + month:price + month:lag0 + 
#     month:lag3 + hrstr:lag0 + hrstr:lag3 + price:lag0 + lag0:lag2
r24mille.alt39.glm <- glm(formula = "kwh ~ billing_active + month + hrstr + weekend + price + lag0 + lag1 + lag2 + lag3 + hrstr:price + billing_active:month + billing_active:price + billing_active:lag0 + billing_active:lag1 + billing_active:lag2 + billing_active:lag3 + month:hrstr + month:price + month:lag0 + month:lag1 + month:lag2 + month:lag3 + hrstr:lag0 + hrstr:lag1 + hrstr:lag2 + hrstr:lag3 + price:lag0 + lag0:lag1 + lag0:lag2 + lag1:lag2", 
                          data = readings.with.cdh.trimmed, 
                          weights = wghts, 
                          family = Gamma(link="log"))

# Start by trying to reduce temperature interacting with itself.
r24mille.alt40.glm <- update(r24mille.alt39.glm, ~ . - lag1:lag2)
anova(r24mille.alt39.glm, r24mille.alt40.glm, test = "LRT") # Accept reduction

r24mille.alt41.glm <- update(r24mille.alt40.glm, ~ . - lag0:lag2)
anova(r24mille.alt40.glm, r24mille.alt41.glm, test = "LRT") # Accept reduction

r24mille.alt42.glm <- update(r24mille.alt41.glm, ~ . - lag0:lag1)
anova(r24mille.alt41.glm, r24mille.alt42.glm, test = "LRT") # Reject reduction

# Reduce price interacting with temperature
r24mille.alt43.glm <- update(r24mille.alt41.glm, ~ . - price:lag0)
anova(r24mille.alt41.glm, r24mille.alt43.glm, test = "LRT") # Reject reduction

# Hour-of-day interacting with temperature
r24mille.alt44.glm <- update(r24mille.alt41.glm, ~ . - hrstr:lag3)
anova(r24mille.alt41.glm, r24mille.alt44.glm, test = "LRT") # Reject reduction

# Month interacting with temperature
r24mille.alt45.glm <- update(r24mille.alt41.glm, ~ . - month:lag3)
anova(r24mille.alt41.glm, r24mille.alt45.glm, test = "LRT") # Reject reduction

# Billing interacting with temperature
r24mille.alt46.glm <- update(r24mille.alt41.glm, ~ . - billing_active:lag3)
anova(r24mille.alt41.glm, r24mille.alt46.glm, test = "LRT") # Reject reduction

BIC(r24mille.alt41.glm) # -9438.884


# Try fewer hours into the past according to the same criteria, that if a 
# past temperature is included all hours in between must be included as well 
# (unless I am confident I can explain/justify a gap).
#
# Model proposed by iterative reduction at cdhlag=2:
#     kwh ~ billing_active + month + hrstr + weekend + price + lag0 + lag1 + 
#     lag2 + hrstr_price + billing_active:price + month:hrstr + month:lag2 + 
#     hrstr:lag0 + hrstr:lag2 + price:lag0 + lag0:lag1
r24mille.alt47.glm <- glm(formula = "kwh ~ billing_active + month + hrstr + weekend + price + lag0 + lag1 + lag2 + hrstr:price + billing_active:price + month:hrstr + month:lag0 + month:lag1 + month:lag2 + hrstr:lag0 + hrstr:lag1 + hrstr:lag2 + price:lag0 + lag0:lag1", 
                          data = readings.with.cdh.trimmed, 
                          weights = wghts, 
                          family = Gamma(link="log"))

# Start by trying to reduce temperature interacting with itself.
r24mille.alt48.glm <- update(r24mille.alt47.glm, ~ . - lag0:lag1)
anova(r24mille.alt47.glm, r24mille.alt48.glm, test = "LRT") # Reject reduction

# Reduce price interacting with temperature
r24mille.alt49.glm <- update(r24mille.alt47.glm, ~ . - price:lag0)
anova(r24mille.alt47.glm, r24mille.alt49.glm, test = "LRT") # Reject reduction

# Hour-of-day interacting with temperature
r24mille.alt50.glm <- update(r24mille.alt47.glm, ~ . - hrstr:lag2)
anova(r24mille.alt47.glm, r24mille.alt50.glm, test = "LRT") # Reject reduction

# Month interacting with temperature
r24mille.alt51.glm <- update(r24mille.alt47.glm, ~ . - month:lag2)
anova(r24mille.alt47.glm, r24mille.alt51.glm, test = "LRT") # Reject reduction

BIC(r24mille.alt47.glm) # -9138.432


# Try fewer hours into the past according to the same criteria, that if a 
# past temperature is included all hours in between must be included as well 
# (unless I am confident I can explain/justify a gap).
#
# Model proposed by iterative reduction at cdhlag=1:
#     kwh ~ billing_active + month + hrstr + weekend + price + lag0 + lag1 + 
#     hrstr_price + billing_active:price + month:hrstr + month:lag1 + 
#     hrstr:lag0 + price:lag0 + lag0:lag1
r24mille.alt52.glm <- glm(formula = "kwh ~ billing_active + month + hrstr + weekend + price + lag0 + lag1 + hrstr:price + billing_active:price + month:hrstr + month:lag0 + month:lag1 + hrstr:lag0 + price:lag0 + lag0:lag1", 
                          data = readings.with.cdh.trimmed, 
                          weights = wghts, 
                          family = Gamma(link="log"))

alt52.add1.rnd1 <- add1(object = r24mille.alt52.glm,
                        scope = formula("(.)^2"),
                        test = "LRT", 
                        k = log(num.observations))