library(segmented) # Find linear regression breakpoint
library(ggplot) # Better plots!
library(glinternet) # Hierarchical group-LASSO variable selection

# Source functions in other files
source("dataframe_processing.R")

# Load SmartMeterReading data from CSV
home <- Sys.getenv("HOME")
fpath <- file.path(home, 
                   "../Dropbox/ISS4E/R/", 
                   "aggregate_readings_01Mar2011_through_17Oct2012.csv")
readings.aggregate <- read.csv(fpath)
readings.aggregate <- InitAggregateReadings(readings.aggregate)

##
# Use 'segmented' package to find the optimal temperature breakpoint.
# 
# TODO(r24mille: I attempted to use psi c(0,16) to find two temperature 
#                breakpoints (ie. Two-Threshold Regression from Moral-Carcedo 
#                et al.) but it does not converge. Looking at the plot, this 
#                seems sensible that only a Switching Regression applies.
model.readings.lm.presegment <- lm(log(kwh) ~ temperature, 
                                     data = readings.aggregate)
seg <- segmented(obj = model.readings.lm.presegment, 
                 seg.Z = ~ temperature,
                 psi = 16)

# plot(seg, 
#      res=TRUE, 
#      res.col = rgb(0, 0, 0, 10, maxColorValue=100),
#      pch = 16,
#      main = "Piecewise Linear Fit with Temperature Breakpoint",
#      xlab = "Summer Temperatures (Celsius)",
#      ylab = "Effect of Temperature",
#      col = "red",
#      lwd = 2)

temperature.break <- seg$psi[1,2]
readings.aggregate$temp_over_break <- ifelse(readings.aggregate$temperature > temperature.break, 
                                             readings.aggregate$temperature - temperature.break, 
                                             0)
readings.aggregate$temp_under_break <- ifelse(readings.aggregate$temperature < temperature.break, 
                                              readings.aggregate$temperature - temperature.break, 
                                              0)

# Clean up some unneeded variables from R environment
rm(model.readings.lm.presegment, seg)

# Create columns that contain the previous hours' temperature > breakpoint 
# and previous hours temperature < breakpoint.
temp.hrs <- 6
temps <- CreatePastTemperatureMatrix(nlags = temp.hrs, 
                                     df.readings = readings.aggregate)
readings.aggregate <- cbind(readings.aggregate, temps)
rm(temps)
readings.trimmed <- TrimExplanatoryVariables(readings.aggregate)

# Prep data into matrices for use in hierarchical group-lasso
Y.mat <- as.matrix(log(readings.trimmed[,1]))
X.mat <- NumericFactorCodedMatrix(readings.trimmed[,-1])
numlvl.vec <- NumberFactorLevels(readings.trimmed[,-1])
cv.hglasso <- glinternet.cv(X = X.mat, 
                            Y = Y.mat, 
                            numLevels = numlvl.vec, 
                            family = "gaussian")

# Select the column names of main effects and interactions
cat.colnames <- names(readings.trimmed[c(FALSE, numlvl.vec > 1)])
cont.colnames <- names(readings.trimmed[c(FALSE, numlvl.vec == 1)])

# Get coefficients for all lambdas
coef.hglasso <- coef(cv.hglasso$glinternetFit)

# Column names for main effects at chosen lambda
lmda.idx <- which(cv.hglasso$lambda == cv.hglasso$lambdaHat1Std)
me.colnames <- c(cat.colnames[coef.hglasso[[lmda.idx]]$mainEffects$cat],
                 cont.colnames[coef.hglasso[[lmda.idx]]$mainEffects$cont])
# Column names for interactions at chosen lambda
lmda.catcat <- coef.hglasso[[lmda.idx]]$interactions$catcat
lmda.contcont <- coef.hglasso[[lmda.idx]]$interactions$contcont
lmda.catcont <- coef.hglasso[[lmda.idx]]$interactions$catcont
int.colnames <- c(paste(cat.colnames[lmda.catcat[,1]], cat.colnames[lmda.catcat[,2]], sep=":"),
                  paste(cont.colnames[lmda.contcont[,1]], cont.colnames[lmda.contcont[,2]], sep=":"),
                  paste(cat.colnames[lmda.catcont[,1]], cont.colnames[lmda.catcont[,2]], sep=":"))
c(me.colnames, int.colnames)



# Backtransform optimal, fitted values and compute more meaningful RMSE and MAPE
Y.bktr <- exp(Y.mat)
hglasso.fitted.bktr <- exp(cv.hglasso$fitted)
hglasso.rmse.bktr <- sqrt(mean((hglasso.fitted.bktr - Y.bktr)^2, 
                                     na.rm = TRUE))
hglasso.mape.bktr <- (1/nrow(Y.bktr)) * sum(abs((Y.bktr - hglasso.fitted.bktr)/Y.bktr))
hglasso.residuals.bktr <- Y.bktr - hglasso.fitted.bktr
hglasso.residuals <- Y.mat - cv.hglasso$fitted

# QQ Plot of fitted vs. observed values
library(car) # Companion to Applied Regression (better qqPlot)
qqPlot(x = hglasso.residuals,
       distribution = "norm",
       envelope = .95,
       mean = 0,
       main = "Hierarchical Group-LASSO Residuals vs. Normal Distribution",
       ylab = "Residuals of log(kWh)",
       xlab = "Norm Quantiles")

# Plot the histogram of residuals
plot(density(hglasso.residuals),
     main = "Histogram of Residuals",
     xlab = "Residuals of log(kWh)")
abline(v = mean(hglasso.residuals),
       col = "darkgray",
       lty = 2)
axis(3, # top label
     padj = 1, # align to bottom of label bounding box
     at = round(mean(hglasso.residuals), 3))

# Plot residuals to look for patterns
hglasso.df.resids <- data.frame(hglasso.residuals, 
                                readings.aggregate$timestamp_dst)
names(hglasso.df.resids) <- c("res", "dst")
ggresiduals <- (ggplot(hglasso.df.resids, aes(y = res, x = dst))
                + geom_point(alpha = 1/3)
                + labs(x = "Date (hourly)", 
                       y = "Residuals of log(kWh)",
                       title= bquote(Residuals~of~Hierarchical~Group-LASSO~Model~at~hat(lambda)==.(cv.hglasso$lambdaHat, 2)))
                + scale_x_datetime(labels = date_format("%b %Y"), 
                                   breaks = "1 month")
                + theme(axis.text.x = element_text(angle=90, 
                                                   vjust=0.5, 
                                                   size=12),
                        axis.text.y = element_text(size=12))
)
ggresiduals