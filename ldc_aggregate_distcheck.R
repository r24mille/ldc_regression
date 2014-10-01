library(MASS) # Generates the shape and rate parameters
library(car) # Companion to Applied Regression (better qqPlot)

# Load SmartMeterReading data from CSV
home <- Sys.getenv("HOME")
fpath <- file.path(home, 
                   "../Dropbox/ISS4E/R/", 
                   "full_aggregate_readings.csv")
readings.aggregate <- read.csv(fpath)

# Fit a reference Gamma distribution to the data
gammafit <- fitdistr(readings.aggregate$kwh, "gamma")
gfitshape <- gammafit$estimate[1]
gfitrate <- gammafit$estimate[2]

# Generate a sample of random deviates
grand <- rgamma(n = gammafit$n,
                shape = gfitshape,
                rate = gfitrate)


qqplot(y = readings.aggregate$kwh, 
       x = grand,
       main = paste0("Q-Q Plot of Gamma Distribution vs. Average Reading (kWh)"),
       xlab = paste0("Expected Fit of Gamma Distribution (shape=",
                     round(gfitshape, 1),
                     ", rate=",
                     round(gfitrate, 1),
                     ")"),
       ylab = "Average Reading (kWh) for Aggregate")
qqline(y = grand,
       distribution = function(p) qgamma(p, shape = gfitshape, rate = gfitrate))

# Provide a Q-Q Plot to do a visual check. CAR package has a lot of 
# convenient defaults (eg. confidence intervals)
qqPlot(x = readings.aggregate$kwh,
       distribution = "gamma",
       shape = gfitshape,
       rate = gfitrate,
       main = paste0("Q-Q Plot of Gamma Distribution vs. Average Reading (kWh)"),
       xlab = paste0("Gamma Distribution (shape=",
                     round(gfitshape, 1),
                     ", rate=",
                     round(gfitrate, 1),
                     ")"),
       ylab = "Average Reading (kWh) for Aggregate")

# Plot the probability density function of the real and model data
plot(density(readings.aggregate$kwh, bw = .25), 
     xlim = c(0,5), 
     ylim = c(0,1),
     col="red")
par(new=TRUE)
plot(density(grand, bw = .25), 
     xlim = c(0,5), 
     ylim = c(0,1),
     col="blue")