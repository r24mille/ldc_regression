library(MASS) # Generates the shape and rate parameters
library(car) # Companion to Applied Regression (better qqPlot)

# Load SmartMeterReading data from CSV
home <- Sys.getenv("HOME")
fpath <- file.path(home, 
                   "../Dropbox/ISS4E/R/", 
                   "full_aggregate_readings.csv")
readings.aggregate <- read.csv(fpath)

# Provide a Q-Q Plot to do a visual check against normal distribution. CAR package has a lot of 
# convenient defaults (eg. confidence intervals)
qqPlot(x = readings.aggregate$kwh,
       distribution = "norm",
       mean = mean(readings.aggregate$kwh),
       main = paste0("Q-Q Plot of Normal Distribution vs. Observed Readings (kWh)"),
       xlab = paste0("Normal Distribution (Mean=",
                     round(mean(readings.aggregate$kwh), 1),
                     ", SD=",
                     round(sd(readings.aggregate$kwh), 1),
                     ")"),
       ylab = "Average Reading (kWh) for Aggregate")

# Fit a reference Gamma distribution to the data
gammafit <- fitdistr(readings.aggregate$kwh, "gamma")
gfitshape <- gammafit$estimate[1]
gfitrate <- gammafit$estimate[2]

# Provide a Q-Q Plot to do a visual check. CAR package has a lot of 
# convenient defaults (eg. confidence intervals)
qqPlot(x = readings.aggregate$kwh,
       distribution = "gamma",
       shape = gfitshape,
       rate = gfitrate,
       main = paste0("Q-Q Plot of Gamma Distribution vs. Observed Readings (kWh)"),
       xlab = paste0("Gamma Distribution (shape=",
                     round(gfitshape, 1),
                     ", rate=",
                     round(gfitrate, 1),
                     ")"),
       ylab = "Average Reading (kWh) for Aggregate")

# Generate a sample of random normal distribution
rndmean <- round(mean(readings.aggregate$kwh), 1)
rndsd <- round(sd(readings.aggregate$kwh), 1)
nrand <- rnorm(n = length(readings.aggregate$kwh),
               mean = mean(readings.aggregate$kwh),
               sd = sd(readings.aggregate$kwh))

# Generate a sample of random Gamma deviates
grand <- rgamma(n = length(readings.aggregate$kwh),
                shape = gfitshape,
                rate = gfitrate)

# Plot the probability density function of the real and Normal distribution
kdbw = 0.25
plot(density(readings.aggregate$kwh, bw = kdbw), 
     xlim = c(0,5), 
     xlab = paste0("Average Meter Reading (kWh, n=", 
                   length(readings.aggregate$kwh),
                   ")"),
     ylim = c(0,1),
     ylab = paste0("Kernel Density (bandwidth=",
                   kdbw,
                   ")"),
     main = "Density of Observations, Gaussian, and Gamma Distributions",
     col="red")
lines(density(nrand, bw = kdbw), 
      col="blue")
lines(density(grand, bw = kdbw), 
      col="green")
legend("topright",
       col = c("red", "blue", "green"), 
       lty = 1, 
       legend=c( "Observed Meter Readings", 
                 as.expression(bquote(paste("Gaussian Distribution (",
                                            mu == .(rndmean),
                                            ", ",
                                            sigma^2 == .(rndsd),
                                            ")     "))),
                 as.expression(bquote(paste("Gamma Distribution (",
                                            alpha == .(round(gfitshape, 1)),
                                            ", ",
                                            beta == .(round(gfitrate, 1)),
                                            ")     ")))))