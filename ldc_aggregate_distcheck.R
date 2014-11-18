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
       main = paste0("Q-Q Plot of Normal Dist. vs. Readings (kWh)"),
       xlab = paste0("Normal Distribution (Mean=",
                     round(mean(readings.aggregate$kwh), 1),
                     ", SD=",
                     round(sd(readings.aggregate$kwh), 1),
                     ")"),
       ylab = "Average Reading (kWh) for Aggregate")

# Provide a Q-Q Plot to do a visual check of log(response) against normal distribution. CAR package has a lot of 
# convenient defaults (eg. confidence intervals)
qqPlot(x = log(readings.aggregate$kwh),
       distribution = "norm",
       mean = mean(readings.aggregate$kwh),
       main = paste0("Q-Q Plot of Normal Dist. vs. Log-Transformed Readings (kWh)"),
       xlab = paste0("Normal Distribution (Mean=",
                     round(mean(readings.aggregate$kwh), 1),
                     ", SD=",
                     round(sd(readings.aggregate$kwh), 1),
                     ")"),
       ylab = "Average log(kWh) for Aggregate")

# Fit a reference Gamma distribution to the data
gamma.fit <- fitdistr(readings.aggregate$kwh, "gamma")
gamma.shape <- gamma.fit$estimate[1]
gamma.rate <- gamma.fit$estimate[2]

# Fit a reference lognormal distribution to the data
lognorm.fit <- fitdistr(readings.aggregate$kwh, "lognormal")
lnorm.meanlog <- lognorm.fit$estimate[1]
lnorm.sdlog <- lognorm.fit$estimate[2]

# Provide a Q-Q Plot to do a visual check. CAR package has a lot of 
# convenient defaults (eg. confidence intervals)
qqPlot(x = readings.aggregate$kwh,
       distribution = "gamma",
       shape = gamma.shape,
       rate = gamma.rate,
       main = paste0("Q-Q Plot of Gamma Dist. vs. Readings (kWh)"),
       xlab = paste0("Gamma Distribution (shape=",
                     round(gamma.shape, 1),
                     ", rate=",
                     round(gamma.rate, 1),
                     ")"),
       ylab = "Average Reading (kWh) for Aggregate")

# Generate a sample of random normal distribution
mean.rounded <- round(mean(readings.aggregate$kwh), 1)
sd.rounded <- round(sd(readings.aggregate$kwh), 1)
rnorm.sample <- rnorm(n = length(readings.aggregate$kwh),
               mean = mean(readings.aggregate$kwh),
               sd = sd(readings.aggregate$kwh))

# Generate a sample of random Gamma deviates
rgamma.sample <- rgamma(n = length(readings.aggregate$kwh),
                shape = gamma.shape,
                rate = gamma.rate)

# Generate a sample of random log-normal deviates
rlnorm.sample <- rlnorm(n = length(readings.aggregate$kwh),
                       meanlog = lnorm.meanlog,
                       sdlog = lnorm.sdlog)

# Plot the probability density function of the real and Normal distribution
par(mfrow=c(1,1))
kdbw = 0.25
width.line = 2
color.kwh = rgb(228,26,28, 255, maxColorValue=255)
color.norm = rgb(55,126,184, 255, maxColorValue=255)
color.gamma = rgb(77,175,74, 255, maxColorValue=255)
color.lnorm = rgb(152,78,163, 255, maxColorValue=255)
plot(density(readings.aggregate$kwh, bw = kdbw), 
     xlim = c(0,5), 
     xlab = paste0("Average Meter Reading (kWh, n=", 
                   length(readings.aggregate$kwh),
                   ")"),
     ylim = c(0,1),
     ylab = paste0("Kernel Density (bandwidth=",
                   kdbw,
                   ")"),
     main = "Density of Observations and Probability Distributions",
     lwd = width.line,
     col = color.kwh)
lines(density(rnorm.sample, bw = kdbw), 
      lwd = width.line,
      col = color.norm)
lines(density(rgamma.sample, bw = kdbw), 
      lwd = width.line,
      col = color.gamma)
lines(density(rlnorm.sample, bw = kdbw),
      lwd = width.line,
      col = color.lnorm)
legend("topright",
       col = c(color.kwh, color.lnorm, color.gamma, color.norm), 
       lty = 1, 
       lwd = width.line, 
       legend=c( "Observed Meter Readings", 
                 as.expression(bquote(paste("Log-Normal Distribution (",
                                            mu == .(mean.rounded),
                                            ", ",
                                            sigma^2 == .(sd.rounded),
                                            ")     "))),
                 as.expression(bquote(paste("Gamma Distribution (",
                                            alpha == .(round(gamma.shape, 1)),
                                            ", ",
                                            beta == .(round(gamma.rate, 1)),
                                            ")     "))),
                 as.expression(bquote(paste("Normal Distribution (",
                                            mu == .(mean.rounded),
                                            ", ",
                                            sigma^2 == .(sd.rounded),
                                            ")     ")))))