require(ggplot2)        # Better plotting!
require(scales)         # For scaling axes (eg. datetime)
require(gridExtra)      # Place multiple ggplot2 plots into a single image
require(xtable)         # LaTeX 
require(visreg)
require(segmented)
require(dlnm)           # Distributed Lag Nonlinear Models
require(car)            # Companion to Applied Regression
require(splines)
require(zoo)

# Some common functions for manipulating the meter reading data.frame
source("dataframe_processing.R")
# Standardize things like plot titles, labels, sizes, and colours
source("thesis_standardization.R")
# Temperature transformation functions
source("temperature_transformation.R")



#####
# Load aggregate smart meter readings and weighted average temperature data 
# from a CSV. Similarly, load the Environment Canada Weather descriptions for 
# the area from CSV.
#####
home <- Sys.getenv("HOME")
fpath <- file.path(home, 
                   "../Dropbox/ISS4E/R", 
                   "aggregate_readings_01Mar2011_through_17Oct2012.csv")

readings.aggregate <- InitReadingsDataFrame(fpath = fpath, 
                                            is.aggregate = TRUE)

fpath2 <- file.path(home, 
                    "../Dropbox/ISS4E/R", 
                    "weather_desc_01Mar2011_through_17Oct2012.csv")
weather <- read.csv(fpath2, 
                    na.strings = c("NULL", "NA", "NaN"), 
                    stringsAsFactors = FALSE)

# Add weather information to readings.aggregate
weather_desc.original <- unique(weather$weather_desc)
readings.aggregate$weather_reduced <- factor(ReduceWeatherCoarseTerms(weather$weather_desc))
readings.aggregate$severe_weather <- ReduceSevereWeather(weather$weather_desc)
readings.aggregate$working_day <- ((readings.aggregate$holiday == TRUE | 
                                      readings.aggregate$dayname == "Sat" | 
                                      readings.aggregate$dayname == "Sun") == FALSE)

# Add column to identify utility rate seasons
readings.aggregate$rateseason <- rep(NA, nrow(readings.aggregate))
readings.aggregate$rateseason[readings.aggregate$month %in% c("m5", "m6", "m7", 
                                                              "m8", "m9", 
                                                              "m10")] <- "summer"
readings.aggregate$rateseason[readings.aggregate$month %in% c("m11", "m12", 
                                                              "m1", "m2", "m3", 
                                                              "m4")] <- "winter"
readings.aggregate$rateseason <- factor(readings.aggregate$rateseason, 
                                        c("winter", "summer"))

# Add a period of day to each observation for quantifying effects for each.
readings.aggregate$dayperiod <- rep(NA, nrow(readings.aggregate))
readings.aggregate$dayperiod[readings.aggregate$hrstr %in% c("h0", "h1", "h2", 
                                                             "h3", "h4", "h5", 
                                                             "h6", "h19", "h20", 
                                                             "h21", "h22", 
                                                             "h23")] <- "overnight"
readings.aggregate$dayperiod[readings.aggregate$hrstr %in% c("h7", "h8", "h9", 
                                                             "h10", "h17", 
                                                             "h18")] <- "transition"
readings.aggregate$dayperiod[readings.aggregate$hrstr %in% c("h11", "h12", 
                                                             "h13", "h14", 
                                                             "h15", "h16")] <- "afternoon"
readings.aggregate$dayperiod[readings.aggregate$working_day == FALSE] <- "nonwd"

# For clarity, reorder explanatory variables which are factors
readings.aggregate <- OrderFactors(readings.aggregate)
rm(fpath, fpath2, home)



#####
# Show the flat rate structure
#####
flatsummer.df <- data.frame("rate" = rep(x = 6.8, times = 25), "hour" = c(0:24))
flatwinter.df <- data.frame("rate" = rep(x = 7.1, times = 25), "hour" = c(0:24))

flatsummer.plot <- (ggplot() + 
                      geom_area(data = flatsummer.df, aes(y = rate, x = hour)) + 
                      scale_y_continuous(limits = c(0,12), 
                                         expand = c(0,0)) + 
                      scale_x_discrete(labels = HourLabels(), 
                                       expand = c(0, -1)) + 
                      ThesisPlotTheme() + 
                      theme(axis.text.x = element_text(angle = 90, 
                                                       vjust = 0.5)) + 
                      labs(x = "Hour of Day", 
                           y = "Price of Electricity (cents per kWh, CAD)",
                           title= "Flat Electricity Rate, Summer"))
flatsummer.plot

flatwinter.plot <- (ggplot() + 
                      geom_area(data = flatwinter.df, aes(y = rate, x = hour)) + 
                      scale_y_continuous(limits = c(0,12), 
                                         expand = c(0,0)) + 
                      scale_x_discrete(labels = HourLabels(), 
                                       expand = c(0, -1)) + 
                      ThesisPlotTheme() + 
                      theme(axis.text.x = element_text(angle = 90, 
                                                       vjust = 0.5)) + 
                      labs(x = "Hour of Day", 
                           y = "Price of Electricity (cents per kWh, CAD)",
                           title= "Flat Electricity Rate, Winter"))
flatwinter.plot

grid.arrange(flatsummer.plot, flatwinter.plot, ncol = 2, nrow = 1)
dev.print(file="../../Figures/FlatRates.png", device=png, height = 600, width = 1600)



#####
# Show the TOU rate structure
#####
tousummer.overnight.df <- data.frame("Rate" = c(6.5, 6.5, 6.5, 6.5, 6.5, 6.5, 6.5, 6.5), 
                                     "Period" = c("Off-peak", "Off-peak", "Off-peak", "Off-peak", "Off-peak", "Off-peak", "Off-peak", "Off-peak"), 
                                     "Hour" = c(0:7))
tousummer.morning.df <- data.frame("Rate" = c(10, 10, 10, 10, 10), 
                                   "Period" = c("Mid-peak", "Mid-peak", "Mid-peak", "Mid-peak", "Mid-peak"), 
                                   "Hour" = c(7:11))
tousummer.afternoon.df <- data.frame("Rate" = c(11.7, 11.7, 11.7, 11.7, 11.7, 11.7, 11.7), 
                                     "Period" = c("On-peak", "On-peak", "On-peak", "On-peak", "On-peak", "On-peak", "On-peak"), 
                                     "Hour" = c(11:17))
tousummer.evening.df <- data.frame("Rate" = c(10, 10, 10), 
                                   "Period" = c("Mid-peak", "Mid-peak", "Mid-peak"), 
                                   "Hour" = c(17:19))
tousummer.night.df <- data.frame("Rate" = c(6.5, 6.5, 6.5, 6.5, 6.5, 6.5), 
                                 "Period" = c("Off-peak", "Off-peak", "Off-peak", "Off-peak", "Off-peak", "Off-peak"), 
                                 "Hour" = c(19:24))

touwinter.overnight.df <- data.frame("Rate" = c(6.2, 6.2, 6.2, 6.2, 6.2, 6.2, 6.2, 6.2), 
                                     "Period" = c("Off-peak", "Off-peak", "Off-peak", "Off-peak", "Off-peak", "Off-peak", "Off-peak", "Off-peak"), 
                                     "Hour" = c(0:7))
touwinter.morning.df <- data.frame("Rate" = c(10.8, 10.8, 10.8, 10.8, 10.8), 
                                   "Period" = c("On-peak", "On-peak", "On-peak", "On-peak", "On-peak"), 
                                   "Hour" = c(7:11))
touwinter.afternoon.df <- data.frame("Rate" = c(9.2, 9.2, 9.2, 9.2, 9.2, 9.2, 9.2), 
                                     "Period" = c("Mid-peak", "Mid-peak", "Mid-peak", "Mid-peak", "Mid-peak", "Mid-peak", "Mid-peak"), 
                                     "Hour" = c(11:17))
touwinter.evening.df <- data.frame("Rate" = c(10.8, 10.8, 10.8), 
                                   "Period" = c("On-peak", "On-peak", "On-peak"), 
                                   "Hour" = c(17:19))
touwinter.night.df <- data.frame("Rate" = c(6.2, 6.2, 6.2, 6.2, 6.2, 6.2), 
                                 "Period" = c("Off-peak", "Off-peak", "Off-peak", "Off-peak", "Off-peak", "Off-peak"), 
                                 "Hour" = c(19:24))

tousummer.plot <- (ggplot() + 
                     geom_area(data = tousummer.overnight.df, aes(y = Rate, 
                                                                  x = Hour, 
                                                                  fill = Period)) + 
                     geom_area(data = tousummer.morning.df, aes(y = Rate, 
                                                                x = Hour, 
                                                                fill = Period)) + 
                     geom_area(data = tousummer.afternoon.df, aes(y = Rate, 
                                                                  x = Hour, 
                                                                  fill = Period)) + 
                     geom_area(data = tousummer.evening.df, aes(y = Rate, 
                                                                x = Hour, 
                                                                fill = Period)) + 
                     geom_area(data = tousummer.night.df, aes(y = Rate, 
                                                              x = Hour, 
                                                              fill = Period)) + 
                     scale_y_continuous(limits = c(0,12), expand = c(0,0)) + 
                     scale_x_discrete(labels = HourLabels(), expand = c(0, -1)) + 
                     scale_fill_manual(values = TouPricePalette(), 
                                       breaks = c("Off-peak", "Mid-peak", "On-peak")) + 
                     ThesisPlotTheme() + 
                     theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) + 
                     labs(x = "Hour of Day", 
                          y = "Electricity Rate (cents per kWh, CAD)",
                          title= "Time-of-Use Rates, Summer"))
tousummer.plot

touwinter.plot <- (ggplot() + 
                     geom_area(data = touwinter.overnight.df, aes(y = Rate, 
                                                                  x = Hour, 
                                                                  fill = Period)) + 
                     geom_area(data = touwinter.morning.df, aes(y = Rate, 
                                                                x = Hour, 
                                                                fill = Period)) + 
                     geom_area(data = touwinter.afternoon.df, aes(y = Rate, 
                                                                  x = Hour, 
                                                                  fill = Period)) + 
                     geom_area(data = touwinter.evening.df, aes(y = Rate, 
                                                                x = Hour, 
                                                                fill = Period)) + 
                     geom_area(data = touwinter.night.df, aes(y = Rate, 
                                                              x = Hour, 
                                                              fill = Period)) + 
                     scale_y_continuous(limits = c(0,12), expand = c(0,0)) + 
                     scale_x_discrete(labels = HourLabels(), expand = c(0, -1)) + 
                     scale_fill_manual(values = TouPricePalette(), 
                                       breaks = c("Off-peak", "Mid-peak", "On-peak")) + 
                     ThesisPlotTheme() + 
                     theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) + 
                     labs(x = "Hour of Day", 
                          y = "Electricity Rate (cents per kWh, CAD)",
                          title= "Time-of-Use Rates, Winter"))
touwinter.plot

grid.arrange(tousummer.plot, touwinter.plot, ncol = 2, nrow = 1)
dev.print(file="../../Figures/TimeOfUseRates.png", device=png, height = 600, width = 1600)



#####
# Visualize the number of meters aggregated hourly
#####
max_agg <- max(readings.aggregate$agg_count)
agg_count.tsplot <- (ggplot() +
                       geom_area(data = readings.aggregate, 
                                 aes(timestamp_std, agg_count)) + 
                       geom_hline(yintercept = max_agg,
                                  linetype = "dashed",
                                  colour = "blue") + 
                       geom_text(aes(x = as.POSIXct("2011-03-12 00:00:00"),
                                 y = max_agg + 550,
                                 label = paste0("max=", max_agg)),
                                 colour = "blue",
                                 size = 4.1) + 
                       labs(x = "Date (hourly)", 
                            y = "Number of Observations",
                            title= "Number of Households Reporting Each Hour (March 1, 2011 - October 17, 2012)") + 
                       ylim(0, (max_agg + 750)) + 
                       scale_x_datetime(labels = date_format(format = "%b %Y"), 
                                        breaks = "1 month") + 
                       ThesisPlotTheme() + 
                       theme(axis.text.x = element_text(angle=90, 
                                                        vjust=0.5)))
agg_count.tsplot
dev.print(file="../../Figures/AggregateCountTimeseries.png", device=png, height = 600, width = 800)



#####
# Average houshold electricity demand as a function of time
#####
kwh.tsplot <- (ggplot() +
                 geom_point(data = readings.aggregate, 
                            aes(timestamp_dst, kwh),
                            alpha = 0.5) +  
                 labs(x = "Time (hourly)", 
                      y = "Aggregate Electricity Demand (kWh)",
                      title= "Aggregate Electricity Demand (Mar 1, 2011 - Oct 17, 2012)") + 
                 scale_x_datetime(labels = date_format(format = "%b %Y"), 
                                  breaks = "1 month") + 
                 ThesisPlotTheme() + 
                 theme(axis.text.x = element_text(angle=90, 
                                                  vjust=0.5)))
kwh.tsplot
dev.print(file="../../Figures/DemandTimeSeries.png", device=png, height = 600, width = 800)



#####
# Average houshold electricity demand density plot
#####
response.bw = 0.05
plot(density(readings.aggregate$kwh, bw = response.bw),
     main = paste0("Density of Response Variables ", 
                   "(N=", nrow(readings.aggregate), ", Bandwidth=", response.bw, ")"), 
     xlab = "Aggregate Electricity Demand (kWh)", 
     lwd = 3)
dev.print(file = "../../Figures/DemandDensity.png", device = png, height = 600, width = 800)



#####
# Plot electricity consumption by hour of day and demonstrate that it needs to 
# be explicitly handed as the periodicity term. Hopefully it improves serially 
# correlated errors.
#####
kwh.hrstr.plot <- (ggplot(readings.aggregate, aes(hrstr, kwh)) + 
                     geom_boxplot() + 
                     theme(legend.position = "none") + 
                     scale_x_discrete(labels = HourLabels()) + 
                     ThesisPlotTheme() + 
                     theme(axis.text.x = element_text(angle = 90, 
                                                      vjust = 0.5)) + 
                     labs(x = "Hour of Day", 
                          y = "Aggregate Electricity Demand (kWh)",
                          title= "Distribution of Aggregate Electricity Demand by Hour of Day"))
kwh.hrstr.plot
dev.print(file="../../Figures/DemandHourOfDay.png", device=png, height = 600, width = 800)

kwh.hrstr.scatterplot <- (ggplot(readings.aggregate, aes(hrstr, kwh)) + 
                     geom_point(alpha = 0.5) + 
                     theme(legend.position = "none") + 
                     scale_x_discrete(labels = HourLabels()) + 
                     ThesisPlotTheme() + 
                     theme(axis.text.x = element_text(angle = 90, 
                                                      vjust = 0.5)) + 
                     labs(x = "Hour of Day", 
                          y = "Aggregate Electricity Demand (kWh)",
                          title= "Distribution of Aggregate Electricity Demand by Hour of Day"))
kwh.hrstr.scatterplot
dev.print(file="../../Figures/DemandHourOfDayScatterPlot.png", device=png, height = 600, width = 800)



#####
# Introduce working_day term, since distribution of weekdays seems similar.
#####
summary(subset(readings.aggregate, subset = holiday == TRUE)) # Holidays occur on Sun/Mon/Sat
dayname.plot <- (ggplot(readings.aggregate, aes(dayname, kwh)) + 
                   geom_boxplot(aes(fill = dayname)) + 
                   theme(legend.position = "none") + 
                   ThesisPlotTheme() + 
                   labs(x = "Day of the Week", 
                        y = "Aggregate Electricity Demand (kWh)",
                        title= "Distribution of Aggregate Electricity Demand by Day of Week"))
holiday.plot <- (ggplot(readings.aggregate, aes(holiday, kwh)) + 
                   geom_boxplot() + 
                   ThesisPlotTheme() + 
                   labs(x = "Holiday", 
                        y = "Aggregate Electricity Demand (kWh)",
                        title= "Distribution of Aggregate Electricity Demand by Holiday Indicator"))
grid.arrange(dayname.plot, holiday.plot, ncol = 2, nrow = 1)
dev.print(file="../../Figures/DemandDaynameAndHoliday.png", device=png, height = 600, width = 1600)

kwh.holiday.median <- median(subset(readings.aggregate, holiday == TRUE)$kwh)
print(kwh.holiday.median)
kwh.nonholiday.median <- median(subset(readings.aggregate, holiday == FALSE)$kwh)
print(kwh.nonholiday.median)
dayname.holiday.datafr <- data.frame(dayname = character(), 
                                     median = numeric(), 
                                     holidays = integer(),
                                     stringsAsFactors = FALSE)
kwh.Sun.median <- median(subset(readings.aggregate, dayname == "Sun")$kwh)
holiday.Sun.count <- nrow(subset(readings.aggregate, dayname == "Sun" & holiday == TRUE)) / 24
dayname.holiday.datafr[1,] <- c("Sunday", kwh.Sun.median, holiday.Sun.count)
kwh.Mon.median <- median(subset(readings.aggregate, dayname == "Mon")$kwh)
holiday.Mon.count <- nrow(subset(readings.aggregate, dayname == "Mon" & holiday == TRUE)) / 24
dayname.holiday.datafr[2,] <- c("Monday", kwh.Mon.median, holiday.Mon.count)
kwh.Tue.median <- median(subset(readings.aggregate, dayname == "Tue")$kwh)
holiday.Tue.count <- nrow(subset(readings.aggregate, dayname == "Tue" & holiday == TRUE)) / 24
dayname.holiday.datafr[3,] <- c("Tuesday", kwh.Tue.median, holiday.Tue.count)
kwh.Wed.median <- median(subset(readings.aggregate, dayname == "Wed")$kwh)
holiday.Wed.count <- nrow(subset(readings.aggregate, dayname == "Wed" & holiday == TRUE)) / 24
dayname.holiday.datafr[4,] <- c("Wednesday", kwh.Wed.median, holiday.Wed.count)
kwh.Thu.median <- median(subset(readings.aggregate, dayname == "Thu")$kwh)
holiday.Thu.count <- nrow(subset(readings.aggregate, dayname == "Thu" & holiday == TRUE)) / 24
dayname.holiday.datafr[5,] <- c("Thursday", kwh.Thu.median, holiday.Thu.count)
kwh.Fri.median <- median(subset(readings.aggregate, dayname == "Fri")$kwh)
holiday.Fri.count <- nrow(subset(readings.aggregate, dayname == "Fri" & holiday == TRUE)) / 24
dayname.holiday.datafr[6,] <- c("Friday", kwh.Fri.median, holiday.Fri.count)
kwh.Sat.median <- median(subset(readings.aggregate, dayname == "Sat")$kwh)
holiday.Sat.count <- nrow(subset(readings.aggregate, dayname == "Sat" & holiday == TRUE)) / 24
dayname.holiday.datafr[7,] <- c("Saturday", kwh.Sat.median, holiday.Sat.count)
xtable(dayname.holiday.datafr)

nonholiday_dayname.plot <- (ggplot(subset(readings.aggregate, subset = holiday == FALSE), 
                                   aes(dayname, kwh)) + 
                              geom_boxplot(aes(fill = dayname)) + 
                              theme(legend.position = "none") + 
                              ThesisPlotTheme() + 
                              labs(x = "Day of Week", 
                                   y = "Aggregate Electricity Demand (kWh)",
                                   title= "Distribution of Aggregate Electricity Demand by Day of Week (non-holidays)"))
working_day.plot <- (ggplot(readings.aggregate, aes(working_day, kwh)) + 
                       geom_boxplot() + 
                       ThesisPlotTheme() + 
                       labs(x = "Working Day", 
                            y = "Aggregate Electricity Demand (kWh)",
                            title= "Distribution of Aggregate Electricity Demand by Working Day Indicator"))
grid.arrange(nonholiday_dayname.plot, working_day.plot, ncol = 2, nrow = 1)
dev.print(file="../../Figures/DemandNonHolidayAndWorkingDay.png", device=png, height = 600, width = 1600)

kwh.wd.median <- median(subset(readings.aggregate, working_day == TRUE)$kwh)
print(kwh.wd.median)
kwh.nonwd.median <- median(subset(readings.aggregate, working_day == FALSE)$kwh)
print(kwh.nonwd.median)
dayname.nonholiday.datafr <- data.frame(dayname = character(), 
                                     median = numeric(), 
                                     holidays = integer(),
                                     stringsAsFactors = FALSE)
kwh.Sun.nonholiday.median <- median(subset(readings.aggregate, dayname == "Sun" & holiday == FALSE)$kwh)
dayname.nonholiday.datafr[1,] <- c("Sunday", kwh.Sun.nonholiday.median, 0)
kwh.Mon.nonholiday.median <- median(subset(readings.aggregate, dayname == "Mon" & holiday == FALSE)$kwh)
dayname.nonholiday.datafr[2,] <- c("Monday", kwh.Mon.nonholiday.median, 0)
kwh.Tue.nonholiday.median <- median(subset(readings.aggregate, dayname == "Tue" & holiday == FALSE)$kwh)
dayname.nonholiday.datafr[3,] <- c("Tuesday", kwh.Tue.nonholiday.median, 0)
kwh.Wed.nonholiday.median <- median(subset(readings.aggregate, dayname == "Wed" & holiday == FALSE)$kwh)
dayname.nonholiday.datafr[4,] <- c("Wednesday", kwh.Wed.nonholiday.median, 0)
kwh.Thu.nonholiday.median <- median(subset(readings.aggregate, dayname == "Thu" & holiday == FALSE)$kwh)
dayname.nonholiday.datafr[5,] <- c("Thursday", kwh.Thu.nonholiday.median, 0)
kwh.Fri.nonholiday.median <- median(subset(readings.aggregate, dayname == "Fri" & holiday == FALSE)$kwh)
dayname.nonholiday.datafr[6,] <- c("Friday", kwh.Fri.nonholiday.median, 0)
kwh.Sat.nonholiday.median <- median(subset(readings.aggregate, dayname == "Sat" & holiday == FALSE)$kwh)
dayname.nonholiday.datafr[7,] <- c("Saturday", kwh.Sat.nonholiday.median, 0)
xtable(dayname.nonholiday.datafr)

lm.temporal <- lm(kwh ~ hrstr + working_day + hrstr:working_day + rateseason,
                  data = readings.aggregate)
xtable(summary(lm.temporal))



#####
# Simple plot of Temperature vs. kWh
#####
temp.vs.demand.plot <- (ggplot(readings.aggregate, 
                               aes(x = temperature, y = kwh)) + 
                          geom_point(alpha = 0.5) + 
                          ThesisPlotTheme() + 
                          ylim(0, max(readings.aggregate$kwh)) + 
                          labs(x = "Outdoor, Dry-bulb Temperature (Celsius)", 
                               y = "Aggregate Electricity Demand (kWh)", 
                               title= "Aggregate Electricity Demand as a Function of Dry-bulb Temperature"))
temp.vs.demand.plot
dev.print(file="../../Figures/TemperatureVsDemand.png", device=png, height = 600, width = 800)



#####
# Conditioned plot of temperature
#####
lm.lintemp <- lm(kwh ~ temperature + hrstr + working_day + hrstr:working_day,
                 data = readings.aggregate)
temperature.visreg <- visreg(fit = lm.lintemp, xvar = "temperature", type = "conditional",
                             alpha = 1, # no confidence interval
                             line.par = list(col = TransparentColor()),
                             points.par = list(col = ScatterBlack()), 
                             cond = list(hrstr="h0", working_day="FALSE"), 
                             ylim = c(0, max(readings.aggregate$kwh)), 
                             main = "Conditional Aggregate Electricity Demand \n as a Function of Dry-bulb Temperature", 
                             ylab = "Conditional Aggregate Electricity Demand (kWh)",
                             xlab = "Outdoor, Dry-bulb Temperature (Celsius)")
dev.print(file="../../Figures/ConditionedTemperatureVsDemand.png", device=png, height = 600, width = 800)



#####
# Linear fit to untransformed, observed dry bulb temperature
#####
lintemp.visreg <- visreg(fit = lm.lintemp, xvar = "temperature", type = "conditional",
                         alpha = 1, # no confidence interval
                         line.par = list(col = "red", lwd = 3),
                         points.par = list(col = ScatterBlack()), 
                         cond = list(hrstr="h0", working_day="FALSE"), 
                         ylim = c(0, max(readings.aggregate$kwh)), 
                         main = "Conditional Aggregate Electricity Demand as a Function of Temperature", 
                         ylab = "Conditional Aggregate Electricity Demand (kWh)",
                         xlab = "Outdoor, Dry-bulb Temperature (Celsius)")
legend("topleft", 
       bty = "n", 
       lty = c(1), 
       lwd = c(3), 
       col = c("red"), 
       legend = c("Linear Regression Fit"))
dev.print(file="../../Figures/TemperatureTransformLinear.png", device=png, height = 600, width = 800)



#####
# Switching regression plot
#####
seg1 <- segmented(obj = lm.lintemp,
                  seg.Z = ~ temperature,
                  psi = 17)
drybulb.break <- seg1$psi[1,2]
seg1.adj.r2 <- summary(seg1)$adj.r.squared
plot(seg1, col = "red", lwd = 3, 
     res = TRUE, res.col = ScatterBlack(), pch = 20, cex = 0.6, 
     ylim = c(0, max(readings.aggregate$kwh)), 
     main = "Temperature Switching Regression Fit to Conditional Aggregate Electricity Demand",
     xlab = "Outdoor, Dry-bulb Temperature (Celsius)",
     ylab = "Conditional Aggregate Electricity Demand (kWh)",
     rug = FALSE)
abline(v = drybulb.break, col = "blue", lty = 2, lwd = 1)
legend("topleft", 
       bty = "n",
       lty = c(1, 2), 
       lwd = c(3, 1), 
       col = c("red", "blue"), 
       legend = c("Switching Regression Fit", "Knot Position"))
dev.print(file="../../Figures/TemperatureTransformSwitchingRegression.png", device=png, height = 600, width = 800)



#####
# Natural cubic splines
#####
lm.ns <- lm(kwh ~ ns(temperature, knots = c(3, 23, 30)) + hrstr + working_day + hrstr:working_day, 
            data = readings.aggregate)
ns.adj.r2 <- summary(lm.ns)$adj.r.squared
ns.adj.r2
ns.visreg <- visreg(fit = lm.ns, xvar = "temperature", type = "conditional",
                    alpha = 1, # no confidence interval
                    line.par = list(col = "red", lwd = 3),
                    points.par = list(col = ScatterBlack()), 
                    cond = list(hrstr="h0", working_day="FALSE"), 
                    ylim = c(0, max(readings.aggregate$kwh)), 
                    main = "Temperature Natural Cubic Splines Fit to Conditional Aggregate Electricity Demand", 
                    ylab = "Conditional Aggregate Electricity Demand (kWh)",
                    xlab = "Outdoor, Dry-bulb Temperature (Celsius)")
ns.temp <- ns(readings.aggregate$temperature, knots = c(3, 23, 30))
abline(v = attr(ns.temp, "knots"), col = "blue", lty = 2, lwd = 1)
legend("topleft", 
       bty = "n",
       lty = c(1, 2), 
       lwd = c(3, 1), 
       col = c("red", "blue"), 
       legend = c("Natural Cubic Spline Fit", "Knot Positions"))
dev.print(file="../../Figures/TemperatureTransformationNaturalSplines.png", device=png, height = 600, width = 800)



#####
# Dry-bulb moving average
#####
ma.order <- 6
ma.trimsize <- (ma.order - 1)
# rollingmean(...) averages elements [i, ..., i+ma.order] and then trims 
# ma.order-1 from the end of the vector
ma.drybulb <- rollmean(x = readings.aggregate$temperature, k = ma.order)
# To align the rollingmean(...) vector to be elements [i-1-ma.order, ..., i], 
# the first ma.order-1 elements should be trimmed from the head of the 
# observation data.frame
readings.ma.pretrim <- tail(readings.aggregate, -ma.trimsize)
readings.ma.pretrim$rolling_average <- ma.drybulb
lm.rollingaverage <- lm(kwh ~ rolling_average + hrstr + working_day + hrstr:working_day, 
            data = readings.ma.pretrim)
rollingaverage.visreg <- visreg(fit = lm.rollingaverage, xvar = "rolling_average", type = "conditional",
                             alpha = 1, # no confidence interval
                             line.par = list(col = TransparentColor()),
                             points.par = list(col = ScatterBlack()), 
                             cond = list(hrstr="h0", working_day="FALSE"), 
                             ylim = c(0, max(readings.aggregate$kwh)), 
                             main = "Conditional Aggregate Electricity Demand \n as a Function of a Six-hour Moving Average of Dry-bulb Temperature", 
                             ylab = "Conditional Aggregate Electricity Demand (kWh)",
                             xlab = "Six-hour Moving Average of Outdoor, Dry-bulb Temperature (Celsius)")
dev.print(file="../../Figures/ConditionedRollingAverageVsDemand.png", device=png, height = 600, width = 800)



#####
# Dry-bulb degree-hours
#####
dh.window <- 6
dh.trimsize <- (dh.window - 1)

dh.drybulb.basismatrix <- DegreeHourBasisMatrix(temperatures = readings.aggregate$temperature, 
                                                balance = drybulb.break, 
                                                window = dh.window)
readings.dh.trim <- tail(cbind(readings.aggregate, dh.drybulb.basismatrix), -dh.trimsize)
lm.degreehour <- lm(kwh ~ CDH + HDH + hrstr + working_day + hrstr:working_day, 
                        data = readings.dh.trim)
par(mfrow = c(1, 2))
hdh.visreg <- visreg(fit = lm.degreehour, xvar = "HDH", type = "conditional",
                     alpha = 1, # no confidence interval
                     line.par = list(col = "red", lwd = 3),
                     points.par = list(col = ScatterBlack()), 
                     cond = list(hrstr="h0", working_day="FALSE", CDH=0), 
                     ylim = c(0, max(readings.aggregate$kwh)), 
                     main = "Conditional Aggregate Electricity Demand vs. HDH", 
                     ylab = "Conditional Aggregate Electricity Demand (kWh)",
                     xlab = "Heating Degree-Hours (Six-hour Window, Dry-bulb)")
cdh.visreg <- visreg(fit = lm.degreehour, xvar = "CDH", type = "conditional",
                     alpha = 1, # no confidence interval
                     line.par = list(col = "red", lwd = 3),
                     points.par = list(col = ScatterBlack()), 
                     cond = list(hrstr="h0", working_day="FALSE", HDH=0), 
                     ylim = c(0, max(readings.aggregate$kwh)), 
                     main = "Conditional Aggregate Electricity Demand vs. CDH", 
                     ylab = "Conditional Aggregate Electricity Demand (kWh)",
                     xlab = "Cooling Degree-Hours (Six-hour Window, Dry-bulb)")
dev.print(file="../../Figures/ConditionedDegreeHoursRegression.png", device=png, height = 600, width = 1600)
par(mfrow = c(1, 1))



#####
# Exposure-lag-response contour plot
#####
# See if using breakpoints as threshold values improves external-response
cb <- crossbasis(x = readings.aggregate$temperature, 
                 argvar = list(fun = "ns", knots = c(10, 24, 25)),
                 lag = 12, 
                 arglag = list(fun = "poly", degree=3))

# Prediction using the crossbasix matrix
lm.cb <- lm(kwh ~ cb + hrstr + working_day + price + hrstr:working_day,
            data = readings.aggregate)
summary(lm.cb)$adj.r.squared
cb.pred <- crosspred(cb, lm.cb)
plot(cb.pred, ptype="contour", 
     key.title = title("kWh"), 
     plot.title = title("Contour Plot of Exposure-Lag-Response Association\n using Aggregate Electricity Demand", 
                        xlab = "Outdoor, Dry-bulb Temperature (Celsius)", 
                        ylab = "Lag (hours)"))
dev.print(file = "../../Figures/ExposureLagResponseAssociationContourPlot.png", device = png, height = 600, width = 800)



#####
# Conditioned plot of feels like temperature and electricity demand
#####
readings.aggregate$feels_like <- FeelsLike(readings.aggregate)
lm.feelslike <- lm(kwh ~ feels_like + hrstr + working_day + hrstr:working_day,
                 data = readings.aggregate)
temperature.visreg <- visreg(fit = lm.feelslike, xvar = "feels_like", type = "conditional",
                             alpha = 1, # no confidence interval
                             line.par = list(col = TransparentColor()),
                             points.par = list(col = ScatterBlack()), 
                             cond = list(hrstr="h0", working_day="FALSE"), 
                             ylim = c(0, max(readings.aggregate$kwh)), 
                             main = "Conditional Aggregate Electricity Demand as a Function of \"Feels Like\" Temperature", 
                             ylab = "Conditional Aggregate Electricity Demand (kWh)",
                             xlab = "Outdoor, \"Feels Like\" Temperature (Celsius)")
dev.print(file="../../Figures/ConditionedFeelsLikeVsDemand.png", device=png, height = 600, width = 800)




#####
# Variance inflation factor of seasonality terms
#####
lm.ns.noseason <- lm(kwh ~ ns(temperature, knots = c(3, 23, 30)) + hrstr + working_day,
                      data = readings.aggregate)
xtable(vif(lm.ns.noseason))


lm.ns.withmonth <- lm(kwh ~ ns(temperature, knots = c(3, 23, 30)) + month + 
                        hrstr + working_day,
                      data = readings.aggregate)
xtable(vif(lm.ns.withmonth))

lm.ns.rateseason <- lm(kwh ~ ns(temperature, knots = c(3, 23, 30)) + 
                         rateseason + tou_active + hrstr + working_day,
                      data = readings.aggregate)
xtable(vif(lm.ns.rateseason))