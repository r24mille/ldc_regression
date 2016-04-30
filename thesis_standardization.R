HourLabels <- function() {
  # Labels for the hour-of-day, commonly used on the x-axis of plots and 
  # rotated.
  #
  # Return:
  #   A vector of 24 characters representing hour-of-day (ie. 00:00, ..., 
  #   23:00).
  
  hour_labels <- c("00:00", "01:00", "02:00", "03:00", "04:00", "05:00", 
                   "06:00", "07:00", "08:00", "09:00", "10:00", "11:00", 
                   "12:00", "13:00", "14:00", "15:00", "16:00", "17:00", 
                   "18:00", "19:00", "20:00", "21:00", "22:00", "23:00")
  return(hour_labels)
}

HourStrings <- function() {
  hour_strings <- c("h0", "h1", "h2", "h3", "h4", "h5", "h6", "h7", "h8", "h9", 
                    "h10", "h11", "h12", "h13", "h14", "h15", "h16", "h17", 
                    "h18", "h19", "h20", "h21", "h22", "h23")
  return(hour_strings)
}

NonWorkingDayTouSequenceDataframe <- function() {
  seq.datafr <- data.frame(tou.sequence = rep("Off-peak", 24),
                           hrstr = HourStrings())
  seq.datafr$tou.sequence <- factor(seq.datafr$tou.sequence, 
                                    levels = c("Off-peak", "Mid-peak", "On-peak"))
  seq.datafr$hrstr <- factor(seq.datafr$hrstr, levels = seq.datafr$hrstr)
  return(seq.datafr)
}

ScatterBlack <- function() {
  scatter.black = rgb(0, 0, 0, 50, maxColorValue = 100)
  return(scatter.black)
}

SummerWorkingDayTouSequenceDataframe <- function() {
  seq.datafr <- data.frame(tou.sequence = c("Off-peak", "Off-peak", "Off-peak", 
                                            "Off-peak", "Off-peak", "Off-peak", 
                                            "Off-peak",
                                            "Mid-peak", "Mid-peak", "Mid-peak", 
                                            "Mid-peak", 
                                            "On-peak", "On-peak", "On-peak", 
                                            "On-peak", "On-peak", "On-peak", 
                                            "Mid-peak", "Mid-peak", 
                                            "Off-peak", "Off-peak", "Off-peak", 
                                            "Off-peak", "Off-peak"),
                           hrstr = HourStrings())
  seq.datafr$tou.sequence <- factor(seq.datafr$tou.sequence, 
                                    levels = c("Off-peak", "Mid-peak", "On-peak"))
  seq.datafr$hrstr <- factor(seq.datafr$hrstr, levels = seq.datafr$hrstr)
  return(seq.datafr)
}

TouPricePalette <- function() {
  # Creates a vector suitable for use as a ggplot2 colour palette. The TOU rate 
  # colour palette (green, yellow, red):(off, mid, on) matches the HEX colour 
  # codes used by the Ontario Energy Board. This may help trigger some sort of 
  # recognition from readers within Ontario.
  #
  # Return:
  #   A 3-element, named vector of HEX colour characters.
  price.palette <- c("Off-peak" = "#c2da8b", 
                     "Mid-peak" = "#fcdf6f", 
                     "On-peak" = "#df9c7d")
  return(price.palette)
}

WinterWorkingDayTouSequenceDataframe <- function() {
  seq.datafr <- data.frame(tou.sequence = c("Off-peak", "Off-peak", "Off-peak", 
                                            "Off-peak", "Off-peak",  "Off-peak", 
                                            "Off-peak",
                                            "On-peak", "On-peak", "On-peak", 
                                            "On-peak", 
                                            "Mid-peak", "Mid-peak", "Mid-peak", 
                                            "Mid-peak", "Mid-peak", "Mid-peak", 
                                            "On-peak", "On-peak", 
                                            "Off-peak", "Off-peak", "Off-peak", 
                                            "Off-peak", "Off-peak"),
                           hrstr = HourStrings())
  seq.datafr$tou.sequence <- factor(seq.datafr$tou.sequence, 
                                    levels = c("Off-peak", "Mid-peak", "On-peak"))
  seq.datafr$hrstr <- factor(seq.datafr$hrstr, levels = seq.datafr$hrstr)
  return(seq.datafr)
}

ThesisPlotTheme <- function() {
  require(ggplot2)
  standard.theme <- theme(legend.position = "right",
                          legend.text = element_text(size = 14),
                          legend.title = element_text(size = 14),
                          axis.text = element_text(size = 14),
                          axis.title = element_text(size = 16),
                          plot.title = element_text(size = 20))
  return(standard.theme)
}

TransparentColor <- function() {
  transparent.col = rgb(0, 0, 0, 0, maxColorValue = 100)
  return(transparent.col)
}

PlotTouEffectResult <- function(effect.datafr, tou.datafr, ggp.title) {
  require("ggplot2")
  demand.palette <- c("#000000", "#377eb8") # Black (pre) and blue (post)
  
  theplot <- (ggplot(data = effect.datafr) + 
                geom_rect(aes(NULL, NULL, 
                              xmin = factor(hrstr, levels = hrstr),
                              xmax = c(tail(hrstr, -1), hrstr[24]),
                              ymin = -Inf,
                              ymax = Inf,
                              fill = tou.sequence),
                          data = tou.datafr, 
                          alpha = 0.6) + 
                geom_line(aes(x = hrstr, 
                              y = fit, 
                              group = factor(tou_active, labels=c("Flat rate", "TOU rates\n(counterfactual)")),
                              linetype = factor(tou_active, labels=c("Flat rate", "TOU rates\n(counterfactual)")),
                              colour = factor(tou_active, labels=c("Flat rate", "TOU rates\n(counterfactual)"))),
                          size = 1) + 
                geom_point(aes(x = hrstr, 
                               y = fit, 
                               group = factor(tou_active, 
                                              labels=c("Flat rate", "TOU rates\n(counterfactual)")),
                               shape = factor(tou_active, 
                                              labels=c("Flat rate", "TOU rates\n(counterfactual)")),
                               colour = factor(tou_active, 
                                               labels=c("Flat rate", "TOU rates\n(counterfactual)"))),
                           size = 4) + 
                geom_errorbar(aes(x = hrstr,
                                  y = fit,
                                  ymin = lower, 
                                  ymax = upper,
                                  group = factor(tou_active, labels=c("Flat rate", "TOU rates\n(counterfactual)")),
                                  colour = factor(tou_active, labels=c("Flat rate", "TOU rates\n(counterfactual)"))),
                              alpha = 0.6,
                              width = 0.3) + 
                labs(colour = "Pricing Method",
                     shape = "Pricing Method",
                     linetype = "Pricing Method", 
                     fill = "TOU Period") + 
                scale_colour_manual(values = demand.palette) + 
                scale_fill_manual(values = TouPricePalette()) + 
                scale_x_discrete(labels = HourLabels()) + 
                #ylim(0.5, 2) + 
                ylim(0, 3.1) + 
                labs(x = "Hour of Day",
                     y = "Average Hourly Household Electricity Demand (kWh)",
                     title = ggp.title) + 
                ThesisPlotTheme() + 
                theme(axis.text.x = element_text(angle = 90, 
                                                 vjust = 0.5),
                      panel.background = element_rect(fill = "#ffffff"),
                      panel.grid.major = element_line(colour = "#bbbbbb"),
                      panel.grid.minor = element_line(colour = "#dddddd"),
                      axis.ticks = element_line(colour = "#bbbbbb")))
  print(theplot)
  return(theplot)
}

GetConfident <- function(observed.datafr, counterfactual.datafr, prediction,
                         tou.idx, rate.season, working.day = TRUE) {
  series <- data.frame(hrstr = factor(levels = HourStrings()), 
                       tou_active = logical(), 
                       fit = numeric(), 
                       se = numeric(), 
                       lower = numeric(), 
                       upper = numeric())
  
  # Trim observations to the pre-TOU index for this case study
  observed.datafr <- observed.datafr[c(1:(tou.idx-1)),]
  
  # Subset the case study sample based on season and working day
  observed.datafr <- subset(observed.datafr, 
                            subset = rateseason == rate.season 
                            & working_day == working.day)
  
  # Group the mean response by hour of day
  mean.obs.resp <- by(data = observed.datafr,
                      INDICES = observed.datafr$hrstr,
                      FUN = function(x) mean(x$kwh))
  
  # Fill out observed portion of the results data.frame
  for(i in 1:24) {
    series[i, "hrstr"] <- paste0("h", (i-1))
    series[i, "tou_active"] <- FALSE
    series[i, "fit"] <- mean.obs.resp[i]
    series[i, "se"] <- NA
    series[i, "lower"] <- NA
    series[i, "upper"] <- NA
  }
  
  
  # Replace values in the counterfactual data.frame with values from the 
  # prediction so that indexed operations can be run
  counterfactual.datafr$kwh <- counterfactual.predict$fit[,"fit"]
  counterfactual.datafr$lower <- counterfactual.predict$fit[,"lwr"]
  counterfactual.datafr$upper <- counterfactual.predict$fit[,"upr"]
  
  # Trim counterfactual to the case study length
  counterfactual.datafr <- counterfactual.datafr[c(1:(tou.idx-1)),]
  
  # Subset sample based on season and working day
  counterfactual.datafr <- subset(counterfactual.datafr, 
                                  subset = rateseason == rate.season 
                                  & working_day == working.day)
  
  mean.counterfactual.fit <- by(data = counterfactual.datafr,
                                INDICES = counterfactual.datafr$hrstr,
                                FUN = function(x) mean(x$kwh))
  mean.counterfactual.lower <- by(data = counterfactual.datafr,
                                  INDICES = counterfactual.datafr$hrstr,
                                  FUN = function(x) mean(x$lower))
  mean.counterfactual.upper <- by(data = counterfactual.datafr,
                                  INDICES = counterfactual.datafr$hrstr,
                                  FUN = function(x) mean(x$upper))
  
  # Fill out counterfactual portion of the results data.frame
  for(i in 1:24) {
    series[i+24, "hrstr"] <- paste0("h", (i-1))
    series[i+24, "tou_active"] <- TRUE
    series[i+24, "fit"] <- mean.counterfactual.fit[i]
    series[i+24, "se"] <- NA
    series[i+24, "lower"] <- mean.counterfactual.lower[i]
    series[i+24, "upper"] <- mean.counterfactual.upper[i]
  }
  
  return(series)
}

OnToOffPeakRatio <- function(observed.datafr, counterfactual.datafr, 
                             prediction, tou.idx) {
  # Replace values in the counterfactual data.frame with values from the 
  # prediction so that indexed operations can be run
  counterfactual.datafr$kwh <- counterfactual.predict$fit[,"fit"]
  
  # Trim observations and counterfactual to the case study length
  observed.datafr <- observed.datafr[c(1:(tou.idx-1)),]
  counterfactual.datafr <- counterfactual.datafr[c(1:(tou.idx-1)),]
  
  # Subset data.frames to summer
  observed.summer.wd <- subset(observed.datafr, 
                               subset = rateseason == "summer"
                               & working_day == TRUE)
  counterfactual.summer.wd <- subset(counterfactual.datafr, 
                                     subset = rateseason == "summer"
                                     & working_day == TRUE)
  
  # Get average daily on- to off-peak ratio
  observed.onoff.ratios <- by(data = observed.summer.wd,
                              INDICES = observed.summer.wd$daynum,
                              FUN = function(x) {
                                on <- subset(x, subset = dayperiod == "afternoon")
                                off <- subset(x, subset = dayperiod == "overnight")
                                return(sum(on$kwh)/sum(off$kwh))
                              })
  observed.onoff.ratio.mean <- mean(observed.onoff.ratios)
  
  counterfactual.onoff.ratios <- by(data = counterfactual.summer.wd,
                                    INDICES = counterfactual.summer.wd$daynum,
                                    FUN = function(x) {
                                      on <- subset(x, subset = dayperiod == "afternoon")
                                      off <- subset(x, subset = dayperiod == "overnight")
                                      return(sum(on$kwh)/sum(off$kwh))
                                    })
  counterfactual.onoff.ratio.mean <- mean(counterfactual.onoff.ratios)
  onoff.ratio.pctchange <- (counterfactual.onoff.ratio.mean/observed.onoff.ratio.mean * 100) - 100
  
  onoff.ratio.results <- data.frame("flat_onoff_ratio" = observed.onoff.ratio.mean,
                                    "tou_onoff_ratio" = counterfactual.onoff.ratio.mean,
                                    "onoff_ratio_pctchange" = onoff.ratio.pctchange)
  return(onoff.ratio.results)
}

PeakToAverageRatio <- function(observed.datafr, counterfactual.datafr, 
                               prediction, tou.idx) {
  # Replace values in the counterfactual data.frame with values from the 
  # prediction so that indexed operations can be run
  counterfactual.datafr$kwh <- prediction$fit[,"upr"]
  
  # Trim observations and counterfactual to the case study length
  observed.datafr <- observed.datafr[c(1:(tou.idx-1)),]
  counterfactual.datafr <- counterfactual.datafr[c(1:(tou.idx-1)),]
  
  # Subset data.frames to summer
  observed.summer.wd <- subset(observed.datafr, 
                               subset = rateseason == "summer"
                               & working_day == TRUE)
  counterfactual.summer.wd <- subset(counterfactual.datafr, 
                                     subset = rateseason == "summer"
                                     & working_day == TRUE)
  
  # Get average daily on- to off-peak ratio
  observed.peakavg.ratios <- by(data = observed.summer.wd,
                                INDICES = observed.summer.wd$daynum,
                                FUN = function(x) {
                                  return(max(x$kwh)/mean(x$kwh))
                                })
  observed.peakavg.ratio.mean <- mean(observed.peakavg.ratios)
  
  counterfactual.peakavg.ratios <- by(data = counterfactual.summer.wd,
                                      INDICES = counterfactual.summer.wd$daynum,
                                      FUN = function(x) {
                                        return(max(x$kwh)/mean(x$kwh))
                                      })
  counterfactual.peakavg.ratio.mean <- mean(counterfactual.peakavg.ratios)
  peakavg.ratio.pctchange <- (observed.peakavg.ratio.mean/counterfactual.peakavg.ratio.mean * 100) - 100
  
  peakavg.ratio.results <- data.frame("flat_peakavg_ratio" = observed.peakavg.ratio.mean,
                                      "tou_peakavg_ratio" = counterfactual.peakavg.ratio.mean,
                                      "peakavg_ratio_pctchange" = peakavg.ratio.pctchange)
}

EffectsTable <- function(observed.datafr, counterfactual.datafr, prediction,
                         tou.idx) {
  effects.datafr <- data.frame(hourly_impact = numeric(), 
                               hourly_ci_upper = numeric(), 
                               hourly_ci_lower = numeric(), 
                               period_impact = numeric(), 
                               period_ci_upper = numeric(), 
                               period_ci_lower = numeric(),
                               percent_impact = numeric(), 
                               percent_ci_upper =numeric(),
                               percent_ci_lower = numeric())
  
  # Trim observations to the pre-TOU index for this case study
  observed.datafr <- observed.datafr[c(1:(tou.idx-1)),]
  
  # Replace values in the counterfactual data.frame with values from the 
  # prediction so that indexed operations can be run
  counterfactual.datafr$kwh <- counterfactual.predict$fit[,"fit"]
  counterfactual.datafr$lower <- counterfactual.predict$fit[,"lwr"]
  counterfactual.datafr$upper <- counterfactual.predict$fit[,"upr"]
  
  # Trim counterfactual to the case study length
  counterfactual.datafr <- counterfactual.datafr[c(1:(tou.idx-1)),]
  
  # Subset the case study sample based on season and working day
  observed.summer <- subset(observed.datafr, 
                            subset = rateseason == "summer")
  observed.winter <- subset(observed.datafr, 
                            subset = rateseason == "winter")
  counterfactual.summer <- subset(counterfactual.datafr, 
                                  subset = rateseason == "summer")
  counterfactual.winter <- subset(counterfactual.datafr, 
                                  subset = rateseason == "winter")
  
  # Print conservation across all hours
  print(paste0("Summer conservation: ", 100 - (sum(counterfactual.summer$kwh) / sum(observed.summer$kwh)) * 100))
  print(paste0("Winter conservation: ", 100 - (sum(counterfactual.winter$kwh) / sum(observed.winter$kwh)) * 100))
  
  # Group the mean response by hour of day
  mean.obs.summer <- by(data = observed.summer,
                        INDICES = observed.summer$dayperiod,
                        FUN = function(x) mean(x$kwh))
  mean.obs.winter <- by(data = observed.winter,
                        INDICES = observed.winter$dayperiod,
                        FUN = function(x) mean(x$kwh))
  
  mean.counterfactual.summer <- by(data = counterfactual.summer,
                                   INDICES = counterfactual.summer$dayperiod,
                                   FUN = function(x) mean(x$kwh))
  mean.counterfactual.winter <- by(data = counterfactual.winter,
                                   INDICES = counterfactual.winter$dayperiod,
                                   FUN = function(x) mean(x$kwh))
  
  mean.counterfactual.upper.summer <- by(data = counterfactual.summer,
                                         INDICES = counterfactual.summer$dayperiod,
                                         FUN = function(x) mean(x$upper))
  mean.counterfactual.lower.summer <- by(data = counterfactual.summer,
                                         INDICES = counterfactual.summer$dayperiod,
                                         FUN = function(x) mean(x$lower))
  
  mean.counterfactual.upper.winter <- by(data = counterfactual.winter,
                                         INDICES = counterfactual.winter$dayperiod,
                                         FUN = function(x) mean(x$upper))
  mean.counterfactual.lower.winter <- by(data = counterfactual.winter,
                                         INDICES = counterfactual.winter$dayperiod,
                                         FUN = function(x) mean(x$lower))
  
  # Add summer results to the dataframe
  effects.datafr[c(1:4), "hourly_impact"] <- mean.counterfactual.summer-mean.obs.summer
  effects.datafr[c(1:4), "hourly_ci_upper"] <- mean.counterfactual.upper.summer-mean.counterfactual.summer
  effects.datafr[c(1:4), "hourly_ci_lower"] <- mean.counterfactual.lower.summer-mean.counterfactual.summer
  
  rownames(effects.datafr) <- paste0("summer ", rownames(mean.obs.summer))
  
  effects.datafr["summer afternoon", "period_impact"] <- effects.datafr["summer afternoon", "hourly_impact"] * 6
  effects.datafr["summer nonwd", "period_impact"] <- effects.datafr["summer nonwd", "hourly_impact"] * 48
  effects.datafr["summer overnight", "period_impact"] <- effects.datafr["summer afternoon", "hourly_impact"] * 12
  effects.datafr["summer transition", "period_impact"] <- effects.datafr["summer afternoon", "hourly_impact"] * 6
  
  effects.datafr["summer afternoon", "period_ci_upper"] <- effects.datafr["summer afternoon", "hourly_ci_upper"] * 6
  effects.datafr["summer nonwd", "period_ci_upper"] <- effects.datafr["summer nonwd", "hourly_ci_upper"] * 48
  effects.datafr["summer overnight", "period_ci_upper"] <- effects.datafr["summer afternoon", "hourly_ci_upper"] * 12
  effects.datafr["summer transition", "period_ci_upper"] <- effects.datafr["summer afternoon", "hourly_ci_upper"] * 6
  
  effects.datafr["summer afternoon", "period_ci_lower"] <- effects.datafr["summer afternoon", "hourly_ci_lower"] * 6
  effects.datafr["summer nonwd", "period_ci_lower"] <- effects.datafr["summer nonwd", "hourly_ci_lower"] * 48
  effects.datafr["summer overnight", "period_ci_lower"] <- effects.datafr["summer afternoon", "hourly_ci_lower"] * 12
  effects.datafr["summer transition", "period_ci_lower"] <- effects.datafr["summer afternoon", "hourly_ci_lower"] * 6
  
  effects.datafr[c(1:4), "percent_impact"] <- ((mean.counterfactual.summer/mean.obs.summer) - 1) * 100
  effects.datafr[c(1:4), "percent_ci_upper"] <- ((mean.counterfactual.upper.summer/mean.counterfactual.summer) - 1) * 100
  effects.datafr[c(1:4), "percent_ci_lower"] <- ((mean.counterfactual.lower.summer/mean.counterfactual.summer) - 1) * 100
  
  # Add winter results to the dataframe
  effects.datafr[c(5:8), "hourly_impact"] <- mean.counterfactual.winter-mean.obs.winter
  effects.datafr[c(5:8), "hourly_ci_upper"] <- mean.counterfactual.upper.winter-mean.counterfactual.winter
  effects.datafr[c(5:8), "hourly_ci_lower"] <- mean.counterfactual.lower.winter-mean.counterfactual.winter
  
  rownames(effects.datafr)[c(5:8)] <- paste0("winter ", rownames(mean.obs.winter))
  
  effects.datafr["winter afternoon", "period_impact"] <- effects.datafr["winter afternoon", "hourly_impact"] * 6
  effects.datafr["winter nonwd", "period_impact"] <- effects.datafr["winter nonwd", "hourly_impact"] * 48
  effects.datafr["winter overnight", "period_impact"] <- effects.datafr["winter afternoon", "hourly_impact"] * 12
  effects.datafr["winter transition", "period_impact"] <- effects.datafr["winter afternoon", "hourly_impact"] * 6
  
  effects.datafr["winter afternoon", "period_ci_upper"] <- effects.datafr["winter afternoon", "hourly_ci_upper"] * 6
  effects.datafr["winter nonwd", "period_ci_upper"] <- effects.datafr["winter nonwd", "hourly_ci_upper"] * 48
  effects.datafr["winter overnight", "period_ci_upper"] <- effects.datafr["winter afternoon", "hourly_ci_upper"] * 12
  effects.datafr["winter transition", "period_ci_upper"] <- effects.datafr["winter afternoon", "hourly_ci_upper"] * 6
  
  effects.datafr["winter afternoon", "period_ci_lower"] <- effects.datafr["winter afternoon", "hourly_ci_lower"] * 6
  effects.datafr["winter nonwd", "period_ci_lower"] <- effects.datafr["winter nonwd", "hourly_ci_lower"] * 48
  effects.datafr["winter overnight", "period_ci_lower"] <- effects.datafr["winter afternoon", "hourly_ci_lower"] * 12
  effects.datafr["winter transition", "period_ci_lower"] <- effects.datafr["winter afternoon", "hourly_ci_lower"] * 6
  
  effects.datafr[c(5:8), "percent_impact"] <- ((mean.counterfactual.winter/mean.obs.winter) - 1) * 100
  effects.datafr[c(5:8), "percent_ci_upper"] <- ((mean.counterfactual.upper.winter/mean.counterfactual.winter) - 1) * 100
  effects.datafr[c(5:8), "percent_ci_lower"] <- ((mean.counterfactual.lower.winter/mean.counterfactual.winter) - 1) * 100
  
  return(effects.datafr)
}