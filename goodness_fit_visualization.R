PlotRegSubSets <- function (bics, adjr2s, title) {
  # Plot results from reg_subset in Bayesian Model Averaging (BMA) package.
  #
  # Args:
  #   bics: A vector of BICs to be plotted
  #   adjr2s: A vector of adjusted R^2 values corresponding to the BICs
  #   title: The main title of the plot
  
  # Adjust padding for labelling on the y2 axis
  par(mar=c(5,4,4,5)+.1)
  xrng = c(1, length(bics))
  rng.bics = c(min(bics), max(bics))
  rng.adjr2s = c(min(adjr2s), max(adjr2s))
  
  # Plot y1
  plot(bics, 
       type = "l", 
       col = "red", 
       ylab = "Bayesian Information Criterion (BIC)", 
       ylim = rng.bics, 
       xlab = "Number of Variables", 
       xlim = xrng, 
       main = title)
  
  # 2nd plot over the existing plot
  par(new=TRUE)
  
  # Ploy y2
  plot(adjr2s, 
       type = "l", 
       col = "blue", 
       xaxt = "n", 
       xlab = "", 
       xlim = xrng, 
       yaxt = "n", 
       ylab = "", 
       ylim = rng.adjr2s)
  axis(4)
  mtext(expression(paste("Variance Explained (Adj.  ", R^2 ,")")),
        side = 4, 
        line = 3)
  legend("topright", 
         col = c("red","blue"), 
         lty = 1, 
         legend = c("BIC",
                    expression(paste("Adj.  ", R^2))))
  
  # 3rd plot over the existing plot
  par(new=TRUE)
  
  # Create a point and label the lowest BIC value
  minbic <- round(min(bics),
                  digits = 2)
  minbic.index <- which.min(bics)
  
  # Add vertical line showing the number of variables
  abline(v = minbic.index,
         xaxt = "n", 
         xlab = "", 
         xlim = xrng, 
         yaxt = "n", 
         ylab = "", 
         lty = 5, # dashed
         col = "gray")
  
  # 4th plot over the existing plot
  par(new=TRUE)
  
  plot(x = minbic.index, 
       xaxt = "n", 
       xlab = "", 
       xlim = xrng, 
       yaxt = "n", 
       ylab = "", 
       y = minbic, 
       ylim = rng.bics, 
       type = "p", 
       pch = 20, # dot
       col = "red")
  
  text(x = minbic.index, 
       xlim = xrng,
       y = minbic, 
       ylim = rng.bics, 
       labels = c(minbic), 
       pos = 3, # above
       col = "black") 
  
  axis(1, #bottom
       at = c(minbic.index))

  # 5th plot over the existing plot
  par(new=TRUE)
  
  # Create a point and label the lowest BIC value
  adjr2.atbic <- round(adjr2s[minbic.index], 
                       digits = 3)
  plot(x = minbic.index, 
       xaxt = "n", 
       xlab = "", 
       xlim = xrng, 
       yaxt = "n", 
       ylab = "", 
       y = adjr2.atbic, 
       ylim = rng.adjr2s, 
       type = "p", 
       pch = 20, # dot
       col = "blue")
  text(x = minbic.index, 
       xlim = xrng,
       y = adjr2.atbic, 
       ylim = rng.adjr2s, 
       labels = c(adjr2.atbic), 
       pos = 1, # below
       col = "black") 
}

PlotGlmFitMeasures <- function (aiccs, bics, y2vals, xvals, 
                                y2title, xtitle, title) {
  # Plots several measures of model fit on a dual y-axis plot. It shows how 
  # information criteria (IC) and residual deviance change as parameters 
  # within the model change iteratively.
  #
  # Args:
  #   aiccs: A vector of AICc values to be plotted
  #   bics: A vector of BIC values to be plotted
  #   y2vals: A vector of values to be plotted on the y2-axis (eg. residual 
  #           deviations or pseudo R^2.
  #   xvals: The values to be plotted along the x-axis. In the context of this 
  #          function they generally represent iterative changes to a model 
  #          parameter.
  #   y2title: The title of the y2-axis
  #   xtitle: The title of the x-axis
  #   title: The main title of the plot
  
  # Adjust padding for labelling on the y2 axis
  par(mar=c(5,4,4,5)+.1)
  xrng = c(min(xvals), max(xvals))
  y1rng = c(min(c(bics, aiccs)),
            max(c(bics, aiccs)))
  y2rng = c(min(y2vals), max(y2vals))
  
  # Colors, helped by http://colorbrewer2.org/
  colaicc = "#4daf4a" # greenish
  colbic = "#377eb8" # blueish
  coly2 = "#e41a1c" # reddish
  
  # y1, plot AICc values
  plot(x = xvals,
       y = aiccs, 
       type = "l", 
       col = colaicc, 
       xlab = xtitle, 
       xlim = xrng, 
       ylab = "Information Criteria (AICc and BIC)", 
       ylim = y1rng, 
       main = title,
       lwd = 2) # line width
  
  # y1, plot BIC values
  lines(x = xvals,
        y = bics,
        col = colbic,
        lwd = 2) # line width
  
  # 2nd plot over the existing plot
  par(new=TRUE)
  
  plot(x = xvals,
       y = y2vals, 
       type = "l", 
       col = coly2, 
       xaxt = "n", 
       xlab = "", 
       xlim = xrng, 
       yaxt = "n", 
       ylab = "", 
       ylim = y2rng,
       lwd = 2) # line width
  axis(4)
  mtext(y2title,
        side = 4, 
        line = 3)
  
  # 3rd plot over the existing plot
  par(new=TRUE)
  
  # Create a point and label the lowest BIC value
  minbic <- round(min(bics),
                  digits = 1)
  minbic.index <- which.min(bics)
  
  # Add vertical line showing the number of variables
  abline(v = xvals[minbic.index],
         xaxt = "n", 
         xlab = "", 
         xlim = xrng, 
         yaxt = "n", 
         ylab = "", 
         lty = 5, # dashed
         col = "gray")
  
  # 4th plot over the existing plot
  par(new=TRUE)
  
  plot(x = xvals[minbic.index], 
       xaxt = "n", 
       xlab = "", 
       xlim = xrng, 
       yaxt = "n", 
       ylab = "", 
       y = minbic, 
       ylim = y1rng, 
       type = "p", 
       pch = 20, # dot
       col = colbic)
  
  text(x = xvals[minbic.index], 
       xlim = xrng,
       y = minbic, 
       ylim = y1rng, 
       labels = c(minbic), 
       pos = 3, # above
       col = "black") 
  
  axis(1, #bottom
       at = c(xvals[minbic.index]))
  
  # 5th plot over the existing plot
  par(new=TRUE)
  
  # Create a point and label the y2 value coinciding with the lowest BIC
  y2.atbic <- y2vals[minbic.index]
  plot(x = xvals[minbic.index],
       xaxt = "n", 
       xlab = "", 
       xlim = xrng, 
       yaxt = "n", 
       ylab = "", 
       y = y2.atbic,
       ylim = y2rng, 
       type = "p", 
       pch = 20, # dot
       col = coly2)
  text(x = xvals[minbic.index], 
       xlim = xrng,
       y = y2.atbic, 
       ylim = y2rng, 
       labels = c(round(y2.atbic, digits = 3)), 
       pos = 3, # above
       col = "black") 
  
  legend("topright", 
         col = c(colaicc, colbic, coly2), 
         lty = 1, 
         lwd = 2, # line width
         legend = c("AICc",
                    "BIC",
                    y2title))
}

Plot2DFitByExplVarCountWithMultiplePastHrsTemp <- function(df.steps, is.bic, 
                                                           title, subtitle) {
  # Plots 2D representation of a goodness-of-fit measure as a function of the
  # explanatory variable count. Multiple lines are plotted, each creating a 
  # plot for the number of past hours' temperature measurements included.
  #
  # Args:
  #   df.steps: A dataframe of stepwise goodness-of-fit results.
  #   is.bic: A boolean indicating whether the goodness-of-fit measure is BIC.
  #           If FALSE, then Pseudo-R^2 is used.
  #   title: The title of the plot, to be passed to the "main" parameter.
  #   subtitle: The subtitle of the plot, to be passed to the "sub" parameter.
  xrng <- c(min(df.steps$num.explvars), max(df.steps$num.explvars))
  if (is.bic == TRUE) {
    yrng <- c(min(df.steps$BIC - 250), max(df.steps$BIC))
  } else {
    yrng <- c(min(df.steps$mcfadden.r2), max(df.steps$mcfadden.r2))
  }
  
  # Picked up 6 hues, 4 shades of each hue, to cover 24 lines
  # Values are from ColorBrewer2 http://colorbrewer2.org/
  seqcols <- c(rgb(158, 202, 225, 255, maxColorValue=255),
               rgb(161, 217, 155, 255, maxColorValue=255),
               rgb(189, 189, 189, 255, maxColorValue=255),
               rgb(253, 174, 107, 255, maxColorValue=255),
               rgb(188, 189, 220, 255, maxColorValue=255),
               rgb(252, 146, 114, 255, maxColorValue=255),
               rgb(66, 146, 198, 255, maxColorValue=255),
               rgb(65, 171, 93, 255, maxColorValue=255),
               rgb(115, 115, 115, 255, maxColorValue=255),
               rgb(241, 105, 19, 255, maxColorValue=255),
               rgb(128, 125, 186, 255, maxColorValue=255),
               rgb(239, 59, 44, 255, maxColorValue=255),
               rgb(8, 81, 156, 255, maxColorValue=255),
               rgb(0, 109, 44, 255, maxColorValue=255),
               rgb(37, 37, 37, 255, maxColorValue=255),
               rgb(166, 54, 3, 255, maxColorValue=255),
               rgb(84, 39, 143, 255, maxColorValue=255),
               rgb(165, 15, 21, 255, maxColorValue=255),
               rgb(8, 48, 107, 255, maxColorValue=255),
               rgb(0, 68, 27, 255, maxColorValue=255),
               rgb(0, 0, 0, 255, maxColorValue=255),
               rgb(127, 39, 4, 255, maxColorValue=255),
               rgb(63, 0, 125, 255, maxColorValue=255),
               rgb(103, 0, 13, 255, maxColorValue=255),
               rgb(139, 69,19, 255, maxColorValue=255)) # 25th color is brown
  
  lags.unique <- unique(df.steps$num.cdhlags)
    
  # Set some plot attributes
  par(xpd = TRUE, # turn off clipping
      mar = c(6, 4.5, 3, 7.5)) # Add margins to the right (b, l, t, r)
  
  # Plot current temperature information
  df.steps.subset <- subset(df.steps, num.cdhlags == 0)
  if (is.bic == TRUE) {
    yvals <- df.steps.subset$BIC
    ylabel <- "Bayesian Informaion Criterion (BIC)"
  } else {
    yvals <- df.steps.subset$mcfadden.r2
    ylabel <- expression(paste("McFadden's Pseudo-", R^2))
  }
  plot(x = df.steps.subset$num.explvars,
       y = yvals, 
       type = "l", 
       xlab = "Number of Explanatory Variables", 
       xlim = xrng, 
       ylab = ylabel, 
       ylim = yrng, 
       main = title,
       sub = subtitle, 
       bty = "L", 
       col = seqcols[1], 
       lty = 1,
       lwd = 2,) # line width
  
  # Iterate through the other numbers of hours into the past
  maxcdhlags <- max(df.steps$num.cdhlags)
  for (i in lags.unique) {
    df.steps.subset <- subset(df.steps, num.cdhlags == i)
    if (is.bic == TRUE) {
      yvals <- df.steps.subset$BIC
      ypnt <- min(df.steps.subset$BIC)
      df.steps.yset <- subset(df.steps.subset, BIC == ypnt)
      xpnt <- df.steps.yset[1, "num.explvars"]
      critpnt.label <- c(round(ypnt, 1))
      textposition <- 1 # below
    } else {
      yvals <- df.steps.subset$mcfadden.r2
      df.steps.yset <- subset(df.steps.subset, BIC == min(df.steps.subset$BIC))
      ypnt <- df.steps.yset[1, "mcfadden.r2"]
      xpnt <- df.steps.yset[1, "num.explvars"]
      critpnt.label <- c(round(ypnt, 3))
      textposition <- 1 # below
    }
    
    lines(x = df.steps.subset$num.explvars,
         y = yvals,
         col = seqcols[(i+1)], 
         lty = (i + 1),
         lwd = 2) # line width
    
    # Label the critcal y value
    points(x = xpnt,
           y = ypnt,
           pch = 20, # dot
           col = seqcols[(i+1)])
    
    text(x = xpnt, 
         xlim = xrng,
         y = ypnt, 
         ylim = yrng, 
         labels = critpnt.label, 
         pos = textposition, 
         cex = 0.6, # 50% font size
         offset = 0.25,
         col = "black") 
  }
  
  legend(x = (xrng[2] + 2),
         y = yrng[2],
         title = "Past Hrs. Temp.",
         legend = paste(lags.unique, "hrs."), 
         lty = lags.unique + 1,
         lwd = 2,
         col = seqcols[c(lags.unique + 1)])
}