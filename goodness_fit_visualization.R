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