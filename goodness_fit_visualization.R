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

PlotGlmFitMeasures <- function (aiccs, bics, resdevs, xvals, xtitle, title) {
  # Plots several measures of model fit on a dual y-axis plot. It shows how 
  # information criteria (IC) and residual deviance change as parameters 
  # within the model change iteratively.
  #
  # Args:
  #   aiccs: A vector of AICc values to be plotted
  #   bics: A vector of BIC values to be plotted
  #   resdevs: A vector of residual deviates to be plotted
  #   xvals: The values to be plotted along the x-axis. In the context of this 
  #          function they generally represent iterative changes to a model 
  #          parameter.
  #   xtitle: The title of the x-axis
  #   title: The main title of the plot
  
  # Adjust padding for labelling on the y2 axis
  par(mar=c(5,4,4,5)+.1)
  xrng = c(min(xvals), max(xvals))
  y1rng = c(min(c(bics, aiccs)),
            max(c(bics, aiccs)))
  y2rng = c(min(resdevs), max(resdevs))
  
  # Colors, helped by http://colorbrewer2.org/
  colaicc = "#4daf4a" # greenish
  colbic = "#377eb8" # blueish
  colresdev = "#e41a1c" # reddish
  
  # y1, plot AICc values
  plot(x = xvals,
       y = aiccs, 
       type = "l", 
       col = colaicc, 
       xlab = xtitle, 
       xlim = xrng, 
       ylab = "Information Criteria (AICc and BIC)", 
       ylim = y1rng, 
       main = title)
  
  # y1, plot BIC values
  lines(x = xvals,
        y = bics,
        col = colbic)
  
  # 2nd plot over the existing plot
  par(new=TRUE)
  
  # y2, plot the residual deviances
  plot(x = xvals,
       y = resdevs, 
       type = "l", 
       col = colresdev, 
       xaxt = "n", 
       xlab = "", 
       xlim = xrng, 
       yaxt = "n", 
       ylab = "", 
       ylim = y2rng)
  axis(4)
  mtext("Residual Deviance",
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
  
  # Create a point and label Res. Deviance occurring at lowest BIC
  resdev.atbic <- round(resdevs[minbic.index], 
                       digits = 1)
  plot(x = xvals[minbic.index],
       xaxt = "n", 
       xlab = "", 
       xlim = xrng, 
       yaxt = "n", 
       ylab = "", 
       y = resdev.atbic, 
       ylim = y2rng, 
       type = "p", 
       pch = 20, # dot
       col = colresdev)
  text(x = xvals[minbic.index], 
       xlim = xrng,
       y = resdev.atbic, 
       ylim = y2rng, 
       labels = c(resdev.atbic), 
       pos = 3, # above
       col = "black") 
  
  legend("topright", 
         col = c(colaicc, colbic, colresdev), 
         lty = 1, 
         legend = c("Corrected Akaike's Information Criterion (AICc)",
                    "Bayesian Information Criterion (BIC)",
                    "Residual Deviance"))
}