plotRegSubSets <- function (bics, adjr2s, title) {
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