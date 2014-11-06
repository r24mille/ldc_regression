library(rgl) # OpenGL functionality in R for better 3D plots
library(reshape) # Reshape data.frame for the heatmap

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

Plot3DStepwiseResults <- function(df, title) {
  # 3D wireframe plot of stepwise linear regression results.
  #
  # TODO(r24mille): Kind of dumping RGL and heatmap versions into this function
  #                 for now. Clean these up some time in the future.
  # 
  # Args:
  #   df: The stepwise regression results in the format of df.stepresults 
  #       convention established elsewhere.
  #   title: The plot title passed to main
  wireframe(BIC ~ num.explvars * num.cdhlags, 
            data = df,
            xlab = list("Number of Explanatory Variables", 
                        rot=-13), 
            ylab = list("Number of Past Hours of Temp. > CDH_break Included", 
                        rot=65),
            zlab = list("Bayesian Information Criterion (BIC)", 
                        rot=90), 
            main = title,
            drape = TRUE,
            colorkey = TRUE,
            scales = list(arrows = FALSE), # Switches unlabelled arrows to ticks
            screen = list(z = -20, x = -60))
  
  # Decent setting for BIC
  # screen = list(z = -20, x = -60)
  
  # That same information plotted with RGL
  df.heat <- df
  df.heat.reshape <- melt(data = df.heat,
                          id.vars = c("num.cdhlags", "num.explvars"),
                          measure.vars = c("BIC"))
  cast.heat <- cast(data = df.heat.reshape, 
                    formula = num.cdhlags ~ num.explvars)
  stepresults.matrix <- as.matrix(cast.heat)
  bic.lim <- range(df.stepresults.trimmed$BIC)
  bic.len <- bic.lim[2] - bic.lim[1] + 1
  hgt.lookup <- terrain.colors(bic.len) # height color lookup table
  bic.colors <- hgt.lookup[stepresults.matrix - bic.lim[1] + 1] # assign colors to heights
  persp3d(x = seq(0, (nrow(stepresults.matrix) - 1), len = nrow(stepresults.matrix)),
          y = seq(1, ncol(stepresults.matrix), len = ncol(stepresults.matrix)),
          z = stepresults.matrix, 
          xlab = "Number of Past Hours Included",
          ylab = "Number of Explanatory Variables",
          zlab = "Bayesian Information Criterion (BIC)",
          main = title,
          col = bic.colors)
  
  # Heatmap
  #maxbic <- max(df.stepresults.trimmed$BIC)
  #stepresults.matrix[which(is.na(stepresults.matrix))] = maxbic;
  heatmap(stepresults.matrix, 
          Rowv=NA, 
          Colv=NA, 
          col = heat.colors(round(bic.len) * 2), 
          margins=c(5,10),
          na.rm = TRUE)
}

PlotLasso <- function(lars.obj, y, design.mat, xvar = "shrinkage", 
                      term.labels = "path", line.marker = "xsteps", 
                      xvar.parsimonious) {
  # 
  # Args:
  #   lars.obj: The LASSO object returned by the lars(...) function.
  #   y : The response variable used in the LARS LASSO operation
  #   design.mat: The design matrix used in the LARS LASSO operation
  #   xvar: Controls the type of the x-axis and should be one of c("shrinkage", 
  #         "rsq", "degf"). "shrinkage" (default) coefficient as a function of 
  #         the shrinkage factor, s. "rsq" plots coefficient as a function of 
  #         the r-squared value of that step. "degf" plots coefficient as a 
  #         function of the degrees of freedom in the model.
  #   term.labels: Controls how terms should be labeled on the plot and should 
  #         be one of c("path", "catleg"). "path" (default) labels 
  #         the term on its LASSO path line. "catlegend" colors the lines by the 
  #         factor/term they come from and creates an a matching legend.
  #   line.marker: Controls what information will be highlighted by the vertical 
  #         dashed line, must be one of c("xsteps", "1se"). "xsteps" simply 
  #         marks ech step of x (default). "1se" highlights the x value that 
  #         corresponds to the parsimonious fit, within 1 SE of the optimal 
  #         model as chosen during cross-validation.
  #   xvar.parsimonious: If line.marker = "1se" then this value is required. It 
  #         is the index of the parsimonious x value that will be highlighted by 
  #         the dotted abline.
     
  # turn off clipping and add margins (b, l, t, r)
  if (term.labels == "catleg") {
    par(xpd = FALSE, mar = c(5, 4.5, 3, 1))
  } else {
    par(xpd = TRUE, mar = c(5, 4.5, 3, 3))
  }
  
  # Set up line colors
  col.lines <- rep("black", ncol(design.mat))
  leg.labels <- list()
  if (term.labels == "catleg") {
    col.month <- "red"
    col.hrstr <- "blue"
    col.weekend <- "green"
    col.templag <- "orange"
    
    colnames.design.mat <- colnames(design.mat)
    indeces.term.month <- grep("^monthm[[:digit:]]*$", colnames.design.mat, 
                               perl = TRUE, value = FALSE)
    indeces.term.hrstr <- grep("^hrstrh[[:digit:]]*$", colnames.design.mat, 
                               perl = TRUE, value = FALSE)
    indeces.term.weekend <- grep("^weekend[[:alpha:]]*$", colnames.design.mat, 
                                 perl = TRUE, value = FALSE)
    indeces.term.templag <- grep("^templag[[:digit:]]*$", colnames.design.mat, 
                                 perl = TRUE, value = FALSE)
    
    col.lines[indeces.term.month] <- col.month
    if(length(indeces.term.month) > 0) {
      leg.labels["Month"] <- col.month
    }
    col.lines[indeces.term.hrstr] <- col.hrstr
    if(length(indeces.term.hrstr) > 0) {
      leg.labels["Hour of Day"] <- col.hrstr
    }
    col.lines[indeces.term.weekend] <- col.weekend
    if(length(indeces.term.weekend) > 0) {
      leg.labels["Weekday or Weekend/Hol."] <- col.weekend
    }
    col.lines[indeces.term.templag] <- col.templag
    if(length(indeces.term.templag) > 0) {
      leg.labels["Temperature History"] <- col.templag
    }
  }
  
  # Create convenient variables for coefficients to be used on y-axis
  coef.max = max(coef(lars.obj))
  coef.min = min(coef(lars.obj))
  coef.range = abs(coef.max - coef.min)
  
  nterms <- ncol(design.mat)
  xvar.values = rep(0, nterms)
  xvar.label = ""
  xvar.lim <- c(0, 1)
  xvar.at <- seq(0, 1, .1)
  
  if (xvar == "shrinkage") {
    # Examine the ordinary least squares fit of the full model using multiple regression(?)
    OLS <- summary(lm(y ~ ., data = as.data.frame(design.mat)))
    absum <- sum(abs(OLS$coeff[-1,1]))
    
    # Manually create the LASSO step plot from the ground up
    t <- apply(X = abs(coef(lars.obj)), 
               MARGIN = 1, # 1 indicates rows
               FUN = sum)
    s <- t/absum
    xvar.values <- s
    xvar.label <- "Shrinkage Factor"
    xvar.lim <- c(min(s), max(s))
    xvar.at <- seq(min(s), max(s), .2)
  } else if (xvar == "rsq") {
    rsq.min <- min(lars.obj$R2)
    rsq.max <- max(lars.obj$R2)
    rsq.max.rndup = ceiling(rsq.max * 10) / 10;
    
    xvar.values <- lars.obj$R2
    xvar.label <- "R-Squared"
    xvar.lim <- c(0, rsq.max.rndup)
    xvar.at <- seq(0, rsq.max.rndup, round((rsq.max/11), 2))
  } else if (xvar == "degf") {
    degf.min <- min(lars.obj$df)
    degf.max <- max(lars.obj$df)
    degf.rng <- degf.max - degf.min
    
    xvar.values <- lars.obj$df
    xvar.label <- "Degrees of Freedom in Model"
    xvar.lim <- c(degf.min, degf.max)
    xvar.at <- seq(degf.min, degf.max, floor(degf.rng/10))
  }

  plot(x = xvar.values,
       y = coef(lars.obj)[, 1], 
       ylim = c(coef.min, coef.max),
       type = "l",
       lwd = 2,
       xlab = xvar.label,
       main = paste("LASSO Path: Coefficients as", xvar.label, "Change"),
       xlim = xvar.lim,
       cex.lab = 1.1,
       axes = FALSE,
       ylab = "Coefficient",
       col = rgb(0, 0, 0, 0, maxColorValue = 1)) # Transparent line
  
  axis(1, at = xvar.at, cex.axis = 1) # Control x-axis
  
  axis(2, 
       at = round(seq(coef.min, coef.max, (coef.range/10)), 2),
       cex.axis = 1, 
       las = 1) # Control y-axis
  
  if (line.marker == "1se") {
    abline(v = xvar.parsimonious,
           col = "darkgray",
           lty = 3)
  } else {
    abline(v = xvar.values,
           col = "darkgray",
           lty = 3)
  }
  
  for(i in 1:nterms) {
    # Draw the LASSO line
    lines(xvar.values, 
          coef(lars.obj)[,i], 
          col = col.lines[i],
          lwd = 1.3)
    
    if (term.labels == "path") {
      # Label the line at its endpoint
      text(xvar.lim[2] + (xvar.at[2]/3), 
           coef(lars.obj)[nrow(coef(lars.obj)), i],
           colnames(design.mat)[i])
    }
  }
  
  if (term.labels == "catleg") {
    legend("bottomleft",
           title = "Terms",
           legend = names(leg.labels), 
           lty = 1,
           lwd = 1.3,
           col = unlist(leg.labels, use.names = FALSE))
  }
  

}

PlotLassoCrossValidation <- function(design.mat, y.vec, k = 10, 
                                     backtransform.mse = "none",
                                     xvar = "shrinkage") {
  # Args:
  #   design.mat: Design matrix suitable to be passed to LARS for LASSO
  #   y.vec: Vector of of response variables associated with the design matrix
  #   k: The number of folds in k-folds cross-validation
  #   backtransform.mse: The type of backtransformation to apply to MSE (if any)
  #                      must be a value of c("none", "log"). The default is 
  #                      "none".
  #   xvar: The variable to plot on the x-axis, must be a value of 
  #         c("shrinkage", "step"). The default is "shrinkage".
  #
  # Return:
  #   The parsimonious value chosen from the x-axis, within 1 standard error 
  #   of optimal value chosen during cross-validation.
  
  # Generate a vector of holdout labels, vector same length as number of rows 
  # in the dataset. Values of holdout labels will be on the range 1-k
  cvlab <- sample(1:k, length(y.vec), replace = TRUE)
  
  # Create a vector of candidate s values (reasonable number for now)
  if (xvar == "step") {
    svec <- seq(1, ncol(design.mat), 1)
  } else {
    svec <- seq(0, 1, .05)
  }
  J <- length(svec) # How many versions of s did I decide to use?
  
  # Going to perform LASSO k times, create a list of k LARS result objects
  lassolist1 <- list()
  
  # Initialize a list to store prediction from each LASSO set
  predtrain1 <- list()
  
  # Compute MSE from predictions
  MSEstore1 <- matrix(NA, J, k) # J values of s (rows), k hold-out sets (columns)
  
  # Use a for loop to get each lasso fit holding out the ith set
  # Then predict the ith set using the holdout model
  for (i in 1:k) {
    lassolist1[[i]] <- lars(x = design.mat[cvlab!=i,], 
                            y = y.vec[cvlab!=i],
                            type = 'lasso',
                            trace = FALSE,
                            normalize = TRUE,
                            intercept = TRUE)
    if (xvar == "step") {
      pred.mode = "step"
    } else {
      pred.mode = "fraction"
    }
    predtrain1[[i]] <- predict.lars(object = lassolist1[[i]],
                                    newx = design.mat[cvlab == i,], 
                                    s = svec,
                                    mode = pred.mode,
                                    type = "fit")$fit
    
    
    # Start a new loop to get MSE for each combination of the ith holdout set and 
    # jth value of s.
    for(j in 1:J) {
      # This computes MSE
      if (backtransform.mse == "log") {
        predtrain.backtranformed <- exp(predtrain1[[i]][,j])
        train.backtransformed <- exp(y.vec[cvlab==i])
        MSEstore1[j,i] <- mean((predtrain.backtranformed - train.backtransformed)^2)
      } else {
        predtrain.backtranformed <- predtrain1[[i]][,j]
        train.backtransformed <- y.vec[cvlab==i]
        MSEstore1[j,i] <- mean((predtrain.backtranformed - train.backtransformed)^2)
      }
    }
  }
  
  # Compute mean and standard error of the observed MSEs at J values
  meanMSE <- apply(MSEstore1, 1, mean)
  stdMSE <- apply(MSEstore1, 1, sd)/sqrt(k)
  
  # y label can change based on whether its backtransformed
  y.label <- "Mean Square Error (MSE)"
  if (backtransform.mse == "log") {
    y.label <- "Backtransformed Mean Square Error (MSE)"
  }
  
  # x label can change based on how predictions were iterated
  x.label <- "Shrinkage Factor"
  if (xvar == "step") {
    x.label <- "LASSO Step"
  }
  
  # Plot the change in MSE as shrinkage factor increases. Error bars established 
  # by kfolds cross-validation.
  mse.lwr <- min(meanMSE) - max(stdMSE)
  mse.upr <- max(meanMSE) + max(stdMSE)
  plot(x = svec, 
       y = meanMSE, 
       ylim = c(mse.lwr, mse.upr), 
       pch = 16, 
       col = colors()[258], 
       axes = FALSE, 
       xlab = x.label,
       ylab = y.label, 
       cex.lab = 1.1,
       main = paste("Avg. Cross-Validation Prediction Error as", x.label, "Changes"))
  
  # Adjust x-axis based on how LASSO is being iterated
  if (xvar == "step") {
    svec.rng <- (max(svec) - min(svec))
    xvar.at <- seq(min(svec), max(svec), ceiling(svec.rng / 10))
  } else {
    svec.rng <- (max(svec) - min(svec))
    xvar.at <- seq(min(svec), max(svec), round((svec.rng / 10), 2))
  }

  axis(1, at = xvar.at, cex.axis = 1)
  axis(2, 
       las = 1, 
       at = seq(round(mse.lwr, 2), 
                round(mse.upr, 2), 
                round(((mse.upr - mse.lwr)/15), 2)), 
       cex.axis = 1)
  lines(x = svec, 
        y = meanMSE, 
        lty = 1, 
        col = colors()[258])
  
  # Standard error interval
  for (i in 1:J) {
    arrows(x0 = svec[i], 
           y0 = (meanMSE[i] - stdMSE[i]), 
           x1 = svec[i], 
           y1 = (meanMSE[i] + stdMSE[i]), 
           angle = 90, 
           code = 3, # Arrow head at both ends of the bar
           length=0.05, # Shorten up the bar
           col = "black", 
           lty = 1)
  }
  
  mse.min.idx <- which.min(meanMSE)
  mse.min.std <- stdMSE[mse.min.idx]
  mse.1se.idx <- which(meanMSE < (min(meanMSE) + mse.min.std))[1]
  mse.1se <- meanMSE[mse.1se.idx]
  # 1SE threshold line
  abline(h = (min(meanMSE) + mse.min.std), 
         lty = 2)
  # Selected shrinkage value
  points(x = svec[mse.1se.idx], 
         y = meanMSE[mse.1se.idx], 
         col = "red", 
         pch = 15, 
         cex = 1.3) 
  
  y.leg <- "Average MSE"
  # Highlighted point can change based on how predictions were iterated
  highlight.leg <- "Chosen shrinkage value (s)"
  if (xvar == "step") {
    highlight.leg <- "Chosen number of LASSO steps"
  }
  
  legend("topright", 
         legend = c(y.leg, 
                    "Standard Error (SE) range", 
                    "1 SE above lowest avg. MSE",
                    highlight.leg),
         pch = c(16, NA, NA, 15),
         col = c(colors()[258], 1, 1, "red"),
         cex = 1,
         lty = c(1, 1, 2, NA))
  
  return(list(xvar.1se = svec[mse.1se.idx]))
}