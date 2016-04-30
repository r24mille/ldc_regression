require(ggplot2)        # Better plotting!
require(scales)         # For scaling the y-axis on the Adj. R^2 results

# Standardize things like plot titles, labels, sizes, and colours
source("thesis_standardization.R")

# Offset of y-axis used for omitted benchmarks and the 
# y-axis translate function.
translate.val <- 0.4

method.names <- c("Non-Weather Only",
                  "Non-Weather Only", 
                  #"None / Dry-bulb",
                  "None /    Dry-bulb",
                  "None /    Dry-bulb",
                  "None / Feels Like",
                  "None / Feels Like",
                  "Lagged (i-1) /    Dry-bulb",
                  "Lagged (i-1) /    Dry-bulb",
                  "Lagged (i-1) / Feels Like",
                  "Lagged (i-1) / Feels Like",
                  "Lagged (i-2) /    Dry-bulb",
                  "Lagged (i-2) /    Dry-bulb",
                  "Lagged (i-2) / Feels Like",
                  "Lagged (i-2) / Feels Like",
                  "Lagged (i-3) /    Dry-bulb",
                  "Lagged (i-3) /    Dry-bulb",
                  "Lagged (i-3) / Feels Like",
                  "Lagged (i-3) / Feels Like",
                  "Lagged (i-4) /    Dry-bulb",
                  "Lagged (i-4) /    Dry-bulb",
                  "Lagged (i-4) / Feels Like",
                  "Lagged (i-4) / Feels Like",
                  #"Lagged (i-5) / Dry-bulb",
                  #"Lagged (i-5) / Dry-bulb",
                  #"Lagged (i-5) / Feels Like",
                  #"Lagged (i-5) / Feels Like",
                  #"Lagged (i-6) / Dry-bulb",
                  #"Lagged (i-6) / Dry-bulb",
                  #"Lagged (i-6) / Feels Like",
                  #"Lagged (i-6) / Feels Like",
                  "CDH & HDH /    Dry-bulb",
                  "CDH & HDH /    Dry-bulb",
                  "CDH & HDH / Feels Like",
                  "CDH & HDH / Feels Like",
                  "Moving Avg /    Dry-bulb",
                  "Moving Avg /    Dry-bulb",
                  "Moving Avg / Feels Like",
                  "Moving Avg / Feels Like",
                  "Exp-Lag-Res /    Dry-bulb",
                  "Exp-Lag-Res /    Dry-bulb",
                  "Exp-Lag-Res / Feels Like",
                  "Exp-Lag-Res / Feels Like")

nonlinear.trnfms <- c("None",
                      " ", # Placeholder for ggplot2
                      #"None",
                      "Switching Regression",
                      "Natural Splines",
                      "Switching Regression",
                      "Natural Splines",
                      "Switching Regression",
                      "Natural Splines",
                      "Switching Regression",
                      "Natural Splines",
                      "Switching Regression",
                      "Natural Splines",
                      "Switching Regression",
                      "Natural Splines",
                      "Switching Regression",
                      "Natural Splines",
                      "Switching Regression",
                      "Natural Splines",
                      "Switching Regression",
                      "Natural Splines",
                      "Switching Regression",
                      "Natural Splines",
                      #"Switching Regression",
                      #"Natural Splines",
                      #"Switching Regression",
                      #"Natural Splines",
                      #"Switching Regression",
                      #"Natural Splines",
                      #"Switching Regression",
                      #"Natural Splines",
                      "Switching Regression",
                      " ", # Placeholder for ggplot2
                      "Switching Regression",
                      " ", # Placeholder for ggplot2
                      "Switching Regression",
                      "Natural Splines",
                      "Switching Regression",
                      "Natural Splines",
                      "Switching Regression",
                      "Natural Splines",
                      "Switching Regression",
                      "Natural Splines")

adjR2s <- c(0.438,
            translate.val, # Placeholder for ggplot2
            #0.580,
            0.854,
            0.862,
            0.857,
            0.862,
            0.875,
            0.884,
            0.878,
            0.884,
            0.881,
            0.889,
            0.883,
            0.890,
            0.873,
            0.880,
            0.875,
            0.882,
            0.853,
            0.859,
            0.855,
            0.862,
            #0.825,
            #0.830,
            #0.828,
            #0.834,
            #0.790,
            #0.796,
            #0.794,
            #0.801,
            0.895,
            translate.val, # Placeholder for ggplot2
            0.896,
            translate.val, # Placeholder for ggplot2
            0.895,
            0.902,
            0.897,
            0.904,
            0.902,
            0.910,
            0.901,
            0.911)

# Stitch together a data.frame of benchmark results
adjR2.results <- data.frame(
  "method.name" = factor(x = method.names, 
                         levels = unique(method.names),
                         ordered = TRUE),
  "nonlinear.trnfm" = factor(x = nonlinear.trnfms,
                             levels = c("None", 
                                        "Switching Regression", 
                                        "Natural Splines",
                                        " "),
                             ordered = TRUE),
  "adj.r2" = adjR2s)

# Custom colour palette for result bars
results.palette <- c("None" = "#1b9e77", 
                     "Switching Regression" = "#d95f02", 
                     "Natural Splines" = "#7570b3", 
                     " " = "#ffffff")

AdjR2_trans <- function() {
  trans <- function(x) x - translate.val
  inv   <- function(x) x + translate.val
  trans_new("AdjR2_trans", trans, inv)
}    

# Create the plot of Adjusted R^2 results
adjR2.results.plot <- (ggplot(adjR2.results, 
                              aes(x = method.name, 
                                  y = adj.r2,
                                  fill = nonlinear.trnfm)) + 
                         geom_bar(stat = "identity",
                                  position="dodge") + 
                         scale_fill_manual(values = results.palette) + 
                         ThesisPlotTheme() + 
                         theme(axis.text.x = element_text(angle = 90, 
                                                          vjust = 0.5,
                                                          hjust = 1)) + 
                         scale_y_continuous(limits=c(0.4, 1),
                                            trans = "AdjR2") + 
                         labs(x = "Method (Temporal Transform / Weather Transform)", 
                              y = expression(paste("Adjusted ", R^{2})), 
                              title= "Temperature Transformation Comparison Results",
                              fill = "Non-Linear Transform")
                       
)
print(adjR2.results.plot)
dev.print(file="../../Figures/AdjustedR2DodgeBarChart.png", device=png, height = 700, width = 850)
