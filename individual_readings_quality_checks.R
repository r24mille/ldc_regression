library(ggplot2)

expl_var_dist_checks <- function(df){
  # Lists some summary statistics, visualizations, etc. to carry out a data 
  # quality check before modeling.
  #
  # Args:
  #   df: A data.frame of smart meter readings. The data.frame must contain a 
  #       column "meterid" that represents the individual/household.
  
  print("Full Data Frame Summary")
  print(summary(df))
  
  print("Density of response variable (kWh)")
  plot(density(df$kwh))
  
  print("Density of log(response)")
  plot(density(log(df$kwh)))
  
  print("Check that each explanatory variable is adequately represented")
  print(qplot(df$meterid))
  print(qplot(df$hrstr))
  print(qplot(df$month))
  print(qplot(df$weekend))
  print(qplot(df$dayname))
  print(qplot(df$holiday))
  print(qplot(df$price))
  print(qplot(df$weather_reduced))
  
  print("Summary of how terms were reduced, check that everything seems valid")
  print("clear:")
  print(unique(df[df$weather_reduced == "clear",]$weather_desc))
  print("cloudy_fog:")
  print(unique(df[df$weather_reduced == "cloudy_fog",]$weather_desc))
  print("rain_tstorms:")
  print(unique(df[df$weather_reduced == "rain_tstorms",]$weather_desc))
  print("snow_ice:")
  print(unique(df[df$weather_reduced == "snow_ice",]$weather_desc))
  
  print("Density plots of continuous explanatory variables")
  plot(density(df$temperature))
  plot(density(df$dewpnt_temp))
}
