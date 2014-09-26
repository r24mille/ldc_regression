##
# Creates a vector of CDH values for the number of hours lagged, specified by 
# nlags.
createCdhLag <- function(nlags, readingsdf) {
  cdhlag <- rep(0, nrow(readingsdf)) # Prealloate cdhlag vector
  
  # Unfortunately iterating over the dataframe seems to be the best method
  for(i in 1:nrow(readingsdf)) {
    cdhlag[i] <- readingsdf[i, "cdh"]
    tmplg <- 1
    
    # Sum what we can or all the values up to nlags
    while (i-tmplg > 0 & tmplg <= nlags) {
      cdhlag[i] <- cdhlag[i] + readingsdf[i-tmplg, "cdh"]
      tmplg <- tmplg + 1
    }
  }
  
  return(cdhlag)
}