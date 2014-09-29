##
# Creates a vector of values based on the "conventional" definition of CDH.
# It sums all temperature values above the CDH breakpoint from the current time 
# through a number of hours into the past (nlags).
#
# readingsdf is a dataframe of smart meter readings.
createCdhLagSum <- function(nlags, readingsdf) {
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

##
# Creates a matrix of values such that each column looks an additional hour 
# into the past. The results is an [n x (nlags + 1)] matrix such that each 
# hour into the past can have its own coefficient when modelled.
#
# readingsdf is a dataframe of smart meter readings.
createCdhLagMatrix <- function(nlags, readingsdf) {
  # Set up matrix
  lagnames <- paste0("lag", c(0:nlags))
  cdhlags <- matrix(nrow = nrow(readingsdf),
                   ncol = (nlags + 1))
  colnames(cdhlags) <- lagnames
  
  # Unfortunately iterating over the dataframe seems to be the best method
  for(i in 1:nrow(readingsdf)) {
    for(j in 0:nlags) {
      if (i-j > 0) {
        cdhlags[i, (j+1)] <- readingsdf[i-j, "cdh"]
      } else {
        cdhlags[i, (j+1)] <- 0
      }
    }
  }
  
  return(cdhlags)
}