##
# An internet discussion thread was used as the starting point for this code.
# https://stackoverflow.com/questions/15874214/piecewise-function-fitting-with-nls-in-r
# 
# Also, as mentioned in the discussion thread this code owes a reference to 
# "Practical Regression and Anova using R" by Julian J. Faraway.
findSingleSetpoint <- function (x, y) {
  f <- function (Cx) {
    lhs <- function(x) ifelse(x < Cx,Cx-x,0)
    rhs <- function(x) ifelse(x < Cx,0,x-Cx)
    fit <- lm(y ~ lhs(x) + rhs(x))
    c(summary(fit)$r.squared, 
      summary(fit)$coef[1], 
      summary(fit)$coef[2],
      summary(fit)$coef[3])
  }
  
  r2 <- function(x) -(f(x)[1])
  
  res <- optimize(r2,interval=c(min(x),max(x)))
  res <- c(res$minimum, f(res$minimum))
  
  Cx.best <- res[1]
  coef1 <- res[3]
  coef2 <- res[4]
  coef3 <- res[5]  
  
  return(c(coef1, coef2, coef3, Cx.best))
}