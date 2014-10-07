McFaddenPseudoR2 <- function(model) {
  # Many online forums describe McFadden's pseudo r-squared as:
  #
  # Pseudo_R^2 = 1 - (residual_deviance / null_deviance)
  #
  # For example:
  #   * http://www.win-vector.com/blog/2011/09/the-simpler-derivation-of-logistic-regression/
  #   * https://stats.stackexchange.com/questions/11676/pseudo-r-squared-formula-for-glms
  # 
  # In the 1974 publication "Conditional Logit Analysis of Qualitative Choice
  # Behavior" (http://www.econ.berkeley.edu/reprints/mcfadden/zarembka.pdf) on
  # page 122 (PDF page 19), equation 37, he states the formula to be:
  #
  # Pseudo_R^2 = 1 - (G / G^H)
  #
  # Where I believe G is residual deviance and G^H is deviance when there are 
  # no explanatory variables (ie. R's null deviance).
  # 
  # TODO(r24mille): Read the paper more carefully to be sure this section is as 
  #                 the internet describes it.
  #
  # Args:
  #   model: A model object created by the glm(...) function.
  #
  # Returns:
  #   The McFadden pseudo r-squared value.
  
  mcfadden.pseudoR2 <- 1 - (model$deviance / model$null.deviance)
  return(mcfadden.pseudoR2)
}