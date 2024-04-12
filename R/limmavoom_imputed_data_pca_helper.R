#' limmavoom_imputed_data_pca_helper
#'
#' Runs limma-voom pipeline on one of the m imputed datasets and returns the coefficients, 
#' Bayesian moderated standard errors and degrees of freedom.
#' This function used while looping through the m imputed datasets.
#' @return A list with coefficients, standard errors, and total degrees of freedom values from limma-voom gene expression analysis.
#' @param i Number from 1-m
#' @param imp Imputed data from mice (mids object)
#' @param DGE A DGEList object.
#' @param voom_formula Formula for design matrix.
#'
#' @importFrom stats model.matrix
#' @importFrom limma eBayes voom lmFit
#' @importFrom mice complete
#' 
#' @keywords internal
limmavoom_imputed_data_pca_helper <- function(i, imp, DGE, voom_formula) {
  # Get imputed data
  imputed_data_i <- mice::complete(imp, i)
  # Design matrix
  design <- model.matrix(as.formula(voom_formula), imputed_data_i)
  # Linear model fitting and empirical Bayes moderation
  eout <- limma::eBayes(limma::lmFit(limma::voom(DGE, design = design)))
  # Extract coefficients and standard errors
  coef <- as.data.frame(eout$coefficients)
  SE <- as.data.frame(sqrt(eout$s2.post) * eout$stdev.unscaled)
  df <- eout$df.total[[1]]
  # Return the result as a list
  return(list(coef = coef, SE = SE, df = df))
}