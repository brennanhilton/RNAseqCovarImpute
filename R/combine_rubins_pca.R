#' combine_rubins_pca
#'
#' Combines results from each imputed dataset using Rubin's rules.
#' @return Coefficient, standard error, and p-value combined across m imputed datasets.
#' @param gene The number gene in the DGE list.
#' @param results_final List with one element per gene storing coefficients and standard errors from limma-voom analysis.
#' @param avg_df The average df.total from limma-voom analysis across the m imputed datasets.
#' @param m Number of imputed data sets.
#'
#' @importFrom stats pt
#' 
#' @keywords internal
combine_rubins_pca <- function(gene, results_final, avg_df, m) {
  # Extract coefficients and standard errors one gene at a time
  coef <- t(results_final[[gene]][seq(1, nrow(results_final[[gene]]), by = 2), ])
  se <- t(results_final[[gene]][seq(2, nrow(results_final[[gene]]), by = 2), ])
  # Variance within is sum of square standard error divided by m number imputed datasets (m)
  Vw <- apply(se, 1, FUN = function(x) sum(x^2)) / m
  # pooled coefs, just the mean
  coefs_combined <- rowMeans(coef)
  # using pooled coefs we get variance between. subtract pooled coef from each imputed datasets coef
  coef_diff <- sweep(coef, 1, coefs_combined, "-")
  # get sum of squares of the pooled coef minus each imputed coef. Divide by m-1
  Vb <- apply(coef_diff, 1, FUN = function(x) sum(x^2)) / (m - 1)
  # Get total variance using using unscaled within variances
  Vtotal <- Vw + Vb + (Vb / m)
  # This is the Rubin's rules pooled SE, the sqrt of total variance
  SE_p <- sqrt(Vtotal)
  # Now get Rubin's rules degrees of freedom adjusted (Barnard and Rubin 1999)
  lambda <- (Vb + (Vb / m)) / Vtotal
  # RIV <- (Vb + (Vb/m))/Vw
  DFold <- (m - 1) / (lambda^2)
  #avg_df is n-k-1 then bayes moderated then averaged across imputations
  DFobs <- (((avg_df) + 1) / ((avg_df) + 3)) * (avg_df+1) * (1 - lambda)
  # The degrees of freedom from voom and lmfit adjusted for rubins rules.
  df <- (DFold * DFobs) / (DFold + DFobs)
  if(is.nan(df[[1]])){ #if information lost (lambda) is 0,df will be nan
    df <-avg_df # in this case we use the original df unadjusted
  }
  # Calculate Rubin's t statistic
  rubins_t_bayes <- coefs_combined / SE_p
  # Calculate combined p-value and return results
  combined_p_bayes <- as.data.frame(t(2 * pt(-abs(rubins_t_bayes), df = df)))
  colnames(combined_p_bayes) <- paste0(colnames(combined_p_bayes), "_p")
  coefs_combined <- as.data.frame(t(coefs_combined))
  colnames(coefs_combined) <- paste0(colnames(coefs_combined), "_coef")
  SE_p <- as.data.frame(t(SE_p))
  colnames(SE_p) <- paste0(colnames(SE_p), "_se")
  return(cbind(combined_p_bayes, coefs_combined, SE_p))
}