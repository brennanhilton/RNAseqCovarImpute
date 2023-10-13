#' combine_rubins
#'
#' Combines results from each imputed dataset using Rubin's rules.
#' @return Dataframe with one row per gene containing coefficients standard errors, degrees of freedom, t-statistics, P-Values, and adjusted P-values from the limma-voom pipeline.
#' \item{coef_combined}{combined logFCs across the multiple imputed datasets using Rubin's rules}
#' \item{SE_P}{pooled standard error across the multiple imputed datasets using Rubin's rules}
#' \item{SE_P_bayes}{pooled standard error across the multiple imputed datasets using Rubin's rules squeezed to global mean variance trend curve with limma-voom Bayesian procedure}
#' \item{df}{limma-voom residual degrees of freedom adjusted for Rubin's rules}
#' \item{df_bayes}{limma-voom residual degrees of freedom adjusted for Rubin's rules and Bayesian procedure}
#' \item{rubins_t}{t-statistic = coef_combined divided by SE_p}
#' \item{rubins_t_bayes}{t-statistic = coef_combined divided by SE_p_bayes}
#' \item{combined_p}{p-value from two-sided t-distribution alpha = 0.05 using rubins_t}
#' \item{combined_p_bayes}{p-value from two-sided t-distribution alpha = 0.05 using rubins_t_bayes}
#' \item{combined_p_adj}{false discovery rate (FDR) adjusted combined_p}
#' \item{combined_p_adj_bayes}{false discovery rate (FDR) adjusted combined_p_bayes}
#' @param DGE A DGEList object.
#' @param model_results Output from limmavoom_imputed_datalist.
#' @param predictor Independent variable of interest, in the form of a linear model contrast. Must be a variable in voom_formula.
#' @param covariate Arguments passed to limma::squeezeVar. If non-NULL, var.prior will depend on this numeric covariate. Otherwise, var.prior is constant.
#' @param robust Arguments passed to limma::squeezeVar. Llogical, should the estimation of df.prior and var.prior be robustified against outlier sample variances?
#' @param winsor.tail.p Arguments passed to limma::squeezeVar. Numeric vector of length 1 or 2, giving left and right tail proportions of x to Winsorize. Used only when robust=TRUE.
#'
#' @importFrom magrittr %>%
#' @importFrom dplyr select starts_with mutate tibble arrange contains
#' @importFrom stats p.adjust pt
#' @importFrom rlang := .data
#' @importFrom limma squeezeVar
#'
#' @examples
#' data(example_data)
#' data(example_DGE)
#' intervals <- get_gene_bin_intervals(example_DGE, example_data, n = 10)
#' gene_bin_impute <- impute_by_gene_bin(example_data,
#'     intervals,
#'     example_DGE,
#'     m = 2
#' )
#' coef_se <- limmavoom_imputed_data_list(
#'     gene_intervals = intervals,
#'     DGE = example_DGE,
#'     imputed_data_list = gene_bin_impute,
#'     m = 2,
#'     voom_formula = "~x + y + z + a + b"
#' )
#'
#' final_res <- combine_rubins(
#'     DGE = example_DGE,
#'     model_results = coef_se,
#'     predictor = "x"
#' )
#' @export


combine_rubins <- function(DGE, model_results, predictor, covariate = NULL, robust = FALSE, winsor.tail.p = c(0.05, 0.1)) {
    # Validity tests
    if (!class(DGE) %in% "DGEList") {
        stop("Input 'DGE' is not a valid DGEList object.")
    }
    if (!any((class(model_results) %in% c("tbl_df", "tbl", "data.frame")))) {
        stop("Input 'predictor' must be a character")
    }
    
    model_results <- model_results %>% dplyr::select(probe, contains(paste0(".", predictor, ".")))
    colnames(model_results) <- sub(paste0("\\.",predictor, "\\."), ".", colnames(model_results))
    # All residual dfs are the same (n-k-1)
    df_residual <- model_results$df_residual.1

    all_coefs <- model_results %>% dplyr::select(starts_with("coef"))
    all_SE <- model_results %>% dplyr::select(starts_with("SE_unscaled"))
    # Number of imputed datasets m
    m <- as.numeric(ncol(all_coefs))

    # Variance within is sum of square standard error divided by m number imputed datasets (m)
    # Calc variance within using unscaled SEs
    Vw <- apply(all_SE, 1, FUN = function(x) sum(x^2)) / m

    # pooled coefs, just the mean
    coefs_combined <- all_coefs %>% mutate(coef_pooled = rowMeans(all_coefs))
    # using pooled coefs we get variance between. subtract pooled coef from each imputed datasets coef
    coef_diff <- sweep(as.matrix(all_coefs), 1, as.matrix(coefs_combined$coef_pooled), "-")
    # get sum of squares of the pooled coef minus each imputed coef. Divide by m-1
    Vb <- apply(coef_diff, 1, FUN = function(x) sum(x^2)) / (m - 1)
    # Get total variance using using unscaled within variances
    Vtotal <- Vw + Vb + (Vb / m)

    # This is the Rubin's rules pooled SE
    SE_p <- sqrt(Vtotal)

    # Now get Rubin's rules degrees of freedom
    lambda <- (Vb + (Vb / m)) / Vtotal
    # RIV <- (Vb + (Vb/m))/Vw
    n <- as.numeric(ncol(DGE))
    DFold <- (m - 1) / (lambda^2)
    DFobs <- (((df_residual) + 1) / ((df_residual) + 3)) * (df_residual) * (1 - lambda)
    # The equation below is the residual degrees of freedom from voom and lmfit adjusted for rubins rules.
    df <- (DFold * DFobs) / (DFold + DFobs)

    #######################
    #######################
    ## eBayes part ########
    #######################
    #######################

    # Get all sigma values across imputed datasets
    all_sigma <- model_results %>% dplyr::select(starts_with("sigma"))
    # We can get std dev by dividing SE by sigma.
    all_std_dev <- all_SE / all_sigma


    # Now we get bayes prior df and bayes SE within each set of imputed datasets.
    # This way we get one set of priors rather than different prior per block of genes.
    priors <- foreach(i = seq(ncol(all_sigma)), .combine = "cbind") %do% {
        # Get sigma and std dev for ith imputed dataset
        sigma_i <- all_sigma[, i]
        std_dev_i <- all_std_dev[, i]

        # Get the priors!
        out <- squeezeVar(sigma_i^2, df_residual, covariate = covariate, robust = robust, winsor.tail.p = winsor.tail.p)
        s2.post <- out$var.post
        df_prior <- out$df.prior
        # Use s2.post with std dev to get bayes standard errors
        SE_bayes <- sqrt(s2.post) * std_dev_i
        SE_name <- paste0("SE_bayes.", i)
        df_name <- paste0("df_prior.", i)
        tibble(df_prior = df_prior, SE_bayes = SE_bayes) %>%
            dplyr::rename(
                !!SE_name := SE_bayes,
                !!df_name := df_prior
            )
    }
    # there is one prior df per set of imputed dataset. I am averaging them here. Each df_prior is
    # based on all genes being tested, but they come from different sets of imputed data.
    df_prior <- priors %>%
        dplyr::select(starts_with("df")) %>%
        rowMeans()
    df_prior <- df_prior[1]

    # Now apply rubins rules to the bayes standard errors
    all_SE_bayes <- priors %>% dplyr::select(starts_with("SE_bayes"))
    # Calc variance within using bayses SEs
    Vw_bayes <- apply(all_SE_bayes, 1, FUN = function(x) sum(x^2)) / m
    # Get total variance using using Bayes shrunk within variances
    Vtotal_bayes <- Vw_bayes + Vb + (Vb / m)
    # This is the Rubin's rules pooled bayses SE
    SE_p_bayes <- sqrt(Vtotal_bayes)

    # Bayes df. Replaces residual degrees of freedom with df.total, which is df.residual + df.prior.

    lambda_bayes <- (Vb + (Vb / m)) / Vtotal_bayes
    # RIV = (Vb + (Vb/m))/Vw
    n <- as.numeric(ncol(DGE))
    DFold <- (m - 1) / (lambda_bayes^2)
    DFobs <- (((df_residual + df_prior) + 1) / ((df_residual + df_prior) + 3)) * (df_residual + df_prior) * (1 - lambda_bayes)
    # The equation below is the residual degrees of freedom from voom and lmfit adjusted for rubins rules.
    df_bayes <- (DFold * DFobs) / (DFold + DFobs)

    # Get the final results!
    rubins_t <- tibble(probe = model_results$probe, coef_combined = coefs_combined$coef_pooled, SE_p = SE_p, SE_p_bayes = SE_p_bayes, df = df, df_bayes = df_bayes) %>%
        mutate(
            rubins_t = coef_combined / SE_p,
            rubins_t_bayes = coef_combined / SE_p_bayes
        ) %>% # new t statistics
        mutate(
            combined_p = 2 * pt(-abs(rubins_t), df = df),
            combined_p_bayes = 2 * pt(-abs(rubins_t_bayes), df = df_bayes)
        ) %>% # get p values with t distribution two sided
        arrange(combined_p) %>% #
        mutate(
            combined_p_adj = p.adjust(combined_p, method = "fdr"),
            combined_p_adj_bayes = p.adjust(combined_p_bayes, method = "fdr")
        )
}
