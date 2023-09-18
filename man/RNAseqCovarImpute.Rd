\name{RNAseqCovarImpute}
\alias{RNAseqCovarImpute}
\title{Impute Covariate Data in RNA-sequencing Studies}
\description{
  The RNAseqCovarImpute package implements multiple imputation of missing covariates and differential gene expression analysis by: 1) Randomly binning genes into smaller groups, 2) Creating M imputed datasets separately within each bin, where the imputation predictor matrix includes all covariates and the log counts per million (CPM) for the genes within each bin, 3) Estimating gene expression changes using voom followed by lmFit functions, separately on each M imputed dataset within each gene bin, 4) Un-binning the gene sets and stacking the M sets of model results before applying the squeezeVar function to apply a variance shrinking Bayesian procedure to each M set of model results, 5) Pooling the results with Rubins’ rules to produce combined coefficients, standard errors, and P-values, and 6) Adjusting P-values for multiplicity to account for false discovery rate (FDR).
}
\author{
  Brennan H Baker
}