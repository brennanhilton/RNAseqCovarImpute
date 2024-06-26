% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/limmavoom_imputed_data_pca.R
\name{limmavoom_imputed_data_pca}
\alias{limmavoom_imputed_data_pca}
\title{limmavoom_imputed_data_pca}
\usage{
limmavoom_imputed_data_pca(imp, DGE, voom_formula, BPPARAM = bpparam())
}
\arguments{
\item{imp}{Imputed data from mice (mids object)}

\item{DGE}{A DGEList object.}

\item{voom_formula}{Formula for design matrix.}

\item{BPPARAM}{A BiocParallelParam object}
}
\value{
Dataframe with one row per gene. Columns contain coefficients, standard errors, and p-values from the limma-voom pipeline.
}
\description{
Combines results from each imputed dataset using Rubin's rules.
}
\examples{
data(example_data)
data(example_DGE)
pca_data = limma::voom(example_DGE)$E
p = PCAtools::pca(pca_data)
pcs = p$rotated[,1:78]
example_data = cbind(example_data, pcs)
imp = mice::mice(example_data, m=3)
mi_pca_res = limmavoom_imputed_data_pca(
   imp = imp,
   DGE = example_DGE,
   voom_formula = "~x + y + z + a + b"
   )
}
