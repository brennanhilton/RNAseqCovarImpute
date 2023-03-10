% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/voom2.R
\name{voom2}
\alias{voom2}
\title{voom2}
\usage{
voom2(
  counts,
  design = NULL,
  lib.size = NULL,
  normalize.method = "none",
  block = NULL,
  correlation = NULL,
  weights = NULL,
  span = 0.5,
  plot = FALSE,
  save.plot = FALSE,
  lib.size.all
)
}
\arguments{
\item{counts}{see voom from limma package}

\item{design}{see voom from limma package}

\item{lib.size}{see voom from limma package}

\item{normalize.method}{see voom from limma package}

\item{block}{see voom from limma package}

\item{correlation}{see voom from limma package}

\item{weights}{see voom from limma package}

\item{span}{see voom from limma package}

\item{plot}{see voom from limma package}

\item{save.plot}{see voom from limma package}

\item{lib.size.all}{argument used by limma_voom_imputed_data_list and limma_voom_imputed_data_list_parallel functions. See description.}
}
\value{
see voom from limma package
}
\description{
Modified voom function used by limma_voom_imputed_data_list and limma_voom_imputed_data_list_parallel functions.
allows input of bins of outcome genes while still accounting for the total library size of all outcome genes, as the total library size is needed to calculate log-cpm values.
}
