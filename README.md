Impute Covariate Data in RNA-sequencing Studies
================

# Introduction

Missing data is a common problem in observational studies, as modeling
techniques such as linear regression cannot be fit to data with missing
points. Missing data is frequently handled using complete case analyses
in which any individuals with missing data are dropped from the study.
Dropping participants can reduce statistical power and, in some cases,
result in biased model estimates. A common technique to address these
problems is to replace or ‘impute’ missing data points with substituted
values. Typically, for a given covariate, missing data points are
imputed using a prediction model including other relevant covariates as
independent variables. In single imputation, a missing value is replaced
with the most likely value based on the predictive model. However, by
ignoring the uncertainty inherent with predicting missing data, single
imputation methods can result in biased coefficients and over-confident
standard errors. Multiple imputation addresses this problem by
generating several predictions, thereby allowing for uncertainty about
the missing data. In a typical multiple imputation procedure: 1) M
imputed data sets are created, 2) each data set is analyzed separately
(e.g., using linear regression), and 3) estimates and standard errors
across the M analyses are pooled using Rubin’s rules.

The RNAseqCovarImpute package implements multiple imputation of missing
covariates and differential gene expression analysis by 1) randomly
binning genes into smaller groups, 2) creating M imputed datasets
separately within each bin, where the imputation predictor matrix
includes all covariates and the log counts per million (CPM) for the
genes within each bin, 3) estimating gene expression changes using
`voom` followed by `lmFit` functions, separately on each M imputed
dataset within each gene bin, 4) un-binning the gene sets and stacking
the M sets of model results before applying the `squeezeVar` function to
apply a variance shrinking Bayesian procedure to each M set of model
results, 5) pooling the results with Rubins’ rules to produce combined
coefficients, standard errors, and P-values, and 6) adjusting P-values
for multiplicity to account for false discovery rate (FDR).

# Installation

``` r
# Install from github
library(devtools)
install_github("brennanhilton/RNAseqCovarImpute")

# Install from Bioconductor (not yet on Bioconductor)

if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("RNAseqCovarImpute")
```

# Generate random data with missing covariate data

Normally you would have your own covariate and RNA-sequencing data. We
generated random data for the purpose of this demonstration. The exact
code used to generate these data are found in the [Example Data for
RNAseqCovarImpute](Example_Data_for_RNAseqCovarImpute.html) vignette. In
short, `example_data` contains 500 rows with data for variables x, y,
and z, which are continuous normally distributed, and a and b, which are
binary variables. Missigness was simulated for all variables other than
x such that a complete case analysis would drop 24.2% of participants.
`example_DGE` contains random count data from the Poisson distribution
for 500 made up genes, ENS1-ENS500

``` r
library(RNAseqCovarImpute)
library(dplyr)
data(RNAseqCovarImpute_data)
```

# RNAseqCovarImpute Demonstration

## 1: Bin the genes into smaller groups

The default is approximately 1 gene per 10 individuals in the study, but
the user can specify a different ratio. For example, in a study with 500
participants and 10,000 genes, 200 bins of 50 genes would be created
using the default ratio. When the total number of genes is not divisible
by the bin size, the method flexibly creates bins of two different
sizes. For example, if the same hypothetical study included 10,001
genes, 199 bins of 50 and 1 bin of 51 genes would be created. The order
of the features (e.g., ENSEMBL gene identifiers) should be randomized
before binning.

``` r
intervals <- get_gene_bin_intervals(example_DGE, example_data, n = 10)
```

| start | end | number |
|------:|----:|-------:|
|     1 |  50 |     50 |
|    51 | 100 |     50 |
|   101 | 150 |     50 |
|   151 | 200 |     50 |
|   201 | 250 |     50 |
|   251 | 300 |     50 |
|   301 | 350 |     50 |
|   351 | 400 |     50 |
|   401 | 450 |     50 |
|   451 | 500 |     50 |

The first 10 gene bins. Start and end columns indicate row numbers for
the beginning and end of each bin. Number indicates the number of genes
in each bin.

Our goal is to bin genes randomly, so we must randomize the order of the
genes in our DGE list. Without this step, genes would be binned together
based on their sequential order within the chosen gene annotation (e.g.,
ENSEMBL or ENTREZ).

``` r
# Randomize the order of gene identifiers
annot <- example_DGE$genes
annot <- annot[sample(seq_len(nrow(annot))), ]
# Match order of the genes in the DGE to the randomized order of genes in the annotation
example_DGE <- example_DGE[annot, ]
```

## 2: Make imputed data sets for each bin of genes and conduct differential expression analysis

Data are imputed using the mice R package with its default predictive
modeling methods, which are predictive mean matching, logistic
regression, polytomous regression, and proportional odds modeling for
continuous, binary, categorical, and unordered variables, respectively.
The user may specify “m”, the number of imputed datasets, and “maxit”,
the number of iterations for each imputation (default = 10). M imputed
datasets are created separately for each gene bin, where the imputation
predictor matrix includes all covariates along with the log-CPM for all
the genes in a particular bin. Thus, each gene bin contains M sets of
imputed data.

The `impute_by_gene_bin` function loops through a DGE list using the
gene bin intervals from the `get_gene_bin_intervals` function. It
returns a list of sets of m imputed datasets, one per gene bin. For
instance, if m = 100 and intervals contains 200 gene bin intervals,
output will be a list of 200 sets of 100 imputed datasets. Each of the
200 sets are imputed using only the genes in one gene bin. This same
procedure can be run more quickly in parallel with the
`impute_by_gene_bin_parallel` function, which uses the ‘PSOCK’ backend
implemented via the doParallel package.

``` r
gene_bin_impute <- impute_by_gene_bin(example_data,
    intervals,
    example_DGE,
    m = 3
)
```

## 3: Estimate gene expression changes using voom followed by lmFit functions, separately on each M imputed dataset within each gene bin

The `limmavoom_imputed_data_list` function loops through the imputed
data list (output from `impute_by_gene_bin` function) and runs RNA-seq
analysis with the limma-voom pipeline. Users specify the formula for the
RNA-seq design matrix and the main predictor of interest for which log
fold-changes will be estimated. This procedure can also be run in
parallel with the `limmavoom_imputed_data_list_parallel` function.

``` r
coef_se <- limmavoom_imputed_data_list(
    gene_intervals = intervals,
    DGE = example_DGE,
    imputed_data_list = gene_bin_impute,
    m = 3,
    voom_formula = "~x + y + z + a + b",
    predictor = "x"
)
```

## 4-6: Apply variance shrinking Bayesian procedure, pooling results with Rubins’ rules, and FDR-adjust P-values

The final step is to combine the results from each imputed dataset using
Rubin’s rules. The argument “model_results” is the output from the
`limmavoom_imputed_data_list` function above. The `combine_rubins`
function applies the `squeezeVar` function before pooling results. The
result is a table with one row per gene. The table includes coefficients
(e.g., logFC values) standard errors, degrees of freedom, t-statistics,
P-Values, and adjusted P-values from the limma-voom pipeline. Both the
raw and empirical Bayes moderated statistics are reported (see limma
package for details).

``` r
final_res <- combine_rubins(
    DGE = example_DGE,
    model_results = coef_se,
    voom_formula = "~x + y + z + a + b"
)
```

| probe  | coef_combined | combined_p_bayes | combined_p_adj_bayes |
|:-------|--------------:|-----------------:|---------------------:|
| ENS432 |        -0.021 |            0.000 |                0.076 |
| ENS65  |        -0.023 |            0.002 |                0.488 |
| ENS327 |        -0.010 |            0.006 |                0.925 |
| ENS239 |         0.015 |            0.007 |                0.925 |
| ENS411 |         0.021 |            0.012 |                0.969 |
| ENS260 |         0.011 |            0.017 |                0.969 |
| ENS291 |         0.017 |            0.021 |                0.969 |
| ENS458 |         0.015 |            0.023 |                0.969 |
| ENS219 |         0.012 |            0.026 |                0.969 |
| ENS107 |         0.007 |            0.029 |                0.969 |

The top 10 genes associated with predictor x sorted by lowest P-value

# Downstream gene set / pathway analysis

The compatibility of pathway and gene set enrichment methods with
multiple imputation depends on their inputs. The RNAseqCovarImpute
multiple imputation method produces one final list of genes with their
associated t-statistics, log fold changes, and P-values for differential
expression. Thus, the method is compatible with gene set enrichment
analyses that utilize gene rankings such as ora, or gene level
statistics such as `camera`, and `gage`. However, RNAseqCovarImpute is
not compatible with gene set enrichment analyses that require as input a
gene expression matrix or data at the individual sample level, as the
nature of multiple imputation requires the creation of multiple gene
expression matrices across each imputed dataset.

# Session info

    ## R Under development (unstable) (2023-03-08 r83956 ucrt)
    ## Platform: x86_64-w64-mingw32/x64 (64-bit)
    ## Running under: Windows 10 x64 (build 19045)
    ## 
    ## Matrix products: default
    ## 
    ## 
    ## locale:
    ## [1] LC_COLLATE=English_United States.utf8 
    ## [2] LC_CTYPE=English_United States.utf8   
    ## [3] LC_MONETARY=English_United States.utf8
    ## [4] LC_NUMERIC=C                          
    ## [5] LC_TIME=English_United States.utf8    
    ## 
    ## time zone: America/Los_Angeles
    ## tzcode source: internal
    ## 
    ## attached base packages:
    ## [1] stats     graphics  grDevices utils     datasets  methods   base     
    ## 
    ## other attached packages:
    ## [1] dplyr_1.1.0              RNAseqCovarImpute_0.99.7 BiocParallel_1.33.9     
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] limma_3.55.5     compiler_4.3.0   tidyselect_1.2.0 Rcpp_1.0.10     
    ##  [5] parallel_4.3.0   mice_3.15.0      tidyr_1.3.0      yaml_2.3.7      
    ##  [9] fastmap_1.1.1    lattice_0.20-45  R6_2.5.1         generics_0.1.3  
    ## [13] knitr_1.43       iterators_1.0.14 backports_1.4.1  tibble_3.2.0    
    ## [17] snow_0.4-4       pillar_1.9.0     rlang_1.1.1      utf8_1.2.3      
    ## [21] broom_1.0.5      xfun_0.39        cli_3.6.0        withr_2.5.0     
    ## [25] magrittr_2.0.3   digest_0.6.31    foreach_1.5.2    grid_4.3.0      
    ## [29] locfit_1.5-9.7   rstudioapi_0.14  edgeR_3.41.6     lifecycle_1.0.3 
    ## [33] vctrs_0.5.2      evaluate_0.21    glue_1.6.2       codetools_0.2-19
    ## [37] fansi_1.0.4      rmarkdown_2.20   purrr_1.0.1      tools_4.3.0     
    ## [41] pkgconfig_2.0.3  htmltools_0.5.4
