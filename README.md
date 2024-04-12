Impute Covariate Data in RNA-sequencing Studies
================

# Introduction

The RNAseqCovarImpute package makes linear model analysis for RNA-seq
read counts compatible with multiple imputation (MI) of missing
covariates. Relying on the Bioconductor `limma` package,
RNAseqCovarImpute is included in Bioconductor as an extension of the
[variance modeling at the observational level (voom)
method](https://doi.org/10.1186/gb-2014-15-2-r29) which can be applied
in circumstances with missing covariate data.

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
standard errors. MI addresses this problem by generating several
predictions, thereby allowing for uncertainty about the missing data. In
a typical MI procedure: 1) M imputed data sets are created, 2) each data
set is analyzed separately (e.g., using linear regression), and 3)
estimates and standard errors across the M analyses are pooled using
Rubin’s rules. A major problem with implementing MI in RNA sequencing
studies is that the outcome data must be included in the imputation
prediction models to avoid bias. This is difficult in omics studies with
high-dimensional data.

The first method we developed in the RNAseqCovarImpute package surmounts
the problem of high-dimensional outcome data by binning genes into
smaller groups to analyze pseudo-independently. This method implements
covariate MI in gene expression studies by 1) randomly binning genes
into smaller groups, 2) creating M imputed datasets separately within
each bin, where the imputation predictor matrix includes all covariates
and the log counts per million (CPM) for the genes within each bin, 3)
estimating gene expression changes using `limma::voom` followed by
`limma::lmFit` functions, separately on each M imputed dataset within
each gene bin, 4) un-binning the gene sets and stacking the M sets of
model results before applying the `limma::squeezeVar` function to apply
a variance shrinking Bayesian procedure to each M set of model results,
5) pooling the results with Rubins’ rules to produce combined
coefficients, standard errors, and P-values, and 6) adjusting P-values
for multiplicity to account for false discovery rate (FDR).

A faster method uses principal component analysis (PCA) to avoid binning
genes while still retaining outcome information in the MI models.
Binning genes into smaller groups requires that the MI and limma-voom
analysis is run many times (typically hundreds). The more
computationally efficient MI PCA method implements covariate MI in gene
expression studies by 1) performing PCA on the log CPM values for all
genes using the Bioconductor `PCAtools` package, 2) creating M imputed
datasets where the imputation predictor matrix includes all covariates
and the optimum number of PCs to retain (e.g., based on Horn’s parallel
analysis or the number of PCs that account for \>80% explained
variation), 3) conducting the standard limma-voom pipeline with the
`voom` followed by `lmFit` followed by `eBayes` functions on each M
imputed dataset, 4) pooling the results with Rubins’ rules to produce
combined coefficients, standard errors, and P-values, and 5) adjusting
P-values for multiplicity to account for false discovery rate (FDR).

# Installation

``` r
# Install from github
library(devtools)
install_github("brennanhilton/RNAseqCovarImpute")

# Install from Bioconductor (not yet on Bioconductor)

if (!require("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
}

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
```

    ## 
    ## Attaching package: 'dplyr'

    ## The following objects are masked from 'package:stats':
    ## 
    ##     filter, lag

    ## The following objects are masked from 'package:base':
    ## 
    ##     intersect, setdiff, setequal, union

``` r
library(BiocParallel)
library(PCAtools)
```

    ## Loading required package: ggplot2

    ## Loading required package: ggrepel

    ## 
    ## Attaching package: 'PCAtools'

    ## The following objects are masked from 'package:stats':
    ## 
    ##     biplot, screeplot

``` r
library(limma)
library(mice)
```

    ## 
    ## Attaching package: 'mice'

    ## The following object is masked from 'package:stats':
    ## 
    ##     filter

    ## The following objects are masked from 'package:base':
    ## 
    ##     cbind, rbind

``` r
data(example_data)
data(example_DGE)
```

# RNAseqCovarImpute Demonstration

## MI PCA method

We demonstrate the MI PCA method first, as it is far faster and performs
just as well, or even better in many cases, compared with the original
gene binning method.

### Conduct PCA

``` r
# We use voom to convert counts to logCPM values, adding 0.5 to all the counts to avoid taking the logarithm of zero, and normalized for library size.
start.pca = Sys.time() # To calculate runtime
pca_data = limma::voom(example_DGE)$E
# Conduct pca
p = PCAtools::pca(pca_data)
# Determine the number of PCs that account for >80% explained variation
which(cumsum(p$variance) > 80)[1]
```

    ## PC78 
    ##   78

``` r
# Extract the PCs and append to our data
pcs = p$rotated[,1:78]
example_data = cbind(example_data, pcs)
```

### Conduct MI with mice

``` r
# Conduct mice. In practice, m should be much larger (e.g., 10-100)
imp = mice::mice(example_data, m=3)
```

    ## 
    ##  iter imp variable
    ##   1   1  y  z  a  b
    ##   1   2  y  z  a  b
    ##   1   3  y  z  a  b
    ##   2   1  y  z  a  b
    ##   2   2  y  z  a  b
    ##   2   3  y  z  a  b
    ##   3   1  y  z  a  b
    ##   3   2  y  z  a  b
    ##   3   3  y  z  a  b
    ##   4   1  y  z  a  b
    ##   4   2  y  z  a  b
    ##   4   3  y  z  a  b
    ##   5   1  y  z  a  b
    ##   5   2  y  z  a  b
    ##   5   3  y  z  a  b

### Conduct limma-voom analysis

Now we use our `limmavoom_imputed_data_pca` function to run the standard
limma-voom pipeline on each MI dataset and then pool the results with
Rubins’ rules. This procedure is run in parallel using the BiocParallel
package with the default back-end. Users can change the back-end using
the `BPPARAM` argument. This argument is passed to
`BiocParallel::bplapply`. For instance, below we specify to run
`limmavoom_imputed_data_pca` in serial. Users specify the formula for
the RNA-seq design matrix for which log fold-changes will be estimated.
For each contrast specified in the voom_formula, p-values, coefficients,
and standard errors are stored under the “\_p”, “\_coef” and “\_se”
columns, respectively.

``` r
mi_pca_res = RNAseqCovarImpute::limmavoom_imputed_data_pca(imp = imp,
                                                           DGE = example_DGE,
                                                           voom_formula = "~x + y + z + a + b",
                                                           BPPARAM = SerialParam())

# Display the results for the first 5 genes for the x variable in the model.
mi_pca_res[1:5, grep("^x", colnames(mi_pca_res))]
```

    ##          x_p        x_coef        x_se
    ## 1 0.06037941  0.0172257239 0.009159629
    ## 2 0.97214805 -0.0001254919 0.003593184
    ## 3 0.16183199 -0.0043427912 0.003101565
    ## 4 0.35100300 -0.0043740371 0.004687267
    ## 5 0.80060273 -0.0017828965 0.007056803

### Adjust for FDR

For a given contrast based on the variables in the design matrix
(voom_formula), we arrange the genes from lowest to highest p-value and
adjust for FDR.

``` r
mi_pca_res_x = mi_pca_res %>% 
  arrange(x_p) %>% 
  mutate(x_p_adj = p.adjust(x_p, method = "fdr")) %>% 
  dplyr::select(probe, x_coef, x_p, x_p_adj)

end.pca = Sys.time() # To calculate runtime
time.pca = end.pca - start.pca # To calculate runtime
# Display the results for the first 5 genes
mi_pca_res_x[1:5,]
```

    ##    probe      x_coef          x_p    x_p_adj
    ## 1 ENS432 -0.02124923 0.0001616095 0.08080476
    ## 2  ENS65 -0.02252621 0.0019819761 0.49549403
    ## 3 ENS327 -0.01044210 0.0062146233 0.91013120
    ## 4 ENS239  0.01454005 0.0072810496 0.91013120
    ## 5 ENS411  0.02146120 0.0111108560 0.97182262

## Gene binning MI method

### Bin the genes into smaller groups

The default is approximately 1 gene per 10 individuals in the study, but
the user can specify a different ratio. For example, in a study with 500
participants and 10,000 genes, 200 bins of 50 genes would be created
using the default ratio. This example demonstrates why this method can
be computationally intensive: the MI and limma-voom analysis will be
conducted separately 200 times. When the total number of genes is not
divisible by the bin size, the method flexibly creates bins of different
sizes. The order of the features (e.g., ENSEMBL gene identifiers) should
be randomized before binning.

``` r
# Get back the original example_data without the PCs appended
data(example_data)
start.old.method = Sys.time() # To calculate runtime
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

### Make imputed data sets for each bin of genes and conduct differential expression analysis

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
200 sets are imputed using only the genes in one gene bin.

``` r
gene_bin_impute <- impute_by_gene_bin(example_data,
    intervals,
    example_DGE,
    m = 3
)
```

    ## 
    ##  iter imp variable
    ##   1   1  y  z  a  b
    ##   1   2  y  z  a  b
    ##   1   3  y  z  a  b
    ##   2   1  y  z  a  b
    ##   2   2  y  z  a  b
    ##   2   3  y  z  a  b
    ##   3   1  y  z  a  b
    ##   3   2  y  z  a  b
    ##   3   3  y  z  a  b
    ##   4   1  y  z  a  b
    ##   4   2  y  z  a  b
    ##   4   3  y  z  a  b
    ##   5   1  y  z  a  b
    ##   5   2  y  z  a  b
    ##   5   3  y  z  a  b
    ##   6   1  y  z  a  b
    ##   6   2  y  z  a  b
    ##   6   3  y  z  a  b
    ##   7   1  y  z  a  b
    ##   7   2  y  z  a  b
    ##   7   3  y  z  a  b
    ##   8   1  y  z  a  b
    ##   8   2  y  z  a  b
    ##   8   3  y  z  a  b
    ##   9   1  y  z  a  b
    ##   9   2  y  z  a  b
    ##   9   3  y  z  a  b
    ##   10   1  y  z  a  b
    ##   10   2  y  z  a  b
    ##   10   3  y  z  a  b
    ## 
    ## 
    ##  iter imp variable
    ##   1   1  y  z  a  b
    ##   1   2  y  z  a  b
    ##   1   3  y  z  a  b
    ##   2   1  y  z  a  b
    ##   2   2  y  z  a  b
    ##   2   3  y  z  a  b
    ##   3   1  y  z  a  b
    ##   3   2  y  z  a  b
    ##   3   3  y  z  a  b
    ##   4   1  y  z  a  b
    ##   4   2  y  z  a  b
    ##   4   3  y  z  a  b
    ##   5   1  y  z  a  b
    ##   5   2  y  z  a  b
    ##   5   3  y  z  a  b
    ##   6   1  y  z  a  b
    ##   6   2  y  z  a  b
    ##   6   3  y  z  a  b
    ##   7   1  y  z  a  b
    ##   7   2  y  z  a  b
    ##   7   3  y  z  a  b
    ##   8   1  y  z  a  b
    ##   8   2  y  z  a  b
    ##   8   3  y  z  a  b
    ##   9   1  y  z  a  b
    ##   9   2  y  z  a  b
    ##   9   3  y  z  a  b
    ##   10   1  y  z  a  b
    ##   10   2  y  z  a  b
    ##   10   3  y  z  a  b
    ## 
    ## 
    ##  iter imp variable
    ##   1   1  y  z  a  b
    ##   1   2  y  z  a  b
    ##   1   3  y  z  a  b
    ##   2   1  y  z  a  b
    ##   2   2  y  z  a  b
    ##   2   3  y  z  a  b
    ##   3   1  y  z  a  b
    ##   3   2  y  z  a  b
    ##   3   3  y  z  a  b
    ##   4   1  y  z  a  b
    ##   4   2  y  z  a  b
    ##   4   3  y  z  a  b
    ##   5   1  y  z  a  b
    ##   5   2  y  z  a  b
    ##   5   3  y  z  a  b
    ##   6   1  y  z  a  b
    ##   6   2  y  z  a  b
    ##   6   3  y  z  a  b
    ##   7   1  y  z  a  b
    ##   7   2  y  z  a  b
    ##   7   3  y  z  a  b
    ##   8   1  y  z  a  b
    ##   8   2  y  z  a  b
    ##   8   3  y  z  a  b
    ##   9   1  y  z  a  b
    ##   9   2  y  z  a  b
    ##   9   3  y  z  a  b
    ##   10   1  y  z  a  b
    ##   10   2  y  z  a  b
    ##   10   3  y  z  a  b
    ## 
    ## 
    ##  iter imp variable
    ##   1   1  y  z  a  b
    ##   1   2  y  z  a  b
    ##   1   3  y  z  a  b
    ##   2   1  y  z  a  b
    ##   2   2  y  z  a  b
    ##   2   3  y  z  a  b
    ##   3   1  y  z  a  b
    ##   3   2  y  z  a  b
    ##   3   3  y  z  a  b
    ##   4   1  y  z  a  b
    ##   4   2  y  z  a  b
    ##   4   3  y  z  a  b
    ##   5   1  y  z  a  b
    ##   5   2  y  z  a  b
    ##   5   3  y  z  a  b
    ##   6   1  y  z  a  b
    ##   6   2  y  z  a  b
    ##   6   3  y  z  a  b
    ##   7   1  y  z  a  b
    ##   7   2  y  z  a  b
    ##   7   3  y  z  a  b
    ##   8   1  y  z  a  b
    ##   8   2  y  z  a  b
    ##   8   3  y  z  a  b
    ##   9   1  y  z  a  b
    ##   9   2  y  z  a  b
    ##   9   3  y  z  a  b
    ##   10   1  y  z  a  b
    ##   10   2  y  z  a  b
    ##   10   3  y  z  a  b
    ## 
    ## 
    ##  iter imp variable
    ##   1   1  y  z  a  b
    ##   1   2  y  z  a  b
    ##   1   3  y  z  a  b
    ##   2   1  y  z  a  b
    ##   2   2  y  z  a  b
    ##   2   3  y  z  a  b
    ##   3   1  y  z  a  b
    ##   3   2  y  z  a  b
    ##   3   3  y  z  a  b
    ##   4   1  y  z  a  b
    ##   4   2  y  z  a  b
    ##   4   3  y  z  a  b
    ##   5   1  y  z  a  b
    ##   5   2  y  z  a  b
    ##   5   3  y  z  a  b
    ##   6   1  y  z  a  b
    ##   6   2  y  z  a  b
    ##   6   3  y  z  a  b
    ##   7   1  y  z  a  b
    ##   7   2  y  z  a  b
    ##   7   3  y  z  a  b
    ##   8   1  y  z  a  b
    ##   8   2  y  z  a  b
    ##   8   3  y  z  a  b
    ##   9   1  y  z  a  b
    ##   9   2  y  z  a  b
    ##   9   3  y  z  a  b
    ##   10   1  y  z  a  b
    ##   10   2  y  z  a  b
    ##   10   3  y  z  a  b
    ## 
    ## 
    ##  iter imp variable
    ##   1   1  y  z  a  b
    ##   1   2  y  z  a  b
    ##   1   3  y  z  a  b
    ##   2   1  y  z  a  b
    ##   2   2  y  z  a  b
    ##   2   3  y  z  a  b
    ##   3   1  y  z  a  b
    ##   3   2  y  z  a  b
    ##   3   3  y  z  a  b
    ##   4   1  y  z  a  b
    ##   4   2  y  z  a  b
    ##   4   3  y  z  a  b
    ##   5   1  y  z  a  b
    ##   5   2  y  z  a  b
    ##   5   3  y  z  a  b
    ##   6   1  y  z  a  b
    ##   6   2  y  z  a  b
    ##   6   3  y  z  a  b
    ##   7   1  y  z  a  b
    ##   7   2  y  z  a  b
    ##   7   3  y  z  a  b
    ##   8   1  y  z  a  b
    ##   8   2  y  z  a  b
    ##   8   3  y  z  a  b
    ##   9   1  y  z  a  b
    ##   9   2  y  z  a  b
    ##   9   3  y  z  a  b
    ##   10   1  y  z  a  b
    ##   10   2  y  z  a  b
    ##   10   3  y  z  a  b
    ## 
    ## 
    ##  iter imp variable
    ##   1   1  y  z  a  b
    ##   1   2  y  z  a  b
    ##   1   3  y  z  a  b
    ##   2   1  y  z  a  b
    ##   2   2  y  z  a  b
    ##   2   3  y  z  a  b
    ##   3   1  y  z  a  b
    ##   3   2  y  z  a  b
    ##   3   3  y  z  a  b
    ##   4   1  y  z  a  b
    ##   4   2  y  z  a  b
    ##   4   3  y  z  a  b
    ##   5   1  y  z  a  b
    ##   5   2  y  z  a  b
    ##   5   3  y  z  a  b
    ##   6   1  y  z  a  b
    ##   6   2  y  z  a  b
    ##   6   3  y  z  a  b
    ##   7   1  y  z  a  b
    ##   7   2  y  z  a  b
    ##   7   3  y  z  a  b
    ##   8   1  y  z  a  b
    ##   8   2  y  z  a  b
    ##   8   3  y  z  a  b
    ##   9   1  y  z  a  b
    ##   9   2  y  z  a  b
    ##   9   3  y  z  a  b
    ##   10   1  y  z  a  b
    ##   10   2  y  z  a  b
    ##   10   3  y  z  a  b
    ## 
    ## 
    ##  iter imp variable
    ##   1   1  y  z  a  b
    ##   1   2  y  z  a  b
    ##   1   3  y  z  a  b
    ##   2   1  y  z  a  b
    ##   2   2  y  z  a  b
    ##   2   3  y  z  a  b
    ##   3   1  y  z  a  b
    ##   3   2  y  z  a  b
    ##   3   3  y  z  a  b
    ##   4   1  y  z  a  b
    ##   4   2  y  z  a  b
    ##   4   3  y  z  a  b
    ##   5   1  y  z  a  b
    ##   5   2  y  z  a  b
    ##   5   3  y  z  a  b
    ##   6   1  y  z  a  b
    ##   6   2  y  z  a  b
    ##   6   3  y  z  a  b
    ##   7   1  y  z  a  b
    ##   7   2  y  z  a  b
    ##   7   3  y  z  a  b
    ##   8   1  y  z  a  b
    ##   8   2  y  z  a  b
    ##   8   3  y  z  a  b
    ##   9   1  y  z  a  b
    ##   9   2  y  z  a  b
    ##   9   3  y  z  a  b
    ##   10   1  y  z  a  b
    ##   10   2  y  z  a  b
    ##   10   3  y  z  a  b
    ## 
    ## 
    ##  iter imp variable
    ##   1   1  y  z  a  b
    ##   1   2  y  z  a  b
    ##   1   3  y  z  a  b
    ##   2   1  y  z  a  b
    ##   2   2  y  z  a  b
    ##   2   3  y  z  a  b
    ##   3   1  y  z  a  b
    ##   3   2  y  z  a  b
    ##   3   3  y  z  a  b
    ##   4   1  y  z  a  b
    ##   4   2  y  z  a  b
    ##   4   3  y  z  a  b
    ##   5   1  y  z  a  b
    ##   5   2  y  z  a  b
    ##   5   3  y  z  a  b
    ##   6   1  y  z  a  b
    ##   6   2  y  z  a  b
    ##   6   3  y  z  a  b
    ##   7   1  y  z  a  b
    ##   7   2  y  z  a  b
    ##   7   3  y  z  a  b
    ##   8   1  y  z  a  b
    ##   8   2  y  z  a  b
    ##   8   3  y  z  a  b
    ##   9   1  y  z  a  b
    ##   9   2  y  z  a  b
    ##   9   3  y  z  a  b
    ##   10   1  y  z  a  b
    ##   10   2  y  z  a  b
    ##   10   3  y  z  a  b
    ## 
    ## 
    ##  iter imp variable
    ##   1   1  y  z  a  b
    ##   1   2  y  z  a  b
    ##   1   3  y  z  a  b
    ##   2   1  y  z  a  b
    ##   2   2  y  z  a  b
    ##   2   3  y  z  a  b
    ##   3   1  y  z  a  b
    ##   3   2  y  z  a  b
    ##   3   3  y  z  a  b
    ##   4   1  y  z  a  b
    ##   4   2  y  z  a  b
    ##   4   3  y  z  a  b
    ##   5   1  y  z  a  b
    ##   5   2  y  z  a  b
    ##   5   3  y  z  a  b
    ##   6   1  y  z  a  b
    ##   6   2  y  z  a  b
    ##   6   3  y  z  a  b
    ##   7   1  y  z  a  b
    ##   7   2  y  z  a  b
    ##   7   3  y  z  a  b
    ##   8   1  y  z  a  b
    ##   8   2  y  z  a  b
    ##   8   3  y  z  a  b
    ##   9   1  y  z  a  b
    ##   9   2  y  z  a  b
    ##   9   3  y  z  a  b
    ##   10   1  y  z  a  b
    ##   10   2  y  z  a  b
    ##   10   3  y  z  a  b

This procedure is run in parallel using the BiocParallel package with
the default back-end. Users can change the back-end using the `BPPARAM`
argument. This argument is passed to `BiocParallel::bplapply`. For
instance, to run `gene_bin_impute` in serial:

``` r
gene_bin_impute <- impute_by_gene_bin(example_data,
    intervals,
    example_DGE,
    m = 3,
    BPPARAM = SerialParam()
)
```

    ## 
    ##  iter imp variable
    ##   1   1  y  z  a  b
    ##   1   2  y  z  a  b
    ##   1   3  y  z  a  b
    ##   2   1  y  z  a  b
    ##   2   2  y  z  a  b
    ##   2   3  y  z  a  b
    ##   3   1  y  z  a  b
    ##   3   2  y  z  a  b
    ##   3   3  y  z  a  b
    ##   4   1  y  z  a  b
    ##   4   2  y  z  a  b
    ##   4   3  y  z  a  b
    ##   5   1  y  z  a  b
    ##   5   2  y  z  a  b
    ##   5   3  y  z  a  b
    ##   6   1  y  z  a  b
    ##   6   2  y  z  a  b
    ##   6   3  y  z  a  b
    ##   7   1  y  z  a  b
    ##   7   2  y  z  a  b
    ##   7   3  y  z  a  b
    ##   8   1  y  z  a  b
    ##   8   2  y  z  a  b
    ##   8   3  y  z  a  b
    ##   9   1  y  z  a  b
    ##   9   2  y  z  a  b
    ##   9   3  y  z  a  b
    ##   10   1  y  z  a  b
    ##   10   2  y  z  a  b
    ##   10   3  y  z  a  b
    ## 
    ##  iter imp variable
    ##   1   1  y  z  a  b
    ##   1   2  y  z  a  b
    ##   1   3  y  z  a  b
    ##   2   1  y  z  a  b
    ##   2   2  y  z  a  b
    ##   2   3  y  z  a  b
    ##   3   1  y  z  a  b
    ##   3   2  y  z  a  b
    ##   3   3  y  z  a  b
    ##   4   1  y  z  a  b
    ##   4   2  y  z  a  b
    ##   4   3  y  z  a  b
    ##   5   1  y  z  a  b
    ##   5   2  y  z  a  b
    ##   5   3  y  z  a  b
    ##   6   1  y  z  a  b
    ##   6   2  y  z  a  b
    ##   6   3  y  z  a  b
    ##   7   1  y  z  a  b
    ##   7   2  y  z  a  b
    ##   7   3  y  z  a  b
    ##   8   1  y  z  a  b
    ##   8   2  y  z  a  b
    ##   8   3  y  z  a  b
    ##   9   1  y  z  a  b
    ##   9   2  y  z  a  b
    ##   9   3  y  z  a  b
    ##   10   1  y  z  a  b
    ##   10   2  y  z  a  b
    ##   10   3  y  z  a  b
    ## 
    ##  iter imp variable
    ##   1   1  y  z  a  b
    ##   1   2  y  z  a  b
    ##   1   3  y  z  a  b
    ##   2   1  y  z  a  b
    ##   2   2  y  z  a  b
    ##   2   3  y  z  a  b
    ##   3   1  y  z  a  b
    ##   3   2  y  z  a  b
    ##   3   3  y  z  a  b
    ##   4   1  y  z  a  b
    ##   4   2  y  z  a  b
    ##   4   3  y  z  a  b
    ##   5   1  y  z  a  b
    ##   5   2  y  z  a  b
    ##   5   3  y  z  a  b
    ##   6   1  y  z  a  b
    ##   6   2  y  z  a  b
    ##   6   3  y  z  a  b
    ##   7   1  y  z  a  b
    ##   7   2  y  z  a  b
    ##   7   3  y  z  a  b
    ##   8   1  y  z  a  b
    ##   8   2  y  z  a  b
    ##   8   3  y  z  a  b
    ##   9   1  y  z  a  b
    ##   9   2  y  z  a  b
    ##   9   3  y  z  a  b
    ##   10   1  y  z  a  b
    ##   10   2  y  z  a  b
    ##   10   3  y  z  a  b
    ## 
    ##  iter imp variable
    ##   1   1  y  z  a  b
    ##   1   2  y  z  a  b
    ##   1   3  y  z  a  b
    ##   2   1  y  z  a  b
    ##   2   2  y  z  a  b
    ##   2   3  y  z  a  b
    ##   3   1  y  z  a  b
    ##   3   2  y  z  a  b
    ##   3   3  y  z  a  b
    ##   4   1  y  z  a  b
    ##   4   2  y  z  a  b
    ##   4   3  y  z  a  b
    ##   5   1  y  z  a  b
    ##   5   2  y  z  a  b
    ##   5   3  y  z  a  b
    ##   6   1  y  z  a  b
    ##   6   2  y  z  a  b
    ##   6   3  y  z  a  b
    ##   7   1  y  z  a  b
    ##   7   2  y  z  a  b
    ##   7   3  y  z  a  b
    ##   8   1  y  z  a  b
    ##   8   2  y  z  a  b
    ##   8   3  y  z  a  b
    ##   9   1  y  z  a  b
    ##   9   2  y  z  a  b
    ##   9   3  y  z  a  b
    ##   10   1  y  z  a  b
    ##   10   2  y  z  a  b
    ##   10   3  y  z  a  b
    ## 
    ##  iter imp variable
    ##   1   1  y  z  a  b
    ##   1   2  y  z  a  b
    ##   1   3  y  z  a  b
    ##   2   1  y  z  a  b
    ##   2   2  y  z  a  b
    ##   2   3  y  z  a  b
    ##   3   1  y  z  a  b
    ##   3   2  y  z  a  b
    ##   3   3  y  z  a  b
    ##   4   1  y  z  a  b
    ##   4   2  y  z  a  b
    ##   4   3  y  z  a  b
    ##   5   1  y  z  a  b
    ##   5   2  y  z  a  b
    ##   5   3  y  z  a  b
    ##   6   1  y  z  a  b
    ##   6   2  y  z  a  b
    ##   6   3  y  z  a  b
    ##   7   1  y  z  a  b
    ##   7   2  y  z  a  b
    ##   7   3  y  z  a  b
    ##   8   1  y  z  a  b
    ##   8   2  y  z  a  b
    ##   8   3  y  z  a  b
    ##   9   1  y  z  a  b
    ##   9   2  y  z  a  b
    ##   9   3  y  z  a  b
    ##   10   1  y  z  a  b
    ##   10   2  y  z  a  b
    ##   10   3  y  z  a  b
    ## 
    ##  iter imp variable
    ##   1   1  y  z  a  b
    ##   1   2  y  z  a  b
    ##   1   3  y  z  a  b
    ##   2   1  y  z  a  b
    ##   2   2  y  z  a  b
    ##   2   3  y  z  a  b
    ##   3   1  y  z  a  b
    ##   3   2  y  z  a  b
    ##   3   3  y  z  a  b
    ##   4   1  y  z  a  b
    ##   4   2  y  z  a  b
    ##   4   3  y  z  a  b
    ##   5   1  y  z  a  b
    ##   5   2  y  z  a  b
    ##   5   3  y  z  a  b
    ##   6   1  y  z  a  b
    ##   6   2  y  z  a  b
    ##   6   3  y  z  a  b
    ##   7   1  y  z  a  b
    ##   7   2  y  z  a  b
    ##   7   3  y  z  a  b
    ##   8   1  y  z  a  b
    ##   8   2  y  z  a  b
    ##   8   3  y  z  a  b
    ##   9   1  y  z  a  b
    ##   9   2  y  z  a  b
    ##   9   3  y  z  a  b
    ##   10   1  y  z  a  b
    ##   10   2  y  z  a  b
    ##   10   3  y  z  a  b
    ## 
    ##  iter imp variable
    ##   1   1  y  z  a  b
    ##   1   2  y  z  a  b
    ##   1   3  y  z  a  b
    ##   2   1  y  z  a  b
    ##   2   2  y  z  a  b
    ##   2   3  y  z  a  b
    ##   3   1  y  z  a  b
    ##   3   2  y  z  a  b
    ##   3   3  y  z  a  b
    ##   4   1  y  z  a  b
    ##   4   2  y  z  a  b
    ##   4   3  y  z  a  b
    ##   5   1  y  z  a  b
    ##   5   2  y  z  a  b
    ##   5   3  y  z  a  b
    ##   6   1  y  z  a  b
    ##   6   2  y  z  a  b
    ##   6   3  y  z  a  b
    ##   7   1  y  z  a  b
    ##   7   2  y  z  a  b
    ##   7   3  y  z  a  b
    ##   8   1  y  z  a  b
    ##   8   2  y  z  a  b
    ##   8   3  y  z  a  b
    ##   9   1  y  z  a  b
    ##   9   2  y  z  a  b
    ##   9   3  y  z  a  b
    ##   10   1  y  z  a  b
    ##   10   2  y  z  a  b
    ##   10   3  y  z  a  b
    ## 
    ##  iter imp variable
    ##   1   1  y  z  a  b
    ##   1   2  y  z  a  b
    ##   1   3  y  z  a  b
    ##   2   1  y  z  a  b
    ##   2   2  y  z  a  b
    ##   2   3  y  z  a  b
    ##   3   1  y  z  a  b
    ##   3   2  y  z  a  b
    ##   3   3  y  z  a  b
    ##   4   1  y  z  a  b
    ##   4   2  y  z  a  b
    ##   4   3  y  z  a  b
    ##   5   1  y  z  a  b
    ##   5   2  y  z  a  b
    ##   5   3  y  z  a  b
    ##   6   1  y  z  a  b
    ##   6   2  y  z  a  b
    ##   6   3  y  z  a  b
    ##   7   1  y  z  a  b
    ##   7   2  y  z  a  b
    ##   7   3  y  z  a  b
    ##   8   1  y  z  a  b
    ##   8   2  y  z  a  b
    ##   8   3  y  z  a  b
    ##   9   1  y  z  a  b
    ##   9   2  y  z  a  b
    ##   9   3  y  z  a  b
    ##   10   1  y  z  a  b
    ##   10   2  y  z  a  b
    ##   10   3  y  z  a  b
    ## 
    ##  iter imp variable
    ##   1   1  y  z  a  b
    ##   1   2  y  z  a  b
    ##   1   3  y  z  a  b
    ##   2   1  y  z  a  b
    ##   2   2  y  z  a  b
    ##   2   3  y  z  a  b
    ##   3   1  y  z  a  b
    ##   3   2  y  z  a  b
    ##   3   3  y  z  a  b
    ##   4   1  y  z  a  b
    ##   4   2  y  z  a  b
    ##   4   3  y  z  a  b
    ##   5   1  y  z  a  b
    ##   5   2  y  z  a  b
    ##   5   3  y  z  a  b
    ##   6   1  y  z  a  b
    ##   6   2  y  z  a  b
    ##   6   3  y  z  a  b
    ##   7   1  y  z  a  b
    ##   7   2  y  z  a  b
    ##   7   3  y  z  a  b
    ##   8   1  y  z  a  b
    ##   8   2  y  z  a  b
    ##   8   3  y  z  a  b
    ##   9   1  y  z  a  b
    ##   9   2  y  z  a  b
    ##   9   3  y  z  a  b
    ##   10   1  y  z  a  b
    ##   10   2  y  z  a  b
    ##   10   3  y  z  a  b
    ## 
    ##  iter imp variable
    ##   1   1  y  z  a  b
    ##   1   2  y  z  a  b
    ##   1   3  y  z  a  b
    ##   2   1  y  z  a  b
    ##   2   2  y  z  a  b
    ##   2   3  y  z  a  b
    ##   3   1  y  z  a  b
    ##   3   2  y  z  a  b
    ##   3   3  y  z  a  b
    ##   4   1  y  z  a  b
    ##   4   2  y  z  a  b
    ##   4   3  y  z  a  b
    ##   5   1  y  z  a  b
    ##   5   2  y  z  a  b
    ##   5   3  y  z  a  b
    ##   6   1  y  z  a  b
    ##   6   2  y  z  a  b
    ##   6   3  y  z  a  b
    ##   7   1  y  z  a  b
    ##   7   2  y  z  a  b
    ##   7   3  y  z  a  b
    ##   8   1  y  z  a  b
    ##   8   2  y  z  a  b
    ##   8   3  y  z  a  b
    ##   9   1  y  z  a  b
    ##   9   2  y  z  a  b
    ##   9   3  y  z  a  b
    ##   10   1  y  z  a  b
    ##   10   2  y  z  a  b
    ##   10   3  y  z  a  b

### Estimate gene expression changes using voom followed by lmFit functions, separately on each M imputed dataset within each gene bin

The `limmavoom_imputed_data_list` function loops through the imputed
data list (output from `impute_by_gene_bin` function) and runs RNA-seq
analysis with the limma-voom pipeline. Users specify the formula for the
RNA-seq design matrix for which log fold-changes will be estimated. This
procedure can also be run with a different parallel back-end or in
serial using the `BPPARAM` argument as above.

``` r
coef_se <- limmavoom_imputed_data_list(
    gene_intervals = intervals,
    DGE = example_DGE,
    imputed_data_list = gene_bin_impute,
    m = 3,
    voom_formula = "~x + y + z + a + b"
)
```

    ## Joining with `by = join_by(probe)`
    ## Joining with `by = join_by(probe)`
    ## Joining with `by = join_by(probe)`
    ## Joining with `by = join_by(probe)`
    ## Joining with `by = join_by(probe)`
    ## Joining with `by = join_by(probe)`
    ## Joining with `by = join_by(probe)`
    ## Joining with `by = join_by(probe)`
    ## Joining with `by = join_by(probe)`
    ## Joining with `by = join_by(probe)`
    ## Joining with `by = join_by(probe)`
    ## Joining with `by = join_by(probe)`
    ## Joining with `by = join_by(probe)`
    ## Joining with `by = join_by(probe)`
    ## Joining with `by = join_by(probe)`
    ## Joining with `by = join_by(probe)`
    ## Joining with `by = join_by(probe)`
    ## 
    ## Joining with `by = join_by(probe)`
    ## Joining with `by = join_by(probe)`
    ## Joining with `by = join_by(probe)`
    ## Joining with `by = join_by(probe)`
    ## Joining with `by = join_by(probe)`
    ## Joining with `by = join_by(probe)`
    ## Joining with `by = join_by(probe)`
    ## Joining with `by = join_by(probe)`
    ## Joining with `by = join_by(probe)`
    ## Joining with `by = join_by(probe)`
    ## Joining with `by = join_by(probe)`
    ## Joining with `by = join_by(probe)`
    ## Joining with `by = join_by(probe)`
    ## Joining with `by = join_by(probe)`
    ## Joining with `by = join_by(probe)`
    ## Joining with `by = join_by(probe)`
    ## Joining with `by = join_by(probe)`
    ## 
    ## Joining with `by = join_by(probe)`
    ## Joining with `by = join_by(probe)`
    ## Joining with `by = join_by(probe)`
    ## Joining with `by = join_by(probe)`
    ## Joining with `by = join_by(probe)`
    ## Joining with `by = join_by(probe)`
    ## Joining with `by = join_by(probe)`
    ## Joining with `by = join_by(probe)`
    ## Joining with `by = join_by(probe)`
    ## Joining with `by = join_by(probe)`
    ## Joining with `by = join_by(probe)`
    ## Joining with `by = join_by(probe)`
    ## Joining with `by = join_by(probe)`
    ## Joining with `by = join_by(probe)`
    ## Joining with `by = join_by(probe)`
    ## Joining with `by = join_by(probe)`
    ## Joining with `by = join_by(probe)`
    ## 
    ## Joining with `by = join_by(probe)`
    ## Joining with `by = join_by(probe)`
    ## Joining with `by = join_by(probe)`
    ## Joining with `by = join_by(probe)`
    ## Joining with `by = join_by(probe)`
    ## Joining with `by = join_by(probe)`
    ## Joining with `by = join_by(probe)`
    ## Joining with `by = join_by(probe)`
    ## Joining with `by = join_by(probe)`
    ## Joining with `by = join_by(probe)`
    ## Joining with `by = join_by(probe)`
    ## Joining with `by = join_by(probe)`
    ## Joining with `by = join_by(probe)`
    ## Joining with `by = join_by(probe)`
    ## Joining with `by = join_by(probe)`
    ## Joining with `by = join_by(probe)`
    ## Joining with `by = join_by(probe)`
    ## 
    ## Joining with `by = join_by(probe)`
    ## Joining with `by = join_by(probe)`
    ## Joining with `by = join_by(probe)`
    ## Joining with `by = join_by(probe)`
    ## Joining with `by = join_by(probe)`
    ## Joining with `by = join_by(probe)`
    ## Joining with `by = join_by(probe)`
    ## Joining with `by = join_by(probe)`
    ## Joining with `by = join_by(probe)`
    ## Joining with `by = join_by(probe)`
    ## Joining with `by = join_by(probe)`
    ## Joining with `by = join_by(probe)`
    ## Joining with `by = join_by(probe)`
    ## Joining with `by = join_by(probe)`
    ## Joining with `by = join_by(probe)`
    ## Joining with `by = join_by(probe)`
    ## Joining with `by = join_by(probe)`
    ## 
    ## Joining with `by = join_by(probe)`
    ## Joining with `by = join_by(probe)`
    ## Joining with `by = join_by(probe)`
    ## Joining with `by = join_by(probe)`
    ## Joining with `by = join_by(probe)`
    ## Joining with `by = join_by(probe)`
    ## Joining with `by = join_by(probe)`
    ## Joining with `by = join_by(probe)`
    ## Joining with `by = join_by(probe)`
    ## Joining with `by = join_by(probe)`
    ## Joining with `by = join_by(probe)`
    ## Joining with `by = join_by(probe)`
    ## Joining with `by = join_by(probe)`
    ## Joining with `by = join_by(probe)`
    ## Joining with `by = join_by(probe)`
    ## Joining with `by = join_by(probe)`
    ## Joining with `by = join_by(probe)`
    ## 
    ## Joining with `by = join_by(probe)`
    ## Joining with `by = join_by(probe)`
    ## Joining with `by = join_by(probe)`
    ## Joining with `by = join_by(probe)`
    ## Joining with `by = join_by(probe)`
    ## Joining with `by = join_by(probe)`
    ## Joining with `by = join_by(probe)`
    ## Joining with `by = join_by(probe)`
    ## Joining with `by = join_by(probe)`
    ## Joining with `by = join_by(probe)`
    ## Joining with `by = join_by(probe)`
    ## Joining with `by = join_by(probe)`
    ## Joining with `by = join_by(probe)`
    ## Joining with `by = join_by(probe)`
    ## Joining with `by = join_by(probe)`
    ## Joining with `by = join_by(probe)`
    ## Joining with `by = join_by(probe)`
    ## 
    ## Joining with `by = join_by(probe)`
    ## Joining with `by = join_by(probe)`
    ## Joining with `by = join_by(probe)`
    ## Joining with `by = join_by(probe)`
    ## Joining with `by = join_by(probe)`
    ## Joining with `by = join_by(probe)`
    ## Joining with `by = join_by(probe)`
    ## Joining with `by = join_by(probe)`
    ## Joining with `by = join_by(probe)`
    ## Joining with `by = join_by(probe)`
    ## Joining with `by = join_by(probe)`
    ## Joining with `by = join_by(probe)`
    ## Joining with `by = join_by(probe)`
    ## Joining with `by = join_by(probe)`
    ## Joining with `by = join_by(probe)`
    ## Joining with `by = join_by(probe)`
    ## Joining with `by = join_by(probe)`
    ## 
    ## Joining with `by = join_by(probe)`
    ## Joining with `by = join_by(probe)`
    ## Joining with `by = join_by(probe)`
    ## Joining with `by = join_by(probe)`
    ## Joining with `by = join_by(probe)`
    ## Joining with `by = join_by(probe)`
    ## Joining with `by = join_by(probe)`
    ## Joining with `by = join_by(probe)`
    ## Joining with `by = join_by(probe)`
    ## Joining with `by = join_by(probe)`
    ## Joining with `by = join_by(probe)`
    ## Joining with `by = join_by(probe)`
    ## Joining with `by = join_by(probe)`
    ## Joining with `by = join_by(probe)`
    ## Joining with `by = join_by(probe)`
    ## Joining with `by = join_by(probe)`
    ## Joining with `by = join_by(probe)`
    ## 
    ## Joining with `by = join_by(probe)`
    ## Joining with `by = join_by(probe)`
    ## Joining with `by = join_by(probe)`
    ## Joining with `by = join_by(probe)`
    ## Joining with `by = join_by(probe)`
    ## Joining with `by = join_by(probe)`
    ## Joining with `by = join_by(probe)`
    ## Joining with `by = join_by(probe)`
    ## Joining with `by = join_by(probe)`
    ## Joining with `by = join_by(probe)`
    ## Joining with `by = join_by(probe)`
    ## Joining with `by = join_by(probe)`
    ## Joining with `by = join_by(probe)`
    ## Joining with `by = join_by(probe)`
    ## Joining with `by = join_by(probe)`
    ## Joining with `by = join_by(probe)`
    ## Joining with `by = join_by(probe)`

### Apply variance shrinking Bayesian procedure, pooling results with Rubins’ rules, and FDR-adjust P-values

The final step is to combine the results from each imputed dataset using
Rubin’s rules. The argument “model_results” is the output from the
`limmavoom_imputed_data_list` function above. The `combine_rubins`
function applies the `squeezeVar` function before pooling results. The
result is a table with one row per gene. The table includes coefficients
(e.g., logFC values) standard errors, degrees of freedom, t-statistics,
P-Values, and adjusted P-values from the limma-voom pipeline. Both the
raw and empirical Bayes moderated statistics are reported. The user
selects the predictor of interest in the form of a linear model contrast
for which model results will be extracted. For a continuous variable
this is just the predictor name. For a categorical variable like `b` in
`example_data` we could specify `predictor = b1` to get the effect of
being in the b = 1 versus the b = 0 group.

``` r
final_res <- combine_rubins(
    DGE = example_DGE,
    model_results = coef_se,
    predictor = "x"
)

end.old.method = Sys.time() # To calculate runtime
time.old.method = end.old.method - start.old.method
```

| probe  | coef_combined | combined_p_bayes | combined_p_adj_bayes |
|:-------|--------------:|-----------------:|---------------------:|
| ENS432 |        -0.021 |            0.000 |                0.092 |
| ENS65  |        -0.022 |            0.002 |                0.537 |
| ENS327 |        -0.011 |            0.006 |                0.987 |
| ENS239 |         0.014 |            0.009 |                0.987 |
| ENS411 |         0.022 |            0.010 |                0.987 |
| ENS260 |         0.011 |            0.012 |                0.987 |
| ENS107 |         0.007 |            0.022 |                0.987 |
| ENS291 |         0.016 |            0.024 |                0.987 |
| ENS219 |         0.012 |            0.026 |                0.987 |
| ENS458 |         0.014 |            0.025 |                0.987 |

The top 10 genes associated with predictor x sorted by lowest P-value

The old method took more time(`time.old.method`) compared with the MI
PCA method (`time.pca`).

# Session info

    ## R version 4.3.3 (2024-02-29 ucrt)
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
    ## [1] mice_3.16.0             limma_3.56.2            PCAtools_2.12.0        
    ## [4] ggrepel_0.9.3           ggplot2_3.4.4           BiocParallel_1.34.2    
    ## [7] dplyr_1.1.0             RNAseqCovarImpute_1.1.2
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] tidyselect_1.2.0          fastmap_1.1.1            
    ##  [3] digest_0.6.31             rsvd_1.0.5               
    ##  [5] rpart_4.1.23              lifecycle_1.0.3          
    ##  [7] survival_3.5-8            magrittr_2.0.3           
    ##  [9] compiler_4.3.3            rlang_1.1.1              
    ## [11] tools_4.3.3               utf8_1.2.3               
    ## [13] yaml_2.3.7                knitr_1.43               
    ## [15] S4Arrays_1.0.6            dqrng_0.3.0              
    ## [17] DelayedArray_0.26.7       plyr_1.8.8               
    ## [19] abind_1.4-5               withr_2.5.0              
    ## [21] purrr_1.0.1               BiocGenerics_0.46.0      
    ## [23] nnet_7.3-19               grid_4.3.3               
    ## [25] stats4_4.3.3              fansi_1.0.4              
    ## [27] jomo_2.7-6                beachmat_2.16.0          
    ## [29] colorspace_2.1-0          edgeR_3.42.4             
    ## [31] scales_1.2.1              iterators_1.0.14         
    ## [33] MASS_7.3-60.0.1           cli_3.6.0                
    ## [35] rmarkdown_2.20            crayon_1.5.2             
    ## [37] generics_0.1.3            rstudioapi_0.14          
    ## [39] reshape2_1.4.4            DelayedMatrixStats_1.22.6
    ## [41] minqa_1.2.5               stringr_1.5.0            
    ## [43] splines_4.3.3             parallel_4.3.3           
    ## [45] matrixStats_1.0.0         vctrs_0.5.2              
    ## [47] boot_1.3-29               glmnet_4.1-7             
    ## [49] Matrix_1.6-5              BiocSingular_1.16.0      
    ## [51] IRanges_2.33.0            S4Vectors_0.37.4         
    ## [53] mitml_0.4-5               irlba_2.3.5.1            
    ## [55] locfit_1.5-9.7            foreach_1.5.2            
    ## [57] tidyr_1.3.0               snow_0.4-4               
    ## [59] glue_1.6.2                nloptr_2.0.3             
    ## [61] pan_1.8                   codetools_0.2-19         
    ## [63] cowplot_1.1.1             stringi_1.7.12           
    ## [65] shape_1.4.6               gtable_0.3.1             
    ## [67] ScaledMatrix_1.8.1        lme4_1.1-34              
    ## [69] munsell_0.5.0             tibble_3.2.0             
    ## [71] pillar_1.9.0              htmltools_0.5.4          
    ## [73] R6_2.5.1                  sparseMatrixStats_1.12.2 
    ## [75] evaluate_0.21             lattice_0.22-5           
    ## [77] backports_1.4.1           broom_1.0.5              
    ## [79] Rcpp_1.0.10               nlme_3.1-164             
    ## [81] xfun_0.39                 MatrixGenerics_1.12.3    
    ## [83] pkgconfig_2.0.3
