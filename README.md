---
title: "Multi-trait Analysis of Gene-by-Environment Interactions Using Summary Statistics of Genetic Marginal Effects"
author: Lan Luo
output: rmarkdown::html_vignette
bibliography: reference.bib
vignette: >
  %\VignetteIndexEntry{MAGENTA}
  %\VignetteEngine{knitr::rmarkdown}
  \usepackage[utf8]{inputenc}
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
```

```{r load, eval = FALSE, echo=FALSE}
library(MAGENTA)
```

## Contents

- Overview
- System Requirements
- Installation Guide
- Demo
- URLs and References


## Overview 
Studying genotype-by-environment interaction (GEI) is fundamental in understanding complex trait variations. Identifying genetic variants with GEI effects is challenging because the GEI analysis generally has low power. MAGENTA (Multi-trait Analysis of Gene-ENvironmenT-wide Association) is a powerful, robust, and computationally efficient method to test the interaction between a gene and environmental groups on multiple traits in large-scale datasets, such as the UK Biobank. More details about MAGENTA can be found in @luo2022magenta.

Specifically, `MAGENTA` package has functions to 

- compute the summary statistics with different types of data input adjusting for the potential overlapping samples under different assumption of genetic marginal effects; and 
- perform summary-statistics-based multi-trait analysis of gene-environment interaction (GEI) tests or genetic main effect and GEI joint effect for both common and rare variants.

```{r, out.width = "650px", echo=FALSE, fig.cap="Fig 1: An overview of MAGENTA workflow. Light blue rectangle represents necessary input. Dark blue rectangle denotes the final output of MAGENTA function. Gray rectangle denotes the intermediate parameters."}
knitr::include_graphics("workflow.png")
```

## System Requirements

The package development version is tested on the following systems:

Mac OSX: Mojave version 10.14.6 (R version 3.6.0)  

Windows 10 (R version 3.6.1)

The CRAN package should be compatible with Windows and Mac operating systems.

## Installing Guide
`MAGENTA` package requires R with version 3.6.1 or higher, which can be downloaded and installed from [CRAN](https://cran.r-project.org/). 

```
install.packages("MAGENTA")
```

### Package dependencies

`MAGENTA` package depends on several R packages, which will be downloaded before installing `MAGENTA`. `MAGENTA` also uses non-exported R code from R packages `ACAT` and `SKAT`. The `MAGENTA` package functions with all packages in their  versions as they appear on `CRAN` or `github` on January 23, 2020 and October 28, 2021, respectively. The versions of software are, specifically:
```
MASS (>= 7.3-51.4),
Matrix (>= 1.2-17),
caret (>= 6.0-84),
stats,
utils,
gtools (>= 3.8.1),
SPAtest (>= 3.0.0),
survival (>= 3.2-3),
SimCorMultRes (>= 1.7.0),
SKAT (>= 1.3.2.1),
expm (>= 0.999-4),
CompQuadForm (>= 1.4.3),
caret (>=6.0-84)
```

## Demo

### Step 1: Prepare Summary Statistics across environmental groups
Suppose we are interested in testing the GEI effect on $K$ traits. Under each environmental group, we estimate the genetic marginal effects of $M$ variants on the $K$ traits and the covariances of these effect estimates. Let $\mathbf{\beta}^{(d)}_k$ denote the effects of the $M$ variants on trait $k$ in environmental group $d$ and $\widehat{\mathbf{\beta}}^{(d)}_k$ denote the corresponding effect estimates. We let $\text{cov}(\widehat{\mathbf{\beta}}_k^{(d)}, \widehat{\mathbf{\beta}}_{k'}^{(d)})$ denote the covariance of  $\widehat{\mathbf{\beta}}_k^{(d)}$ and $\widehat{\mathbf{\beta}}_{k'}^{(d)}$. For the same trait $k=k'$, the $\text{cov}(\widehat{\mathbf{\beta}}_k^{(d)}, \widehat{\mathbf{\beta}}_{k'}^{(d)})$ is the variance-covariance of $\widehat{\mathbf{\beta}}^{(d)}_k$. For different traits $k\neq k'$, $\text{cov}(\widehat{\mathbf{\beta}}_k^{(d)}, \widehat{\mathbf{\beta}}_{k'}^{(d)})$ quantifies the covariance between effect estimates on the two traits when they are measured on the overlapping set of subjects and $\text{cov}(\widehat{\mathbf{\beta}}_k^{(d)}, \widehat{\mathbf{\beta}}_{k'}^{(d)})=0$ if there are no overlapping subjects. For the simplicity of notation, we drop the environmental group index $(d)$ in the remaining part of this section

The variant-level score statistics $\mathbf{U}_k$ for testing the genetic effect $\mathbf{\beta}_k = 0$ and its variance-covariance estimate $\mathbf{V}_k$ are routinely used to construct the marginal effect estimates $\widehat{\mathbf{\beta}}_k = \mathbf{V}^{-1}_k\mathbf{U}_k$.

MAGENTA introduces below two approaches of estimating $\text{cov}(\widehat{\mathbf{\beta}}_k, \widehat{\mathbf{\beta}}_{k'})$. One approach needs individual-level data of all traits but can accurately estimate the correlation induced by the overlapping subjects across trait. The other approach utilizes only the trait-specific summary statistics $\mathbf{U}_k$ and $\mathbf{V}_k$ and is accurate when there are no genetic main effects, which is true under the null hypothesis of main test or joint test, but not necessarily the case under the null hypothesis of GEI test. This apparoch can still provide a reasonable approximation when the main effects are weak but may be inaccurate when the genetic main effects are unusually strong

#### 1. Input data is individual-level data
Assume the users provide the individual-level data $\mathbf{Y}$ for traits, $\mathbf{X}$ for covariates, $\mathbf{G}$ for genotype, respectively. Let $\widehat{\mathbf{\theta}}_k = (\widehat{\mathbf{\alpha}}^{\rm T}_k, \widehat{\mathbf{\beta}}^{\rm T}_k)^{\rm T}$ and $\mathbf{Z}_{ik} = (\mathbf{X}^{\rm T}_{ik}, \mathbf{G}^{\rm T}_{ik})^{\rm T}$ with $\widehat{\mathbf{\alpha}}^{\rm T}_k, \widehat{\mathbf{\beta}}^{\rm T}_k$ denoting the estimated regression coefficients for covariates and genotype. Without loss of generality, we organize the data such that the first $n_{kk'}$ subjects are those with measurements on both traits $k$ and $k'$. Then $\text{cov}(\widehat{\mathbf{\theta}}_{k}, \widehat{\mathbf{\theta}}_{k'})$ can be estimated as
\[
\text{cov}(\widehat{\mathbf{\theta}}_{k}, \widehat{\mathbf{\theta}}_{k'}) = \left(\sum_{i = 1}^{n_k} \widehat{b}^{''}_{ik}\mathbf{Z}_{ik}\mathbf{Z}^{\rm T}_{ik}\right)^{-1}
\left(\sum_{i = 1}^{n_{kk'}} \mathbf{Z}_{ik} \widehat{v}^{1/2}_{ik} \widehat{\omega}_{kk'} \widehat{v}^{1/2}_{ik'} \mathbf{Z}^{\rm T}_{ik'}\right) 
\left(\sum_{i = 1}^{n_{k'}} \widehat{b}^{''}_{ik'}\mathbf{Z}_{ik'}\mathbf{Z}^{\rm T}_{ik'}\right)^{-1},
\]

where $\widehat{b}^{''}_{ik}=1$ and $\widehat{v}_{ik} = \widehat{\sigma}^2_k$  for the linear regression model with residual variance estimate $\widehat{\sigma}^2_k$, and $\widehat{b}^{''}_{ik}=\widehat{v}_{ik}=e^{\widehat{\mathbf{\alpha}}^{\rm T}_k\mathbf{X}_{ik}}/(1+e^{\widehat{\mathbf{\alpha}}^{\rm T}_k\mathbf{X}_{ik}})^2$ for the logistic regression model, 
and $\widehat{\omega}_{kk'}$ is the empirical correlation of the Pearson residuals of the two traits using their overlapping samples. 
The $\widehat{b}^{''}_{ik}$, $\widehat{b}^{''}_{ik'}$, $\widehat{v}_{ik}$, $\widehat{v}_{ik'}$ and $\widehat{\omega}_{kk'}$ need to compute only once over all the genes in the genome-wide analysis.
We then obtain the $\text{cov}(\widehat{\mathbf{\beta}}_{k}, \widehat{\mathbf{\beta}}_{k'})$ as the submatrix of $\text{cov}(\widehat{\mathbf{\theta}}_{k}, \widehat{\mathbf{\theta}}_{k'})$ corresponding to $\widehat{\mathbf{\beta}}_k$ on the rows and $\widehat{\mathbf{\beta}}_{k'}$ on the columns.

Here we use a simulated dataset rawdata to show how to apply *Get_beta_cov_data* function. Assume there are 3 studies, 3 continuous traits and 86 variants. We assume a completely overlapping pattern and there are 1000 subjects with 3 traits measured. The object *dat* and *geno.dat* should have same format of subject ID as the first column of the dataset. 
```{r loaddata, eval = FALSE, echo = TRUE}
data("rawdata")
attach(rawdata)
head(traits.dat$Study1)
head(dat)
head(geno.dat)
```

The function *Get_beta_cov_data* additionally requires user to specify the names of covariates, traits, SNPs and treatment or environmental variable. Currently, it only allows $K$ same type of traits.

```{r get_beta_cov_data, eval = FALSE, echo = TRUE}
## start from the individual-level data
K <- 3
beta.sumstats <- Get_beta_cov_data(geno.dat = geno.dat, dat = dat, 
                                   cov.list = c("cov1", "cov2"),
                                   env = "treatment", test = "GEI",
                                   ID = "ID", type = "continuous",
                                   trait.list = paste0("Trait", 1:K),
                                   SNP.list = colnames(geno.dat)[-1])
```

The returning object contains the marginal genetic effect estimates $\mathbf{\beta}_k^{(d)}$ as well as the variance-covariance estimate $\text{cov}(\widehat{\mathbf{\beta}}_k, \widehat{\mathbf{\beta}}_{k'})$ for $k$th trait in $d$th environmental group.
```{r betasumstats, eval = FALSE, echo = TRUE}
names(beta.sumstats)
```

#### 2. Input data is trait-specific summary statistics
When there are no genetic main effects or the effects are weak,
we can estimate $\text{cov}(\widehat{\mathbf{\beta}}_k, \widehat{\mathbf{\beta}}_{k'})$ only using the trait-specific summary statistics:

\[
  \text{cov}(\widehat{\mathbf{\beta}}_{k}, \widehat{\mathbf{\beta}}_{k'}) = \widehat\zeta_{kk'}\mathbf{V}_{k}^{-1}{\rm diag}(\mathbf{V}_{k})^{1/2}\mathbf{R}{\rm diag}(\mathbf{V}_{k'})^{1/2}\mathbf{V}_{k'}^{-1} %~~~~~~~~~~~~~~~~~~~~~ (1)
\]

When $k=k'$, $\widehat{\zeta}_{kk'}=1$; When $k\neq k'$,  $\widehat{\zeta}_{kk'}$ is estimated using the sample correlation of Z-scores between traits $k$ and $k'$ over a large number of independent common variants not associated with the traits (@luo2020multi). The $\mathbf{R}$ is the correlation matrix among variants calculated using the Pearson correlation coefficient among the genotypes of the variants, which can be from an external reference genotype data. 

Here we show an example of using *Get_beta_cov_UV* function. Given the score summary statistics for $k$th trait in $d$th environmental group, we first need to estimate the covariance matrix of Z-scores between traits *zeta* in each environmental group, which can be calculated using *Get_zeta* function with input of score summary statistics ($U_{kj}^{(d)}$ and its variance $V_{kj}^{(d)}$) of a large number (>=500) common independent variants that are not associated with the traits.

```{r get_zeta, eval = FALSE, echo = TRUE}
data(sumstats.dat)
names(sumstats.dat)
attach(sumstats.dat)
str(common.sumstats)
zeta <- Get_zeta(common.sumstats = common.sumstats, trait.list = paste("Trait", 1:K))
zeta
```
In the returning object *zeta*, it contains three sublists *zeta.main*, *zeta.joint* and *zeta.GEI*, which captures the correlation between traits for genetic main effect test, joint test and GEI test respectively. In the following analysis, *MAGENTA* function will use the corresponding *zeta* matrix for different type of tests.

Given the variant-level score summary statistics and *zeta* matrix, MAGENTA will construct the marginal genetic effect estimates and its variance-covariance matrix. The users need to additionally provide LD correlation matrix among variants *R.C*, the minor allele frequency *MAF*, as well as the genetic correlation matrix *KA*. If users don't provide such a correlation matrix among variants, MAGENTA will try to recover $\mathbf{R}$ from the variance-covariance matrix $V_k^{(d)}, k = 1, \dots, K, d = 1, \dots, D$. 

If the genetic main effect p value is less than the pre-specified threshold (default is $5\times 10^{-3}$), it will prompt a warning that the estimated GEI p value may be biased. 

```{r get_beta_cov_UV, eval = FALSE, echo = TRUE}
str(sumstats)
beta.sumstats <- Get_beta_cov_UV(sumstats = sumstats, MAF = MAF, R.C = R.C, 
                                 zeta.ret = zeta.ret, KA = KA, 
                                 test = "joint",
                                 trait.list = paste0("Trait", 1:K))
names(beta.sumstats)
```
#### 3. Saddle point approximation for binary traits
For binary traits, the normal approximation for the summary statistics can be inaccurate even for variants with MAC $>$ 10 (@dey2017fast,@bi2020fast). We apply a saddlepoint approximation (SPA) and use the SPA-adjusted summary statistics to improve the accuracy of the variant-level p values.

1. Use R package SPAtest to obtain the SPA p value $P_{\rm SPA}$ for testing the association between the variant and the trait. \cite{dey2017fast}
2. Convert $P_{\rm SPA}$ to Z-score by using $Z_{\rm SPA} = {\rm sign}(U) \Phi^{-1}(P_{\rm SPA}/2)$, where  $U$ is the score statistic for the single-variant test and $\Phi^{-1}(\cdot)$ is the quantile function of the standard normal distribution. 
3. The SPA-adjusted variance of the score statistic is calculated as $V_{\rm SPA}=(U/Z_{\rm SPA})^2$.

MAGENTA will authomatically apply SPA approach once the input traits are binary.

### Step 2. Perform joint and GEI tests for common and rare variants (MAGENTA)
To facilitate the meta-analysis of GEI in large consortium-based studies, MAGENTA translates testing for GEI to testing for heterogeneity of genetic marginal effects across environmental groups and generate summary statistics for marginal effects of multiple variants in the gene on multiple traits of interest. MAGENTA combines these summary statistics and produces a GEI test and a joint test of genetic main and GEI effects. 

To run MAGENTA function, the summary statistics *beta*, their covariance matrix *beta.cov* for each trait in each environmental group, and minor allele frequency *MAF* are required. Other inputs are optional. Here an example dataset (*MAGENTA.example*) offers all the information needed in *MAGENTA* function.

```{r MAGENTA, eval = FALSE}
data("MAGENTA.example")
names(MAGENTA.example)
```

Specifically, *genetic_cor.trait* represents the genetic correlation information (@Bulik2015). The genetic correlation can be estimated using *ldsc* software (@bulik2015ld and URLs). Some websites also archive the genetic correlation among many complex traits (URLs). If no genetic correlation is available, *MAGENTA* will use an exchangeable correlation structure among traits. In the default setting, MAGENTA will return both joint and GEI p values. Users can also choose to print out these intermediate joint/GEI p values for multi-/single-traits $\times$ multi-/single-variant analysis.

```{r cMAGENTA, warning=FALSE, message=FALSE, eval = FALSE}
attach(MAGENTA.example)
pval <-  MAGENTA(U = U, V = V, MAF = MAF, genetic_cor.trait = genetic_cor.trait, zeta = zeta)
pval
detach(MAGENTA.example)
```
#### MAGENTA joint effect test

#### MAGENTA GEI effect test

## URLs
ldsc website: https://github.com/bulik/ldsc/wiki/Heritability-and-Genetic-Correlation
Genetic correlation: https://media.nature.com/original/nature-assets/ng/journal/v47/n11/extref/ng.3406-S2.csv

## References
