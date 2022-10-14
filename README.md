Contents
--------

-   Overview
-   System Requirements
-   Installation Guide
-   Demo
-   URLs and References

Overview
--------

Studying genotype-by-environment interaction (GEI) is fundamental in
understanding complex trait variations. Identifying genetic variants
with GEI effects is challenging because the GEI analysis generally has
low power. MAGENTA (Multi-trait Analysis of Gene-ENvironmenT-wide
Association) is a powerful, robust, and computationally efficient method
to test the interaction between a gene and environmental groups on
multiple traits in large-scale datasets, such as the UK Biobank. More
details about MAGENTA can be found in Luo et al (2022).

Specifically, `MAGENTA` package has functions to

-   compute the summary statistics with different types of data input
    adjusting for the potential overlapping samples under different
    assumption of genetic marginal effects; and
-   perform summary-statistics-based multi-trait analysis of
    gene-environment interaction (GEI) tests or genetic main effect and
    GEI joint effect for both common and rare variants.

![Fig 1: An overview of MAGENTA workflow. Light blue rectangle
represents necessary input. Dark blue rectangle denotes the final output
of MAGENTA function. Gray rectangle denotes the intermediate
parameters.](%22/Users/luolan2/Project/MAGENTA/package/MAGENTA/vignettes/workflow.png%22)

System Requirements
-------------------

The package development version is tested on the following systems:

Mac OSX: Mojave version 10.14.6 (R version 3.6.0)

Windows 10 (R version 3.6.1)

The CRAN package should be compatible with Windows and Mac operating
systems.

Installing Guide
----------------

`MAGENTA` package requires R with version 3.6.1 or higher, which can be
downloaded and installed from [CRAN](https://cran.r-project.org/).

    install.packages("MAGENTA")

### Package dependencies

`MAGENTA` package depends on several R packages, which will be
downloaded before installing `MAGENTA`. `MAGENTA` also uses non-exported
R code from R packages `ACAT` and `SKAT`. The `MAGENTA` package
functions with all packages in their versions as they appear on `CRAN`
or `github` on January 23, 2020 and October 28, 2021, respectively. The
versions of software are, specifically:

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

Demo
----

Detailed instructions can be found in Vignetts/MAGENTA.html file.

References
----------

Luo, Lan, Devan V Mehrotra, Judong Shen, and Zheng-Zheng Tang. 2022.
“Multi-Trait Analysis of Gene-by-Environment Interactions Using Summary
Statistics of Genetic Marginal Effects.” Submitted.