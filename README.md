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
low power. MTAGEI (Multi-Trait Analysis of Gene-Environment Interactions) is a powerful, robust, and computationally efficient method
to test the interaction between a gene and environmental groups on
multiple traits in large-scale datasets, such as the UK Biobank. More
details about MTAGEI can be found in Luo et al (2023).

Specifically, `MTAGEI` package has functions to

-   compute the summary statistics with different types of data input
    adjusting for the potential overlapping samples under different
    assumption of genetic marginal effects; and
-   perform summary-statistics-based multi-trait analysis of
    gene-environment interaction (GEI) tests or genetic main effect and
    GEI joint effect for both common and rare variants.

![](%22workflow.png%22) <img
  src="workflow.png"
  alt="Fig 1"
  title="An overview of MTAGEI workflow. Light blue rectangle represents necessary input. Dark blue rectangle denotes the final output of MTAGEI function. Gray rectangle denotes the intermediate parameters."
  style="display: inline-block; margin: 0 auto; max-width: 300px">

System Requirements
-------------------

The package development version is tested on the following systems:

Mac OSX: Mojave version 10.14.6 (R version 3.6.0)

Windows 10 (R version 3.6.1)

The CRAN package should be compatible with Windows and Mac operating
systems.

Installing Guide
----------------

`MTAGEI` package requires R with version 3.6.1 or higher, which can be
downloaded and installed from [Github](https://github.com/lan/MTAGEI).

    install.packages("MTAGEI.tar.gz", repos = NULL, type = "source")

### Package dependencies

`MTAGEI` package depends on several R packages, which will be
downloaded before installing `MTAGEI`. `MTAGEI` also uses non-exported
R code from R packages `ACAT` and `SKAT`. The `MTAGEI` package
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

Detailed instructions can be found in Vignetts/MTAGEI.html file. [Preview](https://htmlpreview.github.io/?)

URLs
----

ldsc website:
<a href="https://github.com/bulik/ldsc/wiki/Heritability-and-Genetic-Correlation" class="uri">https://github.com/bulik/ldsc/wiki/Heritability-and-Genetic-Correlation</a>

Genetic correlation:
<a href="https://media.nature.com/original/nature-assets/ng/journal/v47/n11/extref/ng.3406-S2.csv" class="uri">https://media.nature.com/original/nature-assets/ng/journal/v47/n11/extref/ng.3406-S2.csv</a>

References
----------

Bulik-Sullivan, Brendan, Hilary K Finucane, Verneri Anttila, Alexander
Gusev, Felix R Day, John RB Perry, Nick Patterson, et al. 2015. “An
Atlas of Genetic Correlations Across Human Diseases and Traits.” Nat.
Genet. 47 (11): 1236–41.

Bulik-Sullivan, Brendan K, Po-Ru Loh, Hilary K Finucane, Stephan Ripke,
Jian Yang, Nick Patterson, Mark J Daly, et al. 2015. “LD Score
Regression Distinguishes Confounding from Polygenicity in Genome-Wide
Association Studies.” Nature Genetics 47 (3): 291.

Luo, Lan, Devan V Mehrotra, Judong Shen, and Zheng-Zheng Tang. 2023.
“Multi-Trait Analysis of Gene-by-Environment Interactions in Large-scale Genetic Studies.” Submitted to Biostatistics.

Luo, Lan, Judong Shen, Hong Zhang, Aparna Chhibber, Devan V Mehrotra,
and Zheng-Zheng Tang. 2020. “Multi-Trait Analysis of Rare-Variant
Association Summary Statistics Using MTAR.” Nat. Commun. 11 (1): 1–11.
