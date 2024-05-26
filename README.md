# SpaCOAP
High-Dimensional Spatial Covariate-Augmented Overdispersed Poisson Factor Model

=========================================================================
<!-- badges: start -->

[![](https://www.r-pkg.org/badges/version-ago/SpaCOAP)](https://cran.r-project.org/package=SpaCOAP)
[![](https://cranlogs.r-pkg.org/badges/SpaCOAP?color=orange)](https://cran.r-project.org/package=SpaCOAP)
[![](https://cranlogs.r-pkg.org/badges/grand-total/SpaCOAP?color=orange)](https://cran.r-project.org/package=SpaCOAP)
<!-- badges: end -->


We introduces an efficient latent representation learning approach tailored specifically for high-dimensional, large-scale spatial count data, incorporating additional covariates for enhanced performance.
 To model correlations among variables measured at a shared spatial location, we introduce a covariate-augmented overdispersed Poisson factor model. We distinguish between high-dimensional covariates sharing similar attributes and those serving as control variables to enrich the representation learning process. To capture the spatial dependency of each variable across different locations, we apply a conditional autoregressive model to the latent factors. Furthermore, we propose a variational expectation-maximization algorithm to estimate the model parameters and latent factors, imposing a low-rank constraint on the high-dimensional regression coefficient matrix.



Check out  [Package Website](https://feiyoung.github.io/SpaCOAP/index.html) for a more complete description of the methods and analyses. 

# Installation
"SpaCOAP" depends on the 'Rcpp' and 'RcppArmadillo' package, which requires appropriate setup of computer. For the users that have set up system properly for compiling C++ files, the following installation command will work.
```{Rmd}
## Method 1:

if (!require("remotes", quietly = TRUE))
    install.packages("remotes")
remotes::install_github("feiyoung/SpaCOAP")

## Method 2: install from CRAN
install.packages("SpaCOAP")

```



## Usage
For usage examples and guided walkthroughs, check the `vignettes` directory of the repo. 

* [Simulated data](https://feiyoung.github.io/SpaCOAP/articles/simu.html)

* [Real data](https://feiyoung.github.io/SpaCOAP/articles/simu.html)

## Simulated codes
For the codes in simulation study, check the `simu_code` directory of the repo.


## News

SpaCOAP version 1.2 released! (2024-05-25) 


