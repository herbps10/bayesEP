# bayesEP

<!-- badges: start -->
[![R-CMD-check](https://github.com/herbps10/bayesEP/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/herbps10/bayesEP/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

## Overview
This package implements a flexible framework for Bayesian Expectation Propagation (EP) based on the algorithm 
described in [Vehtari et al. 2020](https://jmlr.org/papers/v21/18-817.html).

## Installation

You can install the development version of `bayesEP` from [Github](https://github.com/herbps10/bayesEP) with:

```r
# install.packages("pak")
pak::pak("herbps10/bayesEP")
```

## References
- Vehtari, A., Gelman, A., Sivula, T., Jylänki, P., Tran, D., Sahai, S., Blomstedt, P., Cunningham, J. P., Schiminovich, D., and Robert, C. P. (2020).
  Expectation Propagation as a Way of Life: A Framework for Bayesian Inference
  on Partitioned Data. Journal of Machine Learning Research, 21(17), 1–53. [http://jmlr.org/papers/v21/18-817.html](http://jmlr.org/papers/v21/18-817.html)
