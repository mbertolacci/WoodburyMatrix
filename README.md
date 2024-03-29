# WoodburyMatrix

[![CRAN
status](https://www.r-pkg.org/badges/version/WoodburyMatrix)](https://cran.r-project.org/package=WoodburyMatrix)
[![R-CMD-check](https://github.com/mbertolacci/WoodburyMatrix/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/mbertolacci/WoodburyMatrix/actions/workflows/R-CMD-check.yaml)
[![Coverage status](https://codecov.io/gh/mbertolacci/WoodburyMatrix/branch/master/graph/badge.svg)](https://codecov.io/github/mbertolacci/WoodburyMatrix?branch=master)

## Overview

WoodburyMatrix is an R package that provides a hierarchy of classes and methods for manipulating matrices formed implicitly from the sums of the inverses of other matrices, a situation commonly encountered in spatial statistics and related fields. It enables easy use of the Woodbury matrix identity and the matrix determinant lemma to allow computation (e.g., solving linear systems) without having to form the actual matrix.

## Installation

The package can be installed from [CRAN](https://cran.r-project.org/package=WoodburyMatrix) using

```{r}
install.packages('WoodburyMatrix')
```

Alternatively, you can get the latest development version from GitHub with

```{r}
devtools::install_github('mbertolacci/WoodburyMatrix')
```

## Using the package

An overview of the related mathematics and instructions on using the package are provided on [the package website](https://mbertolacci.github.io/WoodburyMatrix/articles/WoodburyMatrix.html).

## Getting help

If you find a bug or have suggestions for features, please [file an issue on GitHub](https://github.com/mbertolacci/WoodburyMatrix/issues). Otherwise, if you need help, feel free to drop the [author an email](mailto:m.bertolacci@gmail.com).
