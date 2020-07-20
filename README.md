# WoodburyMatrix

[![Travis build status](https://travis-ci.org/mbertolacci/WoodburyMatrix.svg?branch=master)](https://travis-ci.org/mbertolacci/WoodburyMatrix) [![Coverage status](https://codecov.io/gh/mbertolacci/WoodburyMatrix/branch/master/graph/badge.svg)](https://codecov.io/github/mbertolacci/WoodburyMatrix?branch=master)

## Overview

WoodburyMatrix is an R package that provides a hierarchy of classes and methods for manipulating matrices formed implicitly from the sums of the inverses of other matrices, a situation commonly encountered in spatial statistics and related fields. It enables easy use of the Woodbury matrix identity and the matrix determinant lemma to allow computation (e.g., solving linear systems) without having to form the actual matrix.

## Installation

The package is currently available only on GitHub, and can be installed with:

```{r}
devtools::install_github('mbertolacci/WoodburyMatrix')
```

Submission to CRAN is anticipated shortly.

## Using the package

An overview of the related mathematics and instructions on using the package are provided on [the package website](https://mbertolacci.github.io/WoodburyMatrix/articles/WoodburyMatrix.html).

## Getting help

If you find a bug or have suggestions for features, please [file an issue on GitHub](https://github.com/mbertolacci/WoodburyMatrix/issues). Otherwise, if you need help, feel free to drop the [author an email](mailto:m.bertolacci@gmail.com).
