
<!-- README.md is generated from README.Rmd. Please edit that file -->

# ipcalg

<!-- badges: start -->

<!-- badges: end -->

The goal of ipcalg is to estimate an DAG/IDAG pair for a dataset with
unknown SCM.

## Installation

You can install development of ipcalg from
[GitHub](https://github.com/dhessing/ipcalg) with:

``` r
devtools::install_github("dhessing/ipcalg")
```

## Example

``` r
library(tibble)
library(pcalg)
library(ipcalg)
n <- 1000
data <- tibble::tibble(
  A = rnorm(n),
  B = rnorm(n),
  Y = A + B + A * B + rnorm(n)
)
pc.fit <- pc(
  suffStat = list(C = cor(data), n = nrow(data)),
  indepTest = gaussCItest,
  alpha = 0.01,
  labels = colnames(data)
)
ipc.fit <- ipc(pc.fit, data, x = 'A', y = 'Y', alpha = 0.01)
```
