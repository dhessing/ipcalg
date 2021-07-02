# Title     : Models for simulation study
# Created by: Dzjon Hessing
# Created on: 27/06/2021

library(tibble)

# No moderation
no_modification <- list(
  sample = function(n) {
    tibble(
      A = rnorm(n),
      Q = rnorm(n),
      Y = A + Q + rnorm(n)
    )
  },
  DAG = c(0, 0, 0,
          0, 0, 0,
          1, 1, 0),
  IDAG = c(0, 0,
           0, 0)
)

# Direct effect modification
effect_modification <- list(
  sample = function(n) {
    tibble(
      A = rnorm(n),
      Q = rnorm(n),
      Y = A + Q + A * Q + rnorm(n)
    )
  },
  DAG = c(0, 0, 0,
          0, 0, 0,
          1, 1, 0),
  IDAG = c(0, 0,
           1, 0)
)

total_effect_modification <- list(
  sample = function(n) {
    data <- tibble(
      A = rnorm(n),
      Q = rnorm(n),
      X = A + Q + A*Q + rnorm(n),
      C = rnorm(n),
      Y = X + C + rnorm(n),
    )
  },
  DAG = c(0, 0, 0, 0, 0,
          0, 0, 0, 0, 0,
          1, 1, 0, 0, 0,
          0, 0, 0, 0, 0,
          0, 0, 1, 1, 0),
  IDAG = c(0, 0, 0, 0,
           1, 0, 0, 0,
           0, 0, 0, 0,
           1, 0, 0, 0)
)

# Effect modification by common cause
# Not used as of now
common_cause <- list(
  sample = function(n) {
    data <- tibble(
      A = rnorm(n),
      C = rnorm(n),
      X = C + rnorm(n),
      Q = C + rnorm(n),
      Y = A + Q + A * Q + rnorm(n),
    )
  },
  DAG = c(0, 0, 0, 0, 0,
          0, 0, 0, 0, 0,
          0, 1, 0, 0, 0,
          0, 1, 0, 0, 0,
          1, 0, 0, 1, 0),
  IDAG = c(0, 0, 0, 0,
           1, 0, 0, 0,
           1, 0, 0, 0,
           0, 0, 1, 0)
)

models <- list(
  'no effect modification' = no_modification,
  'direct effect modification' = effect_modification,
  'total effect modification' = total_effect_modification
)