# Title     : Simulation study for IPC function
# Objective : Compare performance of IPC algoithm to PC algorithm
# Created by: Dzjon Hessing
# Created on: 27/06/2021

library(tidyverse)
source('R/simulation_study/models.R')

#' One iteration of estimates for the PC algorithm.
simulation.pc <- function(data, DAG) {
  # Runs PC
  pc.fit <- pcalg::pc(suffStat = list(C = cor(data),
                                      n = nrow(data)),
                      indepTest = gaussCItest,
                      alpha = 0.01,
                      labels = colnames(data))

  # Fail if PC generates an invalid CPDAG
  amat <- wgtMatrix(pc.fit@graph)
  amat[which(amat != 0)] <- 1
  if (!isValidGraph(amat, 'cpdag'))
    return()

  # Runs ida
  ida.results <- ida(1, ncol(data), cov(data), pc.fit@graph, method = 'global')

  # Gets all dags from the cpdag
  ad <- pdag2allDags(as(pc.fit, 'amat'))$dags

  if (nrow(ad) == 1) {
    # If there's just one: picks that one
    ida.result <- ida.results
    idx <- 1
  } else {
    # Multiple: Find if the correct DAG in in the CPDAG
    idx <- which(apply(ad, 1, function(row) { all(row == DAG) }))

    if (length(idx) == 0) {
      # If it is not present: fail
      return()
    }

    # Else it picks that one
    ida.result <- ida.results[idx]
  }

  # Returns estimates for low, average and high, but these are the same for PC.
  list(
    estimates = tibble_row(
      low = ida.result,
      average = ida.result,
      high = ida.result
    ),
    fit = pc.fit,
    idx = idx
  )
}

#' One iteration of getting estimates for the IPC algorithm
simulation.ipc <- function(data, idx, IDAG, pc.fit) {
  # Runs ipc
  ipc.fit <- ipcalg::ipc(pc.fit, data, x = 'A', y = 'Y', 0.01)

  # Do we recover the true IDAG?
  recovered <- FALSE
  if (identical(ipc.fit$idag[idx,], IDAG)) {
    recovered <- TRUE
  }

  # Gets estimates for low, average and high values of Q
  cols.mean <- setdiff(colnames(data), c('A', 'Y', 'Q'))
  means <- as.list(data[cols.mean] %>% summarise_all(mean))
  lang <- str2lang(ipc.fit$equations[idx])
  list(
    estimates = tibble_row(
      low = eval(substituteDirect(lang, c(list(A = 1, Q = -1), means))) -
        eval(substituteDirect(lang, c(list(A = 0, Q = -1), means))),
      average = eval(substituteDirect(lang, c(list(A = 1, Q = 0), means))) -
        eval(substituteDirect(lang, c(list(A = 0, Q = 0), means))),
      high = eval(substituteDirect(lang, c(list(A = 1, Q = 1), means))) -
        eval(substituteDirect(lang, c(list(A = 0, Q = 1), means)))),
    recovered = recovered
  )
}

iterate <- function(model, n.sample) {
  DAG.incorrect <- 0
  IDAG.incorrect <- 0
  while (TRUE) {
    data <- model$sample(n.sample)
    pc.result <- simulation.pc(data, model$DAG)
    if (is.null(pc.result)) {
      DAG.incorrect <- DAG.incorrect + 1
    } else {
      ipc.result <- simulation.ipc(data, pc.result$idx, model$IDAG, pc.result$fit)
      if (!ipc.result$recovered) {
        IDAG.incorrect <- IDAG.incorrect + 1
      }
      return(tibble_row(!!!c(
        pc = pc.result$estimates,
        ipc = ipc.result$estimates,
        DAG.incorrect = DAG.incorrect,
        IDAG.incorrect = IDAG.incorrect)))
    }
  }
}



set.seed(132)
cat('\nRunning simulation...')
tidyr::expand_grid(
  model = c('no effect modification', 'direct effect modification', 'total effect modification'),
  sample.size = c(200, 500, 1000, 2000),
  iteration = 1:1000
) %>%
  rowwise() %>%
  mutate(iterate(models[[model]], sample.size)) %>%
  write_rds('data/simulation_study.RDS')

cat('\nWrote simulation results to data/simulation_study.RDS')
