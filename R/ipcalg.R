library(pcalg)
library(Deriv)
library(igraph)
library(BiocGenerics)
library(graph)
library(Rgraphviz)

#' Estimate an IDAG
#'
#' Estimate an IDAG for the causal effect of x on y from a fitted PC function.
#'
#' @param pc.fit fitted PC object
#' @param data dataset
#' @param x intervention
#' @param y outcome
#' @param alpha significance level for regressions
#'
#' @return fitted IPC object
#' @export
#' @examples
#' library(tibble)
#' library(pcalg)
#' n <- 1000
#' data <- tibble::tibble(
#'   A = rnorm(n),
#'   B = rnorm(n),
#'   Y = A + B + A * B + rnorm(n)
#' )
#'
#' pc.fit <- pc(
#'   suffStat = list(C = cor(data), n = nrow(data)),
#'   indepTest = gaussCItest,
#'   alpha = 0.01,
#'   labels = colnames(data)
#' )
#'
#' ipc.fit <- ipc(pc.fit, data, x = 'A', y = 'Y', alpha = 0.01)
ipc <- function(pc.fit, data, x, y, alpha = 0.01) {
  cpdag <- pcalg::pdag2allDags(methods::as(pc.fit, 'amat'))

  dags <- cpdag$dags
  nodeNms <- cpdag$nodeNms
  idags <- t(apply(cpdag$dags, 1, dag2idag, data = data, x = x, y = y, alpha = alpha, nodeNms = cpdag$nodeNms))
  inodeNms <- inodeNms(cpdag$nodeNms, x, y)
  formulas <- sapply(1:nrow(dags), function(i) {
    toFormula(dag2amat(dags[i,], nodeNms), dag2amat(idags[i,], inodeNms), x, y)
  })
  equations <- unlist(lapply(formulas, function(f) { formula2equation(stats::lm(f, data)) }))

  ret <- list(
    pc.fit = pc.fit,
    data = data,
    x = x,
    y = y,
    dags = dags,
    nodeNms = nodeNms,
    idags = idags,
    inodeNms = inodeNms,
    formulas = formulas,
    equations = equations
  )
  class(ret) <- 'ipcAlgo'
  ret
}



ACE <- function(ipc.fit, idx, keep = list()) {
  do1 <- keep
  do2 <- keep
  x <- ipc.fit[['x']]
  do1[x] <- 1
  do2[x] <- 0
  equation <- ipc.fit$equations[idx]
  Deriv::Simplify(paste0(do(equation, do1), ' - (', do(equation, do2), ')'))
}

do <- function(string, do) {
  deparse1(Deriv::Simplify(methods::substituteDirect(str2lang(string), do)))
}

delta <- function(y, x) {
  paste0('Δ', y, x)
}

formula2equation <- function(model) {
  if (stats::is.empty.model(model)) {
    '0'
  } else {
    cc <- stats::coef(model)
    string <- paste(round(cc[1], 2), paste(round(cc[-1], 2), names(cc[-1]), sep = " * ", collapse = " + "), sep = " + ")
    string <- gsub(':', '*', string)
    string <- gsub('\\+ -', '- ', gsub(' \\* ', '*', string))
  }
}

dag2amat <- function(dag, nodeNms) {
  matrix(dag, length(nodeNms), dimnames = list(nodeNms, nodeNms))
}

inodeNms <- function(nodeNms, x, y) {
  c(setdiff(nodeNms, c(x, y)), delta(y, x))
}

interactionsWith <- function(x, y, amat) {
  parents <- names(which(amat[, y] == 1))

  ret <- NULL
  if (x %in% parents) {
    ret <- c(ret, paste0('(', paste(parents, collapse = '+'), ')^2'))
  }

  for (i in seq_along(parents)) {
    if (parents[i] != x) {
      ret <- c(ret, interactionsWith(x, parents[i], amat))
    }
  }

  return(ret)
}

adjustmentSetDAG <- function(x, y, amat) {
  inner <- function(node) {
    parents <- names(which(amat[, node] == 1))

    if (node == x) {
      return(node)
    }

    unlist(sapply(parents, function(parent) {
      rec <- inner(parent)

      if (!is.null(rec)) {
        rec
      }
      else {
        parent
      }
    }))
  }

  return(unique(inner(y)))
}

adjustmentSetIDAG <- function(y, amat) {
  names(which(amat[, y] == 1))
}


dag2idag <- function(dag, data, x, y, alpha, nodeNms) {
  amat <- dag2amat(dag, nodeNms)

  # Gets parents of y
  interactions <- interactionsWith(x, y, amat)

  # Adds Δy to adjacency matrix
  amat <- cbind(amat, 0)
  amat <- rbind(amat, 0)
  rownames(amat)[nrow(amat)] <- delta(y, x)
  colnames(amat)[ncol(amat)] <- delta(y, x)

  # Estimates significant interaction effects with x on y, then and adds these edges to y.
  if (length(interactions) > 0) {
    data.fit <- stats::lm(paste0(y, ' ~ ', paste(interactions, collapse = ' + ')), data)
    parents <- setdiff(names(which(stats::coef(summary(data.fit))[-1, 4] < alpha)), colnames(data))
    parents <- unlist(sapply(parents, function(p) { strsplit(p, ':') }), use.names = F)
    parents <- parents[parents != 'A']
    amat[parents, ncol(amat)] <- 1
  }

  # Removes columns x and y
  amat <- amat[setdiff(colnames(amat), c(x, y)), setdiff(rownames(amat), c(x, y)), drop = F]

  # Removes disconnected nodes
  clusters <- igraph::clusters(igraph::graph_from_graphnel(methods::as(amat, 'graphNEL')))
  membership <- clusters$membership[delta(y, x)]
  del <- names(which(clusters$membership != membership))
  if (length(del) > 0) {
    amat[del,] <- 0
    amat[, del] <- 0
  }
  amat
}

toFormula <- function(amat, iamat, x, y) {
  parents <- adjustmentSetDAG(x, y, amat)
  iparents <- adjustmentSetIDAG(delta(y, x), iamat)

  f <- paste0(y, ' ~ ')
  if (length(parents) == 0) {
    f <- paste0(f, '0')
  }
  else {
    f <- paste0(f, paste(parents, collapse = ' + '))
    if (length(iparents) == 1) {
      f <- paste0(f, ' + ', x, ':', paste(iparents, collapse = '+'))
    }
    if (length(iparents) > 1) {
      f <- paste0(f, ' + ', x, ':(', paste(iparents, collapse = '+'), ')')
    }
  }
  f
}

#' @export
plot.ipcAlgo <- function(obj, idx = NULL) {
  requireNamespace("Rgraphviz")
  if (is.null(idx)) {
    n <- 1:nrow(obj$dags)
  } else {
    n <- c(idx)
  }
  graphics::par(mfrow = c(1, 2))
  for (i in n) {
    Rgraphviz::plot(methods::as(dag2amat(obj$dags[i,], obj$nodeNms), 'graphNEL'))
    graphics::title('DAG', line = -2)
    Rgraphviz::plot(methods::as(dag2amat(obj$idags[i,], obj$inodeNms), 'graphNEL'))
    graphics::title('IDAG', line = -2)
  }
  graphics::par(mfrow = c(1, 1))
}


#' @export
summary.ipcAlgo <- function(obj) {
  dags <- obj$dags
  n <- nrow(dags)
  x <- obj$x
  y <- obj$y
  cat(sprintf('Number of DAGS: \t%s\n\n', n))
  cat('Formulas:\n')
  for (i in 1:n) {
    cat(sprintf('[%s] %s\n', i, obj$formulas[i]))
  }

  cat('\nCausal function:\n')
  for (i in 1:n) {
    cat(sprintf('[%s] %s = %s\n', i, y, obj$equations[i]))
  }

  cat(sprintf('\nACE = E[%s|do(%s=1)] - E[%s|do(%s=0)]:\n', y, x, y, x))
  for (i in 1:n) {
    cat(sprintf('[%s] %s\n', i, ACE(obj, i)))
  }

  # cat('\nStratified interaction effect:\n')
  # for (i in 1:n) {
  #   cat(sprintf('[%s] %s\n', i, obj$stratifications[i]))
  # }
}