ecos_fit_binomial_topic_model <- function(X,
                                          k,
                                          numiter.main = 100,
                                          numiter.refine = 20,
                                          method.main = "em",
                                          init.method = c("topicscore", "random"),
                                          control.init = list(),
                                          control.main = list(numiter = 4),
                                          control.refine = list(numiter = 4, extrapolate = TRUE),
                                          verbose = c(
                                            "progressbar",
                                            "detailed", "none"
                                          )) {
  verify.count.matrix(X)
  verbose <- match.arg(verbose)
  if (any.allzero.cols(X)) {
    X <- remove.allzero.cols(X)
    warning(sprintf(paste(
      "One or more columns of X are all zero; after",
      "removing all-zero columns, %d columns will be",
      "used for model fitting"
    ), ncol(X)))
  }
  fit_pois <- fastTopics::init_poisson_nmf(X,
    k = k, init.method = init.method,
    control = control.init,
    verbose = ifelse(verbose == "none", "none", "detailed")
  )
  fit_pois <- fastTopics::fit_poisson_nmf(X,
    fit0 = fit_pois, numiter = numiter.main,
    method = method.main, control = list(extrapolate = TRUE),
    verbose = verbose
  )
  if (numiter.refine > 0) {
    if (verbose != "none") {
      cat("Refining binomial model fit with EM updates.\n")
    }
    fit_binom_em <- poisson2binom(X, fit_pois, numem = numiter.refine, verbose = FALSE)
  }
  return(fit_binom_em)
}

#' @importFrom Matrix colSums
any.allzero.cols <- function (X)
  any(colSums(X > 0) == 0)

# Filter out all-zero columns from the matrix.
#
#' @importFrom Matrix colSums
remove.allzero.cols <- function (X)
  X[,colSums(X > 0) >= 1]

# Add back all-zero columns to the theta matrix.
#
#' @importFrom Matrix colSums
addback.allzero.cols <- function(X, K, Z){
  ind <- c(1:dim(X)[2])[colSums(X > 0) >= 1]
  theta <- matrix(data = NA, nrow = dim(X)[2], ncol = 2) 
  theta[,] <- c(0,0)
  theta[ind,] <- Z
  return(theta)
}
  
  

