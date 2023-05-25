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
