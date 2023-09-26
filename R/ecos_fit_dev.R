#' @title Fitting Grade of Membership model for clustering into communities
#'
#' @description A grade of membership model for clustering a site by features
#' data - which could be presence-absence or counts abundance data of species
#' in the sites.
#'
#' @param dat input data matrix of samples along rows and features along sites,
#'            with each entry a 0/1 signifying presence/absence or counts of
#'            abundances.
#' @param K The number of clusters to fit
#' @param tol The tolerance level of the model.
#' @param num_trials Number of EM runs from different starting points. This is'
#'                   key for picking the best fit model across multiple runs.
#' @param fit_control The control parameters for the model.
#'
#' @return Returns a model fit with \code{omega} as cluster membership
#'         probabilities and \code{theta} as cluster feature matrix.
#'
#' @importFrom fastTopics fit_topic_model
#'
#' @examples
#'
#' data("himalayan_birds")
#' species_abundance_counts <- t(exprs(himalayan_birds))
#' fit <- ecostructure_fit(species_abundance_counts, K = 2, tol = 0.1)
#' species_pa_counts <- species_abundance_counts
#' species_pa_counts[species_pa_counts >= 1] <- 1
#' fi2 <- ecos_fit(species_pa_counts, K = 2, tol = 0.1)
#'
#' @export



ecos_fit_dev <- function(dat,
                         max_dat = NULL,
                         K,
                         sparse = F,
                         num_iter_main = 100,
                         num_iter_refine = 100,
                         fit_control_init = list(),
                         fit_control_main = list(numiter = 4),
                         fit_control_refine = list(numiter = 4, extrapolate = TRUE)) {
  row_names <- rownames(dat)
  if (length(row_names) != dim(dat)[1]) {
    warning("row names not provided, or not proper, using fake rownames")
    rownames(dat) <- 1:dim(dat)[1]
  }
  if (all(dat - floor(dat) != 0)) {
    stop("The matrix dat must be a matrix of integers - a binary or a counts matrix")
  }

  ids_na <- which(is.na(dat))
  if (length(ids_na) > 0) {
    warning("NAs in dat matrix: replaced by 0")
    dat[is.na(dat)] <- 0
  }

  dat_by_2 <- dat %% 2
  if (all(dat %% 2 - dat == 0)) {
    if (max(dat) == 1 & min(dat) == 0) {
      message("Binary matrix input: Fitting the Binomial Grade of Membership
              model.")

      fit_control_init_default <- list()
      fit_control_init <- modifyList(fit_control_init_default, fit_control_init)

      fit_control_main_default <- list()
      fit_control_main <- modifyList(fit_control_main_default, fit_control_main)

      fit_control_refine_default <- list()
      fit_control_refine <- modifyList(fit_control_refine_default, fit_control_refine)
      
      if (sparse == T) {
        dat <- as(dat, "sparseMatrix")
      }

      topic_clus <- do.call(
        ecos_fit_binomial_topic_model,
        list(
          X = dat,
          k = K,
          numiter.main = num_iter_main,
          numiter.refine = num_iter_refine,
          control.init = fit_control_init,
          control.main = fit_control_main,
          control.refine = fit_control_refine
        )
      )

      ll <- list(
        "omega" = topic_clus$L,
        "theta" = topic_clus$F
      )
    }
  }

  if (max(dat) > 1) {
    fit_control_init_default <- list()
    fit_control_init <- modifyList(fit_control_init_default, fit_control_init)

    fit_control_main_default <- list()
    fit_control_main <- modifyList(fit_control_main_default, fit_control_main)

    fit_control_refine_default <- list()
    fit_control_refine <- modifyList(fit_control_refine_default, fit_control_refine)

    if (sparse == T) {
      dat <- as(dat, "sparseMatrix")
    }
    topic_clus <- do.call(
      fastTopics::fit_topic_model,
      list(
        X = dat,
        k = K,
        numiter.main = num_iter_main,
        numiter.refine = num_iter_refine,
        control.init = fit_control_refine,
        control.main = fit_control_main,
        control.refine = fit_control_refine
      )
    )
    if (any.allzero.cols(dat)) {
      topic_clus$F <- addback.allzero.cols(X = dat, K = K, Z = topic_clus$F)
      warning(sprintf(paste(
        "One or more columns of X are all zero; Columns that", 
        "were removed from the model fit are being added back",
        "to the theta matrix with zero values for each motif k."
        )))
    }
    ll <- list(
      "omega" = topic_clus$L,
      "theta" = topic_clus$F,
      "BIC" = topic_clus$BIC
    )
  }
  return(ll)
}

# Compute two-tailed p-value from z-score.
#
#' @importFrom stats pnorm
pfromz <- function (z)
  2*pnorm(-abs(z))

# Compute log10 two-tailed p-value from z-score.
#
#' @importFrom stats pnorm
lpfromz <- function (z)
  (log(2) + pnorm(-abs(z),log.p = TRUE))/log(10)

# Set all entries of x less than a to a, and set alll entries of x
# greater than b to b.
clamp <- function (x, a, b)
  pmax(pmin(x,b),a)

# Return true if x is a compressed, sparse, column-oriented numeric
# matrix.
is.sparse.matrix <- function (x)
  inherits(x,"dgCMatrix")

# Efficiently extract the nonzero elements from column j of sparse
# matrix A (a member of class "dgCMatrix"). Output "x" contains the
# nonzero values, and output "i" contains the
get.nonzeros <- function (A, j)
  list(x = A[,j,drop = FALSE]@x,i = A[,j,drop = FALSE]@i + 1)

# Check if the matrix contains one or more all-zero columns.
#
#' @importFrom Matrix colSums
any.allzero.cols <- function (X)
  any(colSums(X > 0) == 0)

# Filter out all-zero columns from the matrix.
#
#' @importFrom Matrix colSums
remove.allzero.cols <- function (X)
  X[,colSums(X > 0) >= 1]

# Apply operation f to all nonzeros of a sparse matrix.
#
#' @importFrom Matrix sparseMatrix
#' 
apply.nonzeros <- function (X, f) {
  d <- summary(X)
  return(sparseMatrix(i = d$i,j = d$j,x = f(d$x),dims = dim(X)))
}

# Compute X/(crossprod(A,B) + e) efficiently when X is a sparse
# matrix.
#
#' @importFrom Matrix sparseMatrix
#' @importFrom Rcpp evalCpp
#' 
x_over_tcrossprod <- function (X, A, B, e) {
  d <- summary(X)
  y <- drop(x_over_crossprod_rcpp(d$i - 1,d$j - 1,d$x,A,B,e))
  return(sparseMatrix(i = d$i,j = d$j,x = y,dims = dim(X)))
}

# Return an m x n matrix rbind(x,...,x), in which length(x) = m.
repmat <- function (x, n)
  matrix(x,n,length(x),byrow = TRUE)

# scale.cols(A,b) scales each column A[,i] by b[i].
scale.cols <- function (A, b)
  t(t(A) * b)

# Scale each row of A so that the entries of each row sum to 1.
#
#' @importFrom Matrix rowSums
#'
normalize.rows <- function (A)
  A / rowSums(A)

# Scale each column of A so that the entries of each column sum to 1.
#
#' @importFrom Matrix colSums
#'
normalize.cols <- function (A)
  t(t(A) / colSums(A))

# Scale each row of A so that the large entry in each row is 1.
normalize.rows.by.max <- function (A) {
  if (!is.matrix(A))
    stop("Input argument \"A\" should be a matrix")
  return(A / apply(A,1,max))
}

# For vector x, return a vector of the same length y containing the
# "least extreme" differences y(i) = x(i) - x(j), in which j is the
# index not equal to i such that abs(x(i) - x(j)) is the smallest
# possible. The length of x should be 2 or more.
le.diff <- function (x) {
  n <- length(x)
  if (n == 2) {
    y <- x[1] - x[2]
    y <- c(y,-y)
  } else {
    y <- rep(0,n)
    for (i in 1:n) {
      d <- x[i] - x
      j <- order(abs(d))[2]
      y[i] <- d[j]
    }
  }
  return(y)
}

# Rescale the factors (F) and loadings (L) with the property that
# tcrossprod(L,F) remains the same after rescaling; specifically,
# rescale the columns of F and L so that, for each k, column k of F
# has the same mean as column k of L.
#
#' @importFrom Matrix colMeans
#'
rescale.factors <- function (F, L) {
  d <- sqrt(colMeans(L)/colMeans(F))
  return(list(F = scale.cols(F,d),
              L = scale.cols(L,1/d)))
}

# This does the same thing as the "rand" function in MATLAB.
#
#' @importFrom stats runif
#'
rand <- function (n, m, min = 0, max = 1) 
  matrix(runif(n*m,min,max),n,m)

# Initialize RcppParallel multithreading using a pre-specified number
# of threads, or using the default number of threads when "n" is NA.
#
#' @importFrom RcppParallel setThreadOptions
#' @importFrom RcppParallel defaultNumThreads
#'
initialize.multithreading <- function (n, verbose = FALSE) {
  if (is.na(n)) {
    setThreadOptions()
    n <- defaultNumThreads()
  } else
    setThreadOptions(numThreads = n)
  if (verbose && n > 1)
    message(sprintf("Using %d RcppParallel threads.",n))
  return(n)
}

# For a Poisson non-negative matrix factorization with rank = 1, the
# maximum-likelihood estimate (MLE) has a closed-form solution (up to
# a scaling factor); this function returns the MLE subject to the
# constraint that mean(F) = mean(L).
#
#' @importFrom Matrix rowMeans
#' @importFrom Matrix colMeans
#'
fit_pnmf_rank1 <- function (X)
  list(F = matrix(colMeans(X)),
       L = matrix(rowMeans(X)))

# Compute the highest posterior density (HPD) interval from a vector
# of random draws from the distribution. See Chen & Shao (1999) for
# background on HPD intervals.
hpd <- function (x, conf.level = 0.68) {
  n <- length(x)
  m <- round(n*(1 - conf.level))
  x <- sort(x)
  y <- x[seq(n-m+1,n)] - x[seq(1,m)]
  i <- which.min(y)
  return(c(x[i],x[n-m+i]))
}

# This replicates the minimum Kullback-Leibler (KL) divergence
# calculation used in ExtractTopFeatures from CountClust, with method
# = "poisson". Input F should be an n x k matrix of frequency
# estimates from the multinomial topic model, where n is the number of
# data columns, and k is the number of topics. The return value is a
# matrix of the same dimension as F containing the minimum
# KL-divergence calculations.
min_kl_poisson <- function (F, e = 1e-15) {
  
  # Get the number of rows (n) and columns (k) of F.
  n <- nrow(F)
  k <- ncol(F)
  
  # Compute the minimum KL-divergence measure for each row and column
  # of F.
  D <- matrix(0,n,k)
  for (i in 1:n) {
    f <- F[i,] + e
    for (j in 1:k) {
      y      <- f[j]*log(f[j]/f) + f - f[j]
      D[i,j] <- min(y[-j])
    }
  }
  
  dimnames(D) <- dimnames(F)
  return(D)
}

# Compute "least extreme" LFC statistics LFC(j) = log2(fj/fk) given
# frequency estimates F. Input F should be an n x k matrix of
# frequency estimates from the multinomial topic model, where n is the
# number of data columns, and k is the number of topics. The return
# value is a matrix of the same dimension as F containing the LFC
# estimates.
le_lfc <- function (F, e = 1e-15) {
  n <- nrow(F)
  k <- ncol(F)
  B <- matrix(0,n,k)
  for (i in 1:n)
    B[i,] <- le.diff(log2(F[i,] + e))
  dimnames(B) <- dimnames(F)
  return(B)
}



poisson2binom <- function (X, fit, numem = 0, umin = 1e-4, verbose = TRUE) {
  if (!requireNamespace("NNLM",quietly = TRUE))
    stop("poisson2binom requires the NNLM package")
  
  # Check input argument "fit".
  if (inherits(fit,"binom_topic_model_fit"))
    return(fit)
  if (!inherits(fit,"poisson_nmf_fit"))
    stop("Input argument \"fit\" should be an object of class ",
         "\"poisson_nmf_fit\"")
  verify.fit(fit)
  if (ncol(fit$F) < 2 | ncol(fit$L) < 2)
    stop("Input matrices \"fit$F\" and \"fit$L\" should have 2 or more",
         "columns")
  
  # Check and process input argument "X".
  if (!((is.numeric(X) & is.matrix(X)) | is.sparse.matrix(X)))
    stop("Input argument \"X\" should be a numeric matrix (a \"matrix\" or ",
         "a \"dgCMatrix\")")
  verify.fit.and.count.matrix(X,fit)
  if (any(X < 0) | any(X > 1))
    warning("Input argument \"X\" should be a \"binary\" matrix ",
            "(that is, all entries should range from 0 and 1)")
  if (is.matrix(X) & is.integer(X))
    storage.mode(X) <- "double"
  
  # Choose U = diag(u) such that L*U is closer to being a matrix of
  # topic proportions.
  L <- fit$L
  F <- fit$F
  n <- nrow(L)
  ones <- matrix(1,n,1)
  if (verbose)
    cat("Rescaling L and F using non-negative linear regression (nnlm).\n")
  u   <- drop(coef(NNLM::nnlm(L,ones)))
  u   <- pmax(u,umin)
  L   <- scale.cols(L,u)
  L   <- normalize.rows(L)
  F   <- scale.cols(F,1/u)
  fit <- list(F = F,L = L,progress = NA)
  
  # Refine the binomial topic model fit by performing several EM updates.
  if (numem > 0) {
    cat("Performing",numem,"EM updates to refine the fit.\n")
    progress <- as.data.frame(cbind(1:numem,0,0))
    names(progress) <- c("iter","delta.f","delta.l")
    if (verbose)
      cat("iter  |F - F'|  |L - L'|\n")
    for (i in 1:numem) {
      fit0 <- fit
      fit  <- fit_binom_topic_model_em(X,fit,numem)
      progress[i,"delta.f"] <- max(abs(fit0$F - fit$F))
      progress[i,"delta.l"] <- max(abs(fit0$L - fit$L))
      if (verbose)
        cat(sprintf("%4d %0.3e %0.3e\n",i,progress[i,"delta.f"],
                    progress[i,"delta.f"]))
    }
    fit$progress <- progress
  }
  
  # Return the Binomial topic model fit.
  fit$s <- rep(1,n)
  names(fit$s) <- rownames(L)
  class(fit) <- c("binom_topic_model_fit","list")
  return(fit)
}

# Perform a single EM udpate for fiitting the binomial topic model to
# binary data matrix X. This code is adapted from the meth_tpxEM
# function in the methClust package by Kushal Dey.
fit_binom_topic_model_em <- function (X, fit, numiter) {
  if (!is.matrix(X))
    X <- as.matrix(X)
  
  # Make sure no parameters are exactly zero or exactly one.
  e <- 1e-8
  L <- fit$L
  F <- fit$F
  F <- clamp(F,e,1 - e)
  L <- clamp(L,e,1 - e)
  L <- normalize.rows(L)
  
  # Perform the E step.
  A  <- X/tcrossprod(L,F)
  M  <- (A %*% F) * L
  Mt <- crossprod(A,L) * F
  A  <- (1 - X)/tcrossprod(L,1 - F)
  U  <- (A %*% (1 - F)) * L
  Ut <- crossprod(A,L) * (1 - F)
  
  # Perform the M step.
  L <- normalize.rows(M + U)
  F <- Mt/(Mt + Ut)
  return(list(F = F,L = L))
}





######## verify args
# Verify that x is a vector with positive entries.
verify.positive.vector <- function (x, arg.name = deparse(substitute(x))) {
  arg.name <- sprintf("\"%s\"",arg.name)
  msg <- paste("Input argument",arg.name,"should be a numeric vector in",
               "which all entries are finite, non-missing and positive")
  if (!is.numeric(x))
    stop(msg)
  else if (any(x <= 0) | any(is.infinite(x)) | anyNA(x))
    stop(msg)
  return(TRUE)
}

# Verify that x is non-negative matrix.
verify.nonnegative.matrix <- function (x, arg.name = deparse(substitute(x))) {
  arg.name <- sprintf("\"%s\"",arg.name)
  msg <- paste("Input argument",arg.name,"should be a non-negative,",
               "numeric matrix (a \"matrix\" or a \"dgCMatrix\"), and",
               "all entries should be finite and non-missing")
  if (!((is.matrix(x) & is.numeric(x)) | is.sparse.matrix(x)))
    stop(msg)
  else if (any(x < 0) | any(is.infinite(x)) | anyNA(x))
    stop(msg)
  return(TRUE)
}

# Verify that x is a valid count matrix.
verify.count.matrix <- function (x, arg.name = deparse(substitute(x))) {
  verify.nonnegative.matrix(x,arg.name)
  arg.name <- sprintf("\"%s\"",arg.name)
  if (!(nrow(x) > 1 & ncol(x) > 1))
    stop(paste("Input matrix",arg.name,"should have at least 2 rows",
               "and 2 columns"))
  return(TRUE)
}

# Verify that x is a valid multinomial topic model fit or Poisson
# non-negative matrix factorization.
verify.fit <- function (x, arg.name = deparse(substitute(x))) {
  arg.name.F <- paste0(arg.name,"$F")
  arg.name.L <- paste0(arg.name,"$L")
  arg.name.s <- paste0(arg.name,"$s")
  arg.name   <- sprintf("\"%s\"",arg.name)
  msg        <- paste("Input argument",arg.name,"should be a list containing",
                      "non-negative matrices \"F\" and \"L\"")
  if (!is.list(x))
    stop(msg)
  else if (!all(is.element(c("F","L"),names(x))))
    stop(msg)
  verify.nonnegative.matrix(x$F,arg.name.F)
  verify.nonnegative.matrix(x$L,arg.name.L)
  arg.name.F <- sprintf("\"%s\"",arg.name.F)
  arg.name.L <- sprintf("\"%s\"",arg.name.L)
  if (ncol(x$F) != ncol(x$L))
    stop(paste("Input matrices",arg.name.F,"and",arg.name.L,"should have",
               "the same number of columns"))
  if (is.element("s",names(x))) {
    
    # Check the vector of "scale factors", s.
    verify.positive.vector(x$s,arg.name.s)
    arg.name.s <- sprintf("\"%s\"",arg.name.s)
    if (length(x$s) != nrow(x$L))
      stop(paste("The length of input vector",arg.name.s,"should equal the",
                 "number of rows in",arg.name.L))
  }
  return(TRUE)
}

# Verify that x is a valid count matrix and "fit" is a valid topic model
# fit or non-negative matrix factorization.
verify.fit.and.count.matrix <-
  function (x, fit,
            arg.name.x   = deparse(substitute(x)),
            arg.name.fit = deparse(substitute(fit))) {
    verify.count.matrix(x,arg.name.x)
    verify.fit(fit,arg.name.fit)
    arg.name.x <- sprintf("\"%s\"",arg.name.x)
    arg.name.F <- sprintf("\"%s$F\"",arg.name.fit)
    arg.name.L <- sprintf("\"%s$L\"",arg.name.fit)
    if (!(nrow(fit$L) == nrow(x) & nrow(fit$F) == ncol(x)))
      stop(paste("Dimensions of input matrices",arg.name.x,",",arg.name.F,
                 "and",arg.name.L,"do not agree"))
    if (!(identical(rownames(fit$L),rownames(x)) &
          identical(rownames(fit$F),colnames(x))))
      stop(paste("Dimnames of input matrices",arg.name.x,",",arg.name.F,
                 "and",arg.name.L,"are not consistent"))
    return(TRUE)
  }

# Return TRUE if x is a finite scalar with no missing entries.
is.scalar <- function (x)
  is.numeric(x) &
  length(x) == 1 &
  all(!is.na(x)) &
  all(is.finite(x))
