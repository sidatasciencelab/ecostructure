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
#' @importFrom maptpx topics
#' @importFrom methClust meth_topics
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

    ll <- list(
      "omega" = topic_clus$L,
      "theta" = topic_clus$F,
      "BIC" = topic_clus$BIC
    )
  }
  return(ll)
}
