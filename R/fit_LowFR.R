# Functions to fit and summarize LowFR models as described in "Low-Rank
# Longitudinal Factor Regression" paper

#library(tidyverse)
#library(rstan)
#source("R/helper.R")
#source("R/sim_data.R")

#' Fit a Low-Rank Longitudinal Factor Regression Model
#'
#' @export
#' @param X Numeric matrix of exposures.
#' @param y Numeric vector of output values.
#' @param ... Arguments passed to `rstan::sampling` (e.g. iter, chains).
#' @return An object of class `stanfit` returned by `rstan::sampling`
#'
fit_LowFR <- function(y, X, Z, p, TT, k=NULL,
                      output="all",
                      burnin=1000, samples=1000, chains=4,
                      random_seed=1234) {

  # Note: output can be either "stan_fit", "all", or a list of parameters
  # for which to return posterior samples. This list can include items from
  #       "alpha_0"
  #       "alpha"
  #       "Gamma"
  #       "zeta"
  #       "sigma2"
  #       "Sigma"
  #       "phi"
  #       "tau"
  # If "all" is specified, samples for the full list above will be returned as
  # a list of arrays.

  # detect n_obs
  n_obs <- length(y)

  # check for missing y or Z values
  missing_y_or_Z <- FALSE
  if (sum(is.na(y)) > 0) {
    print("Please drop or impute missing y values before model fitting. LowFR is currently only able to impute missing X values.")
    missing_y_or_Z <- TRUE
  }
  if (sum(is.na(Z)) > 0) {
    print("Please drop or impute missing Z values before model fitting. LowFR is currently only able to impute missing X values.")
    missing_y_or_Z <- TRUE
  }
  if (missing_y_or_Z) {
    return(NULL)
  }

  # set up machinery for imputation if any missing X values
  X_missing <- is.na(X)
  num_missing <- sum(X_missing)
  if (num_missing > 0) {
    print(paste0(num_missing, " missing exposure values will be imputed during each posterior sample."))
  }
  # replace NAs with zeros so stan will accept it
  X_replace_nas <- X
  for (i in 1:nrow(X)) {
    for (j in 1:ncol(X)) {
      if (is.na(X[i,j])) {
        X_replace_nas[i,j] = 0
      }
    }
  }
  # index missing values to match indices of imputation parameters
  curr_index <- 1
  for (i in 1:N) {
    for (j in 1:(p*TT)) {
      if (X_missing[i,j] > 0) {
        X_missing[i,j] <- curr_index
        curr_index <- curr_index + 1
      }
      else {
        X_missing[i,j] <- 0
      }
    }
  }

  # choose k using SVD if needed
  if (is.null(k)) {
    k <- k_svd_LowFR(X=X, p=p, TT=TT)
  }

  # fit model
  options(mc.cores = parallel::detectCores())
  fit <- sampling(stanmodels$LowFR,
                  list(N=n_obs,
                       p=p,
                       q=ncol(Z),
                       k=k,
                       TT=TT,
                       H=min(p,TT),
                       X=X,
                       y=y,
                       Z=Z),
                  chains=chains,
                  iter=burnin+samples,
                  warmup=burnin,
                  seed=random_seed,
                  init=0)

  # return the stan fit if specified as output
  if (output == "stan_fit") {
    return(fit)
  }

  # otherwise, extract posterior samples and return list of specified parameters
  post_samples <- rstan::extract(fit)
  results <- list()
  result_names <- c()
  if ("alpha_0" %in% output | output == "all") {
    #results <- append(results, post_samples$alpha_0)
    results <- c(results, list(post_samples$alpha_0))
    result_names <- c(result_names, "alpha_0")
  }
  if ("alpha" %in% output | output == "all") {
    #results <- append(results, post_samples$alpha)
    results <- c(results, list(post_samples$alpha))
    result_names <- c(result_names, "alpha")
  }
  if ("Gamma" %in% output | output == "all") {
    #results <- append(results, post_samples$Gamma)
    results <- c(results, list(post_samples$Gamma))
    result_names <- c(result_names, "Gamma")
  }
  if ("zeta" %in% output | output == "all") {
    #results <- append(results, post_samples$Gamma)
    results <- c(results, list(post_samples$zeta))
    result_names <- c(result_names, "zeta")
  }
  if ("sigma2" %in% output | output == "all") {
    #results <- append(results, post_samples$sigma2)
    results <- c(results, list(post_samples$sigma2))
    result_names <- c(result_names, "sigma2")
  }
  if ("Sigma" %in% output | output == "all") {
    #results <- append(results, post_samples$Sigma)
    results <- c(results, list(post_samples$Sigma))
    result_names <- c(result_names, "Sigma")
  }
  if ("phi" %in% output | output == "all") {
    #results <- append(results, post_samples$phi)
    results <- c(results, list(post_samples$phi))
    result_names <- c(result_names, "phi")
  }
  if ("tau" %in% output | output == "all") {
    #results <- append(results, post_samples$tau)
    results <- c(results, list(post_samples$tau))
    result_names <- c(result_names, "tau")
  }
  names(results) <- result_names
  return(results)
}

###############################################################################
######################### Choose k based on SVD ###############################
###############################################################################
k_svd_LowFR <- function(X, p, TT) {
  # create data matrix of x_it's to do SVD and determine number of factors
  Xit_mat <- matrix(nrow=(TT*nrow(X)), ncol=p)

  # loop over data to populate matrix
  for (i in 1:nrow(X)) {
    for (t in 1:TT) {
      Xit_mat[(TT*(i-1) + t),] <- X[i,seq(from=t, to=p*TT, by=TT)]
    }
  }

  # do SVD to get singular values
  sing_vals <- svd(Xit_mat)$d
  frac_sing_vals <- rep(NA, p)
  for (i in 1:p) {
    frac_sing_vals[i] <- sum(sing_vals[1:i]) / sum(sing_vals)
  }

  # find index where frac first exceeds 0.9
  for (i in 1:p) {
    if (frac_sing_vals[i] > 0.9) {
      print(paste0("Setting k = ", i, " based on SVD."))
      return(i)
    }
  }
  print("k selection failed. Setting k = p.")
  return(p)
}
