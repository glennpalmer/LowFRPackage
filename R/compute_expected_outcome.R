# Function to evaluate expected outcome for a given set of exposure values

#' Evaluate expected outcome given exposure and covariate values
#'
#' @export
#' @param fit Model fit returned by fit_lowFR function with default setting "all"
#' @param x_obs Numeric vector of exposure values
#' @param z_obs Numeric vector of covariate values
#'
#' @return A vector of posterior samples of the expected outcome

#'
compute_expected_outcome <- function(fit, x_obs, z_obs) {

  # initialize storage
  n_samps <- dim(fit$alpha_0)
  output_vec <- rep(NA, n_samps)

  # loop over posterior samples and compute expected outcome
  for (i in 1:n_samps) {
    curr_val <- fit$alpha_0[i] +
      as.vector(fit$alpha[i,] %*% x_obs) +
      as.vector(t(x_obs) %*% fit$Gamma[i,,] %*% x_obs) +
      as.vector(t(z_obs) %*% fit$zeta[i,])
    output_vec[i] <- curr_val
  }

  return(output_vec)
}
