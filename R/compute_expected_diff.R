# Function to evaluate change in expected outcome given two sets of exposure values

#' Evaluate expected change in outcome given a shift in exposure values.
#'
#' Evaluate expected change in outcome given a shift from exposures x1 to exposures x2. Optionally can take covariates z1 and z2, but these don't have an effect on the result if z1=z2.
#'
#' @export
#' @param fit Model fit returned by fit_lowFR function with default setting "all"
#' @param x1 Numeric vector of exposure values
#' @param x2 Numeric vector of exposure values
#' @param z1 (Optional) Numeric vector of covariate values
#' @param z2 (Optional) Numeric vector of covariate values
#'
#' @return A vector of posterior samples of the expected outcome given a shift from x1 to x2. (And optionally from z1 to z2.)

#'
compute_expected_diff <- function(fit, x1, x2, z1=NULL, z2=NULL) {

  # initialize z1 and z2 to 0 if not provided
  if (is.null(z1)) {
    if (is.null(z2)) {
      z1 <- rep(0, dim(fit$zeta)[2])
      z2 <- rep(0, dim(fit$zeta)[2])
    }
    else {
      print("z2 provided but z1 is null. Setting z1 equal to z2.")
      z1 <- z2
    }
  }
  else if (is.null(z2)) {
    print("z1 provided but z2 is null. Setting z2 equal to z1.")
    z2 <- z1
  }

  # compute pre and post outcome posterior
  pre_outcome <- compute_expected_outcome(fit=fit, x_obs=x1, z_obs=z1)
  post_outcome <- compute_expected_outcome(fit=fit, x_obs=x2, z_obs=z2)

  # return the difference
  return(post_outcome - pre_outcome)
}
