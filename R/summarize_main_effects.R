#' Summarize main effects for all exposures
#'
#' Compute posterior means and credible intervals for the change in expected outcome if each exposure is shifted from -0.5 to 0.5, with all other exposures held constant at 0.
#'
#' @export
#' @param fit Model fit returned by fit_lowFR function with default setting "all"
#' @param varnames p*TT dimensional vector of names of exposures
#' @param lower (Optional) Lower quantile value between 0 and 1
#' @param upper (Optional) Upper quantile value between 0 and 1
#'
#' @return A dataframe with p*TT rows summarizing each main effect.

#'
summarize_main_effects <- function(fit, varnames, lower=0.025, upper=0.975) {
  # initialize effects
  num_effects <- dim(fit$alpha)[2]
  x1_mat <- diag(num_effects) * -0.5
  x2_mat <- diag(num_effects) * 0.5

  # initialize storage for summaries
  mean_vec <- rep(NA, num_effects)
  lower_vec <- rep(NA, num_effects)
  upper_vec <- rep(NA, num_effects)

  # loop over effects and compute summaries
  for (i in 1:num_effects) {
    curr_diff <- compute_expected_diff(fit=fit, x1=x1_mat[i,], x2=x2_mat[i,])
    curr_summary <- summarize_posterior(samples=curr_diff, lower=lower, upper=upper)
    mean_vec[i] <- curr_summary$mean
    lower_vec[i] <- as.numeric(curr_summary$lower_quantile)
    upper_vec[i] <- as.numeric(curr_summary$upper_quantile)
  }

  # create data frame for output
  output_df <- data.frame(varnames, mean_vec, lower_vec, upper_vec)
  names(output_df) <- c("Exposure", "Mean", "Lower", "Upper")

  # return the df
  return(output_df)
}
