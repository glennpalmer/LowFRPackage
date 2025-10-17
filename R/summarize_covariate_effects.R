#' Summarize covariate effects
#'
#' Compute posterior means and credible intervals for the linear covariate effects.
#'
#' @export
#' @param fit Model fit returned by fit_lowFR function with default setting "all"
#' @param varnames q dimensional vector of names of covariates (i.e., the columns of Z)
#' @param lower (Optional) Lower quantile value between 0 and 1
#' @param upper (Optional) Upper quantile value between 0 and 1
#'
#' @return A dataframe with p*TT rows summarizing each main effect.

#'
summarize_covariate_effects <- function(fit, varnames, lower=0.025, upper=0.975) {
  # initialize storage for summaries
  num_effects <- dim(fit$zeta)[2]
  mean_vec <- rep(NA, num_effects)
  lower_vec <- rep(NA, num_effects)
  upper_vec <- rep(NA, num_effects)

  # loop over effects and compute summaries
  for (i in 1:num_effects) {
    curr_effect <- fit$zeta[,i]
    curr_summary <- summarize_posterior(samples=curr_effect, lower=lower, upper=upper)
    mean_vec[i] <- curr_summary$mean
    lower_vec[i] <- as.numeric(curr_summary$lower_quantile)
    upper_vec[i] <- as.numeric(curr_summary$upper_quantile)
  }

  # create data frame for output
  output_df <- data.frame(varnames, mean_vec, lower_vec, upper_vec)
  names(output_df) <- c("Covariate", "Mean", "Lower", "Upper")

  # return the df
  return(output_df)
}
