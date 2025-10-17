#' Plot covariate effects for all exposures
#'
#' Plot posterior means and credible intervals for the covariate effects.
#'
#' @export
#' @param fit Model fit returned by fit_lowFR function with default setting "all"
#' @param varnames q dimensional vector of names of exposures
#' @param lower (Optional) Lower quantile value between 0 and 1
#' @param upper (Optional) Upper quantile value between 0 and 1
#'
#' @return A ggplot object plotting main effect posterior means and credible intervals.

#'
plot_covariate_effects <- function(fit, varnames, lower=0.025, upper=0.975) {
  # compute covariate effect dataframe
  covariate_effect_df <- summarize_covariate_effects(fit=fit, varnames=varnames, lower=lower, upper=upper)
  covariate_effect_df$Covariate <- factor(covariate_effect_df$Covariate,
                                    levels=rev(covariate_effect_df$Covariate))


  # make plot
  covariate_plot <- ggplot(data=covariate_effect_df, aes(x=Covariate, y=Mean)) +
    geom_errorbar(aes(ymin=Lower, ymax=Upper)) +
    geom_hline(yintercept=0, color="black") +
    geom_point(color="red") +
    coord_flip() +
    labs(title="Covariate effects",
         y="", x="")

  # return ggplot object
  return(covariate_plot)
}
