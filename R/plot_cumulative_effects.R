#' Plot cumulative effects for all chemicals
#'
#' Plot posterior means and credible intervals for the change in expected outcome if each chemical is shifted from -0.5 to 0.5 at all times, with all other chemicals held constant at 0.
#'
#' @export
#' @param fit Model fit returned by fit_lowFR function with default setting "all"
#' @param varnames p dimensional vector of names of chemicals
#' @param lower (Optional) Lower quantile value between 0 and 1
#' @param upper (Optional) Upper quantile value between 0 and 1
#'
#' @return A ggplot object plotting cumulative effect posterior means and credible intervals.

#'
plot_cumulative_effects <- function(fit, varnames, lower=0.025, upper=0.975) {
  # compute main effect dataframe
  cumulative_effect_df <- summarize_cumulative_effects(fit=fit, varnames=varnames, lower=lower, upper=upper)
  cumulative_effect_df$Exposure <- factor(cumulative_effect_df$Exposure,
                                    levels=rev(cumulative_effect_df$Exposure))


  # make plot
  main_plot <- ggplot(data=cumulative_effect_df, aes(x=Exposure, y=Mean)) +
    geom_errorbar(aes(ymin=Lower, ymax=Upper)) +
    geom_hline(yintercept=0, color="black") +
    geom_point(color="red") +
    coord_flip() +
    labs(title="Cumulative effects",
         y="", x="")

  # return ggplot object
  return(main_plot)
}
