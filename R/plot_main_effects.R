#' Plot main effects for all exposures
#'
#' Plot posterior means and credible intervals for the change in expected outcome if each exposure is shifted from -0.5 to 0.5, with all other exposures held constant at 0.
#'
#' @export
#' @param fit Model fit returned by fit_lowFR function with default setting "all"
#' @param varnames p*TT dimensional vector of names of exposures
#' @param lower (Optional) Lower quantile value between 0 and 1
#' @param upper (Optional) Upper quantile value between 0 and 1
#'
#' @return A ggplot object plotting main effect posterior means and credible intervals.

#'
plot_main_effects <- function(fit, varnames, lower=0.025, upper=0.975) {
  # compute main effect dataframe
  main_effect_df <- summarize_main_effects(fit=fit, varnames=varnames, lower=lower, upper=upper)
  main_effect_df$Exposure <- factor(main_effect_df$Exposure,
                                           levels=rev(main_effect_df$Exposure))


  # make plot
  main_plot <- ggplot(data=main_effect_df, aes(x=Exposure, y=Mean)) +
    geom_errorbar(aes(ymin=Lower, ymax=Upper)) +
    geom_hline(yintercept=0, color="black") +
    geom_point(color="red") +
    coord_flip() +
    labs(title="Main effects (LowFR)",
         y="", x="")

  # return ggplot object
  return(main_plot)
}
