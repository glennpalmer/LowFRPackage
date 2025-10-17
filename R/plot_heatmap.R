#' Plot effects for pairs of times
#'
#' Plot a heatmap of posterior mean expected shifts in y when a list of chemicals varies around zero at a pair of observation times, with all other exposures held constant at zero. Regions where the 95% credible interval includes zero are shown as white.
#'
#' @export
#' @param fit Model fit returned by fit_lowFR function with default setting "all"
#' @param chems Vector of indices between 1 and p indicating which chemicals to vary
#' @param t1 Integer between 1 and TT indiating the first time at which to vary exposure
#' @param t2 Integer between 1 and TT indicating the second time at which to vary exposure
#' @param lower (Optional) Lower quantile value between 0 and 1; determines interval width for zeroing out heatmap
#' @param upper (Optional) Upper quantile value between 0 and 1; determines interval width for zeroing out heatmap
#' @param zero_out_insignificant (Optional) Indicates whether to set regions where the specified credible interval includes zero equal to zero.
#' @param grid_resolution (Optional) The resolution at which to make the heatplot
#' @param grid_min (Optional) Minimum exposure level on axes
#' @param grid_max (Optional) Maximum exposure level on axes
#'
#' @return A ggplot object plotting a heatmap of the expected change in the outcome as the specified chemicals vary away from zero, with fixed covariate values and all other chemicals fixed at 0.

#'
plot_heatmap <- function(fit, chems, t1, t2, lower=0.025, upper=0.975,
                         zero_out_insignificant=TRUE,
                         grid_resolution=0.05,
                         grid_min=-2, grid_max=2) {

  # initialize grid values
  val1 <- seq(grid_min, grid_max, grid_resolution)
  val2 <- rep(val1, each=length(val1))
  val1 <- rep(val1, length(val1))
  exp_grid <- data.frame(val1, val2)

  # initialize vectors for mean, lower bound, and upper bound
  mean_diff <- rep(NA, length(val1))
  lower_diff <- rep(NA, length(val1))
  upper_diff <- rep(NA, length(val1))

  #### determine relevant exposure indices
  num_chems <- dim(fit$Sigma)[2]
  num_times <- dim(fit$alpha)[2] / num_chems
  # first time
  indices1 <- rep(0, (num_chems*num_times))
  timeblock <- rep(0, num_times)
  timeblock[t1] <- 1
  for (chem_indices in chems) {
    indices1[((chem_indices-1)*num_times+1):(chem_indices*num_times)] <- timeblock
  }
  # second time
  indices2 <- rep(0, (num_chems*num_times))
  timeblock <- rep(0, num_times)
  timeblock[t2] <- 1
  for (chem_indices in chems) {
    indices2[((chem_indices-1)*num_times+1):(chem_indices*num_times)] <- timeblock
  }

  # populate with expected shifts in outcome
  for (i in 1:length(mean_diff)) {
    curr_diff <- compute_expected_diff(fit=fit,
                                       x1=rep(0, (num_times*num_chems)),
                                       x2=(indices1*val1[i] + indices2*val2[i]))

    mean_diff[i] <- mean(curr_diff)
    lower_diff[i] <- as.numeric(quantile(curr_diff, probs=lower))
    upper_diff[i] <- as.numeric(quantile(curr_diff, probs=upper))

    # zero out mean if credible interval includes zero
    if (zero_out_insignificant) {
      if (lower_diff[i]*upper_diff[i] < 0) {
        mean_diff[i] <- 0
      }
    }
  }

  # add summaries to dataframe
  exp_grid$mean_diff <- mean_diff
  exp_grid$lower_diff <- lower_diff
  exp_grid$upper_diff <- upper_diff

  # make axis titles
  if (length(chems)==1) {
    x_label <- paste0("Chemical ", chems, " exposure level at time ", t1)
    y_label <- paste0("Chemical ", chems, " exposure level at time ", t2)
  }
  else {
    chem_list <- chems[1]
    for (i in 2:length(chems)) {
      if (i==length(chems)) {
        chem_list <- paste0(chem_list, " and ", chems[i])
      }
      else {
        chem_list <- paste0(chem_list, ", ", chems[i])
      }
    }
    x_label <- paste0("Chemical ", chem_list, " exposure levels at time ", t1)
    y_label <- paste0("Chemical ", chem_list, " exposure levels at time ", t2)
  }

  # make graph
  surf_plot <- ggplot(data=exp_grid, aes(x = val1, y = val2, z = mean_diff, fill = mean_diff)) +
    geom_tile() +
    scale_fill_gradient2(low="blue", mid="white", high="red") +
    labs(title="Expected change in outcome relative to exposure level zero",
         subtitle="All other exposures held constant at zero",
         x=x_label,
         y=y_label)

  # return ggplot object
  return(surf_plot)
}
