# Function to summarize mean and quantiles for a vector of posterior samples

#' Summarize a vector of posterior samples in terms of the posterior mean and credible interval
#'
#' @export
#' @param samples Vector of posterior samples of some quantity of interest
#' @param lower Lower quantile value between 0 and 1
#' @param upper Upper quantile value between 0 and 1
#'
#' @return A list containing the mean and lower and upper quantiles (default 2.5th and 97.5th percentile, giving an equal-tailed 95% interval)

#'
summarize_posterior <- function(samples, lower=0.025, upper=0.975) {
  post_mean <- mean(samples)
  lower_quantile <- quantile(samples, probs=c(lower))
  upper_quantile <- quantile(samples, probs=c(upper))

  output_list <- list(post_mean, lower_quantile, upper_quantile)
  names(output_list) <- c("mean", "lower_quantile", "upper_quantile")

  return(output_list)
}
