#' Create a stacked bar chart with customizable parameters
#'
#' This function creates a stacked bar chart using ggplot2, where the bars are
#' stacked based on the proportions of different cell types in each (x, y) spot.
#'
#' @param deconProp A data frame containing the proportions of different cell types.
#' @param pos A data frame containing the positional information for each spot.
#' @param x_scale Scaling factor for the x-axis coordinates (default is 1).
#' @param y_scale Scaling factor for the y-axis coordinates (default is 1).
#' @param width Width of the bars (default is 1).
#'
#' @return A ggplot object representing the stacked bar chart.
#'
#' @examples
#' # Example usage:
#' # Assuming deconProp and pos are defined before this call
#' # create_stacked_bar_chart(deconProp, pos, x_scale = 1, y_scale = 1, width = 1)
#'
#' @import ggplot2
#' @import tidyr
#' @import dplyr
#'
#' @export
#' 
create_scatterbar <- function(deconProp, pos, x_scale = 1, y_scale = 1, width = 1) {
  # Reshape the data
  data_long <- as.data.frame(deconProp)
  data_long$spot <- rownames(data_long)
  tidy_data <- tidyr::pivot_longer(data_long, cols = -spot, names_to = "cell_type", values_to = "proportion")
  colnames(pos)[1] <- "x"
  pos$spot <- rownames(pos)
  rownames(pos) <- NULL
  combined_data <- merge(tidy_data, pos, by = "spot")

  # Calculate cumulative proportions for each (x, y) spot
  combined_data <- combined_data %>%
    group_by(x, y) %>%
    # Ensures that the heights of the bars within a spot add up to 1
    mutate(cumulative_proportion = cumsum(proportion) - proportion)

  # Apply scaling factors to x and y coordinates
  scaled_combined_data <- combined_data %>%
    mutate(x = x * x_scale, y = y * y_scale)

  # Correct positioning of bars within each (x, y) spot
  p <- ggplot(scaled_combined_data, aes(x = x, y = y + cumulative_proportion + proportion/2)) +
    geom_tile(aes(fill = cell_type, height = proportion), width = width, lwd = 0) +
    theme_bw() +
    labs(title = "visual", x = "X", y = "Y")

  return(p)
}
