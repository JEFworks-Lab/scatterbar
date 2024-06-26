#' Create a scattered stacked bar chart plot
#'
#' @description This function creates a scatterbar plot using ggplot2, where the bars are stacked based on the different proportions of groups in each 2-D location/spot. A scatterbar plot is a combination of a scatter plot and a stacked bar chart, allowing for the visualization of proportional data across spatial coordinates.
#' The function allows for customized scaling factors and padding when creating the plot. If no scaling factors are specified, the function automatically determines the optimal scaling factors based on the data.
#'
#' @param proportion_data A data frame containing the proportions of different categories for each location. Each row represents a location, and each column (except the row names) represents a category.
#' @param position_data  A data frame containing the positional information for each location. This data frame includes the x and y coordinates for each location/spot (the respective row names).
#' @param x_scale X-axis scaling factor (default is NULL). If not provided, it will be automatically calculated based on the data.
#' @param y_scale Y-axis scaling factor (default is NULL). If not provided, it will be automatically calculated based on the data.
#' @param padding_x Padding for x-axis (default is 0).
#' @param padding_y Padding for y-axis (default is 0).
#' @param show_legend Boolean indicating whether to display the plot legend (default is TRUE).
#' @param plot_title Title for plot (default is NULL).
#' @param legend_title Custom title for the legend (default is "Group").
#' @param colors Optional vector of colors to use for each category (default is NULL). If not provided, a default palette will be used.
#'
#' @return A ggplot object representing the scattered stacked bar chart plot.
#'
#' @examples
#' data(deconProp)
#' data(pos)
#' create_scatterbar(deconProp, pos, padding_x = 0.3, padding_y = 0.3, legend_title = "Cell Types")
#'
#' @import ggplot2
#' @import tidyr
#' @import dplyr
#'
#' @export
create_scatterbar <- function(proportion_data, position_data, x_scale = NULL, y_scale = NULL, padding_x=0, padding_y=0, show_legend = TRUE, plot_title=NULL, legend_title="Group", colors = NULL) {

  #Check that proportion_data and position_data are either data.frames or matrices
  if( !is.matrix(proportion_data) & !is.data.frame(proportion_data) ){
    stop("`proportion_data` must be a matrix or data.frame.")
  }
  if( !is.matrix(position_data) & !is.data.frame(position_data)){
    stop("`position_data` must be a matrix or data.frame with 2 columns named `x` and `y` and the name of the location/spot.")
  }

  #Check that position_data has exactly 2 columns: x and y
  if( (any(!colnames(pos) %in% c("x", "y")) == TRUE) | (dim(pos)[2] != 2) ){
    stop("`position_data` must have exactly 2 columns named `x` and `y`.")
  }

  #Reshape the data
  data_long <- as.data.frame(proportion_data)
  data_long$spot <- rownames(data_long)
  tidy_data <- tidyr::pivot_longer(data_long, cols = -spot, names_to = "Group", values_to = "proportion")

  #Prepare position_data
  colnames(position_data)[1] <- "x"
  position_data$spot <- rownames(position_data)
  rownames(position_data) <- NULL

  #Filter spots to include only those that exist in both data frames
  common_spots <- intersect(tidy_data$spot, position_data$spot)
  tidy_data <- tidy_data[tidy_data$spot %in% common_spots, ]
  position_data <- position_data[position_data$spot %in% common_spots, ]

  #Merge the data
  combined_data <- merge(tidy_data, position_data, by = "spot")

  #Calculate cumulative proportions for each (x, y) spot
  combined_data <- combined_data %>%
    group_by(x, y) %>%
    # Ensures that the heights of the bars within a spot add up to 1
    dplyr::mutate(cumulative_proportion = cumsum(proportion) - proportion)

  #Determine optimal scaling factors
  if (is.null(x_scale) && is.null(y_scale)) {
  x_range <- range(combined_data$x)
  y_range <- range(combined_data$y)
  x_dist <- (x_range[2] - x_range[1])
  y_dist <- (y_range[2] - y_range[1])
  sq_num_spots <- sqrt(nrow(pos))
  x_scale <-x_dist/sq_num_spots
  y_scale <- y_dist/sq_num_spots

  } else if (is.null(x_scale)){
    x_range <- range(combined_data$x)
    x_dist <- (x_range[2] - x_range[1])
    sq_num_spots <- sqrt(nrow(pos))
    x_scale <-x_dist/sq_num_spots
  } else if(is.null(y_scale)){
    y_range <- range(combined_data$y)
    y_dist <- (y_range[2] - y_range[1])
    sq_num_spots <- sqrt(nrow(pos))
    y_scale <-y_dist/sq_num_spots
  }

  #Apply padding
  x_scale <- x_scale - padding_x
  y_scale <- y_scale - padding_y

  #Plot scatterbar, correcting the position of the bars within each (x, y) spot as they are plotted
  p <- ggplot2::ggplot(combined_data, aes(x = x, y = y + cumulative_proportion*y_scale + proportion*y_scale/2)) +
    geom_tile(aes(fill = Group, height = proportion*y_scale), width = x_scale, lwd = 0) + ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      axis.line = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_blank(),
      axis.text.y = ggplot2::element_blank(),
      axis.ticks = ggplot2::element_blank(),
      axis.title.x = ggplot2::element_blank(),
      axis.title.y = ggplot2::element_blank(),
      panel.background = ggplot2::element_blank(),
      plot.background = ggplot2::element_blank(),
      legend.text = ggplot2::element_text(size = 12, colour = "black"),
      legend.title = ggplot2::element_text(size = 12, colour = "black")
    ) + ggplot2::labs(fill = legend_title)


  #Color of bars customization
  if (is.null(colors)) {
    p <- p + ggplot2::scale_fill_brewer(palette = "Spectral")
  } else {
    p <- p + ggplot2::scale_fill_manual(values = colors)
  }

  #Legend title customization
  if (!show_legend) {
    p <- p + ggplot2::theme(legend.position = "none")
  }

  #Plot title customization
  if (!is.null(plot_title)) {
    p <- p + ggplot2::ggtitle(plot_title)
  }

  return(p)
}
