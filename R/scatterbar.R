# Define global variables
utils::globalVariables(c("x", "y", "spot", "proportion", "cumulative_proportion", "Group"))

#' Create a scattered stacked bar chart plot
#'
#' @description This function creates a scatterbar plot using ggplot2, where the bars are stacked based on the different proportions of groups in each 2-D location/spot. A scatterbar plot is a combination of a scatter plot and a stacked bar chart, allowing for the visualization of proportional data across spatial coordinates.
#' The function allows for customized scaling factors and padding when creating the plot. If no scaling factors are specified, the function automatically determines the optimal scaling factors based on the data.
#'
#' @param data A data frame containing the proportions of different categories for each location. Each row represents a location, and each column (except the row names) represents a category.
#' @param xy  A data frame containing the positional information for each location. This data frame includes the x and y coordinates for each location/spot (the respective row names).
#' @param size_x X-axis scaling factor (default is NULL). If not provided, it will be automatically calculated based on the data.
#' @param size_y Y-axis scaling factor (default is NULL). If not provided, it will be automatically calculated based on the data.
#' @param padding_x Padding for x-axis (default is 0).
#' @param padding_y Padding for y-axis (default is 0).
#' @param show_legend Boolean indicating whether to display the plot legend (default is TRUE).
#' @param legend_title Custom title for the legend (default is "Group").
#' @param colors Optional vector of colors to use for each category (default is NULL). If not provided, a default palette will be used.
#'
#' @return A ggplot object representing the scattered stacked bar chart plot.
#'
#' @examples
#' data(mOB)
#' scatterbar(mOB$data, mOB$xy, padding_x = 0.3, padding_y = 0.3, legend_title = "Cell Types")
#'
#' @importFrom grDevices rainbow
#' @importFrom magrittr %>%
#' @import ggplot2
#' @import tidyr
#' @import dplyr
#'
#' @export
scatterbar <- function(data, xy, size_x = NULL, size_y = NULL, padding_x=0, padding_y=0, show_legend = TRUE, legend_title="Group", colors = NULL) {

  #Check that data and xy are either data.frames or matrices
  if( !is.matrix(data) & !is.data.frame(data) ){
    stop("`data` must be a matrix or data.frame.")
  }
  if( !is.matrix(xy) & !is.data.frame(xy)){
    stop("`xy` must be a matrix or data.frame with 2 columns named `x` and `y` and the name of the location/spot.")
  }

  #Check that xy has exactly 2 columns: x and y
  if( (any(!colnames(xy) %in% c("x", "y")) == TRUE) | (dim(xy)[2] != 2) ){
    stop("`xy` must have exactly 2 columns named `x` and `y`.")
  }

  #Reshape the data
  data_long <- as.data.frame(data)
  data_long$spot <- rownames(data_long)
  tidy_data <- tidyr::pivot_longer(data_long, cols = -spot, names_to = "Group", values_to = "proportion")

  #Set the order of the levels of Group based on the column order of data
  tidy_data$Group <- factor(tidy_data$Group, levels = colnames(data))

  #Prepare xy
  colnames(xy)[1] <- "x"
  xy$spot <- rownames(xy)
  rownames(xy) <- NULL

  #Filter spots to include only those that exist in both data frames
  common_spots <- intersect(tidy_data$spot, xy$spot)
  tidy_data <- tidy_data[tidy_data$spot %in% common_spots, ]
  xy <- xy[xy$spot %in% common_spots, ]

  #Merge the data
  combined_data <- merge(tidy_data, xy, by = "spot")

  #Calculate cumulative proportions for each (x, y) spot
  combined_data <- combined_data %>%
    dplyr::group_by(.data$x, .data$y) %>%
    dplyr::arrange(Group) %>%
    # Ensures that the heights of the bars within a spot add up to 1
    dplyr::mutate(cumulative_proportion = cumsum(proportion) - proportion)

  #Determine optimal scaling factors
  if (is.null(size_x) && is.null(size_y)) {
    x_range <- range(combined_data$x)
    y_range <- range(combined_data$y)
    x_dist <- (x_range[2] - x_range[1])
    y_dist <- (y_range[2] - y_range[1])
    sq_num_spots <- sqrt(nrow(xy))
    size_x <-x_dist/sq_num_spots
    size_y <- y_dist/sq_num_spots

  } else if (is.null(size_x)){
    x_range <- range(combined_data$x)
    x_dist <- (x_range[2] - x_range[1])
    sq_num_spots <- sqrt(nrow(xy))
    size_x <-x_dist/sq_num_spots
  } else if(is.null(size_y)){
    y_range <- range(combined_data$y)
    y_dist <- (y_range[2] - y_range[1])
    sq_num_spots <- sqrt(nrow(xy))
    size_y <-y_dist/sq_num_spots
  }

  #Apply padding
  size_x <- size_x - padding_x
  size_y <- size_y - padding_y

  #Plot scatterbar, correcting the position of the bars within each (x, y) spot as they are plotted
  p <- ggplot2::ggplot(combined_data, ggplot2::aes(x = x, y = y - size_y/2 + cumulative_proportion*size_y + proportion*size_y/2)) +
    ggplot2::geom_tile(ggplot2::aes(fill = Group, height = proportion*size_y), width = size_x, lwd = 0) + ggplot2::theme(
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
    p <- p + ggplot2::scale_fill_manual(values = rev(rainbow(ncol(data), s=0.8)), breaks=rev(colnames(data)))
  } else {
    p <- p + ggplot2::scale_fill_manual(values = colors)
  }

  #Legend title customization
  if (!show_legend) {
    p <- p + ggplot2::theme(legend.position = "none")
  }

  return(p)
}
