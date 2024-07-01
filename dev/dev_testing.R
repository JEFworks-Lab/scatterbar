## Run STdeconvolve
library(STdeconvolve)
## Load built in data
data(mOB)
pos <- mOB$pos
cd <- mOB$counts
## Remove pixels with too few genes
counts <- cleanCounts(cd, min.lib.size = 100)
## Feature select for genes
corpus <- restrictCorpus(counts, removeAbove=1.0, removeBelow = 0.05)
## Choose optimal number of cell-types
ldas <- fitLDA(t(as.matrix(corpus)), Ks = c(8))
## Get best model results
optLDA <- optimalModel(models = ldas, opt = "min")
## Extract deconvolved cell-type proportions (theta) and transcriptional profiles (beta)
results <- getBetaTheta(optLDA, perc.filt = 0.05, betaScale = 1000)
deconProp <- results$theta


#look into later....
colnames(pos)[1] <- "x"
pos$spot <- rownames(pos)
rownames(pos) <- NULL

combined_data <- merge(tidy_data, pos, by = "spot")

create_stacked_bar_chart <- function(combined_data,x_scale=1, y_scale=1, width=1) {
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
    geom_tile(aes(fill = cell_type, height = proportion), width = width, lwd=0) +
    theme_bw() +
    labs(title = "visual", x = "X", y = "Y")

  return(p)
}

# Now updated with scale factors in function
stacked_bar_chart <- create_stacked_bar_chart(combined_data,x_scale=1, y_scale=1, width=1)
print(stacked_bar_chart)

scatterbar::create_scatterbar(deconProp, pos, padding_x = 0.3, padding_y=0.3)

spectral_colors <- brewer.pal(n = 11, name = "Spectral")
print(spectral_colors)


p2 <- ggplot2::ggplot() +
  ggplot2::theme(
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
  ) +
  scatterpie::geom_scatterpie(ggplot2::aes(x=x, y=y, group=spot, r=max(0.4, max(pos)/nrow(pos)*4)),
                              lwd = 0, data = deconProp, cols = my_colors, legend_name = "Groups")
p2


postest <- pos
postest$x <- round(postest$x)
postest$y <- round(postest$y)


p2 <- ggplot2::ggplot() +
  ggplot2::theme(
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
  ) +
  scatterpie::geom_scatterpie(ggplot2::aes(x=x, y=y, group=s, r=max(0.4, max(postest)/nrow(postest)*4), colors = my_colors),
                              lwd = 0.01, data = deconProp, legend_name = "Groups")
p2

vizAllTopics(deconProp, postest, topicCols = my_colors, plotTitle="Rounded")
print(spectral_colors)
print(my_colors)
my_colors <- c("#D53E4F", "#F46D43", "#FDAE61", "#FEE08B","#E6F598", "#ABDDA4", "#66C2A5", "#3288BD")
