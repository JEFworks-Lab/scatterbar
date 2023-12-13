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
deconGexp <- results$beta
## Install & load reshape2 package
install.packages("reshape2")        
library("reshape2")
## Reshape data
data_long <- as.data.frame((deconProp))
data_long$spot <- rownames(data_long)
tidy_data <- tidyr::pivot_longer(data_long, cols = -spot, names_to = "cell_type", values_to = "proportion")
# Graph data
plot <- ggplot(tidy_data, aes(x = spot, y = proportion, fill = cell_type)) +
  geom_bar(stat = "identity") +
  labs(x = "Spot", y = "Proportion") +
  scale_fill_brewer(palette = "Set3") +
  theme_minimal()
plot

# Get unique spot values
unique_spots <- unique(tidy_data$spot)

# Create a list to store individual data frames and plots
plot_list <- list()

# Loop through each unique spot value
for (spot_value in unique_spots) {
  # Subset the data for the current spot value
  subset_data <- tidy_data[tidy_data$spot == spot_value, ]
  
  # Create a plot for the current subset of data
  plot2 <- ggplot(subset_data, aes(x = cell_type, y = proportion, fill = cell_type)) +
    geom_bar(stat = "identity") +
    labs(title = paste("Spot:", spot_value), x = "Cell Type", y = "Proportion") +
    scale_fill_brewer(palette = "Set3") +
    theme_minimal()
  
  # Store the data frame and plot in the list
  plot_list[[as.character(spot_value)]] <- list(data = subset_data, plot = plot2)
}
#Examples
plot_list[["TTTCTAACTCATAAGGAT"]][["plot"]]
plot_list[["ACAACTATGGGTTGGCGG"]][["plot"]]
plot_list[["AGGGACGTTAGTGTGCCA"]][["plot"]]
##Going back to stacked bar
tidy_data2 <- data.frame(
  spot = rep(1:260, each = 8),
  cell_type = rep(1:8, times = 260),
  proportion = tidy_data$proportion 
)
ui <- fluidPage(
  selectInput("spot_select", "Select Spot:", choices = unique(tidy_data2$spot)),
  plotOutput("spot_plot")
)

server <- function(input, output) {
  output$spot_plot <- renderPlot({
    subset_data <- subset(tidy_data2, spot == input$spot_select)
    plot <- ggplot(subset_data, aes(x = factor(cell_type), y = proportion, fill = factor(cell_type))) +
      geom_bar(stat = "identity") +
      labs(x = "Cell Type", y = "Proportion") +
      scale_fill_brewer(palette = "Set3") +
      theme_minimal()
    print(plot)
  })
}

shinyApp(ui = ui, server = server)

# we have 2 cell-types A and B, 3 spots

pos1 <- data.frame(
  x = c(0,0,1),
  y = c(0,1,1)
)
deconProp = t(data.frame(
  'spot1' = c(0,1),
  'spot2' = c(1,0),
  'spot3' = c(0.5, 0.5)
))
colnames(deconProp) = c('CT1', 'CT2')
rownames(pos1) <- rownames(deconProp)
deconProp
pos
STdeconvolve::vizAllTopics(deconProp, pos1, r=0.1) + theme_bw()


## figure out how we get here
df <- data.frame(
  x = c(0,0,0,0,1,1),
  y = c(0,0,1,1,1-0.05/2,1+0.05/2),
  z = c('A','B',
        'A','B',
        'A','B'),
  height = c(0.1,0,
             0,0.1, 
             0.05,0.05)
)
library(ggplot2)
ggplot(df, aes(x, y)) +
  geom_tile(aes(fill = z, height=height), colour = "grey50", width=0.1) +
  theme_bw()
colnames(pos)[1] <- "x"
pos$spot <- rownames(pos)
rownames(pos) <- NULL
combined_data <- merge(tidy_data, pos, by = "spot")
combined_data

create_plot <- function(combined_data) {
  p <- ggplot(combined_data, aes(x = x, y = y)) +
    geom_tile(aes(fill = cell_type, height = proportion), width = 1, color = "grey50") +
    theme_bw() +
    labs(title = "v1 alt. visual", x = "X", y = "Y")
  
  return(p)
}
#Generate the stacked bar chart
stack_plot <- create_plot(combined_data)
print(stack_plot)

#x_scale=1, y_scale=1,
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
    labs(title = "v2 alt. visual", x = "X", y = "Y")
  
  return(p)
}


#Generate the stacked bar chart
foo = combined_data[1:8,]
foo$proportion = foo$proportion
foo$cumulative_proportion = foo$cumulative_proportion
foo$x = 0
foo$y = 0
foo

foo = combined_data[1:2,]
foo$proportion = foo$proportion
foo$cumulative_proportion = foo$cumulative_proportion
foo$x = 0
foo$y = 0
foo

foo = combined_data
foo$x = foo$x * 100
foo$y = foo$y * 100
foo


# Now updated with scale factors in function
stacked_bar_chart <- create_stacked_bar_chart(combined_data,x_scale=100, y_scale=100, width=100)
print(stacked_bar_chart)


##Updated function to include reshaping data in create_stacked_bar_chart
create_stacked_bar_chart <- function(deconProp, pos, x_scale = 1, y_scale = 1, width = 1) {
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
    labs(title = "v2 alt. visual", x = "X", y = "Y")
  
  return(p)
}


