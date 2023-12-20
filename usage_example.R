## modifying dee's function
create_stacked_bar_chart <- function(deconProp, pos, scale = 1, width = 1, ...) {
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
    mutate(x = x * scale, y = y * scale)
  
  # Correct positioning of bars within each (x, y) spot
  p <- ggplot(scaled_combined_data, aes(x = x, y = y + cumulative_proportion + proportion/2)) +
    geom_tile(aes(fill = cell_type, height = proportion), width = width, ...)
  
  return(p)
}

######### https://jef.works/blog/2022/05/03/deconvolution-vs-clustering/
setwd('~/OneDrive - Johns Hopkins/Data_Private/MERFISH_cortex_For_Daniel/Visium2/')
## read in pixel position information
pos.info <- read.csv('adult-mouse-brain-ffpe-1-standard-1-3-0/spatial/tissue_positions_list.csv', header=FALSE)
pos <- pos.info[,c(6,5)]
rownames(pos) <- pos.info[,1]
colnames(pos) <- c('x', 'y')
par(mfrow=c(1,1))
plot(pos)

## read in gene expression counts
cd <- Matrix::readMM('adult-mouse-brain-ffpe-1-standard-1-3-0/filtered_feature_bc_matrix/matrix.mtx.gz')
barcodes <- read.csv('adult-mouse-brain-ffpe-1-standard-1-3-0/filtered_feature_bc_matrix/barcodes.tsv.gz', sep='\t', header=FALSE)
features <- read.csv('adult-mouse-brain-ffpe-1-standard-1-3-0/filtered_feature_bc_matrix/features.tsv.gz', sep='\t', header=FALSE)
rownames(cd) <- features[,2]
colnames(cd) <- barcodes[,1]
cd[1:5,1:5]

## restrict to only pixels with genes
## and flip right side up
pos <- pos[colnames(cd),]
pos[,2] <- -pos[,2]

## visualize depth
par(mfrow=c(1,1), mar=rep(1,4))
MERINGUE::plotEmbedding(pos, col=colSums(cd), cex=1)

library(STdeconvolve)
## remove pixels with too few genes
counts <- STdeconvolve::cleanCounts(cd, 
                                    min.lib.size = 1, 
                                    max.lib.size = Inf,
                                    min.reads = 2000, 
                                    verbose = TRUE)

## feature select for genes
corpus <- STdeconvolve::restrictCorpus(counts, 
                                       removeAbove=1.0, 
                                       removeBelow = 0.05,
                                       alpha = 0.01,
                                       nTopOD = NA)
dim(corpus)

## choose number of cell types
ldas <- STdeconvolve::fitLDA(t(as.matrix(corpus)), Ks = c(12))

## get best model results
optLDA <- STdeconvolve::optimalModel(models = ldas, opt = 12)

## extract deconvolved cell type proportions (theta) and transcriptional profiles (beta)
results <- STdeconvolve::getBetaTheta(optLDA, 
                                      perc.filt = 0.05, 
                                      betaScale = 1000)
deconProp <- results$theta
deconGexp <- results$beta

cols <- rainbow(12)
names(cols) <- 1:12
create_stacked_bar_chart(deconProp, pos, width=1, scale = 0.005, color='grey') + scale_fill_manual(values=cols) + theme_void()
STdeconvolve::vizAllTopics(deconProp, pos, lwd=0.1, r=100)

## issues: how to automate height to be same as width instead of just hacking on scale?
