## read in pixel position information
pos.info <- read.csv('/Users/jeanfan/Desktop/blog/deconvolve_vs_clustering/adult-mouse-brain-ffpe-1-standard-1-3-0/spatial/tissue_positions_list.csv', header=FALSE)
pos <- pos.info[,c(6,5)]
rownames(pos) <- pos.info[,1]
colnames(pos) <- c('x', 'y')
plot(pos)

## read in gene expression counts
cd <- Matrix::readMM('/Users/jeanfan/Desktop/blog/deconvolve_vs_clustering/adult-mouse-brain-ffpe-1-standard-1-3-0/filtered_feature_bc_matrix/matrix.mtx.gz')
barcodes <- read.csv('/Users/jeanfan/Desktop/blog/deconvolve_vs_clustering/adult-mouse-brain-ffpe-1-standard-1-3-0/filtered_feature_bc_matrix/barcodes.tsv.gz', sep='\t', header=FALSE)
features <- read.csv('/Users/jeanfan/Desktop/blog/deconvolve_vs_clustering/adult-mouse-brain-ffpe-1-standard-1-3-0/filtered_feature_bc_matrix/features.tsv.gz', sep='\t', header=FALSE)
rownames(cd) <- features[,2]
colnames(cd) <- barcodes[,1]
cd[1:5,1:5]

## restrict to only pixels with genes
## and flip right side up
pos <- pos[colnames(cd),]
pos[,2] <- -pos[,2]

## visualize depth
library(Matrix)
par(mfrow=c(1,1), mar=rep(1,4))
MERINGUE::plotEmbedding(pos, col=colSums(cd), cex=1)

## stdeconvolve
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


##### scatterbar
library(scatterbar)
library(ggplot2)
# Basic scatterbar plot with default settings
scatterbar(deconProp, pos,
                  size_x = 220, size_y = 220) + ggplot2::coord_fixed()

adult_mouse_brain_ffpe <- list(pos=pos, prop=deconProp)
usethis::use_data(adult_mouse_brain_ffpe)
