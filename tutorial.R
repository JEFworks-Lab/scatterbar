## Tutorial for Dee

## Goal 1: run STdeconvolve
library(STdeconvolve)
## load built in data
data(mOB)
pos <- mOB$pos
cd <- mOB$counts
## remove pixels with too few genes
counts <- cleanCounts(cd, min.lib.size = 100)
## feature select for genes
corpus <- restrictCorpus(counts, removeAbove=1.0, removeBelow = 0.05)
## choose optimal number of cell-types
ldas <- fitLDA(t(as.matrix(corpus)), Ks = c(8))
## get best model results
optLDA <- optimalModel(models = ldas, opt = "min")
## extract deconvolved cell-type proportions (theta) and transcriptional profiles (beta)
results <- getBetaTheta(optLDA, perc.filt = 0.05, betaScale = 1000)
deconProp <- results$theta
deconGexp <- results$beta
## visualize deconvolved cell-type proportions
vizAllTopics(deconProp, pos)

## Goal 2: read the code for the current vizAllTopics() function
vizAllTopics

## Goal 3: given the same input (deconProp and pos for now)
## design and implement an alternative visualization using ggplot2
