### generating data for package
## see https://r-pkgs.org/data.html

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

library(tidyr)
library(dplyr)
library(ggplot2)
pos = round(pos)
create_scatterbar(deconProp, pos)

usethis::use_data(deconProp)
usethis::use_data(pos)

