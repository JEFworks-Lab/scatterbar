## visualize cell-type proportions within aggregated pixels

library(SpatialExperiment)
library(SEraster)
data("merfish_mousePOA")
class(merfish_mousePOA)

## rasterize at 55um resolution
rastCt <- SEraster::rasterizeCellType(merfish_mousePOA,
                                      col_name = "celltype",
                                      resolution = 55,
                                      fun = "sum",
                                      square = TRUE)

SEraster::plotRaster(rastCt, name = "Total cells")

## calculate proportions
cellids_perpixel <- colData(rastCt)$cellID_list
ct <- merfish_mousePOA$celltype; names(ct) <- colnames(merfish_mousePOA)
ct <- as.factor(ct)
prop <- do.call(rbind, lapply(cellids_perpixel, function(x) {
  table(ct[x])/length(x)
}))
rownames(prop) <- rownames(colData(rastCt))
head(prop)

## pull out positions
pos <- spatialCoords(rastCt)
head(pos)

## filter for pseudospots with more than 1 cell
vi <- colData(rastCt)$num_cell > 1
pos <- pos[vi,]
prop <- prop[vi,]
dim(pos)
dim(prop)

library(scatterbar)
library(ggplot2)
scatterbar::scatterbar(prop, data.frame(pos), colors = sample(rainbow(length(levels(ct)))),
                              padding_x = 10, padding_y = 10) + coord_fixed()
