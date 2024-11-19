#' Deconvolved cell-type proportions from STdeconvolve for a Spatial Transcriptomic dataset of the mouse olfactory bulb along with the positions of their respective spots
#'
#' @format A list of 2 objects, each of 260 rows and 2 columns: data, which contains the 8 deconvolved cell-types proportions for each 260 spots/pixels in a spatial transcriptomics experiment performed on mouse olfactory bulb tissue and xy, which contains the x and y-coordinates for each 260 spots/pixels on a mouse olfactory bulb tissue slide.
#'
#' @source \url{https://www.science.org/doi/10.1126/science.aaf2403}
#'
#' @usage data(mOB)
"mOB"

#' Deconvolved cell-type proportions from STdeconvolve for a Visium data from the adult mouse brain along with the positions of their respective spots
#'
#' @format A list of 2 objects, each of 2264 rows and 2 columns: prop, which contains the 12 deconvolved cell-types proportions for each 2264 spots/pixels in a spatial transcriptomics experiment performed on adult mouse brain tissue tissue and pos, which contains the x and y-coordinates for each 2264 spots/pixels on an adult mouse brain tissue slide.
#'
#' @source \url{https://www.10xgenomics.com/datasets/adult-mouse-brain-ffpe-1-standard-1-3-0}
#'
#' @usage data(adult_mouse_brain_ffpe)
"adult_mouse_brain_ffpe"
