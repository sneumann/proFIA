#' Process FIA-HRMS datasets.
#'
#' Process FIA-HRMS datasets passing from raw data (mzMl, CDF,mzXML
#' format) to a peak table suitable for statistical anlysis.
#' 
#' The full workflow is composed of the following chain of function.
#' \code{proFIAset}=>\code{group.FIA}=>\code{fillPeaks.WKNN} and the full process is easy to automate using the
#' analyzeAcquisitionFIA which do all the steps easely. Resulting table may be easely 
#' exported using the 3 exports function (\code{exportDataMatrix},\code{exportVariableMetadata}\code{exportSampleMetadata})
#' Groups may be visualised using \code{plotEICs} to plot all the EICs of a group.
#' 
#' @docType package
#' @importFrom minpack.lm nls.lm nls.lm.control
#' @importFrom FNN get.knn
#' @importFrom maxLik maxLik
#' @importFrom pracma pinv
#' @importFrom grDevices colorRampPalette rainbow as.raster col2rgb rgb
#' @importFrom graphics abline legend lines matplot plot
#'     points polygon layout rasterImage text
#' @importFrom methods new
#' @importFrom stats aggregate coef convolve cor density
#'          dnorm lm loess median nextn pnorm sd
#'          predict quantile t.test
#' @importFrom utils combn object.size write.table
#' @importFrom xcms xcmsRaw plotTIC rawEIC plotRaw peaks
#' @importFrom Biobase AnnotatedDataFrame ExpressionSet
#' @importFrom BiocParallel bplapply bpparam
#' @name proFIA-package
#' @aliases proFIA proFIA-package
NULL
