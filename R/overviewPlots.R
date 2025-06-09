
 makedf = function(x, nm) {
  co = SpatialExperiment::spatialCoords(x)
  colnames(co) = c("x", "y")
  co = data.frame(co, ph=x$phenotype, name=nm)
  co
 }



#' process Breast Cancer data from MCMICRO
#' @importFrom BiocParallel register bplapply 
#' @importFrom SpatialExperiment spatialCoords
#' @param BPPARAM a BiocParallel *Param object
#' @return a list of SpatialFeatureExperiments
#' @export
processBC = function( BPPARAM = BiocParallel::MulticoreParam(4)) {
    pa3 = system.file("h5ad", "Galaxy37-LSP12022-A3.h5ad", package = "GBCC.mcm")
    pa4 = system.file("h5ad", "Galaxy37-LSP12022-A4.h5ad", package = "GBCC.mcm")
    pb3 = system.file("h5ad", "Galaxy37-LSP12022-B3.h5ad", package = "GBCC.mcm")
    pb4 = system.file("h5ad", "Galaxy37-LSP12022-B4.h5ad", package = "GBCC.mcm")
    alld = list(TripNegBC = pa3, `LumBHer2+` = pa4, `LumBHer2-` = pb3, 
        NormJej = pb4)
    BiocParallel::register(BPPARAM)
    bplapply(alld, process_mcmicro)
}

#' provide an overview plot of all breast cancer samples and cell type distributions
#' @importFrom ggplot2 ggplot geom_point facet_grid aes facet_wrap guides guide_legend
#' @importFrom BiocParallel register bplapply 
#' @importFrom SpatialExperiment spatialCoords
#' @note Returns a ggplot object
#' @examples
#' BCoverviewPlots()
#' @export
BCoverviewPlots = function (BPPARAM = BiocParallel::MulticoreParam(4)) 
{
    allpr2 = processBC()
    nms = names(allpr2)
    allsc = lapply(nms, function(x) makedf(allpr2[[x]], nm = x))
    hh = do.call(rbind, allsc)
    ggplot(hh, aes(x = x, y = y, colour = ph)) + geom_point(size = 0.2) + 
        facet_wrap(. ~ name, nrow = 2) + guides(color = guide_legend(override.aes = list(size = 5)))
}
