#' Use Voyager spatial statistics support with MCMICRO proteomics experiments
#' @import zellkonverter
#' @import SpatialFeatureExperiment
#' @import SpatialExperiment
#' @importFrom SummarizedExperiment assayNames colData
#' @import Voyager
#' @param h5ad character(1) path to an MCMICRO output in h5ad format
#' @param coordnames character() names in colData of transformed h5ad used to specify X, Y
#' @param assay4quants character(1) name in assayNames of transformed h5ad used as protein quantification
#' @param verbose logical(1) will message as stages complete if TRUE
#' @examples
#' pa = get_mcm_path("m53.1")
#' mcm53.1 = process_mcmicro(pa)
#' features_use = c("NCAM", "FOXP3", "CD8A", "LDH")
#' Voyager::plotLocalResult(mcm53.1, "localG", features = features_use,
#'                 colGeometryName = "centroids", divergent = TRUE,
#'                 diverge_center = 0)
#' @export
process_mcmicro = function(h5ad,
     coordnames = c("X_centroid", "Y_centroid"), assay4quants="X", verbose=TRUE,
     ...) {
  sce = zellkonverter::readH5AD(h5ad)
  spe = as(sce, "SpatialExperiment")
  stopifnot(coordnames %in% colnames(colData(spe)))
  spatialCoords(spe) = data.matrix(colData(spe)[, coordnames])
  sfe = as(spe, "SpatialFeatureExperiment")
  stopifnot(assay4quants %in% assayNames(sfe))
  allfeat = rownames(sfe)
  if (verbose) message("findSpatialNeighbors")
  g1 = findSpatialNeighbors(sfe) 
  colGraph(sfe, "spatNeigh") = g1 
  if (verbose) message("runUnivariate moran.plot")
  sfe = runUnivariate(sfe, type="moran.plot", features=allfeat, 
         exprs_values=assay4quants, include_self=TRUE)
  if (verbose) message("run localG")
  sfe <- runUnivariate(sfe, type = "localG", features = allfeat, exprs_values="X",
                     colGraphName = "spatNeigh", include_self = TRUE)
  sfe
}
