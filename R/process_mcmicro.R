#' create a default configuration for several runUnivariate tasks
#' @param assay4quants character(1) name of the assay component of SFE on which runUnivariate is used
#' @param colGraphName character(1) name of the colGraph in the SFE
#' @param BPPARAM a BiocParallel *Param object
#' @note Any values supplied to dots in the returned closure must be present for all components for all
#' desired runUnivariate calls.  
#' @return A closure.  The intent is that the value is passed a vector of
#' feature names.  The resulting list has named elements.  The name of
#' each element is the 'type' of univariate statistical analysis, and the
#' elements of the associated element are arguments to runUnivariate used with that
#' analysis type.
#' @param assay4quants character(1) name of assay component to use for univariate analysis,
#; defaults to "X", which comes from ingesting an h5ad 
#' @param colGraphName character(1) defaults to 'spatNeigh' which is the name used within
#' `process_mcmicro` for the colGraph production
#' @param BPPARAM defaults to the `BiocParallel::bpparam()` value for the session at
#' the time this function is called
#' @export
default_uruns = function(assay4quants="X", colGraphName="spatNeigh",
    BPPARAM=BiocParallel::bpparam()) function(features, ...) {
  list(moran.plot=list( features=features, exprs_values=assay4quants, include_self=TRUE,
        BPPARAM=BPPARAM, ...),
       localG=list(features=features, exprs_values=assay4quants, 
                colGraphName=colGraphName, include_self=TRUE, BPPARAM=BPPARAM, ...))
}

#' Use Voyager spatial statistics support with MCMICRO proteomics experiments
#' @import zellkonverter
#' @import SpatialFeatureExperiment
#' @import SpatialExperiment
#' @importFrom SummarizedExperiment assayNames colData
#' @import Voyager
#' @import methods
#' @param h5ad character(1) path to an MCMICRO output in h5ad format
#' @param uconfig this is a closure that builds runUnivariate calls with all available features,
#' see `default_uruns`
#' @param coordnames character() names in colData of transformed h5ad used to specify X, Y
#' @param spneigh_parms list() of arguments (other than x) to Voyager::findSpatialNeighbors, defaults to list()
#' @param assay4quants character(1) name in assayNames of transformed h5ad used as protein quantification
#' @param verbose logical(1) will message as stages complete if TRUE
#' @note Voyager::findSpatialNeighbors uses methods in the spdep package to define neighbors.  In this
#' function, we compute a `colGraph` component, which assumes MARGIN value is 2.  An example
#' spneigh_parms setting is `list(method='knearneigh', k=4)`.
#' @examples
#' pa = get_mcm_path("m53.1")
#' # by default, moran.plot and localG type univariate analyses are conducted
#' mcm53.1 = process_mcmicro(pa)
#' features_use = c("NCAM", "FOXP3", "CD8A", "LDH")
#' Voyager::plotLocalResult(mcm53.1, "localG", features = features_use,
#'                 colGeometryName = "centroids", divergent = TRUE,
#'                 diverge_center = 0, show_axes=TRUE)
#' mcm53.1b = process_mcmicro(pa, spneigh_parms=list(method="knearneigh", k=3))
#' Voyager:::plotColGraph(mcm53.1b)
#' Voyager::plotLocalResult(mcm53.1b, "localG", features = features_use,
#'                 colGeometryName = "centroids", divergent = TRUE,
#'                 diverge_center = 0)
#' @export
process_mcmicro = function(h5ad, uconfig=default_uruns(),
     coordnames = c("X_centroid", "Y_centroid"), spneigh_parms=list(), assay4quants="X", verbose=TRUE) {
  sce = zellkonverter::readH5AD(h5ad)
  spe = as(sce, "SpatialExperiment")
  stopifnot(coordnames %in% colnames(colData(spe)))
  spatialCoords(spe) = data.matrix(colData(spe)[, coordnames])
  sfe = as(spe, "SpatialFeatureExperiment")
  stopifnot(assay4quants %in% assayNames(sfe))
  allfeat = rownames(sfe)
  if (verbose) message("findSpatialNeighbors")
  fsnargs = list(x=sfe)
  g1 = do.call(findSpatialNeighbors, c(fsnargs, spneigh_parms))
  colGraph(sfe, "spatNeigh") = g1 
  uconf = uconfig(rownames(spe))
  for (i in names(uconf)) uconf[[i]]$type = i
  for (i in names(uconf)) {
        if (verbose) message(sprintf("runUnivariate %s", i))
        setup = unlist(list(x = sfe, uconf[[i]]), recursive = FALSE)
        sfe = do.call("runUnivariate", setup)
  }
  sfe
}


get_xylims = function(sfe) {
 }
