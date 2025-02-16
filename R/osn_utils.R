
osn_url = "https://mghp.osn.xsede.org/bir190004-bucket01/BiocGBCCmcm"

mcm.sfe.zips = c( "m62.1.sfe.zip", "m62.2.sfe.zip", 
   "m62.3.sfe.zip", "m62.4.sfe.zip")

# assumes artifacts have sfe.zip as suffix

#' list available mcmicro spatial proteomics for GBCC
#' @param local logical(1), if TRUE, assume available at Sys.getenv("MCM_LOCAL")
#' @return vector of paths
#' @export
available_mcm = function(local=FALSE) {
 if (local) {
  location = Sys.getenv("MCM_LOCAL")
  if (nchar(location)==0) stop("environment variable MCM_LOCAL is empty")
  return(dir(location, full.names=TRUE))
  }
 return(paste(osn_url, "/", mcm.sfe.zips, sep=""))
}

#' retrieve a specified zipped takane-formatted SpatialFeatureExperiment from cache, downloading to cache if necessary
#' @import BiocFileCache
#' @note alabaster.sfe's saveObject was used before zipping these results of `process_mcmicro` for
#' datasets supplied by Jeremy Goecks.
#' @param experiment character(1)
#' @param local logical(1), if TRUE, assume available at Sys.getenv("MCM_LOCAL")
#' @param cache like BiocFileCache() value
#' @examples
#' pa = path_to_zipped_mcm()
#' tf = tempfile()
#' unzip(pa, exdir=tf)
#' requireNamespace("alabaster.sfe")
#' dem = alabaster.base::readObject(file.path(tf, "m62.1.sfe"))
#' dem
#' table(dem$phenotype)
#' @export
path_to_zipped_mcm = function(experiment = "m62.1.sfe.zip", local=FALSE,
     cache = BiocFileCache::BiocFileCache()) {
# check name is available
  avail = available_mcm(local = local)
  stopifnot(length(grep(experiment, avail))==1)
# check if in cache
  chk = BiocFileCache::bfcquery(cache, experiment)
  ind = 1
  if (nrow(chk)>1) {
   message("multiple hits for selected experiment, using last")
   ind = nrow(chk)
   }
  if (nrow(chk)>0) return(chk$rpath[ind])
# not in cache; retrieve, cache, return path
  touse = grep(experiment, avail, value=TRUE)
  stopifnot(length(touse)==1)
  pa = BiocFileCache::bfcadd(cache, experiment, fpath=touse)
  return(pa)
}
  

