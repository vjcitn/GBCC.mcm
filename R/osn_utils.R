
osn_url = "https://mghp.osn.xsede.org/bir190004-bucket01/BiocGBCCmcm"

mcm.sfe.zips = c( "m53.1.sfe.zip", "m53.2.sfe.zip", 
   "m53.3.sfe.zip", "m53.4.sfe.zip", "m62.1.sfe.zip")

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

#' retrieve a specified zipped experiment from cache, downloading to cache if necessary
#' @param experiment character(1)
#' @export
path_to_zipped_mcm = function(experiment = "m53.1.sfe.zip", local=FALSE,
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
  

