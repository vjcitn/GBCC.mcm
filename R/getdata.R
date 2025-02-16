
h5ads = c(m62.1="Galaxy69-[Single_Cell_Phenotyping_on_collection_62__1].h5ad", 
m62.2 = "Galaxy69-[Single_Cell_Phenotyping_on_collection_62__2].h5ad", 
m62.3 = "Galaxy69-[Single_Cell_Phenotyping_on_collection_62__3].h5ad", 
m62.4 = "Galaxy69-[Single_Cell_Phenotyping_on_collection_62__4].h5ad"
)


#' get path to mcmicro h5ad in package
#' @param tag character(1) one of m62.1, ..., m62.4, m62.1
#' @export
get_mcm_path = function(tag = "m62.1") {
  stopifnot(length(tag)==1)
  stopifnot(tag %in% names(h5ads))
  system.file(file.path("h5ad", h5ads[tag]), package="GBCC.mcm")
}

