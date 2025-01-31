h5ads = c(m53.1="Galaxy62-[Convert_McMicro_Output_to_Anndata_on_collection_53__1].h5ad",
  m53.2="Galaxy62-[Convert_McMicro_Output_to_Anndata_on_collection_53__2].h5ad",
  m53.3="Galaxy62-[Convert_McMicro_Output_to_Anndata_on_collection_53__3].h5ad",
  m53.4="Galaxy62-[Convert_McMicro_Output_to_Anndata_on_collection_53__4].h5ad",
  m62.1="Galaxy69-[Single_Cell_Phenotyping_on_collection_62__1].h5ad")

#' get path to mcmicro h5ad in package
#' @param tag character(1) one of m53.1, ..., m53.4, m62.1
#' @export
get_mcm_path = function(tag = "m53.1") {
  stopifnot(length(tag)==1)
  stopifnot(tag %in% names(h5ads))
  system.file(file.path("h5ad", h5ads[tag]), package="GBCC.mcm")
}

