
h5ads.old = c(m62.1="Galaxy69-[Single_Cell_Phenotyping_on_collection_62__1].h5ad", 
m62.2 = "Galaxy69-[Single_Cell_Phenotyping_on_collection_62__2].h5ad", 
m62.3 = "Galaxy69-[Single_Cell_Phenotyping_on_collection_62__3].h5ad", 
m62.4 = "Galaxy69-[Single_Cell_Phenotyping_on_collection_62__4].h5ad"
)

h5ads = c(B4="Galaxy37-LSP12022-B4.h5ad", B3="Galaxy37-LSP12022-B3.h5ad",
  A4="Galaxy37-LSP12022-A4.h5ad", A3="Galaxy37-LSP12022-A3.h5ad")


#' get path to mcmicro h5ad in package
#' @param tag character(1) one of A3, A4, B3, B4
#' @note A3 = Triple Negative Breast Cancer; A4 = Luminal B HER2+ Breast Cancer; 
#' B3 = Luminal B HER- Breast Cancer; B4 = normal jejunum
#' @export
get_mcm_path = function(tag = "A4") {
  stopifnot(length(tag)==1)
  stopifnot(tag %in% names(h5ads))
  system.file(file.path("h5ad", h5ads[tag]), package="GBCC.mcm")
}

