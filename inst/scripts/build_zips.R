library(GBCC.mcm)
library(alabaster.sfe)
library(Voyager)

target = "/Users/vincentcarey/GBCC_ZIPS"
dir.create(target)

nms = names(GBCC.mcm:::h5ads)

for (i in nms) {
  assign(i, process_mcmicro(get_mcm_path(i)))
  saveObject(get(i), file.path(target, paste0(i, ".sfe")))
  }
