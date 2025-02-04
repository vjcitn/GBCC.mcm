#fol = "/Users/vincentcarey/ANVIL/MCMICRO/GBCC_ZIPS"
fol = "."
aa = c("m53.1.sfe", "m53.2.sfe", "m53.3.sfe", "m53.4.sfe", "m62.1.sfe")
fis = paste(fol, aa, sep="/")
names(fis) = c("m53.1", "m53.2", "m53.3", "m53.4", "m62.1")

targets = c("AF488", "AF555", "AF647", "A488_background", "A555_background", 
"A647_background", "FDX1", "CD357", "CD1D", "CD163", "CD3D", 
"CD31", "LDH", "CD66B", "VDAC1", "ELANE", "CD57", "CD45", "CD11B", 
"SMA", "CD16", "ECAD", "FOXP3", "NCAM", "CD4", "KERATIN", "CD14", 
"IBA1", "CD1B", "CD8A")

library(GBCC.mcm)
library(alabaster.sfe)
#fis = GBCC.mcm:::h5ads

library(shiny)

ui = fluidPage(
 sidebarPanel(
  helpText("MCMICRO samples"),
  radioButtons("samples", "samples",
    fis), 
  width=2
  ),
 mainPanel(
  tabsetPanel(
   tabPanel("basic",
    plotOutput("simple"),
    uiOutput("checkboxes")
   ),
   tabPanel("neighbors",
    plotOutput("graph"),
   ),
   tabPanel("about",
    textOutput("iniabout"),
    helpText("sources at github.com/vjcitn/GBCC.mcm"),
    verbatimTextOutput("cursfe"),
    helpText("Session information:"),
    verbatimTextOutput("sessinf")
   )
  )
 )
)
