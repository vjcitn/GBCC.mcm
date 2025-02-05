fol = "."
aa = c("m62.1.sfe", "m62.2.sfe", "m62.3.sfe", "m62.4.sfe")
fis = paste(fol, aa, sep="/")
names(fis) = c("m62.1", "m62.2", "m62.3", "m62.4")

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
    fluidRow( 
      column(checkboxInput("useknn", "use nn", value=FALSE),width=2),
      column(numericInput("kpicked", "kval", value=3, min=2, max=30),width=2)
      ),
    plotOutput("graph"),
   ),
   tabPanel("focus",
    uiOutput("xsliders"),
    uiOutput("ysliders"),
    plotOutput("ellipses")
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
