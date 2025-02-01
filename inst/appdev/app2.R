fol = "/Users/vincentcarey/ANVIL/MCMICRO/GBCC_ZIPS"
aa = c("m53.1.sfe", "m53.2.sfe", "m53.3.sfe", "m53.4.sfe", "m62.1.sfe")
fis = paste(fol, aa, sep="/")
names(fis) = c("m53.1", "m53.2", "m53.3", "m53.4", "m62.1")

library(GBCC.mcm)
library(alabaster.sfe)
#fis = GBCC.mcm:::h5ads

library(shiny)

ui = fluidPage(
 sidebarPanel(
  helpText("MCMICRO samples"),
  radioButtons("samples", "samples",
    fis), width=2
  ),
 mainPanel(
  tabsetPanel(
   tabPanel("basic",
    plotOutput("simple")
   ),
   tabPanel("basic",
    plotOutput("graph")
   ),
   tabPanel("about",
    verbatimTextOutput("iniabout")
   )
  )
 )
)

server = function(input, output) {
 output$iniabout = renderPrint(
  sprintf("GBCC.mcm ver. %s", packageVersion("GBCC.mcm"))
 )

 output$simple = renderPlot({
   cur = readObject(input$samples)
   features_use = c("NCAM", "FOXP3", "CD8A", "LDH")
   Voyager::plotLocalResult(cur, "localG", features = features_use,
                   colGeometryName = "centroids", divergent = TRUE,
                   diverge_center = 0)
   })
 output$graph = renderPlot({
   cur = readObject(input$samples)
   Voyager::plotColGraph(cur)
   })
}

runApp(list(ui=ui, server=server))

