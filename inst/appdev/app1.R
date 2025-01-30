
fis = GBCC.mcm:::h5ads

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
   pa = system.file(file.path("h5ad", input$samples), package="GBCC.mcm")
   mcm53.1 = process_mcmicro(pa)
   features_use = c("NCAM", "FOXP3", "CD8A", "LDH")
   Voyager::plotLocalResult(mcm53.1, "localG", features = features_use,
                   colGeometryName = "centroids", divergent = TRUE,
                   diverge_center = 0)
   })
}

runApp(list(ui=ui, server=server))

