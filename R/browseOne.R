targets_june = function() c("RNA_Pol_II_CTD", "pERK", "p53", "Ki67", "panCK", "CD45", "Ecad",
"aSMA", "Vimentin", "EGFR", "pRB", "p21", "CD20", "aSMA_2", "PD1",
"PCNA", "ER", "HER2", "CK14", "CK19", "CK17", "LaminAC", "AR",
"H2ax", "PR", "HLAA", "CK5")

#' define shiny app to explore one sample
#' @import shiny
#' @param input_processed output of `process_mcmicro`
#' @param targets character() vector
#' @examples
#' pa = system.file("h5ad", "Galaxy37-LSP12022-A3.h5ad", package="GBCC.mcm")
#' a3 = process_mcmicro(pa)
#' # be sure to run example(..., ask=FALSE)
#' if (interactive()) {
#'   browseOne(a3)
#' }
#' @export  
browseOne = function(input_processed, targets = targets_june()) {
  ui = fluidPage(
   sidebarPanel(
    helpText("Single MCMICRO viewer"),
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
      plotOutput("graph")
     ),
     tabPanel("tables",
      uiOutput("feattabsel"),
      DT::dataTableOutput("stattab")
     ),
     tabPanel("focus",
      uiOutput("xsliders"),
      uiOutput("ysliders"),
      actionButton("dofilt", "filter"),
      plotly::plotlyOutput("ellipses")
     ),
     tabPanel("about",
      textOutput("iniabout"),
      textOutput("aboutsamp"),
      helpText("sources at github.com/vjcitn/GBCC.mcm"),
      verbatimTextOutput("cursfe"),
      helpText("Session information:"),
      verbatimTextOutput("sessinf")
     )
    )
   )
  )
  
  
  
  server = function(input, output) {
   getcur = reactive({
    cur = input_processed
    if (input$useknn) {
     spneigh_parms = list(method='knearneigh', k=input$kpicked) 
     fsnargs = list(x=cur)
     g1 = do.call(SpatialFeatureExperiment::findSpatialNeighbors, c(fsnargs, spneigh_parms))
     colGraph(cur, "spatNeigh") = g1
     }
    cur
    })
   output$cursfe = renderPrint({
    getcur()
  ##  cur = readObject(input$samples)
  ##  cur
   })
   output$iniabout = renderText(
    sprintf("GBCC.mcm version. %s.  MCMICRO outputs provided as h5ad were transformed to 
  SpatialFeatureExperiment instances; some statistical computations were performed
  in advance and bound to the instances for viewing with this app.", packageVersion("GBCC.mcm"))
   )
   output$aboutsamp = renderText(
"Samples are from a tissue microarray and are imaged with Cyclic Immunofluorescence
A3 = Triple Negative Breast Cancer; A4 = Luminal B HER2+ Breast Cancer; B3 = Luminal B HER- Breast Cancer; B4 = normal jejunum")
  
   output$simple = renderPlot({
   ##  cur = readObject(input$samples)
     cur = getcur()
     validate(need(isTRUE(nchar(input$genes2[1])>0), "awaiting gene list"))
     features_use = input$genes2
     Voyager::plotLocalResult(cur, "localG", features = features_use,
                     colGeometryName = "centroids", divergent = TRUE,
                     diverge_center = 0, show_axes=TRUE)
     })
   output$graph = renderPlot({
  ##   cur = readObject(input$samples)
     cur = getcur()
     Voyager::plotColGraph(cur)
     })
   output$checkboxes = renderUI({
     checkboxGroupInput("genes2", "Select targets to display:", choices=targets,
         selected=targets[1:4], inline=TRUE)
     })
   output$sessinf = renderPrint({
    sessionInfo()
    })
   output$ellipses = plotly::renderPlotly({
    plotly::ggplotly(build_ellipses()$gg)
    })
   prep_ellipses = reactive({
    x = getcur()
    mm = SummarizedExperiment::colData(x)
    dd = data.matrix(mm[, c("X_centroid", "Y_centroid")])
    cens = lapply(seq_len(nrow(dd)), function(x) sf::st_point(dd[x,]))
    maj = mm$MajorAxisLength/2
    min = mm$MinorAxisLength/2
    angle_in_rad = mm$Orientation
    angle_in_deg = angle_in_rad/0.01745329  # pi/180, assume sign handled well
    tt = lapply(seq_len(length(cens)), function(x)
      sfdep::st_ellipse(cens[[x]], sx=maj[x], sy=min[x], rotation=angle_in_deg[x]))
    dd = data.frame(dd)
    dd$ph = mm$phenotype
    dd$cid = paste(mm$phenotype, ":", seq_len(nrow(mm))) 
    list(ttt=lapply(tt, function(x) matrix(unlist(x), nr=101)), dd=dd)
    })
   build_ellipses = eventReactive(input$dofilt, {   # FIXME -- much of it doesn't need input
#    x = getcur()
#    mm = SummarizedExperiment::colData(x)
#    dd = data.matrix(mm[, c("X_centroid", "Y_centroid")])
#    cens = lapply(seq_len(nrow(dd)), function(x) sf::st_point(dd[x,]))
#    maj = mm$MajorAxisLength/2
#    min = mm$MinorAxisLength/2
#    angle_in_rad = mm$Orientation
#    angle_in_deg = angle_in_rad/0.01745329  # pi/180, assume sign handled well
#    tt = lapply(seq_len(length(cens)), function(x)
#      sfdep::st_ellipse(cens[[x]], sx=maj[x], sy=min[x], rotation=angle_in_deg[x]))
#    ttt = lapply(tt, function(x) matrix(unlist(x), nr=101))
    ttt0 = prep_ellipses()
    ttt = ttt0$ttt
    dd = ttt0$dd
    validate(need(length(input$xlow)>0,"waiting for controls"))
    validate(need(input$xlow<input$xhi,"be sure xlow < xhi"))
    validate(need(input$ylow<input$yhi,"be sure ylow < yhi"))
    tttf = lapply(ttt, function(x) {x[x[,1]>input$xlow & x[,1]<input$xhi & x[,2]>input$ylow & x[,2]<input$yhi,]})
    ok = sapply(tttf, function(x) nrow(x)>0)
    tttt = sf::st_multilinestring(tttf[which(ok)])
    gg = ggplot2::ggplot(tttt) + ggplot2::geom_sf() + ggplot2::geom_point(data=dd[ok,], ggplot2::aes(x=X_centroid,
      y=Y_centroid, colour=ph, text=cid))
    list(gg=gg)
    }  )
   output$ysliders = renderUI({
      x = getcur()
      mm = SummarizedExperiment::colData(x)
      validate(need(length(mm$Y_centroid)>0, "retrieving data"))
      fluidRow(
       column(width=3, sliderInput("ylow", "ylow", min=0, max=max(round(mm$Y_centroid+1,0)), value=250, step=50, animate=TRUE)),
       column(width=3, sliderInput("yhi", "yhi", min=0, max=max(round(mm$Y_centroid+1,0)), 
          value=1350, step=50, animate=TRUE)
       )
      )
    })
   output$xsliders = renderUI({
      x = getcur()
      mm = SummarizedExperiment::colData(x)
      validate(need(length(mm$X_centroid)>0, "retrieving data"))
# value choices for xlow etc based on A3 example
      fluidRow(
       column(width=3, sliderInput("xlow", "xlow", min=0, max=max(round(mm$X_centroid+1,0)), value=500, step=50, animate=TRUE)),
       column(width=3, sliderInput("xhi", "xhi", min=0, max=max(round(mm$X_centroid+1,0)), value=1650, step=50, animate=TRUE))
       )
    })
   output$feattabsel = renderUI({
      x = getcur()
      feats = rownames(x)
      selectInput("feat2tab", "feature", choices=feats, selected=feats[1], multiple=FALSE)
      })
   output$stattab = DT::renderDataTable({
      x = getcur()
      gdat = localResults(x)$localG  # FIXME
      validate(need(nchar(input$feat2tab)>0, "waiting for feature to tabulate"))
      res = as(gdat[[input$feat2tab]], "data.frame")
      res = cbind(res, phenotype=x$phenotype)
      isn = which(sapply(res, is.numeric))
      if (length(isn)==0) return(res)
      for (i in isn) res[[i]] = round(res[[i]], 7)
      res
      })
  }
runApp(list(ui=ui, server=server))
}  
  
