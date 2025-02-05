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


server = function(input, output) {
 getcur = reactive({
  cur = readObject(input$samples)
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
       selected=c("NCAM", "FOXP3", "CD8A", "LDH"), inline=TRUE)
   })
 output$sessinf = renderPrint({
  sessionInfo()
  })
 output$ellipses = renderPlot({
##  x = readObject(input$samples)
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
  ttt = lapply(tt, function(x) matrix(unlist(x), nr=101))
  validate(need(length(input$xlow)>0,"waiting for controls"))
  tttf = lapply(ttt, function(x) {x[x[,1]>input$xlow & x[,1]<input$xhi & x[,2]>input$ylow & x[,2]<input$yhi,]})
  ok = sapply(tttf, function(x) nrow(x)>0)
  tttt = sf::st_multilinestring(tttf[which(ok)])
  dd = data.frame(dd)
  dd$ph = mm$phenotype
  ggplot2::ggplot(tttt) + ggplot2::geom_sf() + ggplot2::geom_point(data=dd[ok,], ggplot2::aes(x=X_centroid,
    y=Y_centroid, colour=ph))
  })
 output$ysliders = renderUI({
##    x = readObject(input$samples)
    x = getcur()
    mm = SummarizedExperiment::colData(x)
    validate(need(length(mm$Y_centroid)>0, "retrieving data"))
    fluidRow(
     column(width=3, sliderInput("ylow", "ylow", min=0, max=max(round(mm$Y_centroid+1,0)), value=0)),
     column(width=3, sliderInput("yhi", "yhi", min=0, max=max(round(mm$Y_centroid+1,0)), 
        value=max(round(mm$Y_centroid+1,0))))
     )
  })
 output$xsliders = renderUI({
##    x = readObject(input$samples)
    x = getcur()
    mm = SummarizedExperiment::colData(x)
    validate(need(length(mm$X_centroid)>0, "retrieving data"))
    fluidRow(
     column(width=3, sliderInput("xlow", "xlow", min=0, max=max(round(mm$X_centroid+1,0)), value=0)),
     column(width=3, sliderInput("xhi", "xhi", min=0, max=max(round(mm$X_centroid+1,0)), value=max(mm$X_centroid)))
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


