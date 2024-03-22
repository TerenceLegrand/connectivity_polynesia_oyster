network_plot <- function(archipelago = NULL) {
  
  library(ggplot2)
  
  matrix_name <- paste0(archipelago,"_back_temporal_mean_2010_2019",".matrix")
  
  matrix.oc <- read.delim(here::here("data","derived-data",matrix_name),
                          sep =" ",
                          header = FALSE)
  colnames(matrix.oc) <- c("from","to","weight")
  
  
  connectivity.plot <- NULL

  for (i in 1:length(matrix.oc$from)) {
    connectivity.i <- data.frame(site_from = site$island[site$matrix_id == matrix.oc$from[i]],
                                 archipelago_from = site$archipelago[site$matrix_id == matrix.oc$from[i]],
                                 lon_from = site$lon[site$matrix_id == matrix.oc$from[i]],
                                 lat_from = site$lat[site$matrix_id == matrix.oc$from[i]],
                                 site_to = site$island[site$matrix_id == matrix.oc$to[i]],
                                 archipelago_to = site$archipelago[site$matrix_id == matrix.oc$to[i]],
                                 lon_to = site$lon[site$matrix_id == matrix.oc$to[i]],
                                 lat_to = site$lat[site$matrix_id == matrix.oc$to[i]],
                                 proba = matrix.oc$weight[i])
    
    connectivity.plot <- rbind(connectivity.plot,connectivity.i)
  }
  
  connectivity.study <- connectivity.plot
  connectivity.study$Connectivity.distance <- -log(connectivity.study$proba)
  
  # connectivity.final <- read.xlsx(paste0(save_dir,"connectivity_",archipelago,".xlsx"))
  # connectivity.study <- connectivity.final[connectivity.final$connection == "explicit" & connectivity.final$nbr_gen == 1,]
  # 
  # connectivity.study <- connectivity.study[! is.na(connectivity.study$Connectivity.distance),]
  
  lineConnections <- list()
  lineStrenght <- numeric(0)
  
  pathCoordinates_lon = numeric(0)
  pathCoordinates_lat = numeric(0)
  pathText <- data.frame()
  
  ##
  
  x.min <- -155 
  x.max <- -134 
  y.min <- -25
  y.max <- -7 
  
  worldMap <- rnaturalearth::ne_countries(scale = 10, returnclass = "sp")
  regionMap <- raster::crop(worldMap,raster::extent(x.min,x.max,y.min,y.max))
  
  for( l.i in sort(connectivity.study$Connectivity.distance, index.return=T, decreasing = TRUE, na.last = FALSE)$ix ){
    lineStrenght <- c(lineStrenght,(connectivity.study$Connectivity.distance[l.i] ) )
    pointFrom <- c(as.numeric(as.character(connectivity.study$lon_from[l.i])),as.numeric(as.character(connectivity.study$lat_from[l.i]))) 
    pointTo <- c(as.numeric(as.character(connectivity.study$lon_to[l.i])),as.numeric(as.character(connectivity.study$lat_to[l.i]))) 
    routes_sl <- geosphere::gcIntermediate(matrix(pointFrom,ncol=2),matrix(pointTo,ncol=2),n = 100, addStartEnd = TRUE, sp = TRUE, breakAtDateLine=TRUE)
    lineConnections = c(lineConnections,sp::SpatialLinesDataFrame(routes_sl, data.frame(ID = l.i), match.ID = F))
  }
  
  
  lineConnectionsSp <- do.call(rbind, lineConnections)
  
  lineStrenght_col <- lineStrenght
  
  cat_distribution <- boxplot(connectivity.study$Connectivity.distance)
  
  
  for( l.s in 1:length(lineStrenght) ) {
    
    if (! is.na(lineStrenght[l.s])) {
      if( lineStrenght[l.s] >= cat_distribution$stats[5] ) { lineStrenght_col[l.s] <- 1 }
      if( lineStrenght[l.s] >= cat_distribution$stats[4] & lineStrenght[l.s] < cat_distribution$stats[5] ) { lineStrenght_col[l.s] <- 2 }
      if( lineStrenght[l.s] >= cat_distribution$stats[3] & lineStrenght[l.s] < cat_distribution$stats[4] ) { lineStrenght_col[l.s] <- 3 }
      if( lineStrenght[l.s] >= cat_distribution$stats[2] & lineStrenght[l.s] < cat_distribution$stats[3] ) { lineStrenght_col[l.s] <- 4 }
      if( lineStrenght[l.s] < cat_distribution$stats[2]) { lineStrenght_col[l.s] <- 5 }
      
    }
  }
  
  myColors <- (c("#6FBBE8","#A1ECD8","#F6F9AB","#FCB46D","#B21414"))
  
  
  # Connected
  map_link <- ggplot() +
    geom_polygon(data = regionMap, aes(x = long, y = lat, group = group), fill="#CDCDCD", colour = "#9E9E9E" , linewidth =0.25 ) +
    geom_sf(data = sf::st_as_sf(lineConnectionsSp) , size= 1 , colour = myColors[lineStrenght_col],alpha =1) +
    geom_point(data = site, aes(x = lon, y = lat),size=1,alpha =0.5, colour="#666666") +
    geom_point(data = pop.coordinates, aes(x = Lon, y = Lat),size=1,alpha =1, colour="black") +
    theme(legend.position='bottom') + xlab("Lon") + ylab("Lat") +
    theme_map
  
  return(map_link)
  
}

legend_network_plot <- function(archipelago = NULL) {
  
  library(ggplot2)
  
  matrix_name <- paste0(archipelago,"_back_temporal_mean_2010_2019",".matrix")
  
  matrix.oc <- read.delim(here::here("data","derived-data",matrix_name),
                          sep =" ",
                          header = FALSE)
  colnames(matrix.oc) <- c("from","to","weight")
  
  connectivity.plot <- NULL
  
  for (i in 1:length(matrix.oc$from)) {
    connectivity.i <- data.frame(site_from = site$island[site$matrix_id == matrix.oc$from[i]],
                                 archipelago_from = site$archipelago[site$matrix_id == matrix.oc$from[i]],
                                 lon_from = site$lon[site$matrix_id == matrix.oc$from[i]],
                                 lat_from = site$lat[site$matrix_id == matrix.oc$from[i]],
                                 site_to = site$island[site$matrix_id == matrix.oc$to[i]],
                                 archipelago_to = site$archipelago[site$matrix_id == matrix.oc$to[i]],
                                 lon_to = site$lon[site$matrix_id == matrix.oc$to[i]],
                                 lat_to = site$lat[site$matrix_id == matrix.oc$to[i]],
                                 proba = matrix.oc$weight[i])
    
    connectivity.plot <- rbind(connectivity.plot,connectivity.i)
  }
  
  connectivity.plot[connectivity.plot$archipelago_from == "marquesas",]
  connectivity.plot[connectivity.plot$archipelago_from != "marquesas" & connectivity.plot$archipelago_from == "marquesas",]
  connectivity.study <- connectivity.plot
  connectivity.study$Connectivity.distance <- -log(connectivity.study$proba)
  
  lineConnections <- list()
  lineStrenght <- numeric(0)
  
  for( l.i in sort(connectivity.study$Connectivity.distance, index.return=T, decreasing = TRUE, na.last = FALSE)$ix ){
    lineStrenght <- c(lineStrenght,(connectivity.study$Connectivity.distance[l.i] ) )
  }
  
  
  
  lineStrenght_col <- lineStrenght
  
  cat_distribution <- boxplot(connectivity.study$Connectivity.distance)
  
  
  for( l.s in 1:length(lineStrenght) ) {
    
    if (! is.na(lineStrenght[l.s])) {
      if( lineStrenght[l.s] >= cat_distribution$stats[5] ) { lineStrenght_col[l.s] <- 1 }
      if( lineStrenght[l.s] >= cat_distribution$stats[4] & lineStrenght[l.s] < cat_distribution$stats[5] ) { lineStrenght_col[l.s] <- 2 }
      if( lineStrenght[l.s] >= cat_distribution$stats[3] & lineStrenght[l.s] < cat_distribution$stats[4] ) { lineStrenght_col[l.s] <- 3 }
      if( lineStrenght[l.s] >= cat_distribution$stats[2] & lineStrenght[l.s] < cat_distribution$stats[3] ) { lineStrenght_col[l.s] <- 4 }
      if( lineStrenght[l.s] < cat_distribution$stats[2]) { lineStrenght_col[l.s] <- 5 }
      
    }
  }
  
  myColors <- (c("#6FBBE8","#A1ECD8","#F6F9AB","#FCB46D","#B21414"))
  
  
  
  legend <- data.frame(strength = lineStrenght,
                       color = lineStrenght_col)
  
  
  map_link <- ggplot(data= legend, aes(x = strength,y=strength, color = color)) +
    geom_point() +
    geom_point(aes(x = legend$strength, y = legend$strength), color = "white", size = 2) + 
    scale_colour_stepsn(colours = myColors, name = "Distribution", breaks = c(0, 0.25, 0.5, 0.75,0.95, 1),labels = c("lower whisker","Q1","Median","Q3","upper whisker","1") ,limits = c(0, 1)) +
    theme_void()
  
  
  
  return(map_link)
  
}


saturation_plot <- function(archipelago = NULL) {
  
  library(ggplot2)
  
  saturation.final <- openxlsx::read.xlsx(here::here("outputs",paste0("saturation_",archipelago,".xlsx")))
  
  saturation.final.df <- rbind(data.frame(connection = "filial",nbr_gen = saturation.final$gen,saturation = saturation.final$filial),
                               data.frame(connection = "coalescent",nbr_gen = saturation.final$gen,saturation = saturation.final$coalescent),
                               data.frame(connection = "shortest_path",nbr_gen = saturation.final$gen[1],saturation = saturation.final$shortest_path[1]))
  
  
  plot_saturation <-  ggplot(saturation.final.df,aes(x=nbr_gen, y=saturation, group = connection)) +
    geom_point(aes(colour = connection),size=1, alpha = 0.5) +
    geom_line(aes(colour = connection),size=1, alpha = 0.5) +
    scale_x_log10(breaks = c(1,10,100,500),limits = c(range(nbr.gen))) + ylim(c(10,100)) +
    ylab(paste0("% island connected")) + xlab(paste0("Nbr generation")) +
    theme_map
  
  # plot_saturation <-  ggplot(saturation.final) +
  #     geom_point(aes(x=gen, y=shortest_path),size=1, alpha = 0.5, color = "grey") +
  #     geom_line(aes(x=gen, y=filial),size=1, alpha = 0.5, color = "#88d8b0") +
  #     geom_point(aes(x=gen, y=filial),size=1, alpha = 0.5, color = "#88d8b0") +
  #     geom_line(aes(x=gen, y=coalescent),size=1, alpha = 0.5, color = "#ff6f69") +
  #     geom_point(aes(x=gen, y=coalescent),size=1, alpha = 0.5, color = "#ff6f69") +
  #     scale_x_log10(breaks = c(1,10,100,500),limits = c(range(nbr.gen))) + ylim(0,100) +
  #     ylab(paste0("% island connected")) + xlab(paste0("Nbr generation")) +
  #     theme_map
  
  return(plot_saturation)
  
  dev.off()
  
  
}


aic_plot <- function(archipelago = NULL) {
  
  results.final <- openxlsx::read.xlsx(here::here("outputs",paste0("results_",archipelago,".xlsx")))
  
  plot_aic <-  ggplot(results.final,aes(x=nbr_gen, y=aic.mean.obs.pred, group = connection)) +
    geom_point(aes(colour = connection),size=1, alpha = 0.5) +
    geom_line(aes(colour = connection),size=1, alpha = 0.5) +
    scale_x_log10(breaks = c(1,10,100,500),limits = c(range(nbr.gen))) +
    ylab(paste0("AIC")) + xlab(paste0("Nbr generation")) +
    theme_map
  
  return(plot_aic)
  
  dev.off()
  
}

scatter_plot <- function(archipelago = NULL, connection = NULL) {
  
  library(ggplot2)
  
  results.final <- openxlsx::read.xlsx(here::here("outputs",paste0("results_",archipelago,".xlsx")))
  connectivity.final <- openxlsx::read.xlsx(here::here("outputs",paste0("connectivity_",archipelago,".xlsx")))
  summary_table <- openxlsx::read.xlsx(here::here("outputs",paste0("summary_table",".xlsx")))
  
  summary_table.i <- summary_table[summary_table$archipelago == archipelago & summary_table$connection == connection,]
  
  connectivity.i <- connectivity.final[connectivity.final$nbr_gen == summary_table.i$generation & connectivity.final$connection == connection,]
  results.i <- results.final[results.final$nbr_gen == summary_table.i$generation & results.final$connection == connection,]
  
  
  range.y <- range(connectivity.i$Differentiation[!is.na(connectivity.i$mean.obs.pred)], na.rm = TRUE)#range(connectivity.study$Differentiation)
  range.x <- range(connectivity.i$mean.obs.pred[!is.na(connectivity.i$mean.obs.pred)], na.rm = TRUE) #range(connectivity.study$predict.mean.hc)
  
  
  if (archipelago == "tuamotuWest_tuamotuEast_gambier_society_marquesas") {
    
    scatter <- ggplot(data =connectivity.i, aes(x=mean.obs.pred, y=Differentiation, group = 1)) + 
      geom_point(size=1,color="#5B5B5B", alpha = 0.5) +
      geom_point(data =connectivity.i[connectivity.i$archipelago_from == "marquesas" | connectivity.i$archipelago_to == "marquesas",], aes(x=mean.obs.pred, y=Differentiation, group = 1), size=1,color= "#6FBBE8", alpha = 0.5) +
      geom_segment(aes(x = max(c(range.x[1],range.y[1])), xend = min(c(range.x[2],range.y[2])), y = max(c(range.x[1],range.y[1])), yend =min(c(range.x[2],range.y[2]))),color="grey",alpha = 0.5, size=0.15) +
      invisible(geom_smooth(method = "lm", color="black", fill="#B5CAE5", se=TRUE ,size=0.35, linetype = "longdash", formula = 'y ~ x')) +
      ylab(paste0("Observed genetic differentiation")) + xlab(paste0("Predicted genetic differentiation")) + 
      theme(aspect.ratio=1) + xlim(range.x) + ylim(range.y) +
      theme_map 
    
  }
  
  if (archipelago == "tuamotuWest_tuamotuEast_gambier_society_marquesas_austral") {
    
    scatter <- ggplot(data =connectivity.i, aes(x=mean.obs.pred, y=Differentiation, group = 1)) + 
      geom_point(size=1,color="#5B5B5B", alpha = 0.5) +
      geom_point(data =connectivity.i[connectivity.i$archipelago_from == "marquesas" | connectivity.i$archipelago_to == "marquesas",], aes(x=mean.obs.pred, y=Differentiation, group = 1), size=1,color= "#6FBBE8", alpha = 0.5) +
      geom_point(data =connectivity.i[connectivity.i$archipelago_from == "austral" | connectivity.i$archipelago_to == "austral",], aes(x=mean.obs.pred, y=Differentiation, group = 1), size=1,color="#A1ECD8", alpha = 0.5) +
      geom_segment(aes(x = max(c(range.x[1],range.y[1])), xend = min(c(range.x[2],range.y[2])), y = max(c(range.x[1],range.y[1])), yend =min(c(range.x[2],range.y[2]))),color="grey",alpha = 0.5, size=0.15) +
      invisible(geom_smooth(method = "lm", color="black", fill="#B5CAE5", se=TRUE ,size=0.35, linetype = "longdash", formula = 'y ~ x')) +
      ylab(paste0("Observed genetic differentiation")) + xlab(paste0("Predicted genetic differentiation")) + 
      theme(aspect.ratio=1) + xlim(range.x) + ylim(range.y) +
      theme_map 
    
  }
  
  if (archipelago == "tuamotuWest_tuamotuEast_gambier_society") {
    
    scatter <- ggplot(data =connectivity.i, aes(x=mean.obs.pred, y=Differentiation, group = 1)) + 
      geom_point(size=1,color="#5B5B5B", alpha = 0.5) +
      geom_point(data =connectivity.i[connectivity.i$island_from == "MaruteaSud" | connectivity.i$island_to == "MaruteaSud",], aes(x=mean.obs.pred, y=Differentiation, group = 1), size=1,color="#FCB46D", alpha = 0.5) +
      geom_point(data =connectivity.i[connectivity.i$island_from == "Mangareva" | connectivity.i$island_to == "Mangareva",], aes(x=mean.obs.pred, y=Differentiation, group = 1), size=1,color="#F6F9AB", alpha = 0.5) +
      geom_segment(aes(x = max(c(range.x[1],range.y[1])), xend = min(c(range.x[2],range.y[2])), y = max(c(range.x[1],range.y[1])), yend =min(c(range.x[2],range.y[2]))),color="grey",alpha = 0.5, size=0.15) +
      invisible(geom_smooth(method = "lm", color="black", fill="#B5CAE5", se=TRUE ,size=0.35, linetype = "longdash", formula = 'y ~ x')) +
      ylab(paste0("Observed genetic differentiation")) + xlab(paste0("Predicted genetic differentiation")) + 
      theme(aspect.ratio=1) + xlim(range.x) + ylim(range.y) +
      theme_map 
    
  }
  
  return(scatter)
  
}

deviation_island <- function(archipelago = NULL, connection = NULL) {
  
  library(ggplot2)
  
  results.final <- openxlsx::read.xlsx(here::here("outputs",paste0("results_",archipelago,".xlsx")))
  connectivity.final <- openxlsx::read.xlsx(here::here("outputs",paste0("connectivity_",archipelago,".xlsx")))
  summary_table <- openxlsx::read.xlsx(here::here("outputs",paste0("summary_table",".xlsx")))
  
  summary_table.i <- summary_table[summary_table$archipelago == archipelago & summary_table$connection == connection,]
  
  connectivity.i <- connectivity.final[connectivity.final$nbr_gen == summary_table.i$generation & connectivity.final$connection == connection,]
  results.i <- results.final[results.final$nbr_gen == summary_table.i$generation & results.final$connection == connection,]
  
  
  fit <- lm(connectivity.i$Differentiation ~ connectivity.i$mean.obs.pred)
  out <- summary(fit)
  fit_curve <- out$coefficients[2,1] * connectivity.i$mean.obs.pred + out$coefficients[1,1] 
  
  
  test_deviation <- data.frame(islands.from = connectivity.i$island_from,
                               islands.to = connectivity.i$island_to,
                               deviation = connectivity.i$Differentiation - fit_curve)
  
  
  islands_i <- unique(c(connectivity.i$island_from,connectivity.i$island_to))
  
  
  deviation_per_island <- data.frame(island = NA,
                                     deviation = NA)
  
  deviation_per_island <- data.frame(NULL)
  
  for (i in 1:length(islands_i)) {
    
    deviation_per_island_i <- data.frame(island = islands_i[i],
                                         deviation = test_deviation$deviation[test_deviation$islands.from == islands_i[i] | test_deviation$islands.to == islands_i[i]])
    
    
    deviation_per_island <- rbind(deviation_per_island,deviation_per_island_i)
  }
  
  
  deviation_per_island <- deviation_per_island[! is.na(deviation_per_island$deviation),]
  
  plot_deviation_island <- ggplot(deviation_per_island, aes(island, deviation)) + 
    geom_boxplot() +
    geom_hline(yintercept=0, linetype="dashed", color = "grey") +
    xlab(paste0("Sampled Island")) +
    ylab(paste0("Deviation from linear fit")) + 
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
    theme_plot
  
  return(plot_deviation_island)
  
}


deviation_island_test <- function(archipelago = NULL, connection = NULL) {
  
  results.final <- openxlsx::read.xlsx(here::here("outputs",paste0("results_",archipelago,".xlsx")))
  connectivity.final <- openxlsx::read.xlsx(here::here("outputs",paste0("connectivity_",archipelago,".xlsx")))
  summary_table <- openxlsx::read.xlsx(here::here("outputs",paste0("summary_table",".xlsx")))
  
  summary_table.i <- summary_table[summary_table$archipelago == archipelago & summary_table$connection == connection,]
  
  connectivity.i <- connectivity.final[connectivity.final$nbr_gen == summary_table.i$generation & connectivity.final$connection == connection,]
  results.i <- results.final[results.final$nbr_gen == summary_table.i$generation & results.final$connection == connection,]
  
  
  fit <- lm(connectivity.i$Differentiation ~ connectivity.i$mean.obs.pred)
  out <- summary(fit)
  fit_curve <- out$coefficients[2,1] * connectivity.i$mean.obs.pred + out$coefficients[1,1] 
  
  
  test_deviation <- data.frame(islands.from = connectivity.i$island_from,
                               islands.to = connectivity.i$island_to,
                               deviation = connectivity.i$Differentiation - fit_curve)
  
  
  islands_i <- unique(c(connectivity.i$island_from,connectivity.i$island_to))
  
  
  deviation_per_island <- data.frame(island = NA,
                                     deviation = NA)
  
  deviation_per_island <- data.frame(NULL)
  
  for (i in 1:length(islands_i)) {
    
    deviation_per_island_i <- data.frame(island = islands_i[i],
                                         deviation = test_deviation$deviation[test_deviation$islands.from == islands_i[i] | test_deviation$islands.to == islands_i[i]])
    
    
    deviation_per_island <- rbind(deviation_per_island,deviation_per_island_i)
  }
  
  
  deviation_per_island <- deviation_per_island[! is.na(deviation_per_island$deviation),]
  
  Kruskal.Wallis.results <- stats::kruskal.test(deviation ~ island, data = deviation_per_island)
  
  summary(Kruskal.Wallis.results)
  
  pairwise.wilcox.test.results <- stats::pairwise.wilcox.test(deviation_per_island$deviation, deviation_per_island$island,
                                                       p.adjust.method = "BH", paired = FALSE)
  
  pairwise.wilcox.test.results <- as.data.frame(pairwise.wilcox.test.results$p.value)
  
  file_name <- here::here("outputs",paste0("pairwise_wilcox_test_table",'.xlsx'))
  openxlsx::write.xlsx(pairwise.wilcox.test.results, file = file_name, rowNames = TRUE, colNames = TRUE)
  
  return(Kruskal.Wallis.results)
  
}


filial_vs_coalescent_connection <- function(archipelago = NULL, nbr_gen = NULL) {
  
  # Filial
  
  matrix_name <- paste0("filial_M_",as.character(nbr_gen),"_",archipelago,"_back_temporal_mean_2010_2019",".matrix")
  
  matrix.oc <- read.delim(here::here("outputs","multistep",matrix_name),
                          sep =" ",
                          header = FALSE)
  
  #matrix.oc <- read.delim(paste0(multi_gen_dir,"filial_M_",as.character(nbr.gen[k]),"_",archipelago,"_back_temporal_mean_2010_2019",".matrix"), sep =" ", header = FALSE)
  colnames(matrix.oc) <- c("from","to","weight")
  
  connectivity.plot <- NULL
  
  
  for (i in 1:length(matrix.oc$from)) {
    connectivity.i <- data.frame(site_from = site$island[site$matrix_id == matrix.oc$from[i]],
                                 archipelago_from = site$archipelago[site$matrix_id == matrix.oc$from[i]],
                                 lon_from = site$lon[site$matrix_id == matrix.oc$from[i]],
                                 lat_from = site$lat[site$matrix_id == matrix.oc$from[i]],
                                 site_to = site$island[site$matrix_id == matrix.oc$to[i]],
                                 archipelago_to = site$archipelago[site$matrix_id == matrix.oc$to[i]],
                                 lon_to = site$lon[site$matrix_id == matrix.oc$to[i]],
                                 lat_to = site$lat[site$matrix_id == matrix.oc$to[i]],
                                 proba = matrix.oc$weight[i])
    
    connectivity.plot <- rbind(connectivity.plot,connectivity.i)
  }
  
  
  connectivity.plot.filial <- subset(connectivity.plot, select =  -proba)
  
  # Coalescent
  
  matrix_name <- paste0("coalescent_M_",as.character(nbr_gen),"_",archipelago,"_back_temporal_mean_2010_2019",".matrix")
  
  matrix.oc <- read.delim(here::here("outputs","multistep",matrix_name),
                          sep =" ",
                          header = FALSE)
  
  #matrix.oc <- read.delim(paste0(multi_gen_dir,"filial_M_",as.character(nbr.gen[k]),"_",archipelago,"_back_temporal_mean_2010_2019",".matrix"), sep =" ", header = FALSE)
  colnames(matrix.oc) <- c("from","to","weight")
  
  connectivity.plot <- NULL
  
  for (i in 1:length(matrix.oc$from)) {
    connectivity.i <- data.frame(site_from = site$island[site$matrix_id == matrix.oc$from[i]],
                                 archipelago_from = site$archipelago[site$matrix_id == matrix.oc$from[i]],
                                 lon_from = site$lon[site$matrix_id == matrix.oc$from[i]],
                                 lat_from = site$lat[site$matrix_id == matrix.oc$from[i]],
                                 site_to = site$island[site$matrix_id == matrix.oc$to[i]],
                                 archipelago_to = site$archipelago[site$matrix_id == matrix.oc$to[i]],
                                 lon_to = site$lon[site$matrix_id == matrix.oc$to[i]],
                                 lat_to = site$lat[site$matrix_id == matrix.oc$to[i]],
                                 proba = matrix.oc$weight[i])
    
    connectivity.plot <- rbind(connectivity.plot,connectivity.i)
  }
  
  connectivity.plot.coalescent <- subset(connectivity.plot, select =  -proba)
  
  length(connectivity.plot.coalescent$site_from)
  length(connectivity.plot.filial$site_from)
  
  compare <- dplyr::anti_join(connectivity.plot.coalescent,connectivity.plot.filial)
  
  return(compare)
  
}

centrality_hm <- function(archipelago = NULL){
  
  library(ggplot2)
  
  matrix_name <- paste0(archipelago,"_back_temporal_mean_2010_2019",".matrix")
  
  matrix.oc <- read.delim(here::here("data","derived-data",matrix_name),
                          sep =" ",
                          header = FALSE)
  
  colnames(matrix.oc) <- c("from","to","weight")
  
  # Create graph object from list
  
  matrix.oc$weight <- -log(matrix.oc$weight)
  matrix.oc <- matrix.oc[matrix.oc$weight > 0,]
  
  #matrix.oc$weight <- matrix.oc$weight + 1
  
  graph.back <- igraph::graph_from_data_frame(matrix.oc,
                                      directed = TRUE,
                                      vertices = NULL)
  
  hc <- igraph::harmonic_centrality(graph.back,
                            vids = igraph::V(graph.back),
                            mode = "out",
                            normalized = TRUE,
                            cutoff = -1)
  
  
  centrality <- data.frame(index_hc = as.numeric(rownames(as.data.frame(hc))),
                           hc = as.numeric(hc))
  site$hc <- NA
  
  for (i in 1:length(centrality$index_hc)) {
    centrality$index_hc[i]
    site$hc[centrality$index_hc[i]+1] <- centrality$hc[i]
  }
  
  
  site_data_filt <- site[!is.na(site$hc),]
  
  site_data_filt$hc_norm <- (site_data_filt$hc-min(site_data_filt$hc))/diff(range(site_data_filt$hc))
  
  centrality_plot <- ggplot() +
    geom_point(data= site_data_filt, aes(x = lon, y = lat, colour = hc_norm)) +
    viridis::scale_color_viridis(option = "viridis") +
    ylab(paste0("Longitude")) + xlab(paste0("Latitude")) + 
    labs(color = "Centrality") +
    ggrepel::geom_label_repel(data= site_data_filt[site_data_filt$hc_norm > quantile(site_data_filt$hc_norm,0.95),], aes(x = lon, y = lat, label = island),
                     box.padding   = 0.35, 
                     point.padding = 0.5,
                     force = 1000,
                     arrow = arrow(angle = 30, length = unit(0.1, "inches"),
                                   ends = "last", type = "open"),
                     segment.color = 'grey50') +
    theme(aspect.ratio=1) +
    #geom_point(data= site_data_filt, aes(x = lon, y = lat, colour = hc_norm)) +
    #theme(legend.position = "bottom") +
    theme_map
  
    return(centrality_plot)
  
}



centrality_btw <- function(archipelago = NULL){
  
  library(ggplot2)
  
  matrix_name <- paste0(archipelago,"_back_temporal_mean_2010_2019",".matrix")
  
  matrix.oc <- read.delim(here::here("data","derived-data",matrix_name),
                          sep =" ",
                          header = FALSE)
  
  colnames(matrix.oc) <- c("from","to","weight")
  # Create graph object from list
  
  matrix.oc$weight <- -log(matrix.oc$weight)
  matrix.oc <- matrix.oc[matrix.oc$weight > 0,]
  
  #matrix.oc$weight <- matrix.oc$weight + 1
  
  graph.back <- igraph::graph_from_data_frame(matrix.oc,
                                      directed = TRUE,
                                      vertices = NULL)
  
  
  btw <- igraph::betweenness(graph.back,
                     v = igraph::V(graph.back),
                     directed = TRUE,
                     weights = NULL,
                     cutoff = -1)
  
  
  
  centrality <- data.frame(index_btw = as.numeric(rownames(as.data.frame(btw))),
                           btw = as.numeric(btw))
  site$btw <- NA
  
  for (i in 1:length(centrality$index_btw)) {
    site$btw[centrality$index_btw[i]+1] <- centrality$btw[i]
    
  }
  
  
  site_data_filt <- site[!is.na(site$btw),]
  
  site_data_filt$btw_norm <- (site_data_filt$btw-min(site_data_filt$btw))/diff(range(site_data_filt$btw))
  
  centrality_plot <- ggplot() +
    geom_point(data= site_data_filt, aes(x = lon, y = lat, colour = btw_norm)) +
    viridis::scale_color_viridis(option = "viridis") +
    ylab(paste0("Longitude")) + xlab(paste0("Latitude")) + 
    labs(color = "Centrality") +
    ggrepel::geom_label_repel(data= site_data_filt[site_data_filt$btw_norm > quantile(site_data_filt$btw_norm,0.95),], aes(x = lon, y = lat, label = island),
                     box.padding   = 0.35, 
                     point.padding = 0.5,
                     force = 1000,
                     arrow = arrow(angle = 30, length = unit(0.1, "inches"),
                                   ends = "last", type = "open"),
                     segment.color = 'grey50') +
    theme(aspect.ratio=1) +
    #geom_point(data= site_data_filt, aes(x = lon, y = lat, colour = btw_norm)) +
    #theme(legend.position = "bottom") +
    theme_map
  
  return(centrality_plot)
  
}
