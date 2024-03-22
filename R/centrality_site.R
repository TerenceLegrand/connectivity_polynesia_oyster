centrality_population <- function(archipelago = NULL, site = NULL, pop.coordinates = NULL){
  
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
  
  
  # Harmonic
  
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
  
  # Betwenness
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
  
  
  
  pop.coordinates$pop_id <- rownames(matrix_genetic)
  
  pop.coordinates$hc <- NA
  pop.coordinates$btw <- NA
  
  for (i in 1:length(pop.coordinates$id)) {
    pop.coordinates$hc[i] <- site$hc[site$matrix_id %in% pop.coordinates$id[i]]
    pop.coordinates$btw[i] <- site$btw[site$matrix_id %in% pop.coordinates$id[i]]
  }

  file_name <- here::here("outputs",paste0("centrality_population_",archipelago,'.xlsx'))
  openxlsx::write.xlsx(pop.coordinates, file = file_name, rowNames = TRUE, colNames = TRUE)
  
  file_name <- here::here("outputs",paste0("centrality_site_",archipelago,'.xlsx'))
  openxlsx::write.xlsx(site, file = file_name, rowNames = TRUE, colNames = TRUE)
  
}