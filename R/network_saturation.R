saturation <- function(archipelago = NULL, nbr.gen = NULL, site_data = NULL) {
  
  
  if (archipelago == "tuamotuWest_tuamotuEast_gambier_society") {
    nbr_island <- length(site_data$island) - sum(site_data$archipelago == "marquesas") - sum(site_data$archipelago == "austral")
  }
  
  if (archipelago == "tuamotuWest_tuamotuEast_gambier_society_marquesas") {
    nbr_island <- length(site_data$island)
  }
  
  if (archipelago == "tuamotuWest_tuamotuEast_gambier_society_marquesas_austral") {
    nbr_island <- length(site_data$island)
  }
  
  
  results_saturation <- data.frame(gen = nbr.gen,
                                   shortest_path = NA,
                                   filial = NA,
                                   coalescent = NA)
  
  for (k in 1:length(nbr.gen)) {
    
    matrix_name <- paste0("filial_M_",as.character(nbr.gen[k]),"_",archipelago,"_back_temporal_mean_2010_2019",".matrix")
    
    matrix.oc <- read.delim(here::here("outputs","multistep",matrix_name),
                            sep =" ",
                            header = FALSE)
    
    #matrix.oc <- read.delim(paste0(multi_gen_dir,"filial_M_",as.character(nbr.gen[k]),"_",archipelago,"_back_temporal_mean_2010_2019",".matrix"), sep =" ", header = FALSE)
    colnames(matrix.oc) <- c("from","to","weight")
    
    results_saturation$filial[k] <- (length(matrix.oc$from)/nbr_island^2)*100
    
    matrix_name <- paste0("coalescent_M_",as.character(nbr.gen[k]),"_",archipelago,"_back_temporal_mean_2010_2019",".matrix")
    
    matrix.oc <- read.delim(here::here("outputs","multistep",matrix_name),
                            sep =" ",
                            header = FALSE)
    
    #matrix.oc <- read.delim(paste0(multi_gen_dir,"coalescent_M_",as.character(nbr.gen[k]),"_",archipelago,"_back_temporal_mean_2010_2019",".matrix"), sep =" ", header = FALSE)
    colnames(matrix.oc) <- c("from","to","weight")
    
    results_saturation$coalescent[k] <- (length(matrix.oc$from)/nbr_island^2)*100
    
  }
  
  
  # test <-  ggplot(results_saturation)+
  #   geom_line(aes(x=nbr.gen, y=filial),size=1, alpha = 0.5, color = "#88d8b0") +
  #   geom_point(aes(x=nbr.gen, y=filial),size=1, alpha = 0.5, color = "#88d8b0") +
  #   geom_line(aes(x=nbr.gen, y=coalescent),size=1, alpha = 0.5, color = "#ff6f69") +
  #   geom_point(aes(x=nbr.gen, y=coalescent),size=1, alpha = 0.5, color = "#ff6f69") +
  #   scale_x_log10(breaks = c(1,10,100,500),limits = c(range(nbr.gen))) + ylim(0,100) +
  #   ylab(paste0("% island connected")) + xlab(paste0("Nbr generation")) +
  #   theme_map
  # 
  # 
  #   ggsave(filename = paste0("test_",archipelago,".pdf"),
  #          plot = test,
  #          device = cairo_pdf(),
  #          path = save_dir
  #   )
  # 
  # 
  
  ## Shortest distances -----
  
  matrix_name <- paste0(archipelago,"_back_temporal_mean_2010_2019",".matrix")
  
  matrix.oc <- read.delim(here::here("data","derived-data",matrix_name),
                          sep =" ",
                          header = FALSE)
  
  #matrix.oc <- read.delim(paste0(dir,"Connectivity_matrices_processed/temp_mean/",archipelago,"_back_temporal_mean_2010_2019",".matrix"), sep =" ", header = FALSE)
  colnames(matrix.oc) <- c("from","to","weight")
  
  # Create graph object from list
  
  matrix.oc$weight <- -log(matrix.oc$weight)
  
  graph.back <- igraph::graph_from_data_frame(matrix.oc,
                                      directed = TRUE,
                                      vertices = NULL)
  
  nbr_path_all <- matrix(NA, nrow = length(V(graph.back)), ncol = length(V(graph.back)))
  
  for (a in 1:length(V(graph.back))) {
    for (b in 1:length(V(graph.back))) {
      
      path <- igraph::shortest_paths(
        graph.back,
        from = V(graph.back)[a],
        to = V(graph.back)[b],
        mode = "out",
        weights = NULL)
      
      nbr_path_all[a,b] <- length(as.numeric(path$vpath[[1]]))-1
      
    }}
  
  # -1 when there are no connection, transform it into NA
  nbr_path_all[nbr_path_all <0] <- NA
  # No tacking self-loops into account 
  #diag(nbr_path_all) <- NA
  
  
  results_saturation$shortest_path[1] <- (sum(!is.na(nbr_path_all))/nbr_island^2)*100
  results_saturation$shortest_path[2] <- mean(as.numeric(nbr_path_all),na.rm = TRUE)
  
  
  # Save
  
  
  file_name <- here::here("outputs",paste0("saturation_",archipelago,'.xlsx'))
  openxlsx::write.xlsx(results_saturation, file = file_name, rowNames = TRUE, colNames = TRUE)
  
}