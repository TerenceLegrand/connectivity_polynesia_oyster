temporal_mean_matrix <- function(site = NULL) {
  
  matrix_site <- matrix(data = 0 ,nrow = length(site$island),ncol = length(site$island),dimnames = list(site$island,site$island))
  
  matrix_temporal_mean <- matrix_site
  
  day_id <- seq(as.Date("2010-01-01"), as.Date("2019-12-31"), by="days")
  count <- 0
  
  issue_file <- data.frame()
  
  for (d in 1:length(day_id)) {
    
    print(day_id[d])
    
    matrix_daily_name <- paste0("extract_count_LPP23d_PLD35d_Pmarg_",day_id[d],"_modif.csv")
    matrix_daily_list <- utils::read.csv2(here::here("data","raw-data","CSV_daily",matrix_daily_name),
                                          sep = ",")
    
    if (length(matrix_daily_list$source) == 0) {issue_file <- rbind(issue_file,matrix_daily_name) }
    
    if (length(matrix_daily_list$source) > 0) {
      
      count = count + 1
      
      matrix_daily <- matrix_site
      
      for (i in 1:length(matrix_daily_list$source)) {
        matrix_daily[rownames(matrix_site) %in% matrix_daily_list$source[i],colnames(matrix_site) %in% matrix_daily_list$destination[i]] <- as.numeric(matrix_daily_list$count[i])
      }
      
      sum.send <- rowSums(matrix_daily,na.rm = TRUE) # TL: sum horizontally, what cells send
      
      sum.receive <- colSums(matrix_daily,na.rm = TRUE) # TL: sum vertically, what cells receive
      
      # Backward
      matrix_daily_backward <- sweep(matrix_daily,2,sum.receive,FUN="/")
      
      matrix_daily_backward[is.nan(matrix_daily_backward)] <- 0 # Remove NaN because 0/0 give NaN (many of colsum is 0)
      matrix_daily_backward[is.infinite(matrix_daily_backward)] <- 0 # Remove Inf
      
      matrix_temporal_mean <- matrix_temporal_mean + matrix_daily_backward
      
    }
  }
  
  
  matrix_temporal_mean <- matrix_temporal_mean/count
  
  # Final normalization
  sum.receive <- colSums(matrix_temporal_mean,na.rm = TRUE) # TL: sum vertically, what cells receive
  
  # Backward
  matrix_temporal_mean_backward <- sweep(matrix_temporal_mean,2,sum.receive,FUN="/")
  
  matrix_temporal_mean_backward[is.nan(matrix_temporal_mean_backward)] <- 0 # Remove NaN because 0/0 give NaN (many of colsum is 0)
  matrix_temporal_mean_backward[is.infinite(matrix_temporal_mean_backward)] <- 0 # Remove Inf
  
  #Test to verify if normalisation went well
  print(paste0("test sum proba, should be 1 ---> ",sum(matrix_temporal_mean_backward[,sample(site$island,1)])))
  
  # Transpose
  matrix_temporal_mean_backward_inv <- t(matrix_temporal_mean_backward)
  
  # Transform into list
  matrix_temporal_mean_backward_inv_list <- reshape2::melt(matrix_temporal_mean_backward_inv)
  colnames(matrix_temporal_mean_backward_inv_list) <- c("from","to","weight")
  # Keep just the connected islands
  matrix_temporal_mean_backward_inv_list <- matrix_temporal_mean_backward_inv_list[matrix_temporal_mean_backward_inv_list$weight != 0,]
  
  
  write.table(matrix_temporal_mean_backward_inv_list,
              file = here::here("data","derived-data","tuamotuWest_tuamotuEast_gambier_society_marquesas_back_temporal_mean_2010_2019.matrix"),
              sep = " ",
              row.names = FALSE,
              col.names = FALSE,
              eol = "\r")
  
  write.table(issue_file,
              file = here::here("data","derived-data","issue_files.txt"),
              sep = " ",
              row.names = FALSE,
              col.names = FALSE,
              eol = "\r")
  
}



filter_matrix <- function(site = NULL) {
  
  # Load big connectivity matrix
  matrix.oc <- read.delim(here::here("data","derived-data","tuamotuWest_tuamotuEast_gambier_society_marquesas_back_temporal_mean_2010_2019.matrix"),
                           sep =" ", header = FALSE)
  colnames(matrix.oc) <- c("from","to","weight")
  
  ## ------------
  
  # All archipelago: add an artificial connection between australes and Toamotu
  # Hereheretue --> Raivavae i.e., Hereheretue receive propagules from Raivavae

  # The weight of this connection is equal to the minimum weight between Tuamotu and Marquesas
  marq_tuam <- matrix.oc[matrix.oc$from %in% site$matrix_id[site$archipelago != "marquesas"] &
                    matrix.oc$to %in% site$matrix_id[site$archipelago == "marquesas"],]
  
  site[site$matrix_id %in% marq_tuam$from,]
  site[site$matrix_id %in% marq_tuam$to,]
  
  # it is Fangatau --> FatuHiva, i.e., Fangatau receive propagules from FatuHiva
  
  matrix.oc.add <- data.frame(from = "Hereheretue",
                              to = "Raivavae",
                              weight = min(matrix.oc$weight[matrix.oc$from %in% site$island[site$archipelago != "marquesas"] &
                                                              matrix.oc$to %in% site$island[site$archipelago == "marquesas"]]))
  
  matrix.oc.add <- rbind(matrix.oc,matrix.oc.add)
  
  # Normalization
  
  for (i in 1:length(site$matrix_id)) {
    matrix.oc.add$weight[matrix.oc.add$from == site$island[i]] <- matrix.oc.add$weight[matrix.oc.add$from == site$island[i]]/sum(matrix.oc.add$weight[matrix.oc.add$from == site$island[i]])
    
  }
  
  write.table(matrix.oc.add,
              file = here::here("data","derived-data","tuamotuWest_tuamotuEast_gambier_society_marquesas_austral_back_temporal_mean_2010_2019.matrix"),
              sep = " ",
              row.names = FALSE,
              col.names = FALSE,
              eol = "\r")
  
  # Remove Marquise and society archipelago from the connections 
  
  matrix.oc.filter <- matrix.oc[! (matrix.oc$from %in% site$island[site$archipelago == "marquesas" | site$archipelago == "austral"] |
                                     matrix.oc$to %in% site$island[site$archipelago == "marquesas" | site$archipelago == "austral"]),]
  
  
  # Normalization
  
  for (i in 1:length(site$matrix_id)) {
    matrix.oc.filter$weight[matrix.oc.filter$from == site$island[i]] <- matrix.oc.filter$weight[matrix.oc.filter$from == site$island[i]]/sum(matrix.oc.filter$weight[matrix.oc.filter$from == site$island[i]])
  }

  write.table(matrix.oc.filter,
              file = here::here("data","derived-data","tuamotuWest_tuamotuEast_gambier_society_back_temporal_mean_2010_2019.matrix"),
              sep = " ",
              row.names = FALSE,
              col.names = FALSE,
              eol = "\r")
  
}

index_matrix <- function(site = NULL, archipelago = NULL) {
  
  matrix_name <- paste0(archipelago,"_back_temporal_mean_2010_2019.matrix")
  
  matrix.oc <- read.delim(here::here("data","derived-data",matrix_name),
                          sep =" ", header = FALSE)
  colnames(matrix.oc) <- c("from","to","weight")
  
  matrix.oc.id <- matrix.oc
  
  for (i in 1:length(matrix.oc$from)) {
    matrix.oc.id$from[i] <- site$matrix_id[site$island %in% matrix.oc$from[i]]
    matrix.oc.id$to[i] <- site$matrix_id[site$island %in% matrix.oc$to[i]]
  }

  matrix.oc.id$from <- as.numeric(matrix.oc.id$from)
  matrix.oc.id$to <- as.numeric(matrix.oc.id$to)
  
  write.table(matrix.oc.id,
              file = here::here("data","derived-data",matrix_name),
              sep = " ",
              row.names = FALSE,
              col.names = FALSE,
              eol = "\r")
  
}
  
  ####
  # 

  # matrix_list <- data.frame(from = as.factor(matrix$source),
  #                           to = as.factor(matrix$destination),
  #                           weight = as.numeric(matrix$count))
  # 
  # graph <- igraph::graph_from_data_frame(matrix_list)
  # 
  # 
  # matrix.ev <- igraph::as_adjacency_matrix(graph, edges = TRUE, names = TRUE, sparse = FALSE)
  # 
  
  

