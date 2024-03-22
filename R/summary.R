table_summary <- function(archipelago = NULL) {
  
  results.final <- openxlsx::read.xlsx(here::here("outputs", paste0("results_",archipelago,".xlsx")))
  connectivity.final <- openxlsx::read.xlsx(here::here("outputs", paste0("connectivity_",archipelago,".xlsx")))
  saturation.final <- openxlsx::read.xlsx(here::here("outputs", paste0("saturation_",archipelago,".xlsx")))
  
  # Shortest path
  connection = "shortest_path"
  
  results.final.i <- results.final[results.final$connection == connection,]
  index <- which.min(results.final.i$aic.mean.obs.pred)  
  
  summary_shortest_path <- data.frame(archipelago = archipelago,
                                      connection = "shortest_path",
                                      percentage = saturation.final$shortest_path[1],
                                      gen = saturation.final$shortest_path[2],
                                      n = results.final.i$oceanographic_pairs[index],
                                      generation = NA,
                                      AIC = results.final.i$aic.mean.obs.pred[index],
                                      R2 = results.final.i$r2.mean.obs.pred[index],
                                      deltaR2 = results.final.i$r2.diff.null[index],
                                      pval = results.final.i$p.mean.obs.pred[index])
  
  # Filial
  connection = "filial"
  
  results.final.i <- results.final[results.final$connection == connection,]
  saturation_gen <- saturation.final$gen[min(which(diff(saturation.final$filial)==0))]
  results.final.i$aic.mean.obs.pred[1:(saturation_gen-1)] <- NA
  index <- which.min(results.final.i$aic.mean.obs.pred[])  
  
  summary_filial <- data.frame(archipelago = archipelago,
                               connection = "filial",
                               percentage = saturation.final$filial[min(which(diff(saturation.final$filial)==0))],
                               gen = saturation.final$gen[min(which(diff(saturation.final$filial)==0))],
                               n = results.final.i$oceanographic_pairs[index],
                               generation = results.final.i$nbr_gen[index],
                               AIC = results.final.i$aic.mean.obs.pred[index],
                               R2 = results.final.i$r2.mean.obs.pred[index],
                               deltaR2 = results.final.i$r2.diff.null[index],
                               pval = results.final.i$p.mean.obs.pred[index])
  
  
  # Coalescent
  connection = "coalescent"
  
  results.final.i <- results.final[results.final$connection == connection,]
  saturation_gen <- saturation.final$gen[min(which(diff(saturation.final$filial)==0))]
  results.final.i$aic.mean.obs.pred[1:(saturation_gen-1)] <- NA
  index <- which.min(results.final.i$aic.mean.obs.pred)
  
  
  summary_coalescent <- data.frame(archipelago = archipelago,
                                   connection = "coalescent",
                                   percentage = saturation.final$coalescent[min(which(diff(saturation.final$coalescent)==0))],
                                   gen = saturation.final$gen[min(which(diff(saturation.final$coalescent)==0))],
                                   n = results.final.i$oceanographic_pairs[index],
                                   generation = results.final.i$nbr_gen[index],
                                   AIC = results.final.i$aic.mean.obs.pred[index],
                                   R2 = results.final.i$r2.mean.obs.pred[index],
                                   deltaR2 = results.final.i$r2.diff.null[index],
                                   pval = results.final.i$p.mean.obs.pred[index])
  
  
  summary_archipelago <- rbind(summary_shortest_path,
                               summary_filial,
                               summary_coalescent)
  
  return(summary_archipelago)
  
}