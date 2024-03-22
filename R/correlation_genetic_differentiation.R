correlation <- function(archipelago = NULL, nbr.gen = NULL, matrix.diff = NULL, pop.coordinates = NULL) {
  
  
  library(igraph)
  library(lme4)
  
  
  # Correlation between connectivity outputs (explict, implicit and shortest path) and genetic differentiation (Fst)
  
  ## ------------
  # Multigen
  
  #nbr.gen <- c(seq(1,10,1),seq(15,500,5))
  
  connectivity.final <- NULL
  results.final <- NULL
  
  ## ------------
  # Filial
  
  for (k in 1:length(nbr.gen)) {
    
    print(paste0(" Filial  ---> ",as.character(round((nbr.gen[k]/max(nbr.gen)*100),3)),"%"))
    
    #archipelago <- "tuamotuWest_tuamotuEast_gambier_society_marquesas_austral"
    
    # Load big connectivity matrix
    matrix_name <- paste0("filial_M_",as.character(nbr.gen[k]),"_",archipelago,"_back_temporal_mean_2010_2019",".matrix")
     
    matrix.oc <- read.delim(here::here("outputs","multistep",matrix_name),
                            sep =" ",
                            header = FALSE)
    colnames(matrix.oc) <- c("from","to","weight")
    
    # Fill study connectivity matrix
    
    matrix.oc.study = matrix(0, nrow = length(matrix.diff[1,]), ncol = length(matrix.diff[,1]))
    for (i in 1:length(pop.coordinates$id)) {
      for (j in 1:length(pop.coordinates$id)) {
        
        if (sum(matrix.oc$from == pop.coordinates$id[i] & matrix.oc$to == pop.coordinates$id[j]) == 1) {
          matrix.oc.study[i,j] <- matrix.oc$weight[matrix.oc$from == pop.coordinates$id[i] &
                                                     matrix.oc$to == pop.coordinates$id[j]]
          
        }
      }
    }
    
    colnames(matrix.oc.study) <- rownames(matrix.oc.study) <- pop.coordinates$names
    
    # Directed to undirected
    
    # Direct probabilities to indirect probabilities (Pab == Pba)
    # Formula from Legrand et al., 2022
    
    # Pab = 1
    # Pba = 1
    # 1 - (1-Pab)*(1-Pba)
    # Pab + Pba - Pab*Pba
    
    matrix.oc.study.ab <- matrix.oc.study.ba <- matrix.oc.study
    
    matrix.oc.study.ab[upper.tri(matrix.oc.study, diag = TRUE)] <- 0
    
    matrix.oc.study.ba[lower.tri(matrix.oc.study, diag = FALSE)] <- 0
    
    #dist.mean.ind <- 1-(1-dist.mean.ab)*(1-t(dist.mean.ba))
    
    matrix.oc.study <- matrix.oc.study.ab + t(matrix.oc.study.ba) - matrix.oc.study.ab*t(matrix.oc.study.ba) # to deal with 1-low proba issues
    
    matrix.oc.study[upper.tri(matrix.oc.study, diag = FALSE)] <- NA
    
    ## ------------
    
    # Find duplicated 
    
    dup <- which(! duplicated(pop.coordinates$names))
    pop.coordinates.dup <- pop.coordinates[dup,]
    
    matrix.oc.study.dup <- matrix.diff.dup <- matrix(NA, nrow = length(dup), ncol = length(dup))
    rownames(matrix.oc.study.dup) <- colnames(matrix.oc.study.dup) <- pop.coordinates.dup$names
    rownames(matrix.diff.dup) <- colnames(matrix.diff.dup) <- pop.coordinates.dup$names
    
    for (da in 1:length(dup)) {
      for (db in 1:length(dup)) {
        if (da > db ) {
          
          matrix.diff.dup[da,db] <- mean(c(as.numeric(as.matrix(matrix.diff[pop.coordinates$names %in% pop.coordinates$names[dup[da]],pop.coordinates$names %in% pop.coordinates$names[dup[db]]])), # ab
                                           as.numeric(as.matrix(matrix.diff[pop.coordinates$names %in% pop.coordinates$names[dup[db]],pop.coordinates$names %in% pop.coordinates$names[dup[da]]]))), na.rm = TRUE) # ba
          
          matrix.oc.study.dup[da,db] <- mean(c(as.numeric(as.matrix(matrix.oc.study[pop.coordinates$names %in% pop.coordinates$names[dup[da]],pop.coordinates$names %in% pop.coordinates$names[dup[db]]])), # ab
                                               as.numeric(as.matrix(matrix.oc.study[pop.coordinates$names %in% pop.coordinates$names[dup[db]],pop.coordinates$names %in% pop.coordinates$names[dup[da]]]))), na.rm = TRUE) # ba
          
        }
      }
    }
    
    
    ## ------------
    
    # Computing Connectivity vs Differentiation
    
    # Constructing dataframe with cell from, cell to, differentiation, connectivity 
    
    size_matrix <- (ncol(matrix.diff.dup)^2-ncol(matrix.diff.dup))/2
    mat_glmm <- matrix(nrow=size_matrix,ncol=11)
    loop <- 0
    
    for (a in 2:ncol(matrix.diff.dup)) { 
      for (b in 1:(a-1)) {
        
        loop <- loop +1
        
        # Island from
        mat_glmm[loop,1] <- pop.coordinates.dup$names[a]
        # Archipelago from
        mat_glmm[loop,2] <- pop.coordinates.dup$archipelago[a]
        # Lon from
        mat_glmm[loop,3] <- pop.coordinates.dup$Lon[a]
        # Lat from
        mat_glmm[loop,4] <- pop.coordinates.dup$Lat[a]
        # Island to
        mat_glmm[loop,5] <- pop.coordinates.dup$names[b]
        # Archipelago to
        mat_glmm[loop,6] <- pop.coordinates.dup$archipelago[b]
        # Lon to
        mat_glmm[loop,7] <- pop.coordinates.dup$Lon[b]
        # Lat to
        mat_glmm[loop,8] <- pop.coordinates.dup$Lat[b]
        # Differentiation
        mat_glmm[loop,9] <- as.numeric(matrix.diff.dup[a,b])
        # Connectivity mean proba
        mat_glmm[loop,10] <- as.numeric(matrix.oc.study.dup[a,b])
        # Connectivity mean distance
        mat_glmm[loop,11] <- as.numeric(-log(matrix.oc.study.dup[a,b]))
      }
    }
    
    colnames(mat_glmm) <- c('pop_A',"archipelago_A",'Lon_pop_A','Lat_pop_A','pop_B',"archipelago_B",'Lon_pop_B','Lat_pop_B','differentiation','connectivity.mean.proba','connectivity.mean.distance')
    mat_glmm <- as.data.frame(mat_glmm)
    
    connectivity.study <- data.frame(island_from = mat_glmm$pop_A,
                                     archipelago_from = mat_glmm$archipelago_A,
                                     lon_from=as.numeric(mat_glmm$Lon_pop_A),
                                     lat_from=as.numeric(mat_glmm$Lat_pop_A),
                                     island_to = mat_glmm$pop_B,
                                     archipelago_to = mat_glmm$archipelago_B,
                                     lon_to=as.numeric(mat_glmm$Lon_pop_B),
                                     lat_to=as.numeric(mat_glmm$Lat_pop_B),
                                     population = as.factor(mat_glmm$pop_A),
                                     connection = "filial",
                                     nbr_gen = nbr.gen[k],
                                     Differentiation = as.numeric(mat_glmm$differentiation),
                                     Connectivity.mean = as.numeric(mat_glmm$connectivity.mean.proba),
                                     Connectivity.distance = as.numeric(mat_glmm$connectivity.mean.distance),
                                     nbr_path = NA,
                                     stringsAsFactors = FALSE )
    
    # Deal with disconnected site
    connectivity.study$Connectivity.distance[is.infinite(connectivity.study$Connectivity.distance)] <- NA
    
    
    ## ------------
    
    # GLMM
    
    # null model
    model.null <- lme4::lmer(Differentiation ~ 1 + (1|population), connectivity.study, REML=F)
    obs.pred.null <- data.frame(observed=(connectivity.study$Differentiation),predicted=predict(model.null),population=connectivity.study$population)
    model <- lm(observed ~ predicted, data = obs.pred.null)
    aic.obs.pred.null <- AIC(model)
    r2.obs.pred.null = summary(model)$adj.r.squared
    p.obs.pred.null = anova(model)$`Pr(>F)`[1]
    
    # Mean connectivity
    
    model.mean <- lme4::lmer(Differentiation ~ Connectivity.distance + (1|population), connectivity.study, REML=F,na.action = na.omit)
    
    # To deal with NA (infinity values)
    predicted.connectivity.distance <- matrix(NA, nrow = length(connectivity.study$Connectivity.distance), ncol = 1)
    predicted.connectivity.distance[as.numeric(rownames(as.data.frame(predict(model.mean)))),] <- predict(model.mean)
    
    obs.pred <- data.frame(observed=(connectivity.study$Differentiation),predicted=predicted.connectivity.distance,population=connectivity.study$population)
    model <- lm(observed ~ predicted, data = obs.pred)
    aic.mean.obs.pred <- AIC(model)
    r2.mean.obs.pred = summary(model)$adj.r.squared
    p.mean.obs.pred =anova(model)$`Pr(>F)`[1]
    
    # Add to connectivity dataframe
    
    connectivity.study <- cbind(connectivity.study,data.frame(null.obs.pred = predict(model.null),
                                                              mean.obs.pred = predicted.connectivity.distance))
    
    results.study <- data.frame(connection = "filial",
                                nbr_gen = nbr.gen[k],
                                differentiation_pairs = length(predict(model.null)),
                                oceanographic_pairs = sum(!is.na(predicted.connectivity.distance)),
                                r2.diff.null = r2.mean.obs.pred-r2.obs.pred.null, # variance explained by connectivity only
                                aic.obs.pred.null = aic.obs.pred.null, # null model
                                r2.obs.pred.null = r2.obs.pred.null,
                                p.obs.pred.null = p.obs.pred.null,
                                aic.mean.obs.pred = aic.mean.obs.pred, # mean connectivity
                                r2.mean.obs.pred = r2.mean.obs.pred,
                                p.mean.obs.pred= p.mean.obs.pred)
    
    connectivity.final <- rbind(connectivity.final,connectivity.study)
    results.final <- rbind(results.final,results.study)
    
    
  }
  
  ## ------------
  # Coalescent
  
  for (k in 1:length(nbr.gen)) {
    
    print(paste0(" Coalescent  ---> ",as.character(round((nbr.gen[k]/max(nbr.gen)*100),3)),"%"))
    
    # Load big connectivity matrix
    
    matrix_name <- paste0("coalescent_M_",as.character(nbr.gen[k]),"_",archipelago,"_back_temporal_mean_2010_2019",".matrix")
    
    matrix.oc <- read.delim(here::here("outputs","multistep",matrix_name),
                            sep =" ",
                            header = FALSE)
    
    colnames(matrix.oc) <- c("from","to","weight")
    
    # Fill study connectivity matrix
    
    matrix.oc.study = matrix(NA, nrow = length(matrix.diff[1,]), ncol = length(matrix.diff[,1]))
    for (i in 1:length(pop.coordinates$id)) {
      for (j in 1:length(pop.coordinates$id)) {
        
        if (sum(matrix.oc$from == pop.coordinates$id[i] & matrix.oc$to == pop.coordinates$id[j]) == 1) {
          matrix.oc.study[i,j] <- matrix.oc$weight[matrix.oc$from == pop.coordinates$id[i] &
                                                     matrix.oc$to == pop.coordinates$id[j]]
          
        }
      }
    }
    
    colnames(matrix.oc.study) <- rownames(matrix.oc.study) <- pop.coordinates$names
    
    # Directed to undirected
    
    # Direct probabilities to indirect probabilities (Pab == Pba)
    # Formula from Legrand et al., 2022
    
    # Pab = 1
    # Pba = 1
    # 1 - (1-Pab)*(1-Pba)
    # Pab + Pba - Pab*Pba
    
    matrix.oc.study.ab <- matrix.oc.study.ba <- matrix.oc.study
    
    matrix.oc.study.ab[upper.tri(matrix.oc.study, diag = TRUE)] <- 0
    
    matrix.oc.study.ba[lower.tri(matrix.oc.study, diag = FALSE)] <- 0
    
    #dist.mean.ind <- 1-(1-dist.mean.ab)*(1-t(dist.mean.ba))
    
    matrix.oc.study <- matrix.oc.study.ab + t(matrix.oc.study.ba) - matrix.oc.study.ab*t(matrix.oc.study.ba) # to deal with 1-low proba issues
    
    matrix.oc.study[upper.tri(matrix.oc.study, diag = FALSE)] <- NA
    
    ## ------------
    
    # Find duplicated 
    
    dup <- which(! duplicated(pop.coordinates$names))
    pop.coordinates.dup <- pop.coordinates[dup,]
    
    matrix.oc.study.dup <- matrix.diff.dup <- matrix(NA, nrow = length(dup), ncol = length(dup))
    rownames(matrix.oc.study.dup) <- colnames(matrix.oc.study.dup) <- pop.coordinates.dup$names
    rownames(matrix.diff.dup) <- colnames(matrix.diff.dup) <- pop.coordinates.dup$names
    
    for (da in 1:length(dup)) {
      for (db in 1:length(dup)) {
        if (da > db ) {
          
          matrix.diff.dup[da,db] <- mean(c(as.numeric(as.matrix(matrix.diff[pop.coordinates$names %in% pop.coordinates$names[dup[da]],pop.coordinates$names %in% pop.coordinates$names[dup[db]]])), # ab
                                           as.numeric(as.matrix(matrix.diff[pop.coordinates$names %in% pop.coordinates$names[dup[db]],pop.coordinates$names %in% pop.coordinates$names[dup[da]]]))), na.rm = TRUE) # ba
          
          matrix.oc.study.dup[da,db] <- mean(c(as.numeric(as.matrix(matrix.oc.study[pop.coordinates$names %in% pop.coordinates$names[dup[da]],pop.coordinates$names %in% pop.coordinates$names[dup[db]]])), # ab
                                               as.numeric(as.matrix(matrix.oc.study[pop.coordinates$names %in% pop.coordinates$names[dup[db]],pop.coordinates$names %in% pop.coordinates$names[dup[da]]]))), na.rm = TRUE) # ba
          
        }
      }
    }
    
    
    ## ------------
    
    # Computing Connectivity vs Differentiation
    
    # Constructing dataframe with cell from, cell to, differentiation, connectivity 
    
    size_matrix <- (ncol(matrix.diff.dup)^2-ncol(matrix.diff.dup))/2
    mat_glmm <- matrix(nrow=size_matrix,ncol=11)
    loop <- 0
    
    for (a in 2:ncol(matrix.diff.dup)) { 
      for (b in 1:(a-1)) {
        
        loop <- loop +1
        
        # Island from
        mat_glmm[loop,1] <- pop.coordinates.dup$names[a]
        # Archipelago from
        mat_glmm[loop,2] <- pop.coordinates.dup$archipelago[a]
        # Lon from
        mat_glmm[loop,3] <- pop.coordinates.dup$Lon[a]
        # Lat from
        mat_glmm[loop,4] <- pop.coordinates.dup$Lat[a]
        # Island to
        mat_glmm[loop,5] <- pop.coordinates.dup$names[b]
        # Archipelago to
        mat_glmm[loop,6] <- pop.coordinates.dup$archipelago[b]
        # Lon to
        mat_glmm[loop,7] <- pop.coordinates.dup$Lon[b]
        # Lat to
        mat_glmm[loop,8] <- pop.coordinates.dup$Lat[b]
        # Differentiation
        mat_glmm[loop,9] <- as.numeric(matrix.diff.dup[a,b])
        # Connectivity mean proba
        mat_glmm[loop,10] <- as.numeric(matrix.oc.study.dup[a,b])
        # Connectivity mean distance
        mat_glmm[loop,11] <- as.numeric(-log(matrix.oc.study.dup[a,b]))
      }
    }
    
    colnames(mat_glmm) <- c('pop_A',"archipelago_A",'Lon_pop_A','Lat_pop_A','pop_B',"archipelago_B",'Lon_pop_B','Lat_pop_B','differentiation','connectivity.mean.proba','connectivity.mean.distance')
    mat_glmm <- as.data.frame(mat_glmm)
    
    connectivity.study <- data.frame(island_from = mat_glmm$pop_A,
                                     archipelago_from = mat_glmm$archipelago_A,
                                     lon_from=as.numeric(mat_glmm$Lon_pop_A),
                                     lat_from=as.numeric(mat_glmm$Lat_pop_A),
                                     island_to = mat_glmm$pop_B,
                                     archipelago_to = mat_glmm$archipelago_B,
                                     lon_to=as.numeric(mat_glmm$Lon_pop_B),
                                     lat_to=as.numeric(mat_glmm$Lat_pop_B),
                                     population = as.factor(mat_glmm$pop_A),
                                     connection = "coalescent",
                                     nbr_gen = nbr.gen[k],
                                     Differentiation = as.numeric(mat_glmm$differentiation),
                                     Connectivity.mean = as.numeric(mat_glmm$connectivity.mean.proba),
                                     Connectivity.distance = as.numeric(mat_glmm$connectivity.mean.distance),
                                     nbr_path = NA,
                                     stringsAsFactors = FALSE )
    
    ## ------------
    # Deal with disconnected site
    connectivity.study$Connectivity.distance[is.infinite(connectivity.study$Connectivity.distance)] <- NA
    
    # GLMM
    
    # null model
    model.null <- lme4::lmer(Differentiation ~ 1 + (1|population), connectivity.study, REML=F)
    obs.pred.null <- data.frame(observed=(connectivity.study$Differentiation),predicted=predict(model.null),population=connectivity.study$population)
    model <- lm(observed ~ predicted, data = obs.pred.null)
    aic.obs.pred.null <- AIC(model)
    r2.obs.pred.null = summary(model)$adj.r.squared
    p.obs.pred.null = anova(model)$`Pr(>F)`[1]
    
    # Mean connectivity
    
    model.mean <- lme4::lmer(Differentiation ~ Connectivity.distance + (1|population), connectivity.study, REML=F,na.action = na.omit)
    
    # To deal with NA (infinity values)
    predicted.connectivity.distance <- matrix(NA, nrow = length(connectivity.study$Connectivity.distance), ncol = 1)
    predicted.connectivity.distance[as.numeric(rownames(as.data.frame(predict(model.mean)))),] <- predict(model.mean)
    
    obs.pred <- data.frame(observed=(connectivity.study$Differentiation),predicted=predicted.connectivity.distance,population=connectivity.study$population)
    model <- lm(observed ~ predicted, data = obs.pred)
    aic.mean.obs.pred <- AIC(model)
    r2.mean.obs.pred = summary(model)$adj.r.squared
    p.mean.obs.pred =anova(model)$`Pr(>F)`[1]
    
    # Add to connectivity dataframe
    
    connectivity.study <- cbind(connectivity.study,data.frame(null.obs.pred = predict(model.null),
                                                              mean.obs.pred = predicted.connectivity.distance))
    
    results.study <- data.frame(connection = "coalescent",
                                nbr_gen = nbr.gen[k],
                                differentiation_pairs = length(predict(model.null)),
                                oceanographic_pairs = sum(!is.na(predicted.connectivity.distance)),
                                r2.diff.null = r2.mean.obs.pred-r2.obs.pred.null, # variance explained by connectivity only
                                aic.obs.pred.null = aic.obs.pred.null, # null model
                                r2.obs.pred.null = r2.obs.pred.null,
                                p.obs.pred.null = p.obs.pred.null,
                                aic.mean.obs.pred = aic.mean.obs.pred, # mean connectivity
                                r2.mean.obs.pred = r2.mean.obs.pred,
                                p.mean.obs.pred= p.mean.obs.pred)
    
    connectivity.final <- rbind(connectivity.final,connectivity.study)
    results.final <- rbind(results.final,results.study)
    
    
  }
  
  
  ## -----------
  # Shortest path
  
  matrix_name <- paste0("filial_M_1_",archipelago,"_back_temporal_mean_2010_2019",".matrix")
  
  matrix.oc <- read.delim(here::here("outputs","multistep",matrix_name),
                          sep =" ",
                          header = FALSE)

  colnames(matrix.oc) <- c("from","to","weight")
  
  # Create graph object from list
  
  graph.back.save <- graph_from_data_frame(matrix.oc,
                                           directed = TRUE,
                                           vertices = NULL)
  
  matrix.oc$weight <- -log(matrix.oc$weight)
  
  graph.back <- graph_from_data_frame(matrix.oc,
                                      directed = TRUE,
                                      vertices = NULL)
  
  dist.mean <- matrix(NA, nrow = length(pop.coordinates$id), ncol = length(pop.coordinates$id))
  path_nbr.mean <- matrix(NA, nrow = length(pop.coordinates$id), ncol = length(pop.coordinates$id))
  
  for (a in 1:length(pop.coordinates$id)) {
    for (b in 1:length(pop.coordinates$id)) {
      
      if (a != b &
          (sum(matrix.oc$from == pop.coordinates[a,"id"]) * sum(matrix.oc$to == pop.coordinates[b,"id"])) > 0 ) { # add an if condition on precense of the vertice in the graph)
        
        path <- shortest_paths(
          graph.back,
          from = as.character(pop.coordinates[a,"id"]),
          to = as.character(pop.coordinates[b,"id"]),
          mode = "out",
          weights = NULL,
          algorithm = c("dijkstra"))
        
        path <- as.numeric(path$vpath[[1]])
        
        if (length(path) > 0) { # it means there is a path
          
          EP = rep(path, each=2)[-1]
          EP = EP[-length(EP)]
          
          E(graph.back.save)[get.edge.ids(graph.back.save,(EP))] # In grpah.back.save is the probabilities
          dist.mean[a,b] <- prod(E(graph.back.save)$weight[get.edge.ids(graph.back.save,(EP))],na.rm = TRUE)
          path_nbr.mean[a,b] <- length(path)-1
          
        } else {
          
          dist.mean[a,b] <- 0
          path_nbr.mean[a,b] <- Inf
          
        }
        
        
      }
      
    }
  }
  
  
  
  dist.mean.ab <- dist.mean.ba <- dist.mean
  
  dist.mean.ab[upper.tri(dist.mean, diag = TRUE)] <- 0
  
  dist.mean.ba[lower.tri(dist.mean, diag = FALSE)] <- 0
  
  #dist.mean.ind <- 1-(1-dist.mean.ab)*(1-t(dist.mean.ba))
  
  dist.mean <- dist.mean.ab + t(dist.mean.ba) - dist.mean.ab*t(dist.mean.ba) # to deal with 1-low proba issues
  
  dist.mean[upper.tri(dist.mean, diag = FALSE)] <- NA
  
  # Find the shortest number of path
  
  path_nbr.mean.min <- path_nbr.mean
  for (a in 1:length(pop.coordinates$id)) {
    for (b in 1:length(pop.coordinates$id)) {
      if (a !=b) {
        path_nbr.mean.min[a,b] <- min(c(path_nbr.mean[a,b],path_nbr.mean[b,a]),na.rm = TRUE)
      }
    }
  }
  
  path_nbr.mean.min[upper.tri(path_nbr.mean.min, diag = TRUE)] <- NA
  path_nbr.mean <- path_nbr.mean.min
  
  ## ------------
  
  # Find duplicated 
  
  dup <- which(! duplicated(pop.coordinates$names))
  pop.coordinates.dup <- pop.coordinates[dup,]
  
  path_nbr.mean.dup <- dist.mean.dup <- matrix.diff.dup <- matrix(NA, nrow = length(dup), ncol = length(dup))
  rownames(dist.mean.dup) <- colnames(dist.mean.dup) <- pop.coordinates.dup$names
  rownames(path_nbr.mean.dup) <- colnames(path_nbr.mean.dup) <- pop.coordinates.dup$names
  rownames(matrix.diff.dup) <- colnames(matrix.diff.dup) <- pop.coordinates.dup$names
  
  for (da in 1:length(dup)) {
    for (db in 1:length(dup)) {
      if (da > db ) {
        
        matrix.diff.dup[da,db] <- mean(c(as.numeric(as.matrix(matrix.diff[pop.coordinates$names %in% pop.coordinates$names[dup[da]],pop.coordinates$names %in% pop.coordinates$names[dup[db]]])), # ab
                                         as.numeric(as.matrix(matrix.diff[pop.coordinates$names %in% pop.coordinates$names[dup[db]],pop.coordinates$names %in% pop.coordinates$names[dup[da]]]))), na.rm = TRUE) # ba
        
        dist.mean.dup[da,db] <- mean(c(as.numeric(as.matrix(dist.mean[pop.coordinates$names %in% pop.coordinates$names[dup[da]],pop.coordinates$names %in% pop.coordinates$names[dup[db]]])), # ab
                                       as.numeric(as.matrix(dist.mean[pop.coordinates$names %in% pop.coordinates$names[dup[db]],pop.coordinates$names %in% pop.coordinates$names[dup[da]]]))), na.rm = TRUE) # ba
        
        path_nbr.mean.dup[da,db] <- mean(c(as.numeric(as.matrix(path_nbr.mean[pop.coordinates$names %in% pop.coordinates$names[dup[da]],pop.coordinates$names %in% pop.coordinates$names[dup[db]]])), # ab
                                           as.numeric(as.matrix(path_nbr.mean[pop.coordinates$names %in% pop.coordinates$names[dup[db]],pop.coordinates$names %in% pop.coordinates$names[dup[da]]]))), na.rm = TRUE) # ba
        
      }
    }
  }
  
  
  ## ------------
  
  # Computing Connectivity vs Differentiation
  
  # Constructing dataframe with cell from, cell to, differentiation, connectivity 
  
  size_matrix <- (ncol(matrix.diff.dup)^2-ncol(matrix.diff.dup))/2
  mat_glmm <- matrix(nrow=size_matrix,ncol=12)
  loop <- 0
  
  for (a in 2:ncol(matrix.diff.dup)) { 
    for (b in 1:(a-1)) {
      
      loop <- loop +1
      
      # Island from
      mat_glmm[loop,1] <- pop.coordinates.dup$names[a]
      # Archipelago from
      mat_glmm[loop,2] <- pop.coordinates.dup$archipelago[a]
      # Lon from
      mat_glmm[loop,3] <- pop.coordinates.dup$Lon[a]
      # Lat from
      mat_glmm[loop,4] <- pop.coordinates.dup$Lat[a]
      # Island to
      mat_glmm[loop,5] <- pop.coordinates.dup$names[b]
      # Archipelago to
      mat_glmm[loop,6] <- pop.coordinates.dup$archipelago[b]
      # Lon to
      mat_glmm[loop,7] <- pop.coordinates.dup$Lon[b]
      # Lat to
      mat_glmm[loop,8] <- pop.coordinates.dup$Lat[b]
      # Differentiation
      mat_glmm[loop,9] <- as.numeric(matrix.diff.dup[a,b])
      # Connectivity mean proba
      mat_glmm[loop,10] <- as.numeric(dist.mean.dup[a,b])
      # Connectivity mean distance
      mat_glmm[loop,11] <- as.numeric(-log(dist.mean.dup[a,b]))
      # Connectivity nbr of step
      mat_glmm[loop,12] <- as.numeric(path_nbr.mean.dup[a,b])
    }
  }
  
  colnames(mat_glmm) <- c('pop_A',"archipelago_A",'Lon_pop_A','Lat_pop_A','pop_B',"archipelago_B",'Lon_pop_B','Lat_pop_B','differentiation','connectivity.mean.proba','connectivity.mean.distance','nbr_step')
  mat_glmm <- as.data.frame(mat_glmm)
  
  connectivity.study <- data.frame(island_from = mat_glmm$pop_A,
                                   archipelago_from = mat_glmm$archipelago_A,
                                   lon_from=as.numeric(mat_glmm$Lon_pop_A),
                                   lat_from=as.numeric(mat_glmm$Lat_pop_A),
                                   island_to = mat_glmm$pop_B,
                                   archipelago_to = mat_glmm$archipelago_B,
                                   lon_to=as.numeric(mat_glmm$Lon_pop_B),
                                   lat_to=as.numeric(mat_glmm$Lat_pop_B),
                                   population = as.factor(mat_glmm$pop_A),
                                   connection = "shortest_path",
                                   nbr_gen = 1,
                                   Differentiation = as.numeric(mat_glmm$differentiation),
                                   Connectivity.mean = as.numeric(mat_glmm$connectivity.mean.proba),
                                   Connectivity.distance = as.numeric(mat_glmm$connectivity.mean.distance),
                                   nbr_path = mat_glmm$nbr_step,
                                   stringsAsFactors = FALSE )
  
  ## ------------
  # Deal with disconnected site
  connectivity.study$Connectivity.distance[is.infinite(connectivity.study$Connectivity.distance)] <- NA
  
  # GLMM
  
  # null model
  model.null <- lme4::lmer(Differentiation ~ 1 + (1|population), connectivity.study, REML=F)
  obs.pred.null <- data.frame(observed=(connectivity.study$Differentiation),predicted=predict(model.null),population=connectivity.study$population)
  model <- lm(observed ~ predicted, data = obs.pred.null)
  aic.obs.pred.null <- AIC(model)
  r2.obs.pred.null = summary(model)$adj.r.squared
  p.obs.pred.null = anova(model)$`Pr(>F)`[1]
  
  # Mean connectivity
  
  model.mean <- lme4::lmer(Differentiation ~ Connectivity.distance + (1|population), connectivity.study, REML=F,na.action = na.omit)
  
  # To deal with NA (infinity values)
  predicted.connectivity.distance <- matrix(NA, nrow = length(connectivity.study$Connectivity.distance), ncol = 1)
  predicted.connectivity.distance[as.numeric(rownames(as.data.frame(predict(model.mean)))),] <- predict(model.mean)
  
  obs.pred <- data.frame(observed=(connectivity.study$Differentiation),predicted=predicted.connectivity.distance,population=connectivity.study$population)
  model <- lm(observed ~ predicted, data = obs.pred)
  aic.mean.obs.pred <- AIC(model)
  r2.mean.obs.pred = summary(model)$adj.r.squared
  p.mean.obs.pred =anova(model)$`Pr(>F)`[1]
  
  # Add to connectivity dataframe
  
  connectivity.study <- cbind(connectivity.study,data.frame(null.obs.pred = predict(model.null),
                                                            mean.obs.pred = predicted.connectivity.distance))
  
  results.study <- data.frame(connection = "shortest_path",
                              nbr_gen = 1,
                              differentiation_pairs = length(predict(model.null)),
                              oceanographic_pairs = sum(!is.na(predicted.connectivity.distance)),
                              r2.diff.null = r2.mean.obs.pred-r2.obs.pred.null, # variance explained by connectivity only
                              aic.obs.pred.null = aic.obs.pred.null, # null model
                              r2.obs.pred.null = r2.obs.pred.null,
                              p.obs.pred.null = p.obs.pred.null,
                              aic.mean.obs.pred = aic.mean.obs.pred, # mean connectivity
                              r2.mean.obs.pred = r2.mean.obs.pred,
                              p.mean.obs.pred= p.mean.obs.pred)
  
  connectivity.final <- rbind(connectivity.final,connectivity.study)
  results.final <- rbind(results.final,results.study)
  results.final$archipelago <- archipelago
  
  ## ------------
  # Save
  
  file_name <- here::here("outputs",paste0("connectivity_",archipelago,'.xlsx'))
  openxlsx::write.xlsx(connectivity.final, file = file_name, rowNames = TRUE, colNames = TRUE)
  
  file_name <- here::here("outputs",paste0("results_",archipelago,'.xlsx'))
  openxlsx::write.xlsx(results.final, file = file_name, rowNames = TRUE, colNames = TRUE)
  
}