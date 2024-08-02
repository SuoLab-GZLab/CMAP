#' Calculate the adjacent spots of central spots
#'
#' @param spatial_location
#' @param spatial_data_type
#' @param n_near_spot
#' @param dis_cut
#'
#' @return
#' @export
#'
#' @examples
spatial_relation_all = function(spatial_location,
                                spatial_data_type=c('honeycomb','square','singlecell_level'),
                                n_near_spot=NULL,
                                dis_cut=NULL){
  spatial_data_type = match.arg(spatial_data_type, choices = c('honeycomb','square','singlecell_level'))

  if(spatial_data_type=='honeycomb'){
    top_n = 7
    cutoff = 30 # normally, 23.9/24.8/25.1
  }else if(spatial_data_type=='square'){
    top_n = 5
    cutoff = 30
  }else if(spatial_data_type=='singlecell_level'){
    top_n = 5
    cutoff = 10
  }

  if(!is.null(dis_cut)){
    cutoff = dis_cut
  }

  if(!is.null(n_near_spot)){
    top_n = n_near_spot
  }

  spot_dis <- rdist::pdist(spatial_location[,c('x','y')])
  rownames(spot_dis) <- rownames(spatial_location)
  colnames(spot_dis) <- rownames(spatial_location)
  # The following results are ordered(from near to far), column is spot_id
  spot_dis_order = apply(spot_dis,1,function(x){sort(x)})
  spot_dis_order_id = apply(spot_dis,1,function(x){names(sort(x))})
  tmp <- spot_dis
  spot_neigh <- list()
  for(i in 1:dim(tmp)[1]){
    spot_neigh[[i]] <- spot_dis_order_id[,i][1:top_n][spot_dis_order[1:top_n,i] < cutoff]
  }
  df_j_all <- spot_neigh # include itself
  names(df_j_all) <- colnames(tmp)

  df_j_up_down=list()
  df_j_left_right=list()

  if(spatial_data_type=='singlecell_level'){
    # such as slide-seq data
    df_j_up_down <- df_j_all
    df_j_left_right <- df_j_all
    spot_neigh_list <- list(all=df_j_all,
                            left_right=df_j_left_right,
                            up_down=df_j_up_down)
    return(spot_neigh_list)
  }


  for(spot_i in names(df_j_all)){
    if(length(df_j_all[[spot_i]])==1){
      df_j_up_down[[spot_i]] <- c()
      df_j_left_right[[spot_i]] <- c()
      next
    }

    near_spot_id <- spatial_location[setdiff(df_j_all[[spot_i]],spot_i),c('x','y'),drop=FALSE]
    central_spot <- spatial_location[spot_i,c('x','y'),drop=FALSE]

    if(any(round(near_spot_id$y) %in% c(round(central_spot$y),ceiling(central_spot$y),floor(central_spot$y)))){
      df_j_left_right[[spot_i]] <- rownames(near_spot_id)[which(round(near_spot_id$y) %in% c(round(central_spot$y),ceiling(central_spot$y),floor(central_spot$y)))]
      df_j_up_down[[spot_i]] <- setdiff(rownames(near_spot_id),df_j_left_right[[spot_i]])
    }else if(any(round(near_spot_id$x) %in% c(round(central_spot$x),ceiling(central_spot$x),floor(central_spot$x)))){
      df_j_up_down[[spot_i]] <- rownames(near_spot_id)[which(round(near_spot_id$x) %in% c(round(central_spot$x),ceiling(central_spot$x),floor(central_spot$x)))]
      df_j_left_right[[spot_i]] <- setdiff(rownames(near_spot_id),df_j_up_down[[spot_i]])
    }else{
      df_j_left_right[[spot_i]] <- rownames(near_spot_id)
      df_j_up_down[[spot_i]] <- rownames(near_spot_id)
    }

  }

  spot_neigh_list <- list(all=df_j_all,
                          left_right=df_j_left_right,
                          up_down=df_j_up_down)

  return(spot_neigh_list)
}

#' neighbor spots used in Spring model
#'
#' @param top_sim_spot
#'
#' @return
#' @export
#'
#' @examples
neighbor_spot <- function(spot_neigh_list,
                          top_sim_spot){
  df_j_up_down <- spot_neigh_list$up_down
  df_j_left_right <- spot_neigh_list$left_right
  df_j_all <- spot_neigh_list$all # include this spot itself

  if(length(df_j_up_down[[top_sim_spot]])!=0 & length(df_j_left_right[[top_sim_spot]])!=0){
    neigh_UD <- df_j_up_down[[top_sim_spot]]
    neigh_LR <- df_j_left_right[[top_sim_spot]]
  }else if(length(df_j_up_down[[top_sim_spot]])==0 & length(df_j_left_right[[top_sim_spot]])!=0){
    neigh_UD <- top_sim_spot
    neigh_LR <- df_j_left_right[[top_sim_spot]]
  }else if(length(df_j_up_down[[top_sim_spot]])!=0 & length(df_j_left_right[[top_sim_spot]])==0){
    neigh_UD <- df_j_up_down[[top_sim_spot]]
    neigh_LR <- top_sim_spot
  }else{
    neigh_LR <- top_sim_spot
    neigh_UD <- top_sim_spot
  }
  neigh <- list(neigh_UD = neigh_UD, neigh_LR = neigh_LR)
  return(neigh)
}

#' Title
#'
#' @param x
#'
#' @return
#' @export
#'
#' @examples
expand_diff <- function(x){
  scale_x <- (x-min(x))/(max(x)-min(x))
  tanh_x <- 2/(1+exp(-2*scale_x))-1
  fx <- tanh_x + var(x)
  return(fx)
}

#' Calculate the correlation between spots and cells
#'
#' @param sc_norm
#' @param st_norm
#' @param genes
#' @param batch
#'
#' @return
#' @export
#'
#' @examples
cor_spot_cell <- function(sc_norm,
                          st_norm,
                          genes,
                          batch = TRUE
){
  matrix <- data_to_transform(sc_norm = sc_norm,
                              st_norm = st_norm,
                              spatial_genes = genes,
                              batch = batch)
  Cl <- cor(as.matrix(matrix),method="spearman")
  cor_sub <- Cl[1:dim(sc_norm)[2],(dim(sc_norm)[2]+1):ncol(Cl),drop=FALSE]

  return(cor_sub)
}


#' Title
#'
#' @param cell_spot_map
#' @param st_meta
#' @param sc_meta
#' @param sc_norm
#' @param st_norm
#' @param batch
#' @param spot_neigh_list
#' @param radius
#'
#' @return
#' @export
#'
#' @examples
calculate_cell_location <- function(cell_spot_map,
                                    st_meta,
                                    sc_meta,
                                    sc_norm,
                                    st_norm,
                                    batch = TRUE,
                                    spot_neigh_list,
                                    radius = 1/2){
  st_gene_list <- map_spot_genes(sc_norm = sc_norm,
                                 sc_meta = sc_meta,
                                 st_norm = st_norm,
                                 spatial_location = st_meta)
  shvgs <- unique(unlist(st_gene_list))
  cell_spot_cor <- cor_spot_cell(sc_norm = sc_norm,
                                 st_norm = st_norm,
                                 genes = shvgs,
                                 batch = batch)
  cell_id  <- cell_spot_map$Single_cell
  sc_meta_coord <- sc_meta[cell_id,,drop=FALSE]
  sc_meta_coord$pred_loc_x <- NA
  sc_meta_coord$pred_loc_y <- NA
  df_j_all <- spot_neigh_list$all

  for(i in 1:length(cell_id)){
    top_sim_spot = as.character(cell_spot_map[cell_spot_map$Single_cell==cell_id[i],,drop=FALSE]$Spot)
    if(length(top_sim_spot)==1){
      top_spot_neighbors = df_j_all[[top_sim_spot]]
      if(length(top_spot_neighbors)>1){
        neighbor <- neighbor_spot(spot_neigh_list,top_sim_spot)
        neigh_LR <- neighbor$neigh_LR
        neigh_UD <- neighbor$neigh_UD
        rela_LR <- reshape2::melt(as.matrix(cell_spot_cor[cell_id[i],neigh_LR,drop=FALSE]), na.rm=TRUE)
        coord_LR = st_meta[as.character(rela_LR$Var2),,drop=FALSE]
        rela_UD <- reshape2::melt(as.matrix(cell_spot_cor[cell_id[i],neigh_UD,drop=FALSE]), na.rm=TRUE)
        coord_UD = st_meta[as.character(rela_UD$Var2),,drop=FALSE]
        n_LR = length(rela_LR$value)
        n_UD = length(rela_UD$value)
        rela_top_spot <- reshape2::melt(as.matrix(cell_spot_cor[cell_id[i],top_sim_spot,drop=FALSE]), na.rm=TRUE)
        coord_top_spot = st_meta[as.character(rela_top_spot$Var2),,drop=FALSE]
        # Do Min-max normalization and non-linear transforamtion(tanh) of k
        coeff_expd <- expand_diff(c(rela_LR$value,rela_UD$value,rela_top_spot$value))
        sc_meta_coord$pred_loc_x[i] <- (sum(coeff_expd[1:n_LR] * coord_LR$x) + coeff_expd[length(coeff_expd)] * coord_top_spot$x) / (sum(coeff_expd[1:n_LR])+coeff_expd[length(coeff_expd)])
        sc_meta_coord$pred_loc_y[i] <- (sum(coeff_expd[(n_LR+1):(n_LR+n_UD)] * coord_UD$y) + coeff_expd[length(coeff_expd)] * coord_top_spot$y) / (sum(coeff_expd[(n_LR+1):(n_LR+n_UD)])+coeff_expd[length(coeff_expd)])
      }else{
        sc_meta_coord$pred_loc_x[i] <- st_meta[top_sim_spot,'x'] + runif(1,min=-radius,max=radius)
        sc_meta_coord$pred_loc_y[i] <- st_meta[top_sim_spot,'y'] + runif(1,min=-radius,max=radius)
      }
    }else if(length(top_sim_spot)>1){
      tmp <- top_sim_spot
      for(k in 1:length(tmp)){
        top_sim_spot  = tmp[k]
        top_spot_neighbors = df_j_all[[top_sim_spot]]
        if(length(top_spot_neighbors)>1){
          neighbor <- neighbor_spot(neighbor_spot,top_sim_spot)
          neigh_LR <- neighbor$neigh_LR
          neigh_UD <- neighbor$neigh_UD
          rela_LR <- reshape2::melt(as.matrix(cell_spot_cor[cell_id[i],neigh_LR,drop=FALSE]), na.rm=TRUE)
          coord_LR = st_meta[as.character(rela_LR$Var2),,drop=FALSE]
          rela_UD <- reshape2::melt(as.matrix(cell_spot_cor[cell_id[i],neigh_UD,drop=FALSE]), na.rm=TRUE)
          coord_UD = st_meta[as.character(rela_UD$Var2),,drop=FALSE]
          n_LR = length(rela_LR$value)
          n_UD = length(rela_UD$value)
          rela_top_spot <- reshape2::melt(as.matrix(cell_spot_cor[cell_id[i],top_sim_spot,drop=FALSE]), na.rm=TRUE)
          coord_top_spot = st_meta[as.character(rela_top_spot$Var2),,drop=FALSE]
          # Do Min-max normalization and non-linear transforamtion(tanh) of k
          coeff_expd <- expand_diff(c(rela_LR$value,rela_UD$value,rela_top_spot$value))

          new_cell <- sc_meta_coord[i,,drop=FALSE]
          new_cell$pred_loc_x <- (sum(coeff_expd[1:n_LR] * coord_LR$x) + coeff_expd[length(coeff_expd)] * coord_top_spot$x) / (sum(coeff_expd[1:n_LR])+coeff_expd[length(coeff_expd)])
          new_cell$pred_loc_y <- (sum(coeff_expd[(n_LR+1):(n_LR+n_UD)] * coord_UD$y) + coeff_expd[length(coeff_expd)] * coord_top_spot$y) / (sum(coeff_expd[(n_LR+1):(n_LR+n_UD)])+coeff_expd[length(coeff_expd)])
          rownames(new_cell) <- paste0(cell_id[i],'_',k)
          sc_meta_coord <- rbind(sc_meta_coord,new_cell)
        }else{
          new_cell <- sc_meta_coord[i,,drop=FALSE]
          new_cell$pred_loc_x <- st_meta[top_sim_spot,'x'] + runif(1,min=-radius,max=radius)
          new_cell$pred_loc_y <- st_meta[top_sim_spot,'y'] + runif(1,min=-radius,max=radius)
          rownames(new_cell) <- paste0(cell_id[i],'_',k)
          sc_meta_coord <- rbind(sc_meta_coord,new_cell)
        }
      }
    }
  }
  sc_meta_coord <- sc_meta_coord[!is.na(sc_meta_coord$pred_loc_x),,drop=FALSE]
  return(sc_meta_coord)
}
