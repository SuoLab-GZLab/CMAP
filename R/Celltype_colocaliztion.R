#' Generate a standard dataframe to store the number of cell-type colocalizations
#'
#' @param df A dataframe of cells, containing position and annotation information.
#' @param cell_type The name of column that stores cell type annotations.
#'
#' @return A empty dataframe designed to store the number of cell-type colocalizations.
#' @export
#'
#' @examples
cell_type_colocalization = function(df,
                                    cell_type='celltype_0916'
){
  cell_type = names(table(df[,cell_type]))
  ct_vec = c()
  for(i in 1:length(cell_type)){
    a <- cell_type[i]
    remain_ct <- cell_type[i:length(cell_type)]
    for(j in 1:length(remain_ct)){
      ct_vec <- c(ct_vec,paste(a,remain_ct[j],sep = '_'))
    }
  }
  col_ct_all <- data.frame(row.names = ct_vec,number=rep(0,length(ct_vec)))
  return(col_ct_all)
}


#' Count the number of cell-type colocalizations
#'
#' @param df A dataframe of cells, containing position and annotation information.
#' @param cut_distance The cutoff distance used to define cell colocaliztion.
#' @param x The name of the column storing the x-coordinate values.
#' @param y The name of the column storing the y-coordinate values.
#' @param cell_type The name of column storing cell type annotations.
#'
#' @return A dataframe storing the number of cell-type colocalizations.
#' @export
#'
#' @examples
celltype_colocalization_count = function(df,
                                         cut_distance=2.5,
                                         x='pred_loc_x',
                                         y='pred_loc_y',
                                         cell_type='celltype_0916'
){
  dist_mat <- as.matrix(proxy::dist(df[,c(x,y)],method = 'Euclidean'))
  dist_mat_bool <- (dist_mat < cut_distance)
  col_ct_tab <- reshape2::melt(dist_mat_bool, na.rm=TRUE)
  col_ct_tab <- col_ct_tab[col_ct_tab$value,]
  col_ct_tab <- data.frame(ct1=df[col_ct_tab$Var1,cell_type],
                           ct2=df[col_ct_tab$Var2,cell_type])
  df_column <- apply(col_ct_tab,1,function(x) {
    if(x['ct1']<=x['ct2']) {
      paste(x['ct1'],x['ct2'],sep = '_')
    }else {
      paste(x['ct2'],x['ct1'],sep = '_')
    }
  })
  # Number of cell-type colocalization
  col_ct_df = cell_type_colocalization(df=df,cell_type=cell_type)
  col_ct_df[names(table(df_column)),'number'] = as.vector(table(df_column))
  return(col_ct_df)
}

#' Randomly shift the spatial localtion of each cell to a random position within 100 μm from its original location
#'
#' @param radius_permute The radius of the random shift distance; should be adjusted according to different datasets.
#' @param df A dataframe of cells, containing position and annotation information.
#' @param x The name of the column storing the x-coordinate values.
#' @param y The name of the column storing the y-coordinate values.
#' @param cut_distance The cutoff distance used to define cell colocaliztion.
#' @param cell_type The name of column storing cell type annotations.
#'
#' @return A dataframe with the randomly shifted cell positions and corresponding colocalization counts.
#'
#' Reference: zhang et al. (2023): https://github.com/ZhuangLab/whole_mouse_brain_MERFISH_atlas_scripts_2023/blob/main/scripts/cell_cell_contacts/randomize_and_count_cell_cell_contacts_15um.ipynb
#'
#' @export
#'
#' @examples
#' In the lung cancer datasets, the coordinates are scaled such that a distance of 25 equals 100μm.
permute_cell_coordinates = function(radius_permute=25,
                                      df,
                                      x='pred_loc_x',
                                      y='pred_loc_y',
                                      cut_distance=2.5,
                                      cell_type='celltype_0916'
){
  df_rand = df
  r = radius_permute * sqrt(runif(dim(df_rand)[1]))
  theta = runif(dim(df_rand)[1]) * 2 * pi
  df_rand[,x] = df_rand[,x] + r * cos(theta)
  df_rand[,y] = df_rand[,y] + r * sin(theta)
  col_ct_df = celltype_colocalization_count(df=df_rand,
                                            cut_distance=cut_distance,
                                            x=x,y=y,cell_type=cell_type
  )
  return(col_ct_df)
}
