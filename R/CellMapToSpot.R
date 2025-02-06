#' Calculate the spatial variable genes
#'
#' @param sc_norm A gene * cell matrix
#' @param sc_meta A dataframe of cells
#' @param st_norm
#' @param spatial_location A dataframe of spots, column `HMRF_cluster` containing the spatial domain label.
#'
#' @return
#' @export
#'
#' @examples
map_spot_genes = function(sc_norm,
                          sc_meta,
                          st_norm,
                          spatial_location){
  # calculate the spatial variable genes of each domain
  st <- Seurat::CreateSeuratObject(counts = st_norm, meta.data = spatial_location,
                                   min.cells = 0, min.features = 0,
                                   assay = "Spatial")
  st@assays$Spatial@data <- as.matrix(st_norm)
  Idents(st) <- st@meta.data$HMRF_cluster
  cluster <- sort(unique(st@meta.data$HMRF_cluster))
  st_gene_list <- list()
  for(i in 1:length(cluster)){
    sub_st <- subset(st, idents = cluster[i])
    sub_st <- sub_st[!apply(sub_st@assays$Spatial@data, 1, var)==0,]
    sub_st <- Seurat::FindVariableFeatures(sub_st, selection.method = "vst", nfeatures = 3000)
    st_gene_list[[i]] <- intersect(rownames(sc_norm),sub_st@assays$Spatial@var.features)
  }

  return(st_gene_list)
}

#' Level 2 mapping (OptimalSpot), globally optimizing the assigned spots of cells within each domain
#'
#' @param sc_norm
#' @param sc_meta
#' @param st_norm
#' @param spatial_location
#' @param batch Whether do batch removal or not.
#' @param st_gene_list The spatial gene lists of each domain
#' @param pca_method Select the PCA function, support `prcomp_irlba` or `prcomp`.
#' @param pred_sc_svm The vector storing the predicted domain label of cells and the vector name is the cell id
#' @param pred_st_svm The vector storing the predicted domain label of spots and the vector name is the spot id
#' @param python_path Provide the path of called python
#' @param num_epochs The number of epochs
#' @param para_distance The weight of expression pattern similarity
#' @param para_density The weight of cellular spatial density
#'
#' @return A dataframe storing the corresponding relationship between cells and spots.
#' @export
#'
#' @examples
map_cell_to_spot <- function(sc_norm,
                             sc_meta,
                             st_norm,
                             spatial_location,
                             batch = TRUE,
                             st_gene_list = NULL,
                             pca_method = "prcomp_irlba",
                             pred_sc_svm,
                             pred_st_svm,
                             python_path,
                             num_epochs = 2000L,
                             para_distance = 1.0,
                             para_density = 1.0){
  if(is.null(st_gene_list)){
    st_gene_list <- map_spot_genes(sc_norm = sc_norm, sc_meta = sc_meta,
                                   st_norm = st_norm, spatial_location = spatial_location)
    shvgs <- unique(unlist(st_gene_list))
  }
  cluster <- sort(unique(pred_sc_svm))
  cell_spot_map <- data.frame()
  for (j in 1:length(cluster)) {
    cluster_id = cluster[j]
    single_cell_id <- names(pred_sc_svm[pred_sc_svm == cluster_id])
    spot_id <- names(pred_st_svm[pred_st_svm == cluster_id])
    sub_embedding <- data_to_transform(sc_norm = sc_norm[, single_cell_id, drop = FALSE],
                                       st_norm = st_norm[, spot_id, drop = FALSE],
                                       spatial_genes = unique(st_gene_list[[j]]),
                                       batch = batch,
                                       pca_method = pca_method)
    cellxgene <- as.data.frame(t(sub_embedding)[single_cell_id, , drop = FALSE])
    spotxgene <- as.data.frame(t(sub_embedding)[spot_id, , drop = FALSE])
    reticulate::use_condaenv(python_path)
    reticulate::py_run_string(paste0("import sys; sys.path.append('", system.file("python", package = "CMAP"), "')"))
    reticulate::source_python(paste(system.file("python", package = "CMAP"), "mapping01.py", sep = "/"))
    cell_spot <- map_cell(cellxgene, spotxgene, num_epochs, para_distance, para_density)
    cell_spot_map_sub <- data.frame(Single_cell = rownames(cellxgene),
                                    Spot = rownames(spotxgene)[apply(cell_spot, 1, which.max)],
                                    Probability = apply(cell_spot, 1, max))
    cell_spot_map <- rbind(cell_spot_map, cell_spot_map_sub)
  }
  return(cell_spot_map)
}
