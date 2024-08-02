#' Title
#'
#' @param X
#' @param MARGIN
#' @param do_safe
#'
#' @return
#' @export
#'
#' @examples
cosine_normalize <- function(X, MARGIN = 1, do_safe = TRUE) {
  # to avoid Inf values, first divide by max
  if (do_safe) {
    X <- sweep(X, MARGIN, apply(X, MARGIN, max), "/")
  }
  sweep(X, MARGIN, apply(X, MARGIN, function(x) sqrt(sum(x^2))), "/")
}

#' Title
#'
#' @param st_norm  gene * spot
#' @param sc_norm
#' @param genes
#' @param batch
#' @param npc
#' @param nclust round(nrow(meta_data)/30)
#' @param max.iter.harmony
#' @param max.iter.cluster
#' @param cosine_norm
#' @param pca_method
#'
#' @return
#' @export
#'
#' @examples
harmony_embedding <- function(st_norm,sc_norm,genes,batch=TRUE,npc = 100,nclust=100,
                              pca_method=c('prcomp','prcomp_irlba'),
                              max.iter.harmony=50,max.iter.cluster=6,cosine_norm=FALSE){

  pca_method = match.arg(pca_method, choices = c('prcomp','prcomp_irlba'))

  st_norm <- st_norm[rowSums(st_norm)>0,,drop=FALSE]
  sc_norm <- sc_norm[rowSums(sc_norm)>0,,drop=FALSE]
  genes <- intersect(intersect(genes,rownames(st_norm)),rownames(sc_norm))
  if(batch){
    tmp_st <- as.matrix(st_norm[genes,,drop=FALSE]) %>% singlecellmethods::scaleData()
    tmp_sc <- as.matrix(sc_norm[genes,,drop=FALSE]) %>% singlecellmethods::scaleData()
    V = cbind(tmp_st,tmp_sc)
  }else{
    tmp_st <- as.matrix(st_norm[genes,,drop=FALSE])
    tmp_sc <- as.matrix(sc_norm[genes,,drop=FALSE])
    V = cbind(tmp_st,tmp_sc) %>% singlecellmethods::scaleData()
  }

  if(cosine_norm){
    V = V %>% cosine_normalize(2)
  }

  meta_data = data.frame(id=c(colnames(st_norm),colnames(sc_norm)),
                         dataset=c(rep('Spatial',dim(st_norm)[2]),rep('SingleCell',dim(sc_norm)[2])))

  # Dimensionality reduction
  # Data has been normalized, so set center=F, scale=F
  exprs_cosine <- V
  if(pca_method=='prcomp_irlba'){
    set.seed(123)
    # retx=TRUE: Returns the data after the rotation: Matrix * Rotated matrix. So the covariance matrix is the diagonal matrix
    pca_res <- irlba::prcomp_irlba(t(exprs_cosine), npc, retx = TRUE, center = FALSE, scale. = FALSE)
  }else{
    set.seed(123)
    # retx=TRUE:
    pca_res <- prcomp(t(exprs_cosine), npc, center = FALSE, scale. = FALSE)
  }

  # Do harmony
  set.seed(123)
  harmonyObj <- harmony::HarmonyMatrix(
    data_mat = pca_res$x, ## PCA embedding matrix of cells
    meta_data = meta_data, ## dataframe with cell labels
    theta = 1, ## cluster diversity enforcement
    vars_use = 'dataset', ## variable to integrate out
    nclust = nclust, ## number of clusters in Harmony model
    max.iter.harmony = max.iter.harmony,
    max.iter.cluster = max.iter.cluster,
    #return_object = TRUE, ## return the full Harmony model object
    do_pca = FALSE, ## don't recompute PCs
    verbose = FALSE
  )

  matrix <- t(rbind(harmonyObj[which(meta_data$dataset=='SingleCell'),,drop=FALSE],
                    harmonyObj[which(meta_data$dataset=='Spatial'),,drop=FALSE]))
  colnames(matrix) <- c(meta_data[meta_data$dataset=='SingleCell','id'],meta_data[meta_data$dataset=='Spatial','id'])

  return(matrix)
}
