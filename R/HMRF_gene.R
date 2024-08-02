suppressPackageStartupMessages(library("smfishHmrf"))
suppressPackageStartupMessages(library("Giotto"))

#' HMRF genes
#'
#' @param clust
#' @param sample_rate
#' @param target
#' @param seed
#'
#' @return
#' @export
#'
#' @examples
sampling_sp_genes = function(clust,
                             sample_rate=2,
                             target=500,
                             seed = 10){
  # clust = spat_cor_netw_DT$cor_clusters$spat_netw_clus
  # sample_rate=2
  # target=500
  tot=0
  num_cluster=length(unique(clust))
  gene_list = list()

  for(i in seq(1, num_cluster)){
    gene_list[[i]] = colnames(t(clust[which(clust==i)]))
  }
  for(i in seq(1, num_cluster)){
    num_g=length(gene_list[[i]])
    tot = tot+num_g/(num_g^(1/sample_rate))
  }
  factor=target/tot
  num_sample=c()
  genes=c()
  for(i in seq(1, num_cluster)){
    num_g=length(gene_list[[i]])
    genes[i] = num_g
    num_sample[i] = round(num_g/(num_g^(1/sample_rate)) * factor)
  }
  set.seed(seed)
  samples=list()
  union_genes = c()
  for(i in seq(1, num_cluster)){
    if(length(gene_list[[i]])<num_sample[i]){
      samples[[i]] = gene_list[[i]]
    }else{
      samples[[i]] = sample(gene_list[[i]], num_sample[i])
    }
    union_genes = union(union_genes, samples[[i]])
  }
  union_genes = unique(union_genes)

  return(list(union_genes = union_genes, num_sample = num_sample, num_gene=genes, gene_list=gene_list))

}

#' Title
#'
#' @param spatial_obj: hmrf_object
#' @param kmtest
#' @param k : the number of spatial domain
#' @param filter_hbb: TRUE or FALSE
#' @param gene_sampling_from_top
#' @param use_spatial_genes
#' @param filter_method
#' @param gene_samples
#' @param spatial_network_name
#'
#' @return
#' @export
#'
#' @examples
hmrf_spatial_gene = function(spatial_obj,
                             kmtest,
                             k=3, # number of clusters to extract
                             filter_hbb=TRUE,
                             gene_sampling_from_top = 2500,
                             use_spatial_genes = "binSpect",
                             filter_method ="elbow",
                             gene_samples = 500,
                             spatial_network_name = 'KNN_network'){
  all_genes = Giotto::fDataDT(spatial_obj)$gene_ID
  if(filter_hbb){
    all_genes = all_genes[!grepl('^mt-|^Hb[ba]-|^Rp[sl]|^Mrp[sl]',all_genes,ignore.case=TRUE,)]
  }
  result_dt = data.table::data.table(genes=kmtest$genes, binSpect.pval=kmtest$adj.p.value)
  spatial_obj = Giotto::addGeneMetadata(spatial_obj, result_dt, by_column=T, column_gene_ID="genes")
  filtered = Giotto::filterSpatialGenes(spatial_obj, all_genes, max=gene_sampling_from_top,
                                        name=use_spatial_genes, method=filter_method)
  if(filtered$num_genes_removed>0){
    cat(paste0("\n Removed ", filtered$num_genes_removed, " from user's input gene list due to being absent or non-spatial genes.\n"))
    cat(paste0("\n Kept ", length(filtered$genes), " spatial genes for the sampling step next\n"))
  }
  spatial_genes = filtered$genes

  if(length(spatial_genes)==0){
    stop(paste0("\n No genes are remaining to do HMRF. Please give a larger gene list.\n"), call.=FALSE)
  }

  n = min(gene_samples,500, length(spatial_genes))

  #@param gene_sampling_rate parameter (1-50) controlling proportion of gene samples from different module when sampling,
  # 1 corresponding to equal gene samples between different modules;
  # 50 corresponding to gene samples proportional to module size.
  gene_sampling_rate = 2
  gene_sampling_seed = 10
  if(n<length(spatial_genes)){
    cat(paste0("\n Computing spatial coexpression modules...\n"))
    spat_cor_netw_DT = Giotto::detectSpatialCorGenes(gobject = spatial_obj ,method = 'network',
                                                     spatial_network_name = spatial_network_name,
                                                     subset_genes = spatial_genes,
                                                     network_smoothing = 0)

    spat_cor_netw_DT = Giotto::clusterSpatialCorGenes(spat_cor_netw_DT,name = 'spat_netw_clus',k = k)
    cat(paste0("\n Sampling spatial genes from coexpression modules...\n"))
    sample_genes = sampling_sp_genes(spat_cor_netw_DT$cor_clusters$spat_netw_clus,
                                     sample_rate=gene_sampling_rate, target=n, seed=gene_sampling_seed)
    spatial_genes_selected = sample_genes$union_genes
    cat(paste0("\n Sampled ", length(spatial_genes_selected), " genes.\n"))
  }else{spatial_genes_selected = spatial_genes}
  cat(paste0("\n Will use ", length(spatial_genes_selected), " genes for init of HMRF.\n"))

  return(spatial_genes_selected)
}
