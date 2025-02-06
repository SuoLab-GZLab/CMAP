suppressMessages(library(CMAP))
suppressMessages(library(Seurat) )
suppressMessages(library(e1071))
suppressMessages(library(purrr)  )
suppressMessages(library(dplyr))
suppressMessages(library(preprocessCore))
suppressMessages(library(reticulate))
suppressMessages(library(smfishHmrf))
suppressMessages(library(Giotto))
suppressMessages(library(ggplot2))
suppressMessages(library(pryr))
suppressMessages(library(cluster))

python_path <- '/home/kejincan/.conda/envs/vers_2copy/bin/python'
use_condaenv(python_path)

save_directory <- "/home/kejincan/Project/SingleCellMapping/Log/MOB/"
if(!file.exists(save_directory)) dir.create(save_directory, recursive = T)

mem_before <- mem_used()
start_time <- Sys.time()

#arg <- commandArgs(T)
#if(length(arg) < 1){
#  cat("Argument: Tune \n")
#  quit('no')
#}
#tune = as.character(arg[1]) 
tune = TRUE
print(paste0("tune the svm parameter:",tune))

load("/home/kejincan/Project/SingleCellMapping/Data/2022_CARD_Simulated_MOB/Data/20221023.CARD_Simulate_Data.Paired_Cell_Resource.sc_count.sc_meta.spatial_count.spatial_location.Rdata")
spatial_count <- spatial_count[rowSums(spatial_count)>0,]
all(colnames(spatial_count)==rownames(spatial_location))
st_norm = log1p(sweep(spatial_count,2,Matrix::colSums(spatial_count),FUN = '/') * 1e4)

#-------------------- Level 1
spatial_location$x_round <- round(spatial_location$x)
spatial_location$y_round <- round(spatial_location$y)

# Create specific instructions for Giotto analysis workflow
instrs <- createGiottoInstructions(save_plot = TRUE,
                                   show_plot = TRUE,
                                   return_plot = TRUE,
                                   python_path = python_path,
                                   save_dir = save_directory)

spatial_obj <- createGiottoObject(raw_exprs = spatial_count, 
                                  spatial_locs = spatial_location[,c('x_round','y_round')], 
                                  instructions = instrs,
                                  cell_metadata = spatial_location)

# Filter genes and cells. If you have filtered some low quality spots before, you can skip this step
spatial_obj <- filterGiotto(gobject = spatial_obj,
                            expression_threshold = 1,
                            gene_det_in_min_cells = 50,
                            min_det_genes_per_cell = 1000,
                            expression_values = c('raw'),
                            verbose = T)
spatial_obj <- normalizeGiotto(gobject = spatial_obj, scalefactor = 6000, verbose = T)
spatial_obj <- createSpatialNetwork(gobject = spatial_obj,
                                    method = 'kNN',
                                    k = 8, # this k represents the number of neighbors
                                    maximum_distance_knn = 1.5, 
                                    minimum_k = 2,
                                    name = 'KNN_network')
kmtest  <- binSpect(spatial_obj, calc_hub = T, hub_min_int = 5,spatial_network_name = 'KNN_network')

cluster_k <- 3
hmrf_folder = paste0(save_directory,'/11_HMRF_k_',cluster_k)
if(!file.exists(hmrf_folder)) dir.create(hmrf_folder, recursive = T)
set.seed(127)
spatial_genes_selected <- hmrf_spatial_gene(spatial_obj,
                                            kmtest,
                                            k=cluster_k) # k: Number of spatial domains; set according to your data.

HMRF_spatial_genes = doHMRF(gobject = spatial_obj,
                            expression_values = 'scaled',
                            spatial_genes = spatial_genes_selected,
                            k = cluster_k, # This value should match the number of spatial domains (k).
                            spatial_network_name="KNN_network",
                            betas = c(0,45,2), 
                            python_path = python_path,
                            output_folder = paste0(hmrf_folder, '/', 'Spatial_genes/SG_topgenes_elbow_k_scaled'))

#@betas_to_add: Results from different betas that you want to add
# Recommendations: Tumor sample: beta=0; Non-tumor: beta=45.
spatial_obj = addHMRF(gobject = spatial_obj,
                      HMRFoutput = HMRF_spatial_genes,
                      k = cluster_k,
                      betas_to_add = c(0,45),  # according to the above beta settings
                      hmrf_name = 'HMRF')

# Add spatial domain to spatial metadata. You can also save the spatial_location as an intermediate file, whichmust include spatial genes and spatial cluster labels.
spatial_location <- as.data.frame(pDataDT(spatial_obj))
rownames(spatial_location) <- spatial_location$cell_ID

#-------------------------------- Level 2
sc_meta <- sc_meta[sc_meta$cellType!='Ependymal',] 
sc_count <- sc_count[,rownames(sc_meta)]
sc_count <- sc_count[rowSums(sc_count)>0,]

all(colnames(spatial_count)==rownames(spatial_location))
spatial_location$HMRF_cluster <- paste0('cluster_',spatial_location$HMRF_k3_b.45)

all(colnames(sc_count)==rownames(sc_meta))
sc_norm = log1p(sweep(sc_count,2,Matrix::colSums(sc_count),FUN = '/') * 1e4)
sc_meta$cellType <- as.character(sc_meta$cellType)

matrix <- data_to_transform(sc_norm,st_norm,spatial_genes_selected,batch=TRUE,pca_method='prcomp_irlba')
train_set <- cbind(as.data.frame(t(matrix[,colnames(st_norm)])),label=spatial_location$HMRF_cluster)
test_set <- as.data.frame(t(matrix[,colnames(sc_norm)]))
train_set$label = as.factor(train_set$label)
# Step1: predict spatial domain of individual cells
if(tune){
	parameters <- tune_parameter(train_set, test_set, kernel = "radial", scale = TRUE, class.weight = TRUE, verbose = TRUE,cross_para=4)
pred_st_svm <- PredictDomain(train_set, test_set, cost=parameters[['cross_4']][['cost']],
                             gamma=parameters[['cross_4']][['gamma']], st_svm=TRUE,verbose = FALSE)
pred_sc_svm <- PredictDomain(train_set, test_set, cost=parameters[['cross_4']][['cost']],
                             gamma=parameters[['cross_4']][['gamma']], scale = TRUE, verbose = TRUE)
}else{
	pred_st_svm <- PredictDomain(train_set, test_set, st_svm=TRUE,verbose = FALSE)
	pred_sc_svm <- PredictDomain(train_set, test_set, scale = TRUE, verbose = TRUE)
}

sc_meta <- sc_meta[apply(attr(pred_sc_svm, "probabilities"),1,max)>0.5,]
pred_sc_svm <- pred_sc_svm[apply(attr(pred_sc_svm, "probabilities"),1,max)>0.5]
 
# Step2: map cell to individual spot
cell_spot_map <- map_cell_to_spot(sc_norm=sc_norm,sc_meta=sc_meta,
                                  st_norm=st_norm,spatial_location=spatial_location,
                                  pred_sc_svm=pred_sc_svm, pred_st_svm=pred_st_svm,
                                  python_path=python_path,
                                  batch=TRUE,
                                  pca_method = "prcomp",
                                  num_epochs=2000L,
                                  para_distance=1.0,
                                  para_density=1.0)

# Step3: give cell an exact position
#source("/home/kejincan/Project/SingleCellMapping/Data/2023_Breast_cancer/Scripts/20241205.CMAP.Function_update.R")
spatial_location_tmp <- spatial_location
spatial_location_tmp$x <- spatial_location_tmp$x_round
spatial_location_tmp$y <- spatial_location_tmp$y_round
spot_neigh_list <- spatial_relation_all(spatial_location_tmp,
                                        spatial_data_type=c('square'),n_near_spot=5,dis_cut=1.2)
sc_meta_coord <- calculate_cell_location(cell_spot_map=cell_spot_map,
                                         st_meta =spatial_location,
                                         sc_meta=sc_meta,
                                         sc_norm=sc_norm,
                                         st_norm=st_norm,
                                         batch = TRUE,
					parallel = TRUE,
                                         spot_neigh_list=spot_neigh_list,
                                         radius = 1/6)
end_time <- Sys.time()
print("Execute time:")
print(end_time - start_time)

mem_after <- mem_used()
print("Execute occpied memory:")
print(mem_after - mem_before)

save(cell_spot_map,sc_meta_coord,file=paste0(save_directory,'/MOB.Tune_',tune,'.CMAP.Rdata'))
