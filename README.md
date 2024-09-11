# CMAP

CMAP is an R package that encapsulates python code, for efficiently mapping large-scale individual cells to their precise spatial locations by integrating single-cell and spatial data through a divide­-and-­conquer strategy.
  
<img width="1140" alt="截屏2024-08-08 下午2 18 18" src="https://github.com/user-attachments/assets/4052a5ea-9ad3-467b-a936-b481ba213273">


This approach allows for progressive refinement of spatial assignments, achieving high accuracy with reduced computational complexity. Unlike existing methods, CMAP achieves sub-spot precision, enabling precise mapping of individual cells within tissues.  

## How to install CMAP
To quickly use CMAP, we recommend setting up a new conda environment:  
```conda env create -f CMAP_EnV.environment.yml```

After activating the created environment, follow the instructions in the `Dependency_packages.txt` file to install the necessary dependency packages. If you encouter any questions or errors during the installation of the dependencies, make sure that all packages listed in the `Dependency_packages.txt` file are installed. This file also provides installation commands for the required dependencies.
```
conda activate CMAP_Env
```bash
After that, the CMAP R package can be easily installed from Github using devtools (few seconds):  
```r
devtools::install_github("SuoLab-GZLab/CMAP")
```

## How to run CMAP
All of these steps take a few hours to complete. 

### 1. Load the packages and set the path of python and saved directory
```
library(CMAP) 
library(Seurat) 
library(e1071)
library(purrr)  
library(dplyr)
library(preprocessCore)
library(reticulate)
library(smfishHmrf)
library(Giotto)

python_path <- '/home/your/conda/envs/cmap/bin/python'
use_condaenv(python_path)
save_directory <- "/home/save/directory"
if(!file.exists(save_directory)) dir.create(save_directory, recursive = T)
```
### 2. Load scRNA-seq and ST data, respectively   
We support seurat object and expression matrix, two formats, as input files.

You can directly load the expression matrix and meta information. 
```
sc_count <- read.csv("sc_count.csv",row.names = 1, check.names=FALSE)
st_count <- read.csv("st_count.csv",row.names = 1, check.names=FALSE)
sc_meta <- read.csv("sc_meta.csv",row.names = 1, check.names=FALSE)
`spatial_location` dataframe must be provided two columns (x and y) which are recorded the coordinates of each spot
spatial_location <- read.csv("st_meta.csv",row.names = 1, check.names=FALSE)
```

You can also load the created Seurat objects. Here, we have provided demo datasets for testing (https://www.dropbox.com/scl/fi/q9axwdl8i5ukoctip12ro/CMAP.Demo.Lung_tumor.Data.Rdata?rlkey=an9l3kchva5yn5i80lqg1v0fv&st=jyvymggp&dl=0).
```
load("CMAP.Demo.Lung_tumor.Data.Rdata")
sc_counts <- sc_object@assays$RNA@counts
sc_meta <- data.frame(sc_object@meta.data,row.names=rownames(sc_object@meta.data))

spatial_count <- as.matrix(st_object@assays$Spatial@counts)
spatial_location <- data.frame(st_object@meta.data,row.names=rownames(st_object@meta.data))
# If you create the spatial object using tissue_hires_image.png, you can follow this command to calculate the spots' coordinates
spatial_location <- cbind(spatial_location,
                          x=st_object@images[[1]]@scale.factors$hires * st_object@images[[1]]@coordinates$imagecol,
                          y=-st_object@images[[1]]@coordinates$imagerow * st_object@images[[1]]@scale.factors$hires)
# Else if you create the spatial object using tissue_lowres_image.png, 
spatial_location <- cbind(spatial_location,
                          x=st_object@images[[1]]@scale.factors$lowres * st_object@images[[1]]@coordinates$imagecol,
                          y=-st_object@images[[1]]@coordinates$imagerow * st_object@images[[1]]@scale.factors$lowres)
```

Normalize the data
```
sc_counts <- sc_counts[rowSums(sc_counts)>0,]
sc_norm = as.matrix(log1p(sweep(sc_counts,2,Matrix::colSums(sc_counts),FUN = '/') * 1e4))

spatial_count <- spatial_count[rowSums(spatial_count)>0,]
st_norm = log1p(sweep(spatial_count,2,Matrix::colSums(spatial_count),FUN = '/') * 1e4)
```

### 3. Use HMRF to do spatial clustering
```
# Create specific instructions for Giotto analysis workflow
instrs <- createGiottoInstructions(save_plot = TRUE,
                                   show_plot = TRUE,
                                   return_plot = TRUE,
                                   python_path = python_path,
                                   save_dir = save_directory)
spatial_obj <- createGiottoObject(raw_exprs = spatial_count,
                                  spatial_locs = spatial_location[,c('x','y')],
                                  instructions = instrs,
                                  cell_metadata = spatial_location)
# Filter genes and cells. If you have filtered some low quality spots before, you can skip this step
spatial_obj <- filterGiotto(gobject = spatial_obj,
                            expression_threshold = 1,
                            gene_det_in_min_cells = 50,
                            min_det_genes_per_cell = 250,
                            expression_values = c('raw'),
                            verbose = T)
spatial_obj <- normalizeGiotto(gobject = spatial_obj, scalefactor = 6000, verbose = T)

# Create spatial network
#@ maximum_distance_knn: Visium data, tissue_hires_scalef set ceiling(24.8/tissue_hires_scalef) or as maximum_distance_knn,  tissue_hires_scalef is saved in scalefactors_json.json; slide-seq/ST: set 1.5
spatial_obj <- createSpatialNetwork(gobject = spatial_obj,
                                    method = 'kNN',
                                    k = 6, # this k represents the number of neighbors
                                    maximum_distance_knn = 370, 
                                    minimum_k = 1,
                                    name = 'KNN_network')
kmtest  <- binSpect(spatial_obj, calc_hub = T, hub_min_int = 5,spatial_network_name = 'KNN_network')

hmrf_folder = paste0(save_directory,'/11_HMRF')
if(!file.exists(hmrf_folder)) dir.create(hmrf_folder, recursive = T)
spatial_genes_selected <- hmrf_spatial_gene(spatial_obj,
                                            kmtest,
                                            k=3) # k: Number of spatial domains; set according to your data.

#@ betas: For detailed settings, see https://search.r-project.org/CRAN/refmans/smfishHmrf/html/smfishHmrf.hmrfem.multi.it.min.html
# For quick results, we recomoned setting betas to 45(non-tumor) or 0(tumor sample).
# If you don't mind taking more time and want the best results, you can iteratively test values between 0 and 100 and select the best one.
HMRF_spatial_genes = doHMRF(gobject = spatial_obj,
                            expression_values = 'scaled',
                            spatial_genes = spatial_genes_selected,
                            k = 3, # This value should match the number of spatial domains (k).
                            spatial_network_name="KNN_network",
                            betas = c(0, 45, 2), 
                            python_path = python_path,
                            output_folder = paste0(hmrf_folder, '/', 'Spatial_genes/SG_topgenes_elbow_k_scaled'))
#@betas_to_add: Results from different betas that you want to add
# Recommendations: Tumor sample: beta=0; Non-tumor: beta=45.
spatial_obj = addHMRF(gobject = spatial_obj,
                      HMRFoutput = HMRF_spatial_genes,
                      k = 3,
                      betas_to_add = 0,  # according to the above beta settings
                      hmrf_name = 'HMRF')
# Add spatial domain to spatial metadata. You can also save the spatial_location as an intermediate file, which must include spatial genes and spatial cluster labels.
spatial_location = spatial_location[as.data.frame(spatial_obj@cell_metadata)[,'cell_ID'],]
spatial_location = cbind(spatial_location,HMRF_cluster = spatial_obj@cell_metadata$HMRF_k3_b.0) # this coloumn needs to be set as described above (the number of domains and beta)
st_norm = st_norm[,rownames(spatial_location)]
```

### 4. Level 1 mapping (DomainDivision), dividing cells into different spatial domains
```
matrix <- data_to_transform(sc_norm,st_norm,spatial_genes_selected,batch=TRUE,pca_method='prcomp_irlba')
train_set <- cbind(as.data.frame(t(matrix[,colnames(st_norm)])),label=spatial_location$HMRF_cluster)
test_set <- as.data.frame(t(matrix[,colnames(sc_norm)]))
train_set$label = as.factor(train_set$label)
# Predict spatial domain of individual cells
# This tuning step requires some time. You can adjust the cross-validation proportion using `cross_para` parameter in the `tune_parameter()` function.
parameters <- tune_parameter(train_set, test_set, kernel = "radial", scale = TRUE, class.weight = TRUE, verbose = TRUE)
pred_st_svm <- PredictDomain(train_set, test_set, cost=parameters[['cross_4']][['cost']],
                             gamma=parameters[['cross_4']][['gamma']], st_svm=TRUE,verbose = FALSE)
pred_sc_svm <- PredictDomain(train_set, test_set, cost=parameters[['cross_4']][['cost']],
                             gamma=parameters[['cross_4']][['gamma']], scale = TRUE, verbose = TRUE)
# If there exists unmatched cells with spatial tissue, you need to set a tunable threshold to filter out cells with low mapping probability
sc_meta <- sc_meta[apply(attr(pred_sc_svm, "probabilities"),1,max)>0.8,] 
pred_sc_svm <- pred_sc_svm[apply(attr(pred_sc_svm, "probabilities"),1,max)>0.8]
sc_norm <- sc_norm[,rownames(sc_meta)]
```

### 5. Level 2 mapping (OptimalSpot), globally optimizing the assigned spots of cells within each domain
```
cell_spot_map <- map_cell_to_spot(sc_norm=sc_norm,sc_meta=sc_meta,
                                  st_norm=st_norm,spatial_location=spatial_location,
                                  pred_sc_svm=pred_sc_svm, pred_st_svm=pred_st_svm,
                                  python_path=python_path,
                                  batch=TRUE,
                                  num_epochs=2000L,
                                  para_distance=1.0,
                                  para_density=1.0)
```

### 6. Level 3 mapping (PreciseLocation), giving each cell an exact location
```
spot_neigh_list <- spatial_relation_all(spatial_location,
                                        spatial_data_type=c('honeycomb'))

sc_meta_coord <- calculate_cell_location(cell_spot_map=cell_spot_map,
                                         st_meta =spatial_location,
                                         sc_meta=sc_meta,
                                         sc_norm=sc_norm,
                                         st_norm=st_norm,
                                         batch = TRUE,
                                         spot_neigh_list=spot_neigh_list,
                                         radius = 1/2)
```
The exact location are saved in Column 'pred_loc_x' and 'pred_loc_y' of sc_meta_scoord dataframe


