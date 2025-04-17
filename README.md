# CMAP

CMAP is an R package that encapsulates python code, for efficiently mapping large-scale individual cells to their precise spatial locations by integrating single-cell and spatial data through a divide­-and-­conquer strategy.
  
<img width="1140" alt="截屏2024-08-08 下午2 18 18" src="https://github.com/user-attachments/assets/4052a5ea-9ad3-467b-a936-b481ba213273">


This approach allows for progressive refinement of spatial assignments, achieving high accuracy with reduced computational complexity. Unlike existing methods, CMAP achieves sub-spot precision, enabling precise mapping of individual cells within tissues.  

## How to install CMAP
To quickly use CMAP, we recommend setting up a new conda environment:  
```conda env create -f CMAP_EnV.environment.yml```

After activating the created environment, follow the instructions in the `Dependency_packages.txt` file to install the necessary dependency packages. If you encouter any questions or errors during the installation of the dependencies, make sure that all packages listed in the `Dependency_packages.txt` file are installed. This file also provides installation commands for the required dependencies.
```bash
conda activate CMAP_Env
```
After that, the CMAP R package can be easily installed from Github using devtools (few seconds):  
```r
devtools::install_github("SuoLab-GZLab/CMAP")
```

## How to run CMAP
The entire process requires a few hours to complete. We have made the MOB dataset available for users to reproduce the mapping results. You can download the simulation data from the following link (https://www.dropbox.com/scl/fo/e06uibevb368ej2v97tpp/APeWNuEbkILSyvS0b0dB2kA?rlkey=nxz6481jurujtm6ecqm4tf7v9&st=i9ho1q34&dl=0). The corresponding code (`MOB.CMAP.R` and `MOB.CMAP.ipynb`) is provided for your convenience.

### 1. Load the packages and set the path of python and saved directory
```r
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
### 2. Load scRNA-seq and ST data
#### Input formats supported:
CMAP supports two input formats: a processed Seurat object and an expression matrix.

You can directly load the expression matrix along with its corresponding meta-information. The expression matrix should have genes as rows and cells as columns.

`spatial_location` A dataframe must contain two columns (x and y), representing the spatial coordinates of each spot. 

`sc_meta` A dataframe may include annotated cell type information; however, if this is not provided, it will not affect the performance of the CMAP mapping prediction.

If you are uncertain about the required input file formats, you can download the provided demo datasets and follow the outlined steps to adjust your data to the appropriate format for your analysis.
```r
sc_count <- read.csv("sc_count.csv",row.names = 1, check.names=FALSE)
st_count <- read.csv("st_count.csv",row.names = 1, check.names=FALSE)
sc_meta <- read.csv("sc_meta.csv",row.names = 1, check.names=FALSE)
spatial_location <- read.csv("st_meta.csv",row.names = 1, check.names=FALSE)
```
You can also load the created Seurat objects. Here, we have provided demo datasets for testing (https://www.dropbox.com/scl/fi/q9axwdl8i5ukoctip12ro/CMAP.Demo.Lung_tumor.Data.Rdata?rlkey=an9l3kchva5yn5i80lqg1v0fv&st=jyvymggp&dl=0).
```r
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
```r
sc_counts <- sc_counts[rowSums(sc_counts)>0,]
sc_norm = as.matrix(log1p(sweep(sc_counts,2,Matrix::colSums(sc_counts),FUN = '/') * 1e4))

spatial_count <- spatial_count[rowSums(spatial_count)>0,]
st_norm = log1p(sweep(spatial_count,2,Matrix::colSums(spatial_count),FUN = '/') * 1e4)
```

### 3. Use HMRF to do spatial clustering
The optimal number of spatial domains is determined based on the anatomical features of the tissue. In case where the number of domains is unknown, we assess different possible values and select the number that yields the highest average Silhouette width.
```r
cluster_k <- 3
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
                                            k = cluster_k) # k: Number of spatial domains; set according to your data.

#@ betas: For detailed settings, see https://search.r-project.org/CRAN/refmans/smfishHmrf/html/smfishHmrf.hmrfem.multi.it.min.html
# For quick results, we recomoned setting betas to 45(non-tumor) or 0(tumor sample).
# If you don't mind taking more time and want the best results, you can iteratively test values between 0 and 100 and select the best one.
HMRF_spatial_genes = doHMRF(gobject = spatial_obj,
                            expression_values = 'scaled',
                            spatial_genes = spatial_genes_selected,
                            k = cluster_k, # This value should match the number of spatial domains (k).
                            spatial_network_name="KNN_network",
                            betas = c(0, 45, 2), 
                            python_path = python_path,
                            output_folder = paste0(hmrf_folder, '/', 'Spatial_genes/SG_topgenes_elbow_k_scaled'))
#@betas_to_add: Results from different betas that you want to add
# Recommendations: Tumor sample: beta=0; Non-tumor: beta=45.
beta = 0
spatial_obj = addHMRF(gobject = spatial_obj,
                      HMRFoutput = HMRF_spatial_genes,
                      k = cluster_k,
                      betas_to_add = beta,  # according to the above beta settings
                      hmrf_name = 'HMRF')
# Add spatial domain to spatial metadata. You can also save the spatial_location as an intermediate file, which must include spatial genes and spatial cluster labels.
spatial_location = spatial_location[as.data.frame(spatial_obj@cell_metadata)[,'cell_ID'],]
column <- paste0('HMRF_k',cluster_k,'_b.',beta)
spatial_location = cbind(spatial_location,HMRF_cluster = spatial_obj@cell_metadata[,column]) # this coloumn needs to be set as described above (the number of domains and beta)
st_norm = st_norm[,rownames(spatial_location)]
```

### 4. Level 1 mapping (DomainDivision), dividing cells into different spatial domains
```r
matrix <- data_to_transform(sc_norm,st_norm,spatial_genes_selected,batch=TRUE,pca_method='prcomp_irlba')
train_set <- cbind(as.data.frame(t(matrix[,colnames(st_norm)])),label=spatial_location$HMRF_cluster)
test_set <- as.data.frame(t(matrix[,colnames(sc_norm)]))
train_set$label = as.factor(train_set$label)
# Predict spatial domain of individual cells
# This tuning step requires some time. You can adjust the cross-validation proportion using `cross_para` parameter in the `tune_parameter()` function.
parameters <- tune_parameter(train_set, test_set, kernel = "radial", scale = TRUE, class.weight = TRUE, verbose = TRUE, cross_para=4)
pred_st_svm <- PredictDomain(train_set, test_set, cost=parameters[['cross_4']][['cost']],
                             gamma=parameters[['cross_4']][['gamma']], st_svm=TRUE,verbose = FALSE)
pred_sc_svm <- PredictDomain(train_set, test_set, cost=parameters[['cross_4']][['cost']],
                             gamma=parameters[['cross_4']][['gamma']], scale = TRUE, verbose = TRUE)
```
If there exists unmatched cells with spatial tissue, you need to set **a tunable threshold** to filter out cells with low mapping probability
```r
sc_meta <- sc_meta[apply(attr(pred_sc_svm, "probabilities"),1,max)>0.8,] 
pred_sc_svm <- pred_sc_svm[apply(attr(pred_sc_svm, "probabilities"),1,max)>0.8]
sc_norm <- sc_norm[,rownames(sc_meta)]
```

### 5. Level 2 mapping (OptimalSpot), globally optimizing the assigned spots of cells within each domain
```r
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
```r
spot_neigh_list <- spatial_relation_all(spatial_location,
                                        spatial_data_type=c('honeycomb'))

sc_meta_coord <- calculate_cell_location(cell_spot_map=cell_spot_map,
                                         st_meta =spatial_location,
                                         sc_meta=sc_meta,
                                         sc_norm=sc_norm,
                                         st_norm=st_norm,
                                         parallel = TRUE,                   
                                         batch = TRUE,
                                         spot_neigh_list=spot_neigh_list,
                                         radius = 1)
```
#### Output: 
The predicted spatial locations are saved in the columns pred_loc_x and pred_loc_y of the sc_meta_scoord dataframe. These coordinates can be utilized for downstream analyses.

Plot the spatial distributions of cells
```r
color_use <- c("T cell" = "#CE4D4C",
               "B cell" = "#EBC948", 
               "Mast cell" = "#DDBEAD",
               "Myeloid cell" = "#8C564B", 
               "Cancer cell" = "#5954A4",
               "Epithelial cell" = "#5279BB")    
ggplot(sc_meta_coord,aes(pred_loc_x,pred_loc_y,color=celltype_0916))+
  geom_point(size=0.01)+ 
  theme_bw()+
  scale_color_manual(values = color_use)+
  theme(panel.grid.major=element_line(colour=NA), 
        panel.background = element_rect(fill = "transparent",colour = NA),
        plot.background = element_rect(fill = "transparent",colour = NA),
        panel.grid.minor = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank())+ 
  theme(plot.title = element_text(hjust = 0.5)) + 
  coord_fixed()+
  labs(x=NULL,y=NULL,color= 'Cell type')+
  theme(legend.position='right',
        legend.text=element_text(size=15),
        legend.title=element_text(size=15))+
  guides(color = guide_legend(override.aes = list(size = 4)))
```
<img width="506" alt="截屏2025-01-31 17 16 04" src="https://github.com/user-attachments/assets/47b59537-6fd2-46d9-98da-4dfc826c5cbc" />

### 7. Cell type co-localization analysis
`sc_meta_coord` A dataframe that included cell type annotations.
`cell_type` The name of column containing cell type annotations in the dataframe.
```r
library(doParallel)
library(foreach)
col_ct_df_1 <- celltype_colocalization_count(df=sc_meta_coord,cell_type = "celltype_0916")
cl <- makeCluster(getOption("cl.cores", 20))  
registerDoParallel(cl)
permutate_time = 1000
merge_permu_col_counts <- foreach(i=1:permutate_time) %dopar% {
  processed_data <- permute_cell_coordinates(df=sc_meta_coord,cell_type = "celltype_0916")
  return(processed_data)
}
final_df <- do.call(cbind, merge_permu_col_counts)
stopCluster(cl)
sub_row <- setdiff(rownames(final_df),paste(unique(sc_meta_coord$celltype_0916),unique(sc_meta_coord$celltype_0916),sep = '_'))
final_df_filt <- final_df[sub_row,]
col_ct_df_sub <- col_ct_df_1[sub_row,,drop=FALSE]
final_df <- final_df_filt
col_ct_df_1 <- col_ct_df_sub
local_null_means = apply(final_df,1,mean)
local_null_sd = apply(final_df,1,sd)
local_null_sd = pmax(local_null_sd,sqrt(1/1000))
local_z_scores=(col_ct_df_1-local_null_means)/local_null_sd
local_p_values = pnorm(local_z_scores$number,lower.tail=FALSE)
adjusted_local_p_values = stats::p.adjust(local_p_values,method='fdr')
fold_changes = col_ct_df_1 / (local_null_means + 1e-4)
contact_result_df = cbind(col_ct_df_1,
                          fold_changes=fold_changes$number,
                          p.adj=adjusted_local_p_values,
                          pval=local_p_values)
cmap_cci <- contact_result_df[contact_result_df$p.adj < 0.05 & contact_result_df$number > 50,]
```
Plot the colocalization results
```r
df <- cmap_cci[order(cmap_cci$fold_changes,decreasing = TRUE),,drop=FALSE]
df$Rank <- seq(1,dim(df)[1],1)
df$label <- ""
df[rownames(df)=='B cell_T cell','label'] <- 'B cell - T cell'
ggplot(df,aes(x=Rank,y=fold_changes))+
  geom_point(aes(color=label != ""),size=8)+
  scale_color_manual(values = c("TRUE" = "red", "FALSE" = "lightblue"))+ # lightblue
  geom_text(aes(label = label),hjust = -0.1, vjust = 0.5, size = 9) + 
  theme_bw()+
  theme(text = element_text(color = "black"),
        panel.grid.major = element_blank(),  
        panel.grid.minor = element_blank(), 
        panel.border = element_blank(), 
        axis.line = element_line(color = "black"), 
        axis.title = element_text(color = "black",size = 20),
        axis.text = element_text(color = "black",size = 20),
        axis.ticks = element_line(color = "black"))+
  labs(x='Rank',y='Co-localization score of cell type pairs')+
  theme(legend.position = "none")
```
<img width="481" alt="截屏2025-01-31 17 31 31" src="https://github.com/user-attachments/assets/d2395fad-531f-4655-b9d0-3a24815468ce" />

