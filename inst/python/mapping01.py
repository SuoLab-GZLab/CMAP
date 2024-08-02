import os,sys
import torch
import scanpy as sc
import numpy as np
import pandas as pd
import anndata as ad
from anndata import AnnData
import pathlib
import matplotlib.pyplot as plt
import matplotlib as mpl
import skimage
import seaborn as sns
sys.path.insert(1,os.path.join(os.path.abspath('.'),"pytorch-msssim"))
import pytorch_msssim_me as pytorch_msssim
sys.path.insert(1,os.path.abspath('.'))
import cmap_utils as mpc
import cmap_optimizer as mpt
import logging
from scipy.sparse import csc_matrix
from scipy.sparse import csr_matrix

print(torch.__version__)
print(torch.cuda.is_available())

device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')

def map_cell(sc_data,st_data,num_epochs,para_distance,para_density,random_seed_set=None):
    ad_st = ad.AnnData(st_data)
    ad_sc = ad.AnnData(sc_data)
    ad_map = mpc.map_cell_to_spot(adata_sc=ad_sc, adata_sp=ad_st, 
    	device=device, num_epochs=num_epochs, learning_rate=0.1, 
    	random_seed_set=random_seed_set,
    	para_distance=para_distance, para_density=para_density)
    return ad_map
    
