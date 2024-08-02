import numpy as np
import pandas as pd
import scanpy as sc
import torch
import logging
from scipy.sparse import csc_matrix
from scipy.sparse import csr_matrix
import os,sys
sys.path.insert(1,os.path.abspath('.'))
import cmap_optimizer as mpt


logging.getLogger().setLevel(logging.INFO)


def map_cell_to_spot(
    adata_sc,
    adata_sp,
    device="cpu",
    learning_rate=0.1,
    num_epochs=1000,
    verbose=True,
    random_seed_set=None,
    para_distance=1.0,
    para_density=1.0,
):
    """
    Function: CMAP-OptimalSpot. Assign cells to optimal spots.

    Parameters:
        adata_sc: AnnData type. Expression matrix of single cells.
        adata_sp: AnnData type. Expression matrix of spots.
        device: 'cpu' or 'gpu'(torch.device).
        random_seed_set: Int(Default: None). Pass an int  to reproduce the results.
        para_distance: The weight for the SSIM term.
        para_density: The weight for the entropy term.

    Return:
        a cell-by-spot corresponding matrix (AnnData type), containing the probability of mapping cell i to spot j.
    """

    if isinstance(adata_sc.X, csc_matrix) or isinstance(adata_sc.X, csr_matrix):
        C = np.array(adata_sc.X.toarray(), dtype="float32")
    elif isinstance(adata_sc.X, np.ndarray):
        C = np.array(adata_sc.X, dtype="float32")
    else:
        X_type = type(adata_sc.X)
        logging.error("AnnData X has unrecognized type: {}".format(X_type))
        raise NotImplementedError

    if isinstance(adata_sp.X, csc_matrix) or isinstance(adata_sp.X, csr_matrix):
        S = np.array(adata_sp.X.toarray(), dtype="float32")
    elif isinstance(adata_sp.X, np.ndarray):
        S = np.array(adata_sp.X, dtype="float32")
    else:
        X_type = type(adata_sp.X)
        logging.error("AnnData X has unrecognized type: {}".format(X_type))
        raise NotImplementedError

    # Set device for running
    device = torch.device(device)

    if verbose:
        print_each = 100
    else:
        print_each = None

    logging.info("Begin...")
    mapobject = mpt.Mapping(C=C, S=S, device=device, random_seed_set=random_seed_set)
    print(mapobject)
        

    # Call `train` function to obtain the mapping matrix
    mapping_matrix = mapobject.train(
        learning_rate=learning_rate, num_epochs=num_epochs, print_each=print_each, para_distance=para_distance, para_density=para_density,
        )

    logging.info("Saving results..")

    return mapping_matrix



