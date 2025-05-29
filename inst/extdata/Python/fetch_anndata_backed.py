# fetch_anndata_backed.py
# This script provides functions to fetch data from AnnData objects in backed mode, handling _CSCDataset and HDF5-backed arrays.

import pandas as pd
import numpy as np
import warnings

def fetch_anndata_backed(obj, fetch_vars, cells=None, layer=None):
    """
    Fetches the specified variables from an AnnData object in backed mode.
    Handles HDF5-backed matrices (_CSCDataset) for .X and layers.
    Arguments
    ----------
    obj: an AnnData object (opened with backed='r').
    fetch_vars: list of variables (genes/features) to fetch from .X.
    cells: Optional, list of cell names to fetch. If None, fetch all.
    layer: Optional, string specifying the layer to fetch from.
    Returns
    ----------
    A pandas DataFrame with the requested data.
    """
    # Only .X and .layers are supported in backed mode
    if isinstance(fetch_vars, str):
        fetch_vars = [fetch_vars]
    if cells is None:
        cells = list(obj.obs_names)
    # Only support gene access from .X or a specified layer
    if layer is None:
        matrix = obj.X
    else:
        matrix = obj.layers[layer]
    # AnnData in backed mode: matrix is _CSCDataset or similar
    # Get indices for requested genes/features
    var_names = list(obj.var_names)
    gene_indices = [var_names.index(g) for g in fetch_vars if g in var_names]
    missing = [g for g in fetch_vars if g not in var_names]
    if missing:
        warnings.warn(f"The following variables were not found: {', '.join(missing)}")
    # Get indices for requested cells
    cell_names = list(obj.obs_names)
    cell_indices = [cell_names.index(c) for c in cells if c in cell_names]
    # Fetch data as numpy array (cells x genes)
    # _CSCDataset supports slicing
    arr = matrix[cell_indices, :][:, gene_indices]
    # Convert to DataFrame
    df = pd.DataFrame(arr, index=[cell_names[i] for i in cell_indices], columns=[var_names[i] for i in gene_indices])
    return df
