import mudata as md
import pandas as pd
import numpy as np

def fetch_cells_mudata(obj, meta_var, meta_levels, modality=None):
    """
    Fetch cell IDs from a MuData object where obs[meta_var] is in meta_levels.
    Only use the provided modality to access the mudata object. meta_var is used as-is.
    """
    if modality is not None:
        obs = obj[modality].obs
        obs_names = obj[modality].obs_names
    else:
        obs = obj.obs
        obs_names = obj.obs_names
    if meta_var not in obs.columns:
        raise ValueError(f"meta_var '{meta_var}' not found in obs.")
    mask = obs[meta_var].isin(meta_levels)
    selected_cells = obs_names[mask].tolist()
    return selected_cells
