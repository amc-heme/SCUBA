import mudata as md
# Numpy, Pandas
import pandas as pd
import numpy as np


def fetch_meta_varnames_targeted(obj, use_mudata_obsm, mod):
    """
    fetch_meta_varnames_targeted
    """
    if use_mudata_obsm == True:
        print("Pulling from main obsm slot of mudata object")
    
        keys = obj.obs_keys()
        
        return keys
    
    else:
        if mod == None:
            raise ValueError(
                "mod must be defined if use_mudata_obsm is False."
                )
    
        print("Pulling from obsm slot of", mod, sep = " ")
        
        keys = obj[mod].obs_keys()
        
        return keys


def fetch_meta_varnames_mudata(obj, modal):
    """
    fetch_meta_varnames_mudata
    
    modality (str)
    """
    # Directs fetch_meta_varnames to pull from the obsm slot of the main mudata object
    # if True
    use_mudata_obsm = False

    # Search for modality in the mudata object first
    # If not found, set to use the main object obsm slot
    if modal in obj.mod_names:
        use_mudata_obsm = False
    elif modal == None:
        # Set use_mudata_obsm to True to pull from the main object obsm slot
        use_mudata_obsm = True
    else:
        # Set use_mudata_obsm to True if modality is not found in the object
        print(modal, " modality not found in mudata object. Pulling metadata variable names from main mudata object obsm slot")
        use_mudata_obsm = True
        
    # Fetch reduction from the location determined above
    if use_mudata_obsm == True:
        keys = fetch_meta_varnames_targeted(
            obj, 
            use_mudata_obsm = True, 
            mod = None
            )
    
        return keys
    else:
        keys = fetch_meta_varnames_targeted(
            obj, 
            use_mudata_obsm = False, 
            mod = modal
            )
    
        return keys
