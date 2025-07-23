import mudata as md
# Numpy, Pandas
import pandas as pd
import numpy as np


def fetch_meta_varnames_targeted(obj, use_mudata_obs, mod):
    """
    fetch_meta_varnames_targeted
    """
    if use_mudata_obs == True:
        print("Pulling from main obs slot of mudata object")
    
        keys = obj.obs_keys()
        
        return keys
    
    else:
        if mod == None:
            raise ValueError(
                "mod must be defined if use_mudata_obs is False."
                )
    
        print("Pulling from obs slot of", mod, sep = " ")
        
        keys = obj[mod].obs_keys()
        
        return keys


def fetch_meta_varnames_mudata(obj, mod):
    """
    fetch_meta_varnames_mudata
    
    modality (str)
    """
    # Directs fetch_meta_varnames to pull from the obsm slot of the main mudata object
    # if True
    use_mudata_obs = False

    # Search for modality in the mudata object first
    # If not found, set to use the main object obsm slot
    if mod in obj.mod_names:
        use_mudata_obs = False
    elif mod == None:
        # Set use_mudata_obs to True to pull from the main object obsm slot
        use_mudata_obs = True
    else:
        # Set use_mudata_obs to True if modality is not found in the object
        warning(
            paste0(
                mod, " modality not found in mudata object. Pulling ",
                "metadata variable names from main mudata object obs slot."
                )
            )
        use_mudata_obs = True
        
    # Fetch reduction from the location determined above
    if use_mudata_obs == True:
        keys = fetch_meta_varnames_targeted(
            obj, 
            use_mudata_obs = True, 
            mod = None
            )
    else:
        keys = fetch_meta_varnames_targeted(
            obj, 
            use_mudata_obs = False, 
            mod = mod
            )
    
    return keys
