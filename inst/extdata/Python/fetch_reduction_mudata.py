import mudata as md
# Numpy, Pandas
import pandas as pd
import numpy as np


def fetch_reduction_targeted(obj, reduction, cells, dims, use_mudata_obsm, mod):
    """
    fetch_reduction_targeted
    """
    if use_mudata_obsm == True:
        print("Pulling from main obsm slot of mudata object")
    
        df = pd.DataFrame(
            obj.obsm[reduction],
            index = obj.obs_names,
            # Construct var names using a 1-based index (for consistency 
            # with Seurat, R objects used in SCUBA)
            columns = [str(x) for x in range(1, obj.obsm[reduction].shape[1] + 1)]
            )
        
        df = df.loc[cells, dims]
        
        return df
    
    else:
        if mod == None:
            raise ValueError(
                "mod must be defined if use_mudata_obsm is False."
                )
    
        print("Pulling from obsm slot of", mod, sep = " ")
        
        df = pd.DataFrame(
            obj[mod].obsm[reduction],
            index = obj[mod].obs_names,
            # Construct var names using a 1-based index (for consistency 
            # with Seurat, R objects used in SCUBA)
            columns = [str(x) for x in range(1, obj[mod].obsm[reduction].shape[1] + 1)]
            )
        
        df = df.loc[cells, dims]
        
        return df


def fetch_reduction_mudata(obj, reduction, cells, dims):
    """
    fetch_reduction_mudata
    
    dims (str)
    """
    # Directs fetch_reduction to pull from the obsm slot of the main mudata object
    # if True
    use_mudata_obsm = False
    
    # Check reduction for invalid inputs ()
    if reduction in obj.mod_names:
        raise ValueError(
            "The reduction entered (" +
            reduction +
            ") is the name of a modality in this object."
            )
    
    # Behavior if modality is not defined
    # Search for reduction in the obsm slot of the mudata object first
    # If not found, search in individual anndata objects
    # If individual anndata objects have a reduction with the same name, use the one
    # in the default modality
    if reduction in obj.obsm_keys():
        # Set use_mudata_obsm to True to pull from the main object obsm slot
        use_mudata_obsm = True
    else:
        # Construct list of modalities for which the reduction is found, if any
        mod_reduction_matches = []
    
        for mod in obj.mod_names:
            if reduction in obj[mod].obsm_keys():
                mod_reduction_matches.append(mod)
        
        # Inspect list of matches
        # No matches: error, reduction not found
        if len(mod_reduction_matches) == 0:
            raise ValueError(
                "The reduction entered (" +
                reduction +
                ") was not found in any modalities or the main obsm slot."
                )
        elif len(mod_reduction_matches) == 1:
            # Tell downstream accession code to use the matching modality
            use_mod = mod_reduction_matches[0]
        
        # Ambiguous matches in more than one modality
        elif len(mod_reduction_matches) > 1:
            # If the reduction was found in the first modality in the object,
            # use that and warn user
            # If not, throw an error and require user to specify the modality
            default_modality = obj.mod_names[0]
        
        if default_modality in mod_reduction_matches:
            use_mod = default_modality
        else:
            raise ValueError(
                "The reduction entered (" +
                reduction +
                ") was found in multiple locations, but not in the default " +
                "modality (the first one in the object). Please specify the " +
                "modality to pull the reduction from."
                # Add a note about syntax to use once that is finalized
                # Mention the list of modalities to pull from?
                )
    
    # Fetch reduction from the location determined above
    if use_mudata_obsm == True:
        df = fetch_reduction_targeted(
            obj, 
            reduction = reduction, 
            cells = cells, 
            dims = dims, 
            use_mudata_obsm = True, 
            mod = None
            )
    
        return df
    else:
        print("Pulling reduction from ", use_mod, sep = "")
        df = fetch_reduction_targeted(
            obj, 
            reduction = reduction, 
            cells = cells, 
            dims = dims, 
            use_mudata_obsm = False, 
            mod = use_mod
            )
    
    return df
