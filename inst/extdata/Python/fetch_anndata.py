# Libraries used
# Numpy, Pandas
import pandas as pd
import numpy as np

# Scipy
from scipy.sparse import csr_matrix, csc_matrix
# Anndata CSCDataset class (for disk-backed anndata objects)
from anndata.abc import CSRDataset, CSCDataset

# Base Python
from collections import Counter
import re
import warnings

def is_match(regex_search, target):
    """
    Returns True if a regex is found in a target string, and 
    False if not.
    
    Arguments
    ----------
    regex_string: a regex expression to search for, formatted 
    as a string.
    
    target (string): the target for which to search regex_string.
    """
    regex = re.compile(regex_search)
    
    if regex.search(target):
        return True
    else: 
        return False

def fetch_vars_from_matrix(matrix, matrix_vars, cells, var_names = None, obs_names = None):
    """
    Returns data for a set of vars from a single matrix (X, 
    obs, obsm, for example). 
    
    Arguments
    ----------
    matrix: a matrix from X, obs, obsm, or var. 
    
    matrix_vars: a list of vars that are present in the matrix passed 
        to `matrix`. All vars passed must be present in the matrix, and 
        must be passed as they appear in the matrix (i.e. without a 
        "key"). **The parent functions calling this function should 
        verify that vars are present in the matrix before running this
        function**.
    
    cells: the cells in `matrix` from which to pull data. To avoid 
        unintended behavior, **the parent functions calling this 
        function should verify that the cells passed here correspond 
        to the cells in the matrix**
    
    var_names: the names of **valid vars present in the matrix**
        which are used to determine which indices to pull data from in
        matrices that do not have column names. This is required for 
        all matrix types that are not pd.DataFrames. In most cases, 
        this will be set to the var_names attribute of the anndata 
        object from which the matrix is being pulled (or for MuData 
        objects, the anndata object corresponding to the modality from 
        which the matrix is pulled). **For obsm matrices corresponding 
        to reductions,** this should be a pandas Index object with 
        string values from **1** to the number of columns in the 
        matrix. A one-based index is used to keep syntax consistent 
        with fetch_data usage in R objects. For example, users will 
        enter "UMAP_1" to pull the first dimension of the UMAP, and 
        fetch_vars_from_matrix will look for "1" in the Pandas Index. 
        This should be index 0 in the Pandas Index, and the numpy 
        matrix will be indexed for column zero).
    
    obs_names: the names of **valid cells in the matrix**, used to 
        determine which indices to pull data from in matrices without 
        row names. This is required for all matrix types that are not 
        pd.DataFrames. In all cases, this will be set to the obs_names 
        attribute of the anndata object from which the matrix is being 
        pulled (or for MuData objects, the anndata object corresponding 
        to the modality from which the matrix is pulled).
    """
    if (
        isinstance(matrix, csr_matrix) | 
        isinstance(matrix, csc_matrix)
        ):
        # Sparse matrices (in-memory)
        # Cells and vars are accessed by index. The indices of cells and 
        # vars are fetched using obs_names and var_names, respectively.
        
        # obs_names and var_names are required for this matrix type. 
        # Throw an error if not present
        if var_names is None:
            raise valueError(
                "If `matrix` is a csc_matrix or a csr_matrix, " + 
                "`var_names` must be defined."
            )
            
        if obs_names is None:
            raise valueError(
                "If `matrix` is a csc_matrix or a csr_matrix, " + 
                "`obs_names` must be defined."
            )
        
        # To avoid an error from get_loc(), filter keyless vars for 
        # only those that are in var_names
        vars_in_X = [var for var in matrix_vars if var in var_names]
        
        vars_idx = [var_names.get_loc(var) for var in vars_in_X]
        cells_idx = [obs_names.get_loc(cell_id) for cell_id in cells]
    
        # Subset matrix for vars and cells 
        if (isinstance(matrix, csc_matrix)):
            # CSC matrix: subsetting by *columns* is computationally cheap,
            # while subsetting by rows is expensive
            # Pull vars, then cells (expensive step of pulling cells is less
            # so when vars have been pulled first)
            matrix_subset = matrix[:,vars_idx][cells_idx,:]
        elif (isinstance(matrix, csr_matrix)):
            # CSC matrix: subsetting by *rows* is computationally cheap,
            # while subsetting by columns is expensive
            # Pull cells, then vars
            matrix_subset = matrix[cells_idx,:][:,vars_idx]
        
        data = pd.DataFrame(
            # Convert subset to numpy array
            # Ensures columns are densified and therefore can be 
            # converted to equivalent R data structures
            matrix_subset.toarray(),
            # Add cell ids and var names back 
            # These should be in the same order as the indices 
            index = cells,
            columns = matrix_vars
            )
    elif (
        isinstance(matrix, CSCDataset) | 
        isinstance(matrix, CSRDataset)
    ):
        # Sparse matrices (disk-backed)
        # Cells and vars are accessed by index. The indices of cells and 
        # vars are fetched using obs_names and var_names, respectively.
        
        # obs_names and var_names are required for this matrix type. 
        # Throw an error if not present
        if var_names is None:
            raise valueError(
                "If `matrix` is a CSCDataset or a CSRDataset, " + 
                "`var_names` must be defined."
            )
            
        if obs_names is None:
            raise valueError(
                "If `matrix` is a CSCDataset or a CSRDataset, " + 
                "`obs_names` must be defined."
            )
        
        # To avoid an error from get_loc(), filter keyless vars for 
        # only those that are in var_names
        vars_in_X = [var for var in matrix_vars if var in var_names]
        
        vars_idx = [var_names.get_loc(var) for var in vars_in_X]
        cells_idx = [obs_names.get_loc(cell_id) for cell_id in cells]
        
        # Subset matrix for vars and cells 
        if (isinstance(matrix, CSCDataset)):
            # CSCDataset: subsetting by *columns* is computationally cheap,
            # while subsetting by rows is expensive
            # Pull vars, then cells (expensive step of pulling cells is less
            # so when vars have been pulled first)
            matrix_subset = matrix[:,vars_idx][cells_idx,:]
        elif (isinstance(matrix, CSRDataset)):
            # CSCDataset: subsetting by *rows* is computationally cheap,
            # while subsetting by columns is expensive
            # Pull cells, then vars
            matrix_subset = matrix[cells_idx,:][:,vars_idx]
        
        data = pd.DataFrame(
            # Convert subset to numpy array
            # Ensures columns are densified and therefore can be 
            # converted to equivalent R data structures
            matrix_subset.toarray(),
            # Add cell ids and var names back 
            # These should be in the same order as the indices 
            index = cells,
            columns = matrix_vars
            )
    elif (isinstance(matrix, pd.DataFrame)):
        # Pandas dataframe: pull data via .loc 
        data = matrix.loc[cells, matrix_vars]
    elif (isinstance(matrix, np.ndarray)):
        # Numpy arrays
        
        # obs_names and var_names are required for this matrix type. 
        # Throw an error if not present
        if var_names is None:
            raise valueError(
                "If `matrix` is a numpy ndarray, " + 
                "`var_names` must be defined."
            )
        
        if obs_names is None:
            raise valueError(
                "If `matrix` is a numpy ndarray, " + 
                "`obs_names` must be defined."
            )
        
        # Convert vars to pull to index
        try: 
          vars_idx = [var_names.get_loc(var) for var in matrix_vars]
        except AttributeError:
          # var_names is defined in the parent function. For reductions, 
          # this is a manually constructed Pandas Index object. If a list is
          # passed instead, get_loc will fail. A more informative error is 
          # returned in this case.
          raise TypeError(
            "input to argument `var_names` is not a Pandas Index object. " +
            "Use `pd.Index()` to convert the input to the correct type."
            )
        
        # Convert cells to pull to index
        cells_idx = [obs_names.get_loc(cell_id) for cell_id in cells]
        
        # Slice numpy array based on the indices above
        # numpy.ix: provides most efficient means of slicing a numpy array
        matrix_subset = matrix[np.ix_(cells_idx, vars_idx)]
        
        # Convert subset to a pandas dataframe
        data = pd.DataFrame(
            matrix_subset,
            # Add cell IDs to dataframe
            index = cells,
            # Add the string representation of the vars pulled as columns
            columns = matrix_vars
            )
    else:
        """
        Strings wrapped according to PEP-8 style guidelines
        https://stackoverflow.com/questions/1874592/
        how-to-write-very-long-string-that-conforms-with-
        pep8-and-prevent-e501
        """
        raise NotImplementedError(
            ("fetch_vars_from_matrix does not know how to handle " + 
            "a matrix of class {0}").format(type(matrix))
            )
            
    return data

# This is used by fetch_keyed_vars_mudata so it must be in this script
def fetch_metadata_mudata(obj, meta_vars):
    """
    Fetches metadata vars in MuData objects.
    
    When axis = 0 (overlap in cells but not in vars):
    If variables are present in the obs matrix of more than one 
    modality, they will be considered "ambiguous" since it can't be 
    verified that they describe the same observation. For example, a 
    QC metric with the same name for RNA and ATAC-seq likely describes 
    something different, and it's not clear which should be returned.
    In this case, users will need to enter the name of the modality 
    with a ":" (for example, "RNA:n_counts").
    
    Variables not found will not throw a warning (this is 
    expected to be caught in parent functions, for example 
    fetch_mudata). "Ambiguous" entries where a variable is present 
    in the obs matrices of multiple modalities will be thrown by this 
    function.
    
    Arguments
    ----------
    obj: a MuData object
    
    meta_vars: a list of vars present in any of the obs matrices in the 
        modalities, or the obs matrix of the 
        MuData object.
    """
    # Initialize dataframe used to pull obs variables 
    obs_data = pd.DataFrame(index=obj.obs.index)
    
    # If any variables entered have a modality prefix with an underscore, 
    # convert to a MuData-style modality prefix (a colon)
    # This allows users to enter either mod1_qc or mod1:qc, and both 
    # will be properly interpereted
    # Downstream code interperets a MuData style prefix
    
    # Variables with an underscore that match existing variable names in 
    # any obs matrix should not be transformed (i.e. a variable in the rna
    # modality named rna_counts should not be transformed to rna:counts)
    # Form a set of existing obs vars in constituent anndata objects, then take 
    # the union with the obs keys in the main object
    obs_vars = (
        {var
         for adata in obj.mod.values()
         for var in adata.obs.keys()}
        | set(obj.obs.keys())
    )
    
    # Form regex to transform feature names
    modalities = obj.mod_names
    # Creates a pattern that matches any of the modalities in the object
    # For example, .join contents evaluate to ^(mod1|mod2|mod3)_(.+)$
    pattern = re.compile(rf"^({'|'.join(map(re.escape, modalities))})_(.+)$")
    
    # Transform names matching the regex if they do not appear in 
    # the set of existing vars
    meta_vars = [
        var if var in obs_vars
        else pattern.sub(r"\1:\2", var)
        for var in meta_vars
    ]
    
    # Pull metadata vars into MuData obs table
    obj.pull_obs(columns = meta_vars)
    
    # Suggested code for aggregation of metadata accross modalities
    if obj.axis == 0:
        for var in meta_vars:
            # When axis = 0, columns pulled will always contain the modality 
            # prefix with a ":" (i.e. "mod1:qc")
            # In the case a user enters a variable with a modality prefix, 
            # it will exactly match columns pulled
            # This will also match cases where a varible is only in the 
            # obs matrix of the MuData object
            if var in obj.obs.columns:
                # In this case, pull the variable that exactly matches the query
                obs_data[var] = obj.obs[var]
            else:
                # Vars entered without a modality prefix
                # Extract columns returned via regex
                # Columns will either consist of just the variable name, or the 
                # variable name with a modality prefix ending in ":"
                # for example, "clusters" or "RNA:clusters" 
                pattern = re.compile(rf'^(?:.*:)?{re.escape(var)}$')
                matches = [col for col in obj.obs.columns if pattern.match(col)]
                
                if len(matches) == 0:
                    # No matches: the variable does not exist
                    # This could be from entering a modality key 
                    # that does not match a variable 
                    
                    # In this case, continue without returning anything
                    # (this is caught outside the function, for example in 
                    # fetch_data)
                    pass
                elif len(matches) == 1:
                    # If one entry is found, the metadata variable appears 
                    # either in just one modality (when axis = 0, obs_names 
                    # are shared)
                    # Copy the matching column to the metadata table output
                    obs_data[var] = obj.obs[matches[0]]
                else:
                    # Multiple entires found: variable is in more than one 
                    # modality, and may not be describing the same property 
                    # (i.e. qc metrics on RNA and qc metrics on ATAC-seq may 
                    # be both named "qc", but don't describe the same property)
                    # Warn the user in this case and require the specification 
                    # of a MuData-style prefix (example: "mod1:qc")
                    r.warning(
                        'The metadata variable "' + var + 
                        '" is present in multiple modalities. This object ' +
                        'has shared obs_names, and variables with the same ' +
                        'name in different modalities do not necesarially ' +
                        'describe the same property. To pull this variable, ' +
                        'please add a prefix to the variable name with a ' +
                        '":", for example "' +
                        matches[0] + '".',
                        **{'call.': False}
                        )
    elif obj.axis == 1:
        raise NotImplementedError(
            "Code to fetch metadata is not yet implemented for MuData " +
            "objects with axis of 1 (shared var_names, concatenated " +
            "obs_names)."
            )
    elif obj.axis == -1:
        raise NotImplementedError(
            "Code to fetch metadata is not yet implemented for MuData " +
            "objects with axis of -1 (shared var_names and " +
            "obs_names)."
            )
    
    # Return constructed table of matching metadata vars
    return obs_data

def fetch_keyed_vars(obj, target_vars, cells, layer):
    """
    Returns expression data for all variables in target_vars that are
    found in the object at obj. Expression data will be returned only 
    for the variables that have an location key associated with them
    (i.e. X_FIS1 for a gene, obs_sample_id for a metadata variable, or
    (obsm_key)_var for a variable in an obsm matrix).
    
    Arguments
    ----------
    obj: an anndata object.
    
    target_vars: a set of variables to search the object for. Variables
    without a location key may be entered, but data will not be 
    returned for these variables.
    
    cells: Optional, a list of cells to return data for. If left as 
    None, data will be returned for all cells in the object. This is 
    passed down from the parent fetch_anndata function.
    
    layer: Optional, a string specifying the layer to return data from, 
    provided the variable in question is from the X matrix. The layer
    argument is ignored for variables from all other locations. This is 
    passed down from the parent fetch_anndata function.
    
    Return
    ----------
    A pandas DataFrame with expression data for each of the requested 
    variables in target_vars.
    """
    # Dataframe will store expression data for all variables found
    data_return = pd.DataFrame()
    
    # Get a list of all "keys" in the object 
    # Currently supported locations: X, obs, and matrices within obsm
    key_names = ["X", "obs"] + list(obj.obsm_keys())
    
    # Construct a list of keys with the indices of vars that match each key
    keyed_var_locations = []
  
    for key in key_names:
        # Conditional that is set to True in 2.2. if there are issues accessing
        # any variables in a given matrix
        key_error = False
        
        ## 2.1. Search for keyed vars for the current location in target_vars ####
        # Create regex object from the current key for searching
        # Capture group will return the text in the var after the key and "_" 
        # (the name of the var in the matrix)
        key_regex = re.compile("^" + key + "_(.*)")
        
        # matches: a dictionary mapping the original keyed var to the "keyless" 
        # var produced from removing the key plus an underscore. The keyless var 
        # is how the var should appear in the matrix it is being pulled from.
        
        # Dictionary comprehension is used to return variables for which the 
        # regex search is "truthy". When a match is found, a a regex match 
        # object will be returned, which will evaluate as true in the if 
        # statement below.
        matches = {var:key_regex.search(var).groups()[0] 
            for var in target_vars 
            if key_regex.search(var)
            }
        
        # For the X matrix, the user may use a key that conflicts with the names 
        # of common obsm matrices, such as "X_umap". To prevent "X_umap_1" from 
        # causing the function to search the gene matrix for "umap_1", keyless 
        # vars that contain an object key before an underscore will be removed.
        
        # This was found to also happen in obsm_key matrices. For example, if
        # an object has two reductions named mofa and mofa_umap, and vars 
        # contains "mofa_umap_1", this function will attempt to search mofa for
        # "umap_1" (unintended), and it will then search mofa_umap for 
        # "mofa_umap_1" (intended). Instead of expanding the solution below to
        # other matrices, which could cause unintended issues with variable not
        # being searched for at all, I decided to use error handling when 
        # searching for a feature, so no errors ocurr in the unintended search
        # example above, and the function will keep iterating through obsm keys
        # to execute the intended example.
        
        # NOTE: if there happens to be a gene that contains an obsm matrix name 
        # plus an underscore, this gene will be unsearchable due to this 
        # workaround. This does not seem likely however, as there are no known 
        # human genes with underscores or reduction names in the gene symbol.
        if key == "X":
            # Define other keys in object 
            exclusion_keys = [key for key in key_names if key != "X"]
            
            # Return regex matches...
            matches = {key:value for (key, value) in matches.items() 
                # ... where none of the other keys in the object are found in 
                # the original var entry (before the key is removed)
                if not any([is_match(regex_search = exclusion_key, target = key) 
                    for exclusion_key in exclusion_keys])
                }
            
        # Extract values from curated matches dictionary (keyless vars to search 
        # X/obs/obsm matrix for)
        keyless_vars = list(matches.values())
        
        ## 2.2. Fetch Data if vars exist in the current location ####
        if (len(keyless_vars) > 0):
            
            ### 2.2.1. Pull expression matrix for the current key location ####
            # Also, subset keyless vars for those that are in the matrix. It 
            # is possible to enter inputs starting with a modality key that are 
            # not in the matrix, which will cause uninformative error messages
            if key == "X":
                # The X matrix alone supports "layer" (via layers)
                if layer == None:
                    matrix = obj.X
                else:
                    matrix = obj.layers[layer]
                    
                # Subset for vars in matrix
                keyless_vars = [
                    var for var in keyless_vars if var in obj.var_names
                    ]
                    
                # Form input to var_names argument
                # For X, the var_names argument is equal to the var_names
                # attribute of the anndata object
                var_names = obj.var_names
            elif key == "obs":
                # Metadata (obs)
                matrix = obj.obs
                
                # Subset for vars in matrix
                # Expectation: if vars are in obs_keys(), they should be in obs.
                # There should not be obs vars that are not in obs_keys()
                keyless_vars = [
                    var for var in keyless_vars if var in obj.obs_keys()
                    ]
                        
                # Form input to var_names parameter
                # Obs matrices are always pd.DataFrames, so var_names is not
                # required as an input to fetch_vars_from_matrix
                var_names = None
                        
            elif key in obj.obsm_keys():
                # For a key in the list of obsm keys, pull matrix for that entry
                matrix = obj.obsm[key]
                
                # Subset for vars in matrix (depending on matrix types)
                if (isinstance(matrix, pd.DataFrame)):
                    # Pandas dataframes: use .columns method to pull valid 
                    # entries
                    keyless_vars = [var for var in keyless_vars if var in matrix.columns]
                    
                    # Form input to var_names argument
                    # var_names is not required as an input to 
                    # fetch_vars_from_matrix when matrix is a pd.DataFrame
                    var_names = None
                elif (
                    isinstance(matrix, np.ndarray) |
                    isinstance(matrix, csr_matrix) | 
                    isinstance(matrix, csc_matrix) |
                    isinstance(matrix, CSCDataset)
                ):
                    # For numpy arrays (most likely, used for reductions), 
                    # and sparse matrices (less likely but possible)
                    # No column names, valid entries are based on index 
                    
                    # Form input to var_names argument 
                    # var_names is constructed for obs matrices, using the 
                    # same type as var_names in anndata objects (Pandas Index)
                    # var_names consists of the valid vars, which are the 
                    # dimensions of the matrix using a 1-based index (for 
                    # consistency with R SCUBA methods, and based on the 
                    # assumption that the user is requesting vars from R)
                    var_names = pd.Index(
                        [str(x) for x in range(1, matrix.shape[1] + 1)]
                        )
                    
                    # Subset vars to valid dimensions as defined above
                    keyless_vars = [
                        var for var in keyless_vars if var in var_names
                        ]
                    
            ### 2.2.2. Pull data for cells, keyless vars from matrix ####
            # Use fetch_vars_from_matrix from above
            try:
                data = fetch_vars_from_matrix(
                    matrix = matrix, 
                    matrix_vars = keyless_vars, 
                    cells = cells, 
                    var_names = var_names, 
                    obs_names = obj.obs_names
                    )
            except NotImplementedError as e:
                # There is a NotImplementedError raised from 
                # fetch_vars_from_matrix that appears when an unsupported 
                # matrix type is entered. However, the context of which matrix 
                # triggered the error (what key) is not available to that 
                # function. The error is caught here and reported
                # with that context
                raise NotImplementedError(
                    ("The SCUBA Python function fetch_anndata does not know "
                    "how to handle a matrix of class {0}. "
                    "This ocurred with the matrix at key '{1}'."
                    ).format(type(matrix), key)
                    )
            
            ### 2.2.3. Concatenate with data already fetched
            # Add location key back in to variable names and concatenate data
            # to data already fetched
            if data.shape[1] > 0:
                # Lambda function prepends the key and "_" to each element
                data = data.rename(
                    lambda x: key + "_" + x, 
                    axis = "columns"
                    )
                    
                # Concatenate data frame with data already fetched
                data_return = pd.concat(
                    [data_return, data],
                    axis = 1
                    )
                    
    # When finished with iteration, return the dataframe
    return data_return

def fetch_keyed_vars_mudata(obj, target_vars, cells, layer):
    """
    Returns expression data for all variables in target_vars that are
    found in the MuData object passed to obj. Expression data will be 
    returned only for the variables that have an location key 
    associated with them (i.e. (modality_key)_FIS1 for a feature, 
    obs_sample_id for a metadata variable, or (obsm_key)_var for a 
    variable in an obsm matrix).
    
    Arguments
    ----------
    obj: a MuData object.
    
    target_vars: a set of variables to search the object for. Variables
    without a location key may be entered, but data will not be 
    returned for these variables.
    
    cells: Optional, a list of cells to return data for. If left as 
    None, data will be returned for all cells in the object. This is 
    passed down from the parent fetch_anndata function.
    
    layer: Optional, a string specifying the layer to return data from, 
    provided the variable in question is from the X matrix. The layer
    argument is ignored for variables from all other locations. This is 
    passed down from the parent fetch_anndata function.
    
    Return
    ----------
    A pandas DataFrame with expression data for each of the requested 
    variables in target_vars.
    """
    # Initialize dataframe to store expression data for all variables found
    data_return = pd.DataFrame()
    
    # 1. Define keys
    # Modality keys, plus "obs", plus obsm keys from each modality
    
    # Compile obsm keys
    # Construct a set of unique keys by looping through the obsm keys of each 
    # constituent anndata object
    obsm_keys = list({key 
        for adata in obj.mod.values() 
        for key in adata.obsm.keys()})
    
    key_names = list(obj.mod_names) + ["obs"] + obsm_keys
    
    # 2. For each key, identify vars and fetch data if found
    # Before iterating through keys, initialize a list of vars found 
    keyed_var_locations = []
    
    for key in key_names:
        # Conditional that is set to True in 2.2. if there are issues accessing
        # any variables in a given matrix
        # See if this is still needed for mudata
        # key_error = False
        
        ## 2.1. Search for keyed vars for the current location in target_vars ####
        # Create regex object from the current key for searching
        # Capture group will return the text in the var after the key and "_" 
        # (the name of the var in the matrix)
        key_regex = re.compile("^" + key + "_(.*)")
        
        # matches: a dictionary mapping the original keyed var to the "keyless" 
        # var produced from removing the key plus an underscore. The keyless var 
        # is how the var should appear in the matrix it is being pulled from.
        
        # Dictionary comprehension is used to return variables for which the 
        # regex search is "truthy". When a match is found, a a regex match 
        # object will be returned, which will evaluate as true in the if 
        # statement below.
        matches = {var:key_regex.search(var).groups()[0] 
            for var in target_vars 
            if key_regex.search(var)
            }
        
        # Extract values from curated matches dictionary (keyless vars to search
        # matrix for)
        keyless_vars = list(matches.values())
        
        ## 2.2. Fetch Data if vars exist in the current location ####
        if (len(keyless_vars) > 0):
            # Means of fetching data based on key
            # Feature data
            if key in list(obj.mod_names):
                # Define matrix associated with key
                if layer == None:
                    matrix = obj[key].X
                else:
                    # Check if layer is in object before pulling
                    if layer in obj.layers.keys:
                        obj[key].layers[layer]
                    else:
                        r.warning(
                            "Requested layer " +
                            layer +
                            " not found in modality " +
                            key +
                            "."
                            )
            
                # Before fetching data, subset for keyless vars 
                # that are in the var names of the anndata object 
                # representing the modality
                keyless_vars = [
                        var for var in keyless_vars if var in obj[key].var_names
                        ]
                        
                # Form input to var_names argument of fetch_matrix,
                # to properly index matrices for features
                # Uses the var names for the modality corresponding to the key
                var_names = obj[key].var_names
                
                ### Pull data for cells, keyless vars from matrix ####
                # Use fetch_vars_from_matrix from above
                try:
                    print("Key = {}".format(key))
                    breakpoint()
                    data = fetch_vars_from_matrix(
                        matrix = matrix, 
                        matrix_vars = keyless_vars, 
                        cells = cells, 
                        var_names = var_names, 
                        obs_names = obj.obs_names
                        )
                except NotImplementedError as e:
                    # There is a NotImplementedError raised from 
                    # fetch_vars_from_matrix that appears when an unsupported 
                    # matrix type is entered. However, the context of which matrix 
                    # triggered the error (what key) is not available to that 
                    # function. The error is caught here and reported
                    # with that context
                    raise NotImplementedError(
                        ("The SCUBA Python function fetch_anndata does not know "
                        "how to handle a matrix of class {0}. "
                        "This ocurred with the matrix at key '{1}'."
                        ).format(type(matrix), key)
                        )
                
            elif key == "obs":
                # Metadata
                # Pass variables to the fetch_metadata_mudata function
                data = fetch_metadata_mudata(
                    obj = obj, 
                    meta_vars = keyless_vars
                    )
            elif key in obsm_keys:
                pass
            
            # Concatenate with data already fetched
            # Add location key back in to variable names and concatenate data
            # to data already fetched
            if data.shape[1] > 0:
                # Lambda function prepends the key and "_" to each element
                data = data.rename(
                    lambda x: key + "_" + x, 
                    axis = "columns"
                    )
                    
                # Concatenate data frame with data already fetched
                data_return = pd.concat(
                    [data_return, data],
                    axis = 1
                    )
                    
    # When finished with iteration through keys, return the dataframe
    return data_return

def remove_key(data, key, vars_modify=None):
    """
    Removes the key prefix from all var entries with the defined
    key. For example, if the key is "obs", all entries with a metadata
    variable will be changed from "obs_<variable_name>" to
    "<variable_name>".
    
    Arguments
    ----------
    data: expression data/metadata returned by fetch_keyed_vars().
    
    key: a string giving a key to remove. The key should not contain
    an underscore. For example, use "obs" instead of "obs_", and "X"
    instead of "X_".
    
    vars_modify (optional): a list giving a subset of vars to remove 
    the key from. For example, if the key is "X" and vars_modify is
    ["X_gene1", "X_gene2"], only "gene1" and "gene2" will have the "X"
    key removed, and all other vars with the "X" key will be unchanged.
    
    Returns
    ----------
    A pandas dataframe with modified column names. The data itself will
    be unchanged.
    """
    # Identify variables with obs key
    key_regex = re.compile("^" + key + "_(.*)")
    
    if vars_modify==None:
        target_vars = data.columns.values
    else:
        target_vars = vars_modify
    
    # Generate names of obs variables without the obs key
    # Also maps var names with obs keys to var names without the key
    remove_key = {var:key_regex.search(var).groups()[0] 
                for var in target_vars
                if key_regex.search(var)
                }
    
    data.rename(remove_key, axis=1, inplace=True)
    
    return data

def fetch_anndata(obj, fetch_vars, cells=None, layer=None):
    """
    Fetches the specified variables from an anndata object.
    
    Arguments
    ----------
    obj: an anndata object.
    
    fetch_vars: a list giving the variables to fetch from the object.
    Variables entered should ideally have a key prepended with the 
    location of the variable in the object, for example, to pull the
    FIS1 gene from the genes matrix (X), specify "X_FIS1" instead of
    "FIS1". To pull metadata, use "obs_", and to pull data from a 
    matrix in obsm, use the name of that matrix, not obsm. For 
    example, if a matrix in obsm is named "protein", use "protein_" 
    to pull data from that matrix. Variables that are entered without
    a key can still be found, as long as there is only one matrix in 
    the object with that variable name. Variables that do not have a 
    valid key (X_, obs_, and a key from obj.obsm_names()) will be 
    ignored, as will duplicate variables.
    
    cells: Optional, a list of cells to return data for. If left as 
    None, data will be returned for all cells in the object.
    
    layer: Optional, a string specifying the layer to return data from, 
    provided the variable in question is from the X matrix. The layer 
    argument is ignored for variables from all other locations, including 
    alternate modalities.
    """
    # If target_vars was passed as a one-element vector from R, it will be
    # a string now. This must be converted to a list to avoid issues during
    # iteration. Multi-element character vectors are properly converted to a 
    # list.
    if isinstance(fetch_vars, str):
        fetch_vars = [fetch_vars]
    
    # 1. Set default values
    # Layer (assay): if null, data will be pulled from object$X

    # Cells: if None (NULL), use all cells in the object
    if cells == None:
        cells = list(obj.obs_names)
    
    # 2. Check fetch_vars for duplicate entries
    # First ensure fetch_vars is a list (otherwise fetch_anndata will attempt
    # to pull data for each letter of a feature entered, instead of the 
    # feature itself)
    # if not isinstance("X_CATSPER4", list):
    #     raise TypeError("fetch_vars must be a list.")
    
    # "fetch_vars" is used instead of "vars" since vars is a base 
    # python function
    
    # Use collections.Counter to find duplicates
    # Must iterate through the keys of the dictionary produced by counter
    # If iterating through fetch_vars, duplicates will display twice (or as
    # many times as they appear in fetch_vars)
    duplicate_vars = [var
        for var in Counter(fetch_vars).keys()
        if Counter(fetch_vars)[var] > 1
        ]
    
    # Report duplicates to the user and notify that only one entry 
    # will be included 
    if len(duplicate_vars) > 0:
        r.warning(
            ("Duplicate entries passed to vars: " +
            ", ".join(duplicate_vars) +
            ". Only one entry for each variable will be returned."),
            # The call context for the warning prints as 
            # do.call(fn, c(args, named_args)) 
            # which is not useful for the end user. The call is suppressed by
            # setting call. to FALSE. 
            # The . in call. causes an error with Python syntax, so it is 
            # specified using dictionary unpacking. "call." is passed as a 
            # string in a dictionary, and tells warning() to accept a call. 
            # argument and set it to False (FALSE in R)
            **{'call.': False}
            )

    # Remove duplicates by creating a dictionary from the list, and converting
    # back to a list
    # https://stackoverflow.com/questions/480214/how-do-i-remove-duplicates-from-a-list-while-preserving-order
    # Dictionaries automatically remove duplicate keys, while preserving the 
    # order keys appear (sets do not)
    fetch_vars = list(dict.fromkeys(fetch_vars))
    
    # 3. Pull data for keyed features in each key location
    # (uses fetch_keyed_vars function)
    fetched_data = fetch_keyed_vars(
        obj = obj, 
        target_vars = fetch_vars,
        cells = cells, 
        layer = layer
        )
    
    # 4. Identify location of remaining vars
    # Feature data can appear in object.obs (metadata), so ambiguous vars 
    # should be checked to avoid returning incorrect data in the case of data 
    # being found in multiple assays
    
    # Identify remaining variables via list comprehension
    remaining_vars = [
        var_i 
        for var_i in fetch_vars 
        if var_i not in list(fetched_data.columns)
        ]
    
    # Determine the location of each remaining variable
    remaining_locations = {}
    
    # Names of variables in X are likely very large, so the variables that 
    # are members of X should be determined first
    X_vars = [
        remaining_var 
        for remaining_var in remaining_vars 
        if remaining_var in obj.var_names
        ]
    
    # Do the same for metadata
    metadata_vars = [
        remaining_var
        for remaining_var in remaining_vars 
        if remaining_var in obj.obs_keys()
        ]
    
    for remaining_var in remaining_vars:
        # For each variable, create a dictionary entry with an empty list
        # Append locations to the list for each variable when it is found
        remaining_locations[remaining_var] = []
        
        # Check X vars computed above
        if remaining_var in X_vars:
            remaining_locations[remaining_var].append("X")
        
        # Check metadata vars
        if remaining_var in metadata_vars:
            remaining_locations[remaining_var].append("obs")
        
        # Look in obsm locations
        for obsm_key in obj.obsm_keys():
            # Must test pandas DataFrames separately to avoid 
            # errors with .columns 
            if isinstance(obj.obsm[obsm_key], pd.DataFrame):
                # Look for var in the columns of the dataframe 
                if remaining_var in obj.obsm[obsm_key].columns:
                    remaining_locations[remaining_var].append(obsm_key)
            else:
                # Ignore matrices that are not pandas dataframes for now
                # (reductions should not make it to this step since a key is 
                # required to distinguish reduction dimensions (X_umap_1 vs. 
                # X_tsne_1), "1" would not be specific enough). Matrices with 
                # searchable variable names should not be numpy arrays
                pass
    
    # 5. Warn the user if remaining vars are found in multiple locations
    # Need a warning that does not interrupt the function (warnings behave 
    # like errors when called with raise)
    
    # Identify remaining vars in multiple locations
    ambiguous_vars = {key:value 
                          for (key, value) in remaining_locations.items() 
                          if len(value) > 1
                          }
                        
    if len(ambiguous_vars) > 0:
        # Display an example of how to remove ambiguous vars to the user
        # Store an ambiguous var, and the possible keys that can be added
        example_var = list(ambiguous_vars.keys())[0]
        example_keys = ambiguous_vars[example_var]
        
        # Generate a list of the features to be en
        example_entries = [key + "_" + example_var for key in example_keys]
        
        # Construct string from examples
        if len(example_entries) == 2:
            # If there are two possible entries, add "or" between the examples 
            example_str = ", or ".join(example_entries)
        elif len(example_entries) > 2:
            # If there are more than two possibilities, add "or" before the 
            # last entry and join all other entries with a comma.
            example_entries[-1] = "or " + example_entries[-1]
            example_str = ", ".join(example_entries)
        
        # Display warning message
        # Warn in R (better user experience, and can be detected by testthat)
        r.warning(
            ("The following variables were found in multiple locations: " +
            ", ".join(list(ambiguous_vars.keys())) +
            ". These variables will not be retrieved due to ambiguity. " +
            "To pull data for these variables, please specify which " +
            "location in the object to pull the variable from using an " + 
            "underscore (for example: " +
            example_str +
            ")."
            ),
            # Remove call context (not useful, see comment above)
            **{'call.': False}
            )
            
    # 6. Fetch data for remaining vars ####
    # Keys will be added, and then data will be fetched for the new keyed vars
    
    ## 6.1. Add keys (for remaining vars where a key was identified) ####
    # Identify vars to add key to (vars with exactly one location identified)
    ambiguous_vars_in_one_location = {key:value 
        for (key, value) in remaining_locations.items() 
        if len(value) == 1
        }
    
    # Create dictionary mapping the vars above to their keyed equivalents
    # Add key with underscore to var name 
    map_keyed_vars = {key:value[0] + "_" + key 
        for (key, value) in ambiguous_vars_in_one_location.items()}
    
    # Identify ambiguous obsm vars (these will result in a warning while 
    # ambiguous vars in X or obs will not)
    ambiguous_obsm_vars = {key:value 
        for (key, value) in ambiguous_vars_in_one_location.items() 
        # if value (checks that value exists, though it always should)
        if value and value[0] not in {"X", "obs"}
        }
    
    new_keyed_vars = list(map_keyed_vars.values())

    ## 6.2. Catch duplicate vars ####
    # Vars without keys may be the same as those entered with keys. For 
    # example, FIS1 and X_FIS1 should describe the same variable. If 
    # duplicates are added to the fetch_data dataframe, unexpected behavior
    # will occur later on when the variables are sorted, and duplicate entries
    # of the same variable do not add any value to the user. If variables 
    # with keys added are the same as the keyed variables already fetched, the 
    # variable should not be added again. 
    
    # Construct a dictionary with variables that match existing columns in 
    # fetched_data after the key is added
    duplicate_vars = {
        var_entered:keyed_var 
        for (var_entered, keyed_var) in map_keyed_vars.items() 
        if keyed_var in fetched_data.columns
        }

    if len(duplicate_vars) > 0:
        # Warn the user that a variable is duplicated and only one 
        # copy will be returned
        # For each key:value pair in duplicate vars, create a set of 
        # human-readable strings describing the feature entered and its 
        # equivalent duplicate
        duplicate_description = [key + ", equivalent to " + value 
            for (key, value) in duplicate_vars.items()]
        duplicate_description = "; ".join(duplicate_description)
        
        # Example of the keyed variable(s) to be returned
        will_be_returned = ", ".join(list(duplicate_vars.values()))
        
        # Warn in R (better user experience, and can be detected by testthat)
        r.warning(
            ("The following entries describe the same variable: " +
            duplicate_description + 
            ". Only the copy of the variable with the key entered will be " +
            "returned (" +
            will_be_returned +
            ")."),
            # Remove call context (not useful, see comment above)
            **{'call.': False}
            )
        
        # Remove the duplicate variable(s) before fetching data
        new_keyed_vars = [var 
            for var in new_keyed_vars 
            if var not in duplicate_vars.values()]
            
        # Also remove duplicate vars from the list of X_vars 
        # (X_vars is used to remove the modality key from features in X entered
        # without a modality key. If this is not done, an entry of "X_GAPDH" 
        # and "GAPDH" will result in a return of the column name "GAPDH" 
        # instead of "X_GAPDH").
        X_vars = [var for var in X_vars if var not in duplicate_vars.keys()]
        
        # Also remove duplicates from fetch_vars (this variable is used later on 
        # to sort the final dataframe, and duplicate entries in fetch_vars will
        # cause duplication in the dataframe).
        fetch_vars = [var 
            for var in fetch_vars
            if var not in duplicate_vars.values()]

    ## 6.3. Fetch data for new keyed variables and append to fetched_data ####
    new_data = fetch_keyed_vars(
        obj = obj, 
        target_vars = new_keyed_vars,
        cells = cells, 
        layer = layer
        )
        
    fetched_data = pd.concat(
        [fetched_data, new_data],
        axis = 1
        )

    # 7. Determine which vars, if any, were not returned
    vars_not_found = [var for var in fetch_vars if var not in fetched_data.columns]
    
    # Record the vars found in 6 (these will show up in vars_not_found incorrectly
    # since the vars in fetch_data contain the new key added)
    new_keyed_vars_found = [var 
        for var in map_keyed_vars.keys() 
        if map_keyed_vars[var] in fetched_data.columns
        ]
        
    # Remove any vars in new_keyed_vars_found from vars_not_found
    vars_not_found = [var 
        for var in vars_not_found 
        if var not in new_keyed_vars_found
        ]
    
    # 8. Warn for ambiguous entries with one match in obsm matrices
    for ambig_var, key in ambiguous_obsm_vars.items():
        if ambig_var not in vars_not_found:
            r.warning(
                (ambig_var +
                " was passed to fetch_data without specifying the key of the " +
                "matrix to pull the feature from, and it is not present in " +
                "X. The feature was found in " +
                "'obsm."+ "".join(key) + "'"+ 
                " and successfully returned. Future returns without a key will " +
                "fail if features are present in multiple obsm matrices."),
                # Remove call context (not useful, see comment above)
                **{'call.': False}
                )
            
    # 9. Warnings/errors for variables not found ####
    # ten_plus_message: added to message if there are more than 10 
    # missing variables
    if (len(vars_not_found) > 10):
        ten_plus_message = "(10 out of {} shown)".format(len(vars_not_found))
    else:
        ten_plus_message = ""
    
    # Show an error if all vars were not found, and a warning if at least 
    # one var was not found
    if len(vars_not_found) == len(fetch_vars):
        raise ValueError(
            "None of the requested variables were found " +
            ten_plus_message +
            ": " +
            ", ".join(vars_not_found)
            )
    elif len(vars_not_found) > 0:
        # Warn in R (better user experience, and can be detected by testthat)
        r.warning(
            "The following requested variables were not found " +
            ten_plus_message +
            ": " +
            ", ".join(vars_not_found),
            # Remove call context (not useful, see comment above)
            **{'call.': False}
            )
    
    # 10. Sort columns ####
    # Order of columns in data should reflect the order entered, not the 
    # order fetched
    # Columns in dataframe use keys for all vars, including those entered 
    # without a key. To sort properly, keys are added to var names using 
    # the mapping produced in 6.
    # Get method: if a mapped keyed variable exists for var, that variable is 
    # returned. If not, the var is returned as-is.
    var_order = [map_keyed_vars.get(var, var) for var in fetch_vars]
    # Remove variables that are not in the column names of the dataframe
    # to avoid errors
    var_order = [var for var in var_order if var in fetched_data.columns]
    
    # Pass list to fetched_data to order columns accordingly
    fetched_data = fetched_data[var_order]

    # 11. Remove keys from obs variables, X variables entered without a key
    # This is done for consistency with Seurat FetchData
    fetched_data = remove_key(
        data = fetched_data,
        key = "obs"
        )

    # Identify X vars entered without a key
    # Add key at the beginning to pass to vars_modify (vars must be 
    # identified as they currently appear in fetched_data)
    X_target_vars = ["X_" + var for var in X_vars]
    
    fetched_data = remove_key(
        data = fetched_data,
        key = "X",
        vars_modify = X_target_vars
        )
    
    return fetched_data

