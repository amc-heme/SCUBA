# Libraries used
# Numpy, Pandas
import pandas as pd
import numpy as np

# Scipy
from scipy.sparse import csr_matrix, csc_matrix

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

def fetch_keyed_vars(obj, target_vars, cells, slot):
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
    
    slot: Optional, a string specifying the layer to return data from, 
    provided the variable in question is from the X matrix. The slot
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
            if key == "X":
                # The X matrix alone supports "slot" (via layers)
                # TODO: pull from layers when slot != None
                matrix = obj.X
            elif key == "obs":
                # Metadata (obs)
                matrix = obj.obs
            elif key in obj.obsm_keys():
                # For a key in the list of obsm keys, pull matrix for that entry
                matrix = obj.obsm[key]
                    
            ### 2.2.2. Pull data for cells, keyless vars from matrix ####
            # Type checking via isinstance is used here due to the variety of
            # data types possible for the matrix
            # This is not considered "pythonic" and may need to be re-thought 
            # https://stackoverflow.com/questions/2225038/determine-the-type-of-an-object
            if (isinstance(matrix, csr_matrix) | isinstance(matrix, csc_matrix)):
                # Sparse matrix format (X): subset *object for genes, then 
                # construct a pandas dataframe (sparse matrices don't subset 
                # easily in python)
                data = pd.DataFrame.sparse.from_spmatrix(
                    obj[:, keyless_vars].X,
                    # csr_matrices don't have row, column names. 
                    # These are added here
                    index = obj.obs_names,
                    columns = keyless_vars
                    )
                    
                # Columns returned will be pandas sparse arrays.
                # Arrays must be densified to be accesssible downstream in R 
                # Densify with np.asarray
                for column in data:
                    data[column] = np.asarray(data[column])
            elif (isinstance(matrix, pd.DataFrame)):
                # Pandas dataframe: pull data via .loc 
                data = matrix.loc[cells, keyless_vars]
            elif (isinstance(matrix, np.ndarray)):
                # Numpy arrays
                
                # The process below was designed with the assumption that only 
                # reductions would have this format. It should be applicable 
                # to other arrays however since arrays are stored without column
                # names, and the user would likely reference variables in these
                # arrays using an index (i.e. "X_umap_1", "(matrix_name)_1")
                
                # Convert numpy array to pandas datafrane, then slice for vars
                matrix = pd.DataFrame(
                    matrix,
                    index = obj.obs_names,
                    # Construct var names using a 1-based index (for consistency 
                    # with Seurat, R objects used in SCUBA)
                    columns = [str(x) for x in range(1, matrix.shape[1] + 1)]
                    )
                    
                # Slice based on keyless vars (should be string format of the 
                # column desired, "1" for the first column)
                data = matrix.loc[cells, keyless_vars]
                
            else:
                """
                Strings wrapped according to PEP-8 style guidelines
                https://stackoverflow.com/questions/1874592/
                how-to-write-very-long-string-that-conforms-with-
                pep8-and-prevent-e501
                """
                raise NotImplementedError(
                    ("FetchData does not know how to handle a matrix of class {0}. "
                    "This ocurred with the matrix at key '{1}'."
                    ).format(type(matrix), key)
                    )
            
            # Add location key back in to variable names (if ncol is greater than 
            # zero)
            if data.shape[1] > 0:
                # Lambda function prepends the key and "_" to each element
                data = data.rename(lambda x: key + "_" + x, axis = "columns")
                
            # Concatenate data frame with data already fetched
            data_return = pd.concat(
                [data_return, data],
                axis = 1
                )
            
    # When finished with iteration, return the dataframe
    return data_return

def fetch_anndata(obj, fetch_vars, cells=None, slot=None):
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
    
    slot: Optional, a string specifying the layer to return data from, 
    provided the variable in question is from the X matrix. The slot
    argument is ignored for variables from all other locations.
    """
    print("Begin python script")
    # If target_vars was passed as a one-element vector from R, it will be
    # a string now. This must be converted to a list to avoid issues during
    # iteration. Multi-element character vectors are properly converted to a 
    # list.
    if isinstance(fetch_vars, str):
        fetch_vars = [fetch_vars]
    
    # 1. Set default values
    # Slot (assay): if null, data will be pulled from object$X

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
    warnings.warn(
        ("Duplicate entries passed to vars: " +
        ", ".join(duplicate_vars) +
        ". Only one entry for each variable will be returned.")
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
        slot = slot
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
        warnings.warn(
            ("The following variables were found in multiple locations: " +
            ", ".join(list(ambiguous_vars.keys())) +
            ". These variables will not be retrieved due to ambiguity. " +
            "To pull data for these variables, please specify which " +
            "location in the object to pull the variable from using an " + 
            "underscore (for example: " +
            example_str +
            ")."
            )
            )
            
    # 6. Fetch data for remaining vars ####
    # Keys will be added, and then data will be fetched for the new keyed vars
    
    ## 6.1. Add keys (for remaining vars where a key was identified) ####
    # Identify vars to add key to (vars with exactly one location identified)
    vars_add_key = {key:value 
        for (key, value) in remaining_locations.items() 
        if len(value) == 1
        }
    
    # Create dictionary mapping the vars above to their keyed equivalents
    # Add key with underscore to var name 
    map_keyed_vars = {key:value[0] + "_" + key 
        for (key, value) in vars_add_key.items()}
    
    
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
    
    warnings.warn(
        ("The following variables entered describe the same variable: " +
        duplicate_description + 
        ". Only the copy of the variable with the key entered will be " +
        "returned (" +
        will_be_returned +
        ").")
        )
    
    # Remove the duplicate variable(s) before fetching data
    new_keyed_vars = [var 
        for var in new_keyed_vars 
        if var not in duplicate_vars.values()]
        
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
        slot = slot
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
    
    # 8. Warnings/errors for variables not found ####
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
        warnings.warn(
            "The following requested variables were not found " +
            ten_plus_message +
            ": " +
            ", ".join(vars_not_found)
            )

    # 9. Sort columns, return fetched data ####
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
    
    print("var_order", var_order, sep = "\n")
    
    # Pass list to fetched_data to order columns accordingly
    fetched_data = fetched_data[var_order]

    return fetched_data
