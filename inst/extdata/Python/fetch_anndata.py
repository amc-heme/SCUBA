# Libraries used
# Numpy, Pandas
import pandas as pd
import numpy as np

# Scipy
from scipy.sparse import csr_matrix, csc_matrix
# Anndata CSCDataset class (for disk-backed anndata objects)
from anndata.abc import CSCDataset, CSRDataset

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
                keyless_vars = [var for var in keyless_vars if var in obj.var_names]
            elif key == "obs":
                # Metadata (obs)
                matrix = obj.obs
                
                # Subset for vars in matrix
                
            elif key in obj.obsm_keys():
                # For a key in the list of obsm keys, pull matrix for that entry
                matrix = obj.obsm[key]
                
                # Subset for vars in matrix (depending on matrix types)
                if (isinstance(matrix, pd.DataFrame)):
                    # Pandas dataframes: use .columns method to pull valid 
                    # entries
                    keyless_vars = [var for var in keyless_vars if var in matrix.columns]
                elif (
                    isinstance(matrix, np.ndarray) |
                    isinstance(matrix, csr_matrix) | 
                    isinstance(matrix, csc_matrix) |
                    isinstance(matrix, CSCDataset)
                ):
                    # For numpy arrays (most likely, used for reductions), 
                    # and sparse matrices (less likely but possible)
                    # No column names, valid entries are based on index 
                    # Valid vars are the index values as strings, using a 
                    # 1-based index (for consistency with R SCUBA methods)
                    valid_idx = [str(x) for x in range(1, matrix.shape[1] + 1)] 
                    
                    keyless_vars = [var for var in keyless_vars if var in valid_idx]
                    
            ### 2.2.2. Pull data for cells, keyless vars from matrix ####
            # Type checking via isinstance is used here due to the variety of
            # data types possible for the matrix
            # This is not considered "pythonic" and may need to be re-thought 
            # https://stackoverflow.com/questions/2225038/determine-the-type-of-an-object
            if (isinstance(matrix, csr_matrix) | isinstance(matrix, csc_matrix)):
                # Sparse matrix format (X and layers) 
                # subset *object for genes, then construct a pandas dataframe 
                # (sparse matrices don't subset easily in python)
                if layer == None:
                    # Subset object and pull from X if layer is underfined
                    data = pd.DataFrame.sparse.from_spmatrix(
                        obj[cells, keyless_vars].X,
                        # csr_matrices don't have row, column names. 
                        # These are added here
                        index = cells,
                        columns = keyless_vars
                        )
                else:
                    # Otherwise, subset and pull from the specified layer
                    data = pd.DataFrame.sparse.from_spmatrix(
                        obj[cells, keyless_vars].layers[layer],
                        # csr_matrices don't have row, column names. 
                        # These are added here
                        index = cells,
                        columns = keyless_vars
                        )
                    
                # Columns returned will be pandas sparse arrays.
                # Arrays must be densified to be accesssible downstream in R 
                # Densify with np.asarray
                for column in data:
                    data[column] = np.asarray(data[column])
            
            elif (isinstance(matrix, CSCDataset) | isinstance(matrix, CSRDataset)):
                # Disk-backed anndata: X matrix uses _CSCDataset class from 
                # anndata. This class is considered internal to the anndata 
                # package and may change in the future
                
                # _CSCDataset matrices only support accessing cells and 
                # features by index. To get the index of each var requested,
                # get_loc is used
                
                # To avoid an error from get_loc(), filter keyless vars for 
                # only those that are in var_names
                vars_in_X = [var for var in keyless_vars if var in obj.var_names]
                
                vars_idx = [obj.var_names.get_loc(var) for var in vars_in_X]
                cells_idx = [obj.obs_names.get_loc(cell_id) for cell_id in cells]

                # Index the CSCDataset. The result will be a scipy csc_matrix
                backed_data = matrix[cells_idx, vars_idx]

                data = pd.DataFrame.sparse.from_spmatrix(
                    backed_data,
                    # Add cell ids and var names back 
                    # These should be in the same order as the indices 
                    index = cells,
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
                
                # Convert numpy array to pandas dataframe, then slice for vars
                matrix = pd.DataFrame(
                    matrix,
                    index = obj.obs_names,
                    # Construct var names using a 1-based index (for consistency 
                    # with Seurat, R objects used in SCUBA)
                    columns = [str(x) for x in range(1, matrix.shape[1] + 1)]
                    )
                    
                # Slice based on keyless vars (should be string format of the 
                # column desired, "1" for the first column)
                try:
                    data = matrix.loc[cells, keyless_vars]
                except KeyError:
                    # If any values are not present, a KeyError is returned.
                    # If this happens, don't pull anything and keep going
                    # (the assumption is that for similarly named matrices
                    # like mofa and mofa_umap, mofa will return an error 
                    # when searched, but mofa_umap will not). This may result
                    # in values that are there not being found, but fixes 
                    # SCUBA #90.
                    # Set key_error to TRUE to prevent attempts to access the
                    # non-existent data variable that would have been created
                    key_error = True
                    pass
                  
            else:
                """
                Strings wrapped according to PEP-8 style guidelines
                https://stackoverflow.com/questions/1874592/
                how-to-write-very-long-string-that-conforms-with-
                pep8-and-prevent-e501
                """
                raise NotImplementedError(
                    ("The SCUBA Python function fetch_anndata does not know "
                    "how to handle a matrix of class {0}. "
                    "This ocurred with the matrix at key '{1}'."
                    ).format(type(matrix), key)
                    )
            
            # Add location key back in to variable names and concatenate data
            # to data already fetched
            # (if there are no errors fetching keyed_vars and ncol of the data
            # frame is greater than zero)
            if key_error == False:
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
