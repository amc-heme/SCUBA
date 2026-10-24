# Pseudobulk aggregation of an anndata object for SCUBA::pseudobulk().
#
# The group (profile) of every cell is decided in R and arrives as 1-based integer codes
# (0 = cell left out). Sums are accumulated chunk by chunk through a sparse indicator product,
# so neither an in-memory nor a backed matrix is ever densified as a whole:
#   - row-major sources (csr_matrix, ndarray, backed CSRDataset, backed dense datasets, and
#     in-memory csc_matrix) are read in blocks of `chunk_cells` consecutive cells;
#   - a backed CSCDataset is read in blocks of consecutive features (its cheap axis).
import numpy as np
import pandas as pd
from scipy.sparse import csr_matrix, issparse

try:
    from anndata.abc import CSRDataset, CSCDataset
except ImportError:  # anndata < 0.11
    from anndata._core.sparse_dataset import CSRDataset, CSCDataset


def _as_index_array(x):
    arr = np.atleast_1d(np.asarray(x))
    return arr.astype(np.int64)


def _matrix(obj, modality, layer):
    """The matrix to aggregate and its feature names."""
    if modality is None or modality == "X":
        if layer is None or layer == "X":
            return obj.X, obj.var_names
        if layer == "raw":
            if obj.raw is None:
                raise ValueError("The object has no raw matrix (adata.raw is None)")
            return obj.raw.X, obj.raw.var_names
        if layer not in obj.layers:
            available = ["X"] + [str(k) for k in obj.layers.keys() if k is not None] + (["raw"] if obj.raw is not None else [])
            raise ValueError(
                "Layer '%s' not found; available: %s. Pass the layer holding raw counts "
                "(layer = \"X\" when X holds counts)." % (layer, ", ".join(available))
            )
        return obj.layers[layer], obj.var_names
    if modality == "raw":
        if obj.raw is None:
            raise ValueError("The object has no raw matrix (adata.raw is None)")
        return obj.raw.X, obj.raw.var_names
    if modality in obj.obsm_keys() and isinstance(obj.obsm[modality], pd.DataFrame):
        df = obj.obsm[modality]
        return df, df.columns
    raise ValueError("Modality '%s' not found; modalities: %s" % (modality, ", ".join(["X"] + list(obj.obsm_keys()))))


def _dense_rows(block):
    """A chunk of rows as something `indicator @ block` accepts (scipy sparse or ndarray)."""
    if hasattr(block, "to_memory"):
        block = block.to_memory()
    if isinstance(block, pd.DataFrame):
        block = block.to_numpy(dtype=np.float64)
    if issparse(block):
        return block
    return np.asarray(block, dtype=np.float64)


def aggregate_anndata(obj, codes, n_groups, modality=None, layer=None, cells=None, features=None, chunk_cells=5000):
    """
    Sum every feature over groups of cells.

    obj: anndata.AnnData (in memory or backed).
    codes: 1-based group of each cell (0 = skip), aligned with `cells`, or with obs_names when
        cells is None.
    n_groups: number of groups.
    modality: None / "X" (the X matrix and its layers), "raw", or the key of an obsm data frame.
    layer: None / "X", a key of adata.layers, or "raw".
    cells: optional cell names (any order); codes must have the same length.
    features: optional feature names to keep (the result keeps this order).
    chunk_cells: cells per block for the row-major path.

    Returns {"sum": n_groups x n_features float64 array, "n": cells per group, "features": names}.
    """
    codes = _as_index_array(codes)
    n_groups = int(np.atleast_1d(np.asarray(n_groups))[0])
    chunk_cells = max(1, int(np.atleast_1d(np.asarray(chunk_cells))[0]))
    matrix, var_names = _matrix(obj, modality, layer)
    n_obs = obj.n_obs
    full = np.zeros(n_obs, dtype=np.int64)
    if cells is None:
        if len(codes) != n_obs:
            raise ValueError("codes has %d entries for %d cells" % (len(codes), n_obs))
        full[:] = codes
    else:
        cells = [str(c) for c in (cells if not isinstance(cells, str) else [cells])]
        if len(codes) != len(cells):
            raise ValueError("codes has %d entries for %d cells" % (len(codes), len(cells)))
        idx = obj.obs_names.get_indexer(cells)
        if (idx < 0).any():
            missing = [c for c, i in zip(cells, idx) if i < 0]
            raise ValueError("Cells not found: %s" % ", ".join(missing[:10]))
        full[idx] = codes
    if codes.size and (codes.min() < 0 or codes.max() > n_groups):
        raise ValueError("group codes must be in 0:%d" % n_groups)

    var_index = pd.Index([str(v) for v in var_names])
    if features is None:
        feat_idx = np.arange(len(var_index))
        feat_names = list(var_index)
    else:
        features = [str(f) for f in (features if not isinstance(features, str) else [features])]
        feat_idx = var_index.get_indexer(features)
        if (feat_idx < 0).any():
            missing = [f for f, i in zip(features, feat_idx) if i < 0]
            raise ValueError("Features not found: %s" % ", ".join(missing[:10]))
        feat_names = features
    subset_features = len(feat_idx) != len(var_index) or not np.array_equal(feat_idx, np.arange(len(var_index)))

    acc = np.zeros((n_groups, len(feat_idx)), dtype=np.float64)
    n = np.bincount(full[full > 0] - 1, minlength=n_groups) if n_groups > 0 else np.zeros(0, dtype=np.int64)
    if n_groups == 0 or len(feat_idx) == 0:
        return {"sum": acc, "n": n, "features": feat_names}

    if isinstance(matrix, CSCDataset):
        # backed CSC: columns are the cheap axis, so walk consecutive features with one
        # indicator over all cells
        rows = np.flatnonzero(full)
        indicator = csr_matrix((np.ones(len(rows)), (full[rows] - 1, rows)), shape=(n_groups, n_obs))
        n_var = matrix.shape[1]
        chunk_feats = max(1, (chunk_cells * 64) // max(1, n_groups))
        for a in range(0, n_var, chunk_feats):
            b = min(a + chunk_feats, n_var)
            inside = np.flatnonzero((feat_idx >= a) & (feat_idx < b))
            if inside.size == 0:
                continue
            block = _dense_rows(matrix[:, a:b])
            prod = indicator @ block
            prod = prod.toarray() if issparse(prod) else np.asarray(prod)
            acc[:, inside] += prod[:, feat_idx[inside] - a]
        return {"sum": acc, "n": n, "features": feat_names}

    # row-major path: consecutive cells
    for a in range(0, n_obs, chunk_cells):
        b = min(a + chunk_cells, n_obs)
        block_codes = full[a:b]
        rows = np.flatnonzero(block_codes)
        if rows.size == 0:
            continue
        block = _dense_rows(matrix[a:b])
        if subset_features:
            block = block[:, feat_idx]
        indicator = csr_matrix((np.ones(len(rows)), (block_codes[rows] - 1, rows)), shape=(n_groups, b - a))
        prod = indicator @ block
        acc += prod.toarray() if issparse(prod) else np.asarray(prod, dtype=np.float64)
    return {"sum": acc, "n": n, "features": feat_names}
