import numpy as np
from anndata import AnnData
import mudata as md

# --- 1. Define sample IDs ---
sample_ids = [f"sample_{i}" for i in range(1, 13)]

# --- 2. Create gene expression modality ---
# simulate expression of 500 genes for 100 cells
gene_expression = AnnData(
    X=np.random.poisson(1.0, size=(100, 500))
)
gene_expression.obs_names = [f"cell_RNA_{i}" for i in range(100)]
# assign each cell to one of 12 samples
gene_expression.obs['sample_id'] = np.random.choice(sample_ids, size=gene_expression.n_obs)

# --- 3. Create peak-count modality ---
# simulate binary peak calls for 10,000 peaks over 80 cells
peak_counts = AnnData(
    X=np.random.binomial(1, 0.1, size=(80, 10000))
)
peak_counts.obs_names = [f"cell_ATAC_{i}" for i in range(80)]
# assign each cell to one of 12 samples
peak_counts.obs['sample_id'] = np.random.choice(sample_ids, size=peak_counts.n_obs)

# --- 4. Combine into a MuData container ---
mdata = md.MuData({
    'rna': gene_expression,
    'atac': peak_counts
})

print(mdata)
