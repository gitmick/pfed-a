"""Step 1: QC. Drop low-count cells, high-mito cells, and likely doublets."""
import numpy as np
import scanpy as sc
sc.settings.verbosity = 1
np.random.seed(0)

adata = sc.read_h5ad("steps/00_raw.h5ad")
n0 = adata.n_obs

# basic gene filter
sc.pp.filter_genes(adata, min_cells=3)

# QC metrics; mito genes by MT- prefix
adata.var["mt"] = adata.var_names.str.startswith("MT-")
sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], inplace=True, percent_top=None)

# low-count / low-gene cells
sc.pp.filter_cells(adata, min_genes=200)
# upper bounds to remove likely multiplets / outliers
adata = adata[adata.obs["n_genes_by_counts"] < 5000, :].copy()
# high mitochondrial fraction
adata = adata[adata.obs["pct_counts_mt"] < 15, :].copy()
n_after_basic = adata.n_obs

# doublet detection (Scrublet), drop predicted doublets
try:
    sc.pp.scrublet(adata, random_state=0)
    n_doublet = int(adata.obs["predicted_doublet"].sum())
    adata = adata[~adata.obs["predicted_doublet"], :].copy()
    doublet_method = "scrublet"
except Exception as e:
    doublet_method = "scrublet-failed:" + type(e).__name__
    n_doublet = -1

print(f"cells: raw={n0} after_basicQC={n_after_basic} doublets_removed={n_doublet} final={adata.n_obs} genes={adata.n_vars} doublet_method={doublet_method}")
adata.write("steps/01_qc.h5ad")
