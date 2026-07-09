"""Step 2: normalize (library-size to 1e4) + log1p."""
import scanpy as sc
sc.settings.verbosity = 1
adata = sc.read_h5ad("steps/01_qc.h5ad")
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
adata.raw = adata  # keep log-normalized full matrix for marker expression later
print("normalized+log1p:", adata.shape)
adata.write("steps/02_norm.h5ad")
