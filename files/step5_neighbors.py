"""Step 5: neighbor graph -> the PREPARED, analysis-ready state (common ground)."""
import scanpy as sc
sc.settings.verbosity = 1
adata = sc.read_h5ad("steps/04_pca.h5ad")
sc.pp.neighbors(adata, n_neighbors=15, n_pcs=30, random_state=0)
print("prepared:", adata.shape, "| neighbors graph in .obsp")
adata.write("steps/05_prepared.h5ad")
