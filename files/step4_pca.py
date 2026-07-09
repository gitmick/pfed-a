"""Step 4: subset to HVG, scale, PCA."""
import scanpy as sc
sc.settings.verbosity = 1
adata = sc.read_h5ad("steps/03_hvg.h5ad")
adata = adata[:, adata.var["highly_variable"]].copy()
sc.pp.scale(adata, max_value=10)
sc.tl.pca(adata, n_comps=50, svd_solver="arpack", random_state=0)
print("PCA:", adata.obsm["X_pca"].shape)
adata.write("steps/04_pca.h5ad")
