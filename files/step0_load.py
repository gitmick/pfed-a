"""Step 0: load raw 10x filtered feature-barcode matrix -> AnnData."""
import scanpy as sc
sc.settings.verbosity = 1
adata = sc.read_10x_h5("data/pbmc5k_filtered.h5")
adata.var_names_make_unique()
print("raw:", adata.shape)
adata.write("steps/00_raw.h5ad")
