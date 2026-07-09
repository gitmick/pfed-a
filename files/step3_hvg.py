"""Step 3: highly-variable-gene selection (Seurat flavor, top 2000)."""
import scanpy as sc
sc.settings.verbosity = 1
adata = sc.read_h5ad("steps/02_norm.h5ad")
sc.pp.highly_variable_genes(adata, n_top_genes=2000, flavor="seurat")
print("HVG flagged:", int(adata.var["highly_variable"].sum()))
adata.write("steps/03_hvg.h5ad")
