"""Step 6: cluster at a coarse and a finer resolution; annotate by canonical markers."""
import json
import numpy as np
import scanpy as sc
sc.settings.verbosity = 1

adata = sc.read_h5ad("steps/05_prepared.h5ad")

# coarse + fine leiden. With scanpy's igraph-flavor leiden, resolution runs
# finer than the classic leidenalg default, so a genuinely COARSE grouping
# (each major lineage = one cluster, monocytes together) needs res~0.05, and a
# finer subtype resolution that splits the monocytes needs res~0.2.
COARSE_RES, FINE_RES = 0.05, 0.2
sc.tl.leiden(adata, resolution=COARSE_RES, key_added="leiden_coarse",
             flavor="igraph", n_iterations=2, directed=False, random_state=0)
sc.tl.leiden(adata, resolution=FINE_RES, key_added="leiden_fine",
             flavor="igraph", n_iterations=2, directed=False, random_state=0)
adata.uns["cluster_resolutions"] = {"coarse": COARSE_RES, "fine": FINE_RES}

print("coarse clusters:", adata.obs["leiden_coarse"].value_counts().to_dict())
print("fine clusters:", adata.obs["leiden_fine"].value_counts().to_dict())

# canonical PBMC markers
markers = {
    "T_CD3": ["CD3D", "CD3E", "CD3G"],
    "CD4T": ["CD4", "IL7R"],
    "CD8T": ["CD8A", "CD8B"],
    "NK": ["GNLY", "NKG7", "KLRD1"],
    "B": ["CD79A", "MS4A1", "CD79B"],
    "Mono_classical": ["CD14", "LYZ", "S100A8", "FCN1"],
    "Mono_nonclassical": ["FCGR3A", "MS4A7"],
    "DC": ["FCER1A", "CST3"],
    "Platelet": ["PPBP"],
}

# mean log-normalized expression per cluster, from .raw (full gene set)
def mean_expr_per_cluster(key):
    raw = adata.raw
    groups = adata.obs[key].cat.categories.tolist()
    out = {}
    for g in groups:
        mask = (adata.obs[key] == g).values
        sub = raw[mask]
        X = sub.X
        X = X.toarray() if hasattr(X, "toarray") else np.asarray(X)
        gene_idx = {gn: i for i, gn in enumerate(raw.var_names)}
        row = {}
        for label, genes in markers.items():
            vals = [float(X[:, gene_idx[gn]].mean()) for gn in genes if gn in gene_idx]
            row[label] = float(np.mean(vals)) if vals else None
        # per-gene detail for the key mono markers
        for gn in ["CD14", "FCN1", "S100A8", "LYZ", "FCGR3A", "MS4A7",
                   "CD3D", "CD4", "IL7R", "CD8A", "CD8B", "MS4A1", "GNLY", "NKG7"]:
            if gn in gene_idx:
                row["g:" + gn] = float(X[:, gene_idx[gn]].mean())
        row["_n"] = int(mask.sum())
        out[g] = row
    return out

summary = {
    "n_cells": int(adata.n_obs),
    "coarse_counts": {k: int(v) for k, v in adata.obs["leiden_coarse"].value_counts().items()},
    "fine_counts": {k: int(v) for k, v in adata.obs["leiden_fine"].value_counts().items()},
    "coarse_marker_means": mean_expr_per_cluster("leiden_coarse"),
    "fine_marker_means": mean_expr_per_cluster("leiden_fine"),
    # crosstab: which fine clusters live inside each coarse cluster
    "coarse_to_fine": {
        cc: adata.obs.loc[adata.obs["leiden_coarse"] == cc, "leiden_fine"].value_counts().to_dict()
        for cc in adata.obs["leiden_coarse"].cat.categories
    },
}
with open("steps/06_cluster_markers.json", "w") as fh:
    json.dump(summary, fh, indent=2, default=str)

adata.write("steps/06_clustered.h5ad")
print("wrote steps/06_clustered.h5ad and steps/06_cluster_markers.json")
