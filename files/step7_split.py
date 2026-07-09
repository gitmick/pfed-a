"""Step 7: confirm whether a coarse cluster genuinely resolves into two populations.

Primary test: does the coarse MONOCYTE cluster split, at fine resolution, into
CD14+ classical (CD14/FCN1/S100A8/LYZ high) and CD16+/FCGR3A+ non-classical
(FCGR3A/MS4A7 high) monocytes? Fallback: T-cell CD4+/CD8+ split.

Reads the clustered AnnData; writes a decision JSON. Runs no clustering itself.
"""
import json
import numpy as np
import scanpy as sc

adata = sc.read_h5ad("steps/06_clustered.h5ad")
raw = adata.raw
gene_idx = {g: i for i, g in enumerate(raw.var_names)}


def m(mask, gene):
    if gene not in gene_idx:
        return None
    X = raw[mask].X
    col = X[:, gene_idx[gene]]
    col = col.toarray().ravel() if hasattr(col, "toarray") else np.asarray(col).ravel()
    return float(col.mean())


def frac_pos(mask, gene):
    if gene not in gene_idx:
        return None
    X = raw[mask].X
    col = X[:, gene_idx[gene]]
    col = col.toarray().ravel() if hasattr(col, "toarray") else np.asarray(col).ravel()
    return float((col > 0).mean())


coarse = adata.obs["leiden_coarse"]
fine = adata.obs["leiden_fine"]

# --- identify the coarse monocyte cluster: highest combined CD14+LYZ+S100A8 mean ---
mono_score = {}
t_score = {}
for cc in coarse.cat.categories:
    mask = (coarse == cc).values
    mono_score[cc] = np.mean([m(mask, g) for g in ["CD14", "LYZ", "S100A8", "FCN1"]])
    t_score[cc] = np.mean([m(mask, g) for g in ["CD3D", "CD3E"]])
mono_coarse = max(mono_score, key=mono_score.get)
t_coarse = max(t_score, key=t_score.get)


def analyze_split(coarse_label, subtype_genes):
    """For the given coarse cluster, look at the fine clusters it contains and
    report per-fine-cluster marker means for the two candidate subpopulations."""
    mask_c = (coarse == coarse_label).values
    fine_here = fine[mask_c].value_counts()
    fine_here = fine_here[fine_here >= 20]  # ignore trace spill
    rows = {}
    for fc in fine_here.index:
        mask = (coarse == coarse_label).values & (fine == fc).values
        row = {"n": int(mask.sum())}
        for grp, genes in subtype_genes.items():
            for g in genes:
                row[f"{grp}:{g}"] = m(mask, g)
                row[f"{grp}:{g}:frac+"] = frac_pos(mask, g)
        rows[str(fc)] = row
    return rows


mono_genes = {
    "classical": ["CD14", "FCN1", "S100A8", "LYZ"],
    "nonclassical": ["FCGR3A", "MS4A7"],
}
t_genes = {
    "CD4": ["CD4", "IL7R"],
    "CD8": ["CD8A", "CD8B"],
}

mono_rows = analyze_split(mono_coarse, mono_genes)
t_rows = analyze_split(t_coarse, t_genes)


def decide_mono(rows):
    """A genuine split = >=1 fine cluster is classical-dominant AND >=1 is
    non-classical-dominant. Comparative, not fixed cutoffs:
      classical    = CD14 > FCGR3A and (FCN1 or S100A8) clearly expressed
      non-classical= FCGR3A > CD14 and FCGR3A high and MS4A7 expressed
    (only fine clusters that actually belong to the monocyte coarse cluster)."""
    classical, nonclassical = [], []
    for fc, r in rows.items():
        cd14 = r.get("classical:CD14") or 0
        fcn1 = r.get("classical:FCN1") or 0
        s100a8 = r.get("classical:S100A8") or 0
        fcgr3a = r.get("nonclassical:FCGR3A") or 0
        ms4a7 = r.get("nonclassical:MS4A7") or 0
        if fcgr3a > cd14 and fcgr3a > 0.8 and ms4a7 > 0.4:
            nonclassical.append(fc)
        elif cd14 > fcgr3a and (fcn1 > 0.8 or s100a8 > 0.8):
            classical.append(fc)
    return classical, nonclassical


mono_classical, mono_nonclassical = decide_mono(mono_rows)
mono_splits = bool(mono_classical) and bool(mono_nonclassical)

out = {
    "coarse_monocyte_cluster": mono_coarse,
    "coarse_tcell_cluster": t_coarse,
    "mono_score_per_coarse": {k: round(v, 3) for k, v in mono_score.items()},
    "monocyte_fine_breakdown": mono_rows,
    "tcell_fine_breakdown": t_rows,
    "decision": {
        "monocyte_splits": mono_splits,
        "classical_fine_clusters": mono_classical,
        "nonclassical_fine_clusters": mono_nonclassical,
    },
}
with open("steps/07_split_decision.json", "w") as fh:
    json.dump(out, fh, indent=2, default=str)

print(json.dumps(out["decision"], indent=2))
print("coarse monocyte cluster =", mono_coarse, "| coarse T cluster =", t_coarse)
print("wrote steps/07_split_decision.json")
