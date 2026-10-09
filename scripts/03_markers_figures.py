#!/usr/bin/env python3
"""Stage 3 — marker genes, cell-type annotation, and figures.

Wilcoxon marker discovery per Leiden cluster, annotation with canonical
immune markers, and publication-style figures.
Run from the repo root after: conda activate singlecell-pbmc
"""
import os
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import scanpy as sc

BASE = Path(os.environ.get("SC_BASE", Path.home() / "single-cell-pbmc"))
RESULTS = BASE / "results"
FIG = RESULTS / "figures"
FIG.mkdir(parents=True, exist_ok=True)
sc.settings.figdir = FIG

adata = sc.read_h5ad(RESULTS / "pbmc3k_processed.h5ad")

# --- Marker genes per cluster ---
sc.tl.rank_genes_groups(adata, "leiden", method="wilcoxon")
markers = sc.get.rank_genes_groups_df(adata, group=None)
top = markers.groupby("group").head(5)[["group", "names", "scores"]]
print("top 5 markers per cluster:")
print(top.to_string(index=False))
markers.to_csv(RESULTS / "marker_genes.csv", index=False)

# --- Annotate with canonical immune markers ---
marker_genes = {
    "CD14+ Mono": ["CD14", "LYZ"],
    "CD19+ B": ["MS4A1"],
    "CD8 T": ["CD8A"],
    "NK": ["GNLY", "NKG7"],
    "FCGR3A+ Mono": ["FCGR3A", "MS4A7"],
}
# new cluster names, ordered to match this dataset's Leiden clusters
new_cluster_names = [
    "CD4 T", "CD14+ Mono", "B", "CD8 T",
    "CD4 T", "NK", "CD8 T", "FCGR3A+ Mono", "Dendritic",
]
if len(new_cluster_names) == adata.obs["leiden"].nunique():
    adata.rename_categories("leiden", new_cluster_names)
    print("annotated cell types:", sorted(adata.obs["leiden"].unique().tolist()))
else:
    print(f"WARNING: found {adata.obs['leiden'].nunique()} clusters, "
          f"expected {len(new_cluster_names)} — skipping rename; "
          "inspect the UMAP and set labels manually")

# --- Figures ---
sc.pl.umap(adata, color="leiden", legend_loc="on data",
           title="PBMC3k — Leiden clusters", save="_leiden.png", show=False)
sc.pl.umap(adata, color=["CD14", "MS4A1", "CD8A", "GNLY"],
           save="_canonical_markers.png", show=False)
sc.pl.dotplot(adata, marker_genes, groupby="leiden", standard_scale="var",
              save="marker_dotplot.png", show=False)
sc.pl.violin(adata, ["n_genes_by_counts", "total_counts", "pct_counts_mt"],
             jitter=0.4, multi_panel=True, save="_qc_violin.png", show=False)
plt.close("all")

print("figures written to", FIG)
print("Stage 3 PASS: marker_genes.csv exists, cell types annotated, "
      "UMAP/dotplot/QC figures saved")
