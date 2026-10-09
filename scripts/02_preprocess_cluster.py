#!/usr/bin/env python3
"""Stage 2 — QC, preprocessing, and clustering of PBMC3k.

QC filters -> normalize_total -> log1p -> highly variable genes ->
regress_out -> scale -> PCA -> neighbors -> UMAP -> Leiden.
Run from the repo root after: conda activate singlecell-pbmc
"""
import os
from pathlib import Path
import scanpy as sc

BASE = Path(os.environ.get("SC_BASE", Path.home() / "single-cell-pbmc"))
DATA = BASE / "data"
RESULTS = BASE / "results"
RESULTS.mkdir(parents=True, exist_ok=True)

adata = sc.read_h5ad(DATA / "pbmc3k_raw.h5ad")
n_start = adata.n_obs
print(f"starting: {adata.n_obs} cells x {adata.n_vars} genes")

# --- QC metrics ---
adata.var["mt"] = adata.var_names.str.startswith("MT-")
sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], percent_top=None,
                           log1p=False, inplace=True)

# --- Filters: drop low-quality cells ---
adata = adata[adata.obs.n_genes_by_counts < 2500, :]
adata = adata[adata.obs.pct_counts_mt < 5, :]
print(f"QC: {n_start} -> {adata.n_obs} cells retained "
      f"({adata.n_obs / n_start:.1%})")
print(f"median genes/cell after QC: {adata.obs.n_genes_by_counts.median():.0f}")

# --- Normalise ---
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
adata.raw = adata  # keep full gene set for marker discovery

# --- Highly variable genes, regress, scale ---
sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
print(f"highly variable genes: {int(adata.var.highly_variable.sum())}")
adata = adata[:, adata.var.highly_variable]
sc.pp.regress_out(adata, ["total_counts", "pct_counts_mt"])
sc.pp.scale(adata, max_value=10)

# --- PCA / neighbours / UMAP / Leiden ---
sc.tl.pca(adata, svd_solver="arpack")
var40 = float(adata.uns["pca"]["variance_ratio"][:40].sum())
print(f"variance captured by first 40 PCs: {var40:.1%}")
sc.pp.neighbors(adata, n_neighbors=10, n_pcs=40)
sc.tl.umap(adata)
sc.tl.leiden(adata, resolution=0.5)
n_clusters = adata.obs["leiden"].nunique()
print(f"Leiden clusters (resolution=0.5): {n_clusters}")

adata.write(RESULTS / "pbmc3k_processed.h5ad")
print(f"saved: {RESULTS / 'pbmc3k_processed.h5ad'}")

# PASS checkpoint: most cells retained, sane cluster count
assert adata.n_obs / n_start > 0.85, "too many cells lost in QC — inspect"
assert 4 <= n_clusters <= 15, "cluster count outside plausible range — inspect"
print("Stage 2 PASS")
