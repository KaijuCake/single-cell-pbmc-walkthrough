#!/usr/bin/env python3
"""Stage 1 — download the PBMC3k dataset (10x Genomics, ~2,700 PBMCs).

Run from the repo root after: conda activate singlecell-pbmc
"""
import os
from pathlib import Path
import scanpy as sc

BASE = Path(os.environ.get("SC_BASE", Path.home() / "single-cell-pbmc"))
DATA = BASE / "data"
DATA.mkdir(parents=True, exist_ok=True)

adata = sc.datasets.pbmc3k()
print(adata)
adata.write(DATA / "pbmc3k_raw.h5ad")

# PASS checkpoint: exactly one file, ~2,700 cells, ~1,400 genes
print(f"saved: {DATA / 'pbmc3k_raw.h5ad'} "
      f"({adata.n_obs} cells x {adata.n_vars} genes)")
assert adata.n_obs > 2000, "unexpectedly few cells — re-download"
print("Stage 1 PASS")
