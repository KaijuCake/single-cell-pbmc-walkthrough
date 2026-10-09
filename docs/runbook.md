# Single-cell Walkthrough Runbook: PBMC3k with Scanpy

A weekend-scale, laptop-friendly single-cell RNA-seq analysis: ~2,700 human
PBMCs from 10x Genomics, processed entirely in Python with Scanpy.

- **Workflow:** download → QC + preprocessing → PCA/UMAP/Leiden clustering → marker genes + cell-type annotation → figures
- **How to use this runbook:** run stages in order; do not start the next stage until the ✅ PASS checkpoint of the current stage holds. All commands assume the repo is cloned to `~/single-cell-pbmc`.
- **Time/disk expectations:** <2 GB disk; the full run takes under an hour on a laptop.

---

## Stage 0 — Environment

```bash
cd ~/single-cell-pbmc
conda env create -f environment.yml
conda activate singlecell-pbmc
python -c "import scanpy; print('scanpy', scanpy.__version__)"
```

**✅ PASS checkpoint:** the environment activates and the scanpy version prints.

---

## Stage 1 — Download PBMC3k

```bash
cd ~/single-cell-pbmc
python scripts/01_download.py
ls -lh data/
```

**✅ PASS checkpoint:** `data/pbmc3k_raw.h5ad` exists, ~2,700 cells × ~1,400 genes.

---

## Stage 2 — QC, preprocessing, clustering

```bash
cd ~/single-cell-pbmc
python scripts/02_preprocess_cluster.py
```

### What the script does, and what to check

1. **QC filters** — drops cells with ≥2,500 detected genes (likely doublets) or ≥5% mitochondrial reads (likely dying cells). Expect to retain >85% of cells.
2. **Normalization** — counts scaled to 10,000 per cell, then log1p.
3. **Feature selection** — highly variable genes only (noise reduction).
4. **Regression + scaling** — removes total-count and mitochondrial effects.
5. **PCA → kNN graph → UMAP → Leiden** — the standard single-cell geometry.

**✅ PASS checkpoint:** the script prints the retained-cell fraction (>85%), the variance captured by 40 PCs, and a cluster count between 4 and 15. A cluster count outside that range means the resolution parameter needs adjusting — inspect before continuing.

---

## Stage 3 — Marker genes, annotation, figures

```bash
cd ~/single-cell-pbmc
python scripts/03_markers_figures.py
ls -lh results/figures/
```

### What to check

- **Marker table** (`results/marker_genes.csv`): each cluster's top markers should include recognizable immune genes (CD14/LYZ for monocytes, MS4A1 for B cells, CD3D/CD8A for T cells, GNLY/NKG7 for NK).
- **Dotplot**: canonical markers should light up in exactly the expected annotated clusters. A marker expressed everywhere is a failed annotation — rename that cluster manually.
- **UMAP**: annotated cell types should form coherent territories; heavy intermixing of two "types" suggests over-clustering.

**✅ PASS checkpoint:** `marker_genes.csv` exists, the annotated types print without a WARNING, and the four figures open cleanly. If the cluster-count WARNING fired, set labels manually against the dotplot before claiming cell types.

---

## GitHub-ready layout

Large `.h5ad` files stay out of git (see `.gitignore`); the processed object is reproducible from the scripts.

```text
single-cell-pbmc-walkthrough/
├── README.md
├── environment.yml
├── .gitignore
├── data/
│   └── pbmc3k_raw.h5ad        # downloaded in Stage 1 (not in git)
├── scripts/
│   ├── 01_download.py
│   ├── 02_preprocess_cluster.py
│   └── 03_markers_figures.py
├── results/
│   ├── pbmc3k_processed.h5ad  # not in git
│   ├── marker_genes.csv
│   └── figures/               # UMAP, dotplot, QC violin
└── docs/
    └── runbook.md
```

Freeze the environment and initialise the repo:

```bash
cd ~/single-cell-pbmc
conda env export --name singlecell-pbmc --file environment.yml
git init
git add README.md environment.yml .gitignore data/README.md scripts docs
git commit -m "Single-cell PBMC walkthrough (Scanpy): QC, Leiden clustering, marker-based annotation"
git status
```

**Final ✅ PASS checkpoint:** every bracketed number in the README came from a file in `results/` that you opened and read — no placeholder brackets remain.
