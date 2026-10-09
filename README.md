# What immune cells are in human blood? A single-cell tour (PBMC3k, Scanpy)

Your blood contains millions of immune cells — T cells, B cells, monocytes,
natural killer cells — all mixed together, all looking fairly similar under a
microscope. **Single-cell RNA sequencing** reads, for each individual cell,
which genes are switched on. Think of it as taking attendance: instead of
asking "what is the blood doing on average?", we ask 2,700 individual cells
"what are *you* doing right now?"

This project takes one public dataset of ~2,700 human blood cells and answers
three questions:

1. **How many distinct cell types are in the mix?** (We let the data decide —
   no labels provided up front.)
2. **What are they?** (We identify each group using genes known to mark
   specific immune cells.)
3. **Which genes define each type?** (The "marker genes" — the molecular ID
   badges of each cell type.)

## Why does it matter?

Almost every modern immunology paper — cancer immunotherapy, autoimmune
disease, vaccine response — uses exactly this workflow to figure out *which*
cells are doing *what*. It is the standard lens for looking at the immune
system one cell at a time.

## What I did, step by step

```mermaid
flowchart LR
    A[Download 2,700 cells] --> B[Quality control]
    B --> C[Normalize]
    C --> D[Find variable genes]
    D --> E[Map cells in 2D]
    E --> F[Group into clusters]
    F --> G[Name each group]
    G --> H[Figures]
```

In plain language:

1. **Download** — fetched a public dataset (10x Genomics PBMC3k) programmatically.
2. **Quality control** — threw out broken or doubled-up cells (too few genes
   detected, or too many mitochondrial reads, both signs of a dying cell).
3. **Normalize** — put every cell on the same scale, so a cell sequenced more
   deeply doesn't look "more active" than its neighbors.
4. **Find the interesting genes** — most genes are boring (on in every cell);
   kept the ~2,000 that vary most between cells, since those carry the signal.
5. **Map and group** — placed cells on a 2D map (UMAP) where similar cells sit
   together, then grouped them with the Leiden algorithm. No labels were used —
   the groups emerged from the data alone.
6. **Name each group** — checked which known immune marker genes each group
   expresses (e.g. CD14 → monocytes, MS4A1 → B cells) and assigned cell types.
7. **Figures** — UMAP maps, a marker-gene dotplot, and QC diagnostics.

*New to the jargon? See [docs/glossary.md](docs/glossary.md).*

## What came out (from my run)

- Cells kept after quality control: **[N] of 2,700 ([P]%)**
- Cell groups found: **[K]**
- Cell types identified: **[e.g. CD4 T cells, CD8 T cells, B cells, CD14+ monocytes, FCGR3A+ monocytes, NK cells, dendritic cells]**
- Example marker genes: **[e.g. LYZ and CD14 mark monocytes; MS4A1 marks B cells; GNLY marks NK cells]**
- Full results: `results/marker_genes.csv` · figures in `results/figures/`

## Reproduce it (technical)

```bash
git clone <this-repo> ~/single-cell-pbmc
cd ~/single-cell-pbmc
conda env create -f environment.yml
conda activate singlecell-pbmc
python scripts/01_download.py        # Stage 1: fetch PBMC3k
python scripts/02_preprocess_cluster.py  # Stage 2: QC → UMAP → Leiden
python scripts/03_markers_figures.py    # Stage 3: markers → annotation → figures
```

Each script ends with a pass/fail checkpoint — see [docs/runbook.md](docs/runbook.md)
for what each checkpoint verifies. Large data files (`.h5ad`) are reproducible
from the scripts and stay out of git.

**Environment:** Python 3.11, Scanpy, Leiden algorithm via `leidenalg`
(see `environment.yml`).

## What this demonstrates

End-to-end single-cell analysis in Python: QC judgement on real data,
dimensionality reduction, graph-based clustering, statistical marker
discovery, and biology-grounded cell-type annotation — the standard workflow
behind most human single-cell papers.

## Limitations

- Reanalysis of a small public dataset; cluster count depends on the
  resolution parameter, which is recorded in the scripts.
- Cell-type labels are annotation by canonical markers, not ground truth —
  a marker expressed everywhere is a failed annotation, and the scripts warn
  rather than mislabel.
