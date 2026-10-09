# Glossary — single-cell RNA-seq in plain language

**Single-cell RNA sequencing (scRNA-seq)**
Reading which genes are switched on, separately for each individual cell —
instead of averaging over a whole tissue sample.

**Gene expression**
Whether a gene is "switched on" and how strongly. A cell's pattern of
switched-on genes largely determines what the cell is and does.

**PBMC (peripheral blood mononuclear cell)**
The immune cells circulating in your blood: T cells, B cells, monocytes,
natural killer (NK) cells, and a few others. "PBMC3k" = ~3,000 of them
(2,700 after the dataset's own filtering).

**Quality control (QC)**
Throwing out bad measurements before analysis. Two main filters here:
cells with suspiciously few detected genes, and cells with a high share of
*mitochondrial* reads.

**Mitochondrial reads**
Reads coming from the mitochondria (the cell's power plants) rather than the
nucleus. A high share usually means the cell was dying and leaking its
contents when captured — so we discard those cells.

**Doublets**
Two cells accidentally captured as one. They look like a cell expressing
two cell types' genes at once; the gene-count filter catches most of them.

**Normalization**
Putting every cell on the same scale (here: 10,000 counts per cell), so a
cell that was simply sequenced more deeply doesn't look more active.

**Highly variable genes**
The ~2,000 genes whose activity differs most between cells. Most genes are
"housekeeping" (on everywhere, uninformative); these are the ones that
distinguish cell types.

**PCA (principal component analysis)**
Compresses thousands of gene measurements into a smaller set of summary
axes ("principal components") that capture the biggest differences
between cells.

**kNN graph**
Each cell is linked to its ~10 most similar neighbors. Clustering then
becomes a graph problem: find densely connected neighborhoods.

**UMAP**
A way to draw thousands of cells on a flat 2D map so that similar cells
land near each other. Distances on the map are approximate — treat
neighborhoods as real, exact gaps with suspicion.

**Leiden clustering**
The algorithm that cuts the kNN graph into groups ("clusters") of similar
cells. A *resolution* parameter controls how fine the groups are; higher
resolution = more, smaller groups.

**Marker genes**
Genes that are strongly switched on in one group and off in the others —
the molecular ID badges used to recognize cell types.

**Wilcoxon rank-sum test**
The statistical test used here to find marker genes: for each gene, it
asks "is this gene consistently more active in this group than in all
other cells?"

**Dotplot**
A figure where each dot's size = share of cells expressing a gene and its
color = expression strength. Lets you verify at a glance that each marker
lights up in exactly the expected cell type.

**Annotation**
Assigning biological names ("CD14+ monocyte", "B cell") to the numbered
clusters, based on their marker genes. The weakest link in the chain —
treated here as a hypothesis to check, not ground truth.
