# EasySC

**Easy Single-Cell RNA Sequencing Analysis Pipeline**

EasySC is a streamlined, user-friendly R package that wraps popular Bioconductor and CRAN tools into a cohesive workflow for scRNA-seq data analysis. It provides one-step functions for quality control, normalisation, batch correction, clustering, doublet removal, and ambient RNA decontamination.

---

## Features

| Step | Function | Purpose |
|------|----------|---------|
| 1 | `Easy_DbFinder` | Doublet detection & filtering (optional) |
| 2 | `Easy_QC` | Compute MT / RP / HB percentages; filter low-quality cells |
| 3 | `Easy_Normal` | Log-normalisation, variable-feature selection, scaling; optional SCTransform |
| 4 | `Easy_harmony` | PCA + Harmony batch correction |
| 5 | `Easy_PC` | Determine optimal number of principal components |
| 6 | `RunUMAP` / `RunTSNE` | Non-linear dimensionality reduction (Seurat native) |
| 7 | `Easy_clustree` | Multi-resolution clustering + cluster-tree visualisation |
| 8 | `Easy_CDC` | Ambient-RNA decontamination with scCDC + decontX (optional) |

---

## Installation

### 1. Install CRAN dependencies

```r
install.packages(c("Seurat", "ggplot2", "cli", "tidyr", "dplyr",
                   "ggpubr", "clustree", "harmony"))
```

### 2. Install Bioconductor dependencies

```r
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install(c("scDblFinder", "BiocParallel", "scCDC",
                       "decontX", "glmGamPoi", "SCTransform"))
```

### 3. Install EasySC

```r
# From GitHub (once published)
# devtools::install_github("cycy12345/EasySC")

# Or install from local source
devtools::install("path/to/EasySC")
```

---

## Standard Analysis Workflow

Below is the recommended end-to-end pipeline. Each function returns either a modified `Seurat` object or a named list containing the object and plots.

```r
library(EasySC)
library(Seurat)

# ------------------------------
# 1. Doublet removal (optional)
# ------------------------------
res_db   <- Easy_DbFinder(sce, group_by = "Patient", dbr = 0.08)
res_db$plot
sce      <- res_db$sce_filter

# ------------------------------
# 2. Quality control
# ------------------------------
res_qc   <- Easy_QC(sce, species = "human", MT = 5, RP = 5, HB = 1,
                    Feature_high = 5000, Feature_low = 200,
                    plot_group = "Patient")
res_qc$plot_Count_Feature
res_qc$plot_percent
sce      <- res_qc$sce_filter

# ------------------------------
# 3. Normalisation
# ------------------------------
sce_norm <- Easy_Normal(sce, use_sct = FALSE)

# ------------------------------
# 4. Batch correction (Harmony)
# ------------------------------
sce_norm <- Easy_harmony(sce_norm, group.by = "Patient")

# ------------------------------
# 5. Choose optimal PC number
# ------------------------------
PC <- Easy_PC(sce_norm, reduction = "pca", cum = 90, var = 5)

# ------------------------------
# 6. Non-linear dimensionality reduction (choose one)
# ------------------------------
sce_norm <- RunUMAP(sce_norm,  reduction = "harmony", dims = 1:PC)
# sce_norm <- RunTSNE(sce_norm, reduction = "harmony", dims = 1:PC)

# ------------------------------
# 7. Multi-resolution clustering & cluster tree
# ------------------------------
res_cl   <- Easy_clustree(sce_norm, assay = "RNA", PC = PC)
sce_cl   <- res_cl$object
clustree_plot <- res_cl$plot

# Inspect the cluster tree and pick a resolution
colnames(sce_cl@meta.data)   # e.g. RNA_snn_res.0.1, RNA_snn_res.0.2, ...

# ------------------------------
# 8. Decontamination (optional)
# ------------------------------
sce_CDC  <- Easy_CDC(sce_cl, cluster_col = "RNA_snn_res.0.1",
                     decontX_threshold = 0.2)

# !! After Easy_CDC you MUST re-run steps 3-7 !!
sce_CDC  <- Easy_Normal(sce_CDC)
sce_CDC  <- Easy_harmony(sce_CDC, group.by = "Patient")
PC2      <- Easy_PC(sce_CDC)
sce_CDC  <- RunUMAP(sce_CDC, reduction = "harmony", dims = 1:PC2)
```

---

## Function Reference

### `Easy_DbFinder(object, dbr = 0.08, group_by = NULL)`
Detects and removes doublets using **scDblFinder**.

- **dbr** -- expected doublet rate (default 8%). Rule of thumb for 10x data: ~1% per 1,000 cells.
- **group_by** -- metadata column for per-sample detection (e.g. `"Patient"`). If `NULL`, all cells are pooled.
- **Returns**: `list(sce_filter = <Seurat>, plot = <ggplot>)`

---

### `Easy_QC(object, species = "human", MT = 5, RP = 5, HB = 1, ...)`
Computes mitochondrial, ribosomal, and haemoglobin percentages and filters cells.

- **species** -- `"human"` or `"mouse"`; automatically sets gene-name patterns.
- **MT / RP / HB** -- maximum allowed percentages.
- **Feature_high / Feature_low** -- upper and lower bounds for `nFeature_RNA`.
- **plot_group** -- metadata column used to group violin plots.
- **Returns**: `list(sce_filter = <Seurat>, plot_Count_Feature = <ggplot>, plot_percent = <ggplot>)`

---

### `Easy_Normal(object, use_sct = FALSE, vars.to.regress = NULL, ...)`
Standard Seurat normalisation pipeline: `NormalizeData` -> `FindVariableFeatures` -> `ScaleData`.

- **use_sct** -- if `TRUE`, additionally runs `SCTransform` and sets the `SCT` assay as default.
- **vars.to.regress** -- variables to regress out (e.g. `c("percent.mt", "nCount_RNA")`).
- **Returns**: modified `Seurat` object.

---

### `Easy_harmony(object, group.by, features = VariableFeatures(object), npcs = 50, ...)`
PCA followed by **Harmony** batch correction.

- **group.by** -- batch variable in metadata (required).
- Internally calls `Easy_PC()` to determine the optimal number of PCs before running Harmony.
- **Returns**: modified `Seurat` object with `reductions$pca` and `reductions$harmony`.

---

### `Easy_PC(seurat, reduction = "pca", cum = 90, var = 5)`
Selects the optimal number of principal components using cumulative-variance and drop-off criteria.

- **cum** -- cumulative variance threshold (%).
- **var** -- individual PC contribution threshold (%).
- Prints an annotated scatter plot to the active graphics device.
- **Returns**: integer (recommended number of PCs).

---

### `Easy_clustree(object, assay = c("RNA", "SCT"), resolutions = ..., PC = 20, ...)`
Multi-resolution clustering and cluster-tree visualisation with **clustree**.

- **resolutions** -- vector of resolution values to test (default `c(0.01, 0.05, 0.1, 0.2, 0.3, 0.5, 0.8, 1)`).
- **PC** -- number of Harmony/PCA dimensions to use for `FindNeighbors`.
- **Returns**: `list(object = <Seurat>, plot = <clustree ggplot>)`

---

### `Easy_CDC(object, cluster_col, assay = "RNA", decontX_threshold = 0.2, ...)`
Ambient-RNA decontamination using **scCDC** (and optionally **decontX**).

- **cluster_col** -- metadata column with cluster labels (e.g. `"RNA_snn_res.0.1"`).
- **decontX_threshold** -- cells with `decontX` contamination >= this value are removed. Set to `NULL` to skip filtering.
- If contamination ratio > 0.0003, only scCDC correction is used.
- If <= 0.0003, decontX is run for deeper cleaning.
- **Important**: after `Easy_CDC`, re-run normalisation, Harmony, UMAP/t-SNE, and clustering.
- **Returns**: decontaminated `Seurat` object.

---

## Dependencies

### Required (Imports)
- `Seurat` (>= 4.0.0)
- `ggplot2` (>= 3.3.0)
- `cli` (>= 3.0.0)
- `tidyr` (>= 1.2.0)
- `dplyr` (>= 1.0.0)
- `ggpubr` (>= 0.4.0)

### Suggested (for full functionality)
| Package | Source | Used by |
|---------|--------|---------|
| `scDblFinder` | Bioconductor | `Easy_DbFinder` |
| `BiocParallel` | Bioconductor | `Easy_DbFinder` |
| `harmony` | CRAN | `Easy_harmony` |
| `clustree` | CRAN | `Easy_clustree` |
| `scCDC` | Bioconductor | `Easy_CDC` |
| `decontX` | Bioconductor | `Easy_CDC` |
| `glmGamPoi` | Bioconductor | `Easy_Normal` (SCT speed-up) |
| `SCTransform` | CRAN / Seurat | `Easy_Normal` |

---

## Notes & Best Practices

1. **Order matters**: DbFinder -> QC -> Normal -> Harmony -> PC -> UMAP/TSNE -> clustree -> CDC.
2. **After `Easy_CDC`**: the count matrix has changed. Always re-run `Easy_Normal`, `Easy_harmony`, `RunUMAP`/`RunTSNE`, and `Easy_clustree`.
3. **`Easy_clustree` prerequisite**: `FindNeighbors` is called internally, but you must already have a Harmony (or PCA) reduction in the object.
4. **Parallelism**: `Easy_DbFinder` uses `BiocParallel::MulticoreParam(3)` by default; adjust if needed.
5. **Species-specific patterns**: `Easy_QC` automatically switches regex patterns for human (`^MT-`, `^RP[sl]`, `^HB[^(p)]`) vs. mouse (`^Mt-`, `^Rp[sl]`, `^Hb[^(p)]`).

---

## Citation

If you use EasySC in your research, please cite the underlying packages that make the analysis possible:

- **Seurat**: Hao et al., 2021 *Cell*
- **Harmony**: Korsunsky et al., 2019 *Nature Methods*
- **scDblFinder**: Germain et al., 2022 *F1000Research*
- **scCDC**: (refer to the scCDC publication)
- **decontX**: Yang et al., 2020 *Bioinformatics*

---

## License

GPL-3

## Contact

For bug reports and feature requests, please open an issue on GitHub.
