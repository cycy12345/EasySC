# EasySC

**简易单细胞 RNA 测序分析流程**

EasySC 是一个简洁、易用的 R 包，将热门的 Bioconductor 和 CRAN 工具整合为一套完整的单细胞 RNA 测序（scRNA-seq）数据分析流程。它提供了一步式函数用于质量控制、标准化、批次校正、聚类、双细胞去除和环境 RNA 去污染，并通过进度条和信息提示让每一步都清晰透明、可复现。

---

## 功能概览

| 步骤 | 函数 | 功能 |
|------|------|------|
| 1 | `Easy_DbFinder` | 双细胞检测与过滤（可选） |
| 2 | `Easy_QC` | 计算 MT / RP / HB 比例；过滤低质量细胞 |
| 3 | `Easy_Normal` | 对数标准化、高变基因筛选、数据缩放；可选 SCTransform |
| 4 | `Easy_harmony` | PCA + Harmony 批次校正 |
| 5 | `Easy_PC` | 确定最优主成分数量 |
| 6 | `RunUMAP` / `RunTSNE` | 非线性降维（Seurat 原生函数） |
| 7 | `Easy_clustree` | 多分辨率聚类 + 聚类树可视化 |
| 8 | `Easy_CDC` | 环境 RNA 去污染（scCDC + decontX，可选） |

---

## 安装

### 1. 安装 CRAN 依赖

```r
install.packages(c("Seurat", "ggplot2", "cli", "tidyr", "dplyr",
                   "ggpubr", "clustree", "harmony"))
```

### 2. 安装 Bioconductor 依赖

```r
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install(c("scDblFinder", "BiocParallel", "scCDC",
                       "decontX", "glmGamPoi", "SCTransform"))
```

### 3. 安装 EasySC

```r
# 从 GitHub 安装（发布后可使用）
# devtools::install_github("cycy12345/EasySC")

# 或从本地源码安装
devtools::install("path/to/EasySC")
```

---

## 标准分析流程

以下是推荐的端到端分析流程。每个函数返回一个修改后的 `Seurat` 对象，或一个包含对象和图表的命名列表。

```r
library(EasySC)
library(Seurat)

# ------------------------------
# 1. 双细胞去除（可选）
# ------------------------------
res_db   <- Easy_DbFinder(sce, group_by = "Patient", dbr = 0.08)
res_db$plot
sce      <- res_db$sce_filter

# ------------------------------
# 2. 质量控制
# ------------------------------
res_qc   <- Easy_QC(sce, species = "human", MT = 5, RP = 5, HB = 1,
                    Feature_high = 5000, Feature_low = 200,
                    plot_group = "Patient")
res_qc$plot_Count_Feature
res_qc$plot_percent
sce      <- res_qc$sce_filter

# ------------------------------
# 3. 标准化
# ------------------------------
sce_norm <- Easy_Normal(sce, use_sct = FALSE)

# ------------------------------
# 4. 批次校正（Harmony）
# ------------------------------
sce_norm <- Easy_harmony(sce_norm, group.by = "Patient")

# ------------------------------
# 5. 选择最优 PC 数量
# ------------------------------
PC <- Easy_PC(sce_norm, reduction = "pca", cum = 90, var = 5)

# ------------------------------
# 6. 非线性降维（二选一）
# ------------------------------
sce_norm <- RunUMAP(sce_norm,  reduction = "harmony", dims = 1:PC)
# sce_norm <- RunTSNE(sce_norm, reduction = "harmony", dims = 1:PC)

# ------------------------------
# 7. 多分辨率聚类与聚类树
# ------------------------------
res_cl   <- Easy_clustree(sce_norm, assay = "RNA", PC = PC)
sce_cl   <- res_cl$object
clustree_plot <- res_cl$plot

# 查看聚类树并选择分辨率
colnames(sce_cl@meta.data)   # 例如 RNA_snn_res.0.1, RNA_snn_res.0.2, ...

# ------------------------------
# 8. 去污染（可选）
# ------------------------------
sce_CDC  <- Easy_CDC(sce_cl, cluster_col = "RNA_snn_res.0.1",
                     decontX_threshold = 0.2)

# !! 运行 Easy_CDC 后必须重新执行步骤 3-7 !!
sce_CDC  <- Easy_Normal(sce_CDC)
sce_CDC  <- Easy_harmony(sce_CDC, group.by = "Patient")
PC2      <- Easy_PC(sce_CDC)
sce_CDC  <- RunUMAP(sce_CDC, reduction = "harmony", dims = 1:PC2)
```

---

## 函数参考

### `Easy_DbFinder(object, dbr = 0.08, group_by = NULL)`
使用 **scDblFinder** 检测并去除双细胞。

- **dbr** -- 预期双细胞率（默认 8%）。10x 数据的估算规则：每 1,000 个细胞约 1%。
- **group_by** -- 元数据中的样本分组列（例如 `"Patient"`）。若为 `NULL`，则所有细胞合并处理。
- **返回**：`list(sce_filter = <Seurat>, plot = <ggplot>)`

---

### `Easy_QC(object, species = "human", MT = 5, RP = 5, HB = 1, ...)`
计算线粒体、核糖体和血红蛋白基因比例并过滤细胞。

- **species** -- `"human"` 或 `"mouse"`；自动设置基因名称匹配模式。
- **MT / RP / HB** -- 允许的最大百分比。
- **Feature_high / Feature_low** -- `nFeature_RNA` 的上下限。
- **plot_group** -- 用于分组绘制小提琴图的元数据列。
- **返回**：`list(sce_filter = <Seurat>, plot_Count_Feature = <ggplot>, plot_percent = <ggplot>)`

---

### `Easy_Normal(object, use_sct = FALSE, vars.to.regress = NULL, ...)`
标准 Seurat 标准化流程：`NormalizeData` -> `FindVariableFeatures` -> `ScaleData`。

- **use_sct** -- 若为 `TRUE`，额外运行 `SCTransform` 并将 `SCT` assay 设为默认。
- **vars.to.regress** -- 要回归去除的变量（例如 `c("percent.mt", "nCount_RNA")`）。
- **返回**：修改后的 `Seurat` 对象。

---

### `Easy_harmony(object, group.by, features = VariableFeatures(object), npcs = 50, ...)`
PCA 后进行 **Harmony** 批次校正。

- **group.by** -- 元数据中的批次变量（必需）。
- 内部调用 `Easy_PC()` 确定最优 PC 数量后再运行 Harmony。
- **返回**：包含 `reductions$pca` 和 `reductions$harmony` 的 `Seurat` 对象。

---

### `Easy_PC(seurat, reduction = "pca", cum = 90, var = 5)`
使用累积方差和衰减标准选择最优主成分数量。

- **cum** -- 累积方差阈值（%）。
- **var** -- 单个 PC 贡献阈值（%）。
- 在当前图形设备上打印带注释的散点图。
- **返回**：整数（推荐 PC 数量）。

---

### `Easy_clustree(object, assay = c("RNA", "SCT"), resolutions = ..., PC = 20, ...)`
使用 **clustree** 进行多分辨率聚类和聚类树可视化。

- **resolutions** -- 要测试的分辨率数值向量（默认 `c(0.01, 0.05, 0.1, 0.2, 0.3, 0.5, 0.8, 1)`）。
- **PC** -- 用于 `FindNeighbors` 的 Harmony/PCA 维度数量。
- **返回**：`list(object = <Seurat>, plot = <clustree ggplot>)`

---

### `Easy_CDC(object, cluster_col, assay = "RNA", decontX_threshold = 0.2, ...)`
使用 **scCDC**（以及可选的 **decontX**）进行环境 RNA 去污染。

- **cluster_col** -- 包含聚类标签的元数据列（例如 `"RNA_snn_res.0.1"`）。
- **decontX_threshold** -- `decontX` 污染分数 >= 此值的细胞将被去除。设为 `NULL` 则跳过过滤。
- 若污染率 > 0.0003，仅使用 scCDC 校正。
- 若 <= 0.0003，额外运行 decontX 进行深度清洁。
- **重要提示**：运行 `Easy_CDC` 后，必须重新执行标准化、Harmony、UMAP/t-SNE 和聚类。
- **返回**：去污染后的 `Seurat` 对象。

---

## 依赖包

### 必需（Imports）
- `Seurat` (>= 4.0.0)
- `ggplot2` (>= 3.3.0)
- `cli` (>= 3.0.0)
- `tidyr` (>= 1.2.0)
- `dplyr` (>= 1.0.0)
- `ggpubr` (>= 0.4.0)

### 建议安装（完整功能）
| 包名 | 来源 | 使用函数 |
|------|------|---------|
| `scDblFinder` | Bioconductor | `Easy_DbFinder` |
| `BiocParallel` | Bioconductor | `Easy_DbFinder` |
| `harmony` | CRAN | `Easy_harmony` |
| `clustree` | CRAN | `Easy_clustree` |
| `scCDC` | Bioconductor | `Easy_CDC` |
| `decontX` | Bioconductor | `Easy_CDC` |
| `glmGamPoi` | Bioconductor | `Easy_Normal`（SCT 加速） |
| `SCTransform` | CRAN / Seurat | `Easy_Normal` |

---

## 注意事项与最佳实践

1. **顺序很重要**：DbFinder -> QC -> Normal -> Harmony -> PC -> UMAP/TSNE -> clustree -> CDC。
2. **`Easy_CDC` 之后**：计数矩阵已改变，必须重新运行 `Easy_Normal`、`Easy_harmony`、`RunUMAP`/`RunTSNE` 和 `Easy_clustree`。
3. **`Easy_clustree` 前提**：函数内部会调用 `FindNeighbors`，但对象中必须已包含 Harmony（或 PCA）降维结果。
4. **并行计算**：`Easy_DbFinder` 默认使用 `BiocParallel::MulticoreParam(3)`；如需调整可修改源码。
5. **物种特异性匹配**：`Easy_QC` 自动切换人类（`^MT-`、`^RP[sl]`、`^HB[^(p)]`）与小鼠（`^Mt-`、`^Rp[sl]`、`^Hb[^(p)]`）的基因匹配模式。

---

## 引用

如果在研究中使用 EasySC，请同时引用其底层依赖包：

- **Seurat**：Hao et al., 2021 *Cell*
- **Harmony**：Korsunsky et al., 2019 *Nature Methods*
- **scDblFinder**：Germain et al., 2022 *F1000Research*
- **scCDC**：（请参见 scCDC 相关文献）
- **decontX**：Yang et al., 2020 *Bioinformatics*

---

## 许可证

GPL-3

## 联系方式

如有 bug 报告或功能建议，请在 GitHub 上提交 Issue。
