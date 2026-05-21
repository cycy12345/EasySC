\name{Easy_CDC}
\alias{Easy_CDC}
\title{使用 scCDC 和 decontX 去除环境 RNA 污染}
\usage{
Easy_CDC(
  object,
  cluster_col,
  assay = "RNA",
  decontX_threshold = 0.2,
  verbose = TRUE
)
}
\arguments{
  \item{object}{一个 \code{Seurat} 对象。}

  \item{cluster_col}{\code{object@meta.data} 中包含聚类标签的
    列名（例如 \code{"RNA_snn_res.0.1"}）。}

  \item{assay}{调用 \pkg{decontX} 时用于提取计数数据的
    assay 名称。默认 \code{"RNA"}。}

  \item{decontX_threshold}{\code{decontX} 污染分数的数值阈值。
    污染分数 \code{>=} 此值的细胞将被去除。
    若为 \code{NULL}，则跳过过滤。默认 \code{0.2}。}

  \item{verbose}{逻辑值。是否打印进度信息。
    默认 \code{TRUE}。}
}
\value{
污染校正后的 \code{Seurat} 对象。
\pkg{scCDC} 的 \code{Corrected} assay 被复制到
\code{RNA} assay slot 中。如果使用了 \pkg{decontX}，
元数据中会新增一列 \code{contamination}。
}
\description{
使用 \pkg{scCDC} 检测并校正环境 RNA 污染。
如果污染率较低（\code{<= 0.0003}），则额外运行
\pkg{decontX} 进行深度去污染，并可选择按污染阈值
过滤细胞。
}
\details{
工作流程：
\enumerate{
  \item 使用 \code{scCDC::ContaminationDetection} 检测污染基因。
  \item 使用 \code{scCDC::ContaminationQuantification} 量化污染率。
  \item 使用 \code{scCDC::ContaminationCorrection} 校正计数矩阵。
  \item 若污染率 \code{> 0.0003}，直接使用 scCDC 结果。
  \item 若 \code{<= 0.0003}，在 scCDC 校正后的计数上运行
    \pkg{decontX}，并可选择过滤细胞。
}

\strong{重要提示：} 运行 \code{Easy_CDC} 后，由于计数矩阵
已发生变化，必须重新运行标准化、PCA/Harmony、
UMAP/t-SNE 和聚类。
}
\examples{
\dontrun{
  sce_CDC <- Easy_CDC(sce_c, cluster_col = "RNA_snn_res.0.1")
  # 去污染后重新运行下游分析步骤
  sce_CDC <- Easy_Normal(sce_CDC)
  sce_CDC <- Easy_harmony(sce_CDC, group.by = "Patient")
}
}
\seealso{
\code{\link[scCDC]{ContaminationDetection}},
\code{\link[scCDC]{ContaminationCorrection}},
\code{\link[decontX]{decontX}}
}
