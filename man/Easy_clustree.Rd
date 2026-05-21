\name{Easy_clustree}
\alias{Easy_clustree}
\title{多分辨率聚类与聚类树可视化}
\usage{
Easy_clustree(
  object,
  assay = c("RNA", "SCT"),
  resolutions = c(0.01, 0.05, 0.1, 0.2, 0.3, 0.5, 0.8, 1),
  algorithm = 1,
  PC = 20,
  graph.name = NULL,
  prefix = NULL
)
}
\arguments{
  \item{object}{一个 \code{Seurat} 对象。必须已包含
    Harmony（或 PCA）降维结果。}

  \item{assay}{用于命名 SNN 图的 assay。
    可选 \code{"RNA"} 或 \code{"SCT"}。
    默认 \code{c("RNA", "SCT")}（取第一个匹配项）。}

  \item{resolutions}{要测试的聚类分辨率数值向量。
    默认 \code{c(0.01, 0.05, 0.1, 0.2, 0.3, 0.5, 0.8, 1)}。}

  \item{algorithm}{传递给 \code{FindClusters} 的聚类算法。
    默认 \code{1}（Louvain 算法）。}

  \item{PC}{用于 \code{FindNeighbors} 的主成分数量
    （Harmony 维度）。通常为 \code{\link{Easy_PC}} 的输出。}

  \item{graph.name}{SNN 图的名称。若为 \code{NULL}，
    则从 \code{assay} 推断为 \code{"{assay}_snn"}。}

  \item{prefix}{元数据中聚类列名的前缀。
    若为 \code{NULL}，则默认为 \code{"{graph.name}_res."}。}
}
\value{
一个命名的 \code{list}：
\itemize{
  \item \code{object}: 包含多分辨率聚类结果的
    \code{Seurat} 对象（结果存储在元数据中）。
  \item \code{plot}: 一张 \code{clustree} 聚类树图
    （\code{ggplot} 对象）。
}
}
\description{
构建共享最近邻（SNN）图（若尚不存在），在多个分辨率下
运行 \code{Seurat::FindClusters}，并使用 \pkg{clustree}
绘制聚类树以帮助选择最佳分辨率。
}
\details{
\code{FindNeighbors} 使用 \code{"harmony"} 降维和
维度 \code{1:PC} 运行。如果 SNN 图已存在，函数会跳过
重建，但会在 \code{FindNeighbors} 之后进行检查以确保兼容。
}
\examples{
\dontrun{
  sce <- Easy_harmony(sce, group.by = "Patient")
  PC <- Easy_PC(sce)
  cl <- Easy_clustree(sce, assay = "RNA", PC = PC)
  cl$plot
  sce <- cl$object
}
}
\seealso{
\code{\link[Seurat]{FindNeighbors}},
\code{\link[Seurat]{FindClusters}},
\code{\link[clustree]{clustree}}
}
