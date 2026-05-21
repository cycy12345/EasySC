\name{Easy_harmony}
\alias{Easy_harmony}
\title{使用 Harmony 进行批次校正}
\usage{
Easy_harmony(
  object,
  group.by,
  features = Seurat::VariableFeatures(object),
  npcs = 50,
  verbose = TRUE,
  ...
)
}
\arguments{
  \item{object}{一个 \code{Seurat} 对象。应已包含高变基因
    （例如在运行 \code{\link{Easy_Normal}} 之后）。}

  \item{group.by}{\code{object@meta.data} 中指示批次变量的
    列名（例如 \code{"Patient"}、\code{"sample"}）。}

  \item{features}{用于 PCA 的基因名称向量。
    默认为 \code{VariableFeatures(object)}。}

  \item{npcs}{要计算的主成分数量。
    默认 \code{50}。}

  \item{verbose}{逻辑值。是否打印进度信息。
    默认 \code{TRUE}。}

  \item{...}{额外参数传递给 \code{harmony::RunHarmony}。}
}
\value{
输入的 \code{Seurat} 对象，PCA 结果存储在
\code{reductions$pca} 中，Harmony 校正后的嵌入存储在
\code{reductions$harmony} 中。
}
\description{
对高变基因运行 PCA，然后进行 Harmony 批次校正。
该函数是 \code{Seurat::RunPCA} 和 \code{harmony::RunHarmony}
的便捷封装。
}
\details{
该函数首先通过 \code{\link{Easy_PC}} 确定最优主成分数，
然后在 PCs \code{1:optimal_PC} 上运行 Harmony。
}
\examples{
\dontrun{
  sce_norm <- Easy_Normal(sce)
  sce_norm <- Easy_harmony(sce_norm, group.by = "Patient")
}
}
\seealso{
\code{\link[Seurat]{RunPCA}},
\code{\link[harmony]{RunHarmony}},
\code{\link{Easy_PC}}
}
