\name{Easy_DbFinder}
\alias{Easy_DbFinder}
\title{使用 scDblFinder 检测并过滤双细胞}
\usage{
Easy_DbFinder(object, dbr = 0.08, group_by = NULL)
}
\arguments{
  \item{object}{一个 \code{Seurat} 对象。}

  \item{dbr}{预期的双细胞率。默认值为 \code{0.08}（8\%）。
    设置为 \code{NULL} 可让 \code{scDblFinder} 自动估计。
    10x 数据的估算规则为：每捕获 1,000 个细胞约 1\%
    （例如 4,000 个细胞对应 4\%）。}

  \item{group_by}{\code{object@meta.data} 中指示样本分组的列名。
    如果提供，则按样本独立进行双细胞检测。
    如果为 \code{NULL}（默认），则所有细胞一起处理。}
}
\value{
一个命名的 \code{list}，包含两个元素：
\itemize{
  \item \code{sce_filter}: 去除双细胞后的 \code{Seurat} 对象
    （仅保留 \code{"singlet"} 细胞）。
  \item \code{plot}: 一个 \code{ggplot2} 柱状图，可视化每个样本
    （或总体，当 \code{group_by = NULL} 时）的单细胞/双细胞分布。
}
}
\description{
使用 \code{scDblFinder} 检测并去除 Seurat 对象中的双细胞（multiplets）。
该函数先将 Seurat 对象转换为 \code{SingleCellExperiment}，
运行 \code{scDblFinder}，然后返回过滤后的对象以及
单细胞/双细胞分布柱状图。
}
\details{
\code{scDblFinder} 使用 \code{BiocParallel::MulticoreParam(3)} 进行
中等程度的并行计算。进度信息通过 \code{cli} 输出。
}
\examples{
\dontrun{
  # 基本用法（所有细胞一起处理）
  res <- Easy_DbFinder(sce, dbr = 0.08)
  res$plot
  sce_clean <- res$sce_filter

  # 按样本独立检测双细胞
  res <- Easy_DbFinder(sce, group_by = "Patient", dbr = 0.08)
  sce_clean <- res$sce_filter
}
}
\seealso{
\code{\link[scDblFinder]{scDblFinder}},
\code{\link[Seurat]{as.SingleCellExperiment}}
}
