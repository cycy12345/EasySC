\name{Easy_Normal}
\alias{Easy_Normal}
\title{数据标准化、高变基因筛选与缩放}
\usage{
Easy_Normal(
  object,
  use_sct = FALSE,
  vars.to.regress = NULL,
  scale.factor = 10000,
  nfeatures = 2000,
  assay = "RNA",
  verbose = TRUE,
  ...
)
}
\arguments{
  \item{object}{一个 \code{Seurat} 对象。}

  \item{use_sct}{逻辑值。若为 \code{TRUE}，则额外运行
    \code{SCTransform} 并将结果存储在 \code{SCT} assay 中。
    默认值为 \code{FALSE}。}

  \item{vars.to.regress}{在 \code{ScaleData} 和/或
    \code{SCTransform} 中要回归去除的变量（例如
    \code{c("percent.mt", "nCount_RNA")}）。默认 \code{NULL}。}

  \item{scale.factor}{\code{NormalizeData} 的缩放因子。
    默认 \code{10000}（CP10K）。}

  \item{nfeatures}{要筛选的高变基因数量。默认 \code{2000}。}

  \item{assay}{使用的 assay 名称。默认 \code{"RNA"}。}

  \item{verbose}{逻辑值。是否打印详细信息。
    默认 \code{TRUE}。}

  \item{...}{额外参数传递给 \code{NormalizeData}、
    \code{FindVariableFeatures}、\code{ScaleData} 或
    \code{SCTransform}。}
}
\value{
输入的 \code{Seurat} 对象，在指定的 assay 中存储了标准化数据、
缩放数据和高变基因。若 \code{use_sct = TRUE}，则额外创建
\code{SCT} assay 并将其设为默认。
}
\description{
执行标准的 Seurat 标准化流程（\code{LogNormalize}）、
筛选高变基因并缩放数据。可选运行 \code{SCTransform}
进行方差稳定变换。
}
\details{
当 \code{use_sct = TRUE} 时，函数会打印下游分析建议：
\itemize{
  \item 差异表达分析 --- 使用 \code{RNA} assay 的
    \code{data} slot。
  \item 可视化 / 降维 --- 使用 \code{SCT} assay 的
    \code{scale.data} slot。
  \item 细胞通讯分析 --- 使用 \code{SCT} assay 的
    \code{data} slot。
}
如果安装了 \pkg{glmGamPoi} 包，\code{SCTransform} 会自动
利用它来加速计算。
}
\examples{
\dontrun{
  # 标准 LogNormalize 流程
  sce_norm <- Easy_Normal(sce, use_sct = FALSE)

  # 启用 SCTransform 并回归线粒体基因比例
  sce_norm <- Easy_Normal(sce, use_sct = TRUE,
                          vars.to.regress = "percent.mt")
}
}
\seealso{
\code{\link[Seurat]{NormalizeData}},
\code{\link[Seurat]{FindVariableFeatures}},
\code{\link[Seurat]{ScaleData}},
\code{\link[Seurat]{SCTransform}}
}
