\name{Easy_PC}
\alias{Easy_PC}
\title{选择最优主成分数量}
\usage{
Easy_PC(seurat, reduction = "pca", cum = 90, var = 5)
}
\arguments{
  \item{seurat}{一个包含 PCA 降维结果的 \code{Seurat} 对象。}

  \item{reduction}{要检查的降维方法名称。
    默认 \code{"pca"}。}

  \item{cum}{累积方差阈值（\%）。累积方差超过此值的
    第一个主成分作为候选 1。默认 \code{90}。}

  \item{var}{单个主成分方差贡献阈值（\%）。只有贡献低于
    此值的主成分才会被纳入 \code{co1} 的计算。
    默认 \code{5}。}
}
\value{
一个整数：推荐使用的主成分数量。
同时会在当前图形设备上打印一张带注释的散点图。
}
\description{
通过结合两个标准来确定要保留的最优主成分（PC）数量：
累积方差 > \code{cum}\% 且单个方差贡献 < \code{var}\%；
再加上最小衰减标准（相邻 PC 之间的差值 > 0.1\%）。
}
\details{
两个判定标准：
\enumerate{
  \item \code{co1}: 累积百分比 > \code{cum} 且
    单个百分比 < \code{var} 的第一个 PC。
  \item \code{co2}: 与前一个 PC 的差值 > 0.1\%
    的最后一个 PC。
}
返回值为 \code{min(co1, co2)}。
}
\examples{
\dontrun{
  sce <- Seurat::RunPCA(sce, features = VariableFeatures(sce))
  npc <- Easy_PC(sce, reduction = "pca", cum = 90, var = 5)
  npc
}
}
\seealso{
\code{\link[Seurat]{RunPCA}}
}
