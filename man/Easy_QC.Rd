\name{Easy_QC}
\alias{Easy_QC}
\title{单细胞 RNA 测序数据质控}
\usage{
Easy_QC(
  object,
  species = "human",
  MT = 5,
  RP = 5,
  HB = 1,
  Feature_high = 5000,
  Feature_low = 200,
  plot_group = "sample",
  color = c("#fbb4ae", "#b3cde3", "#ccebc5", "#decbe4",
            "#fed9a6", "#ffffcc", "#e5d8bd", "#fddaec",
            "#f2f2f2")
)
}
\arguments{
  \item{object}{一个 \code{Seurat} 对象。其元数据中必须包含
    \code{nFeature_RNA} 和 \code{nCount_RNA}。}

  \item{species}{物种标识。可选 \code{"human"}（默认）或
    \code{"mouse"}。决定用于检测 MT、RP 和 HB 的基因名称匹配模式。}

  \item{MT}{线粒体基因百分比阈值（\%）。\code{percent.mt >= MT} 的
    细胞将被去除。默认值为 \code{5}。}

  \item{RP}{核糖体基因百分比阈值（\%）。\code{percent.rp >= RP} 的
    细胞将被去除。默认值为 \code{5}。}

  \item{HB}{血红蛋白基因百分比阈值（\%）。\code{percent.hb >= HB} 的
    细胞将被去除。默认值为 \code{1}。}

  \item{Feature_high}{\code{nFeature_RNA} 的上限。\code{nFeature_RNA > Feature_high}
    的细胞将被去除。默认 \code{5000}。}

  \item{Feature_low}{\code{nFeature_RNA} 的下限。\code{nFeature_RNA < Feature_low}
    的细胞将被去除。默认 \code{200}。}

  \item{plot_group}{\code{object@meta.data} 中用于分组绘制小提琴图的
    列名（例如 \code{"sample"}、\code{"Patient"}）。
    如果为 \code{NULL}，则默认为 \code{"orig.ident"}。}

  \item{color}{用于 \code{ggpubr::ggviolin} 的颜色向量。
    默认为 9 色柔和配色板。}
}
\value{
一个命名的 \code{list}，包含三个元素：
\itemize{
  \item \code{sce_filter}: 过滤后的 \code{Seurat} 对象。
  \item \code{plot_Count_Feature}: \code{nFeature_RNA} 和
    \code{nCount_RNA} 的小提琴图。
  \item \code{plot_percent}: \code{percent.mt}、
    \code{percent.rp} 和 \code{percent.hb} 的小提琴图。
}
}
\description{
计算线粒体（MT）、核糖体（RP）和血红蛋白（HB）基因百分比；
过滤低质量细胞；并返回 \code{nFeature_RNA}、\code{nCount_RNA}
以及上述三个百分比指标的小提琴图。
}
\details{
基因匹配模式：
\itemize{
  \item 人类 --- MT: \code{^MT-}, RP: \code{^RP[sl]}, HB: \code{^HB[^(p)]}
  \item 小鼠 --- MT: \code{^Mt-}, RP: \code{^Rp[sl]}, HB: \code{^Hb[^(p)]}
}
}
\examples{
\dontrun{
  sce_qc <- Easy_QC(sce, species = "human", MT = 5, RP = 5, HB = 1,
                    Feature_high = 5000, Feature_low = 200,
                    plot_group = "Patient")
  sce_qc$plot_Count_Feature
  sce_qc$plot_percent
  sce_clean <- sce_qc$sce_filter
}
}
\seealso{
\code{\link[Seurat]{PercentageFeatureSet}},
\code{\link[Seurat]{FetchData}}, \code{\link[ggpubr]{ggviolin}}
}
