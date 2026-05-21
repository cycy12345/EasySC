\name{Easy_QC}
\alias{Easy_QC}
\title{Quality Control for scRNA-seq Data}
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
  \item{object}{A \code{Seurat} object. Must contain \code{nFeature_RNA}
    and \code{nCount_RNA} in its metadata.}

  \item{species}{Species identifier. Either \code{"human"} (default) or
    \code{"mouse"}. Determines the gene-name patterns used for MT, RP,
    and HB detection.}

  \item{MT}{Threshold for mitochondrial gene percentage (\\%). Cells with
    \code{percent.mt >= MT} are removed. Default is \code{5}.}

  \item{RP}{Threshold for ribosomal gene percentage (\\%). Cells with
    \code{percent.rp >= RP} are removed. Default is \code{5}.}

  \item{HB}{Threshold for haemoglobin gene percentage (\\%). Cells with
    \code{percent.hb >= HB} are removed. Default is \code{1}.}

  \item{Feature_high}{Upper bound for \code{nFeature_RNA}. Cells with
    \code{nFeature_RNA > Feature_high} are removed. Default \code{5000}.}

  \item{Feature_low}{Lower bound for \code{nFeature_RNA}. Cells with
    \code{nFeature_RNA < Feature_low} are removed. Default \code{200}.}

  \item{plot_group}{Column name in \code{object@meta.data} for grouping
    the violin plots (e.g. \code{"sample"}, \code{"Patient"}). If
    \code{NULL}, defaults to \code{"orig.ident"}.}

  \item{color}{A vector of colours for \code{ggpubr::ggviolin}.
    Default is a 9-colour pastel palette.}
}
\value{
A named \code{list} with three elements:
\itemize{
  \item \code{sce_filter}: The \code{Seurat} object after filtering.
  \item \code{plot_Count_Feature}: Violin plot of
    \code{nFeature_RNA} and \code{nCount_RNA}.
  \item \code{plot_percent}: Violin plot of \code{percent.mt},
    \code{percent.rp}, and \code{percent.hb}.
}
}
\description{
Computes mitochondrial (MT), ribosomal (RP), and haemoglobin (HB)
gene percentages; filters low-quality cells; and returns violin plots
for \code{nFeature_RNA}, \code{nCount_RNA}, and the three percentage
metrics.
}
\details{
Gene patterns:
\itemize{
  \item Human --- MT: \code{^MT-}, RP: \code{^RP[sl]}, HB: \code{^HB[^(p)]}
  \item Mouse --- MT: \code{^Mt-}, RP: \code{^Rp[sl]}, HB: \code{^Hb[^(p)]}
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
