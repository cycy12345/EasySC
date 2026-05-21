\name{Easy_Normal}
\alias{Easy_Normal}
\title{Normalisation, Variable Feature Selection, and Scaling}
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
  \item{object}{A \code{Seurat} object.}

  \item{use_sct}{Logical. If \code{TRUE}, additionally
    runs \code{SCTransform} and stores the result in the
    \code{SCT} assay. Default is \code{FALSE}.}

  \item{vars.to.regress}{Variables to regress out during
    \code{ScaleData} and/or \code{SCTransform} (e.g.
    \code{c("percent.mt", "nCount_RNA")}). Default is
    \code{NULL}.}

  \item{scale.factor}{Scale factor for
    \code{NormalizeData}. Default \code{10000} (CP10K).}

  \item{nfeatures}{Number of highly variable features to
    select. Default \code{2000}.}

  \item{assay}{Assay to use. Default \code{"RNA"}.}

  \item{verbose}{Logical. Print detailed messages.
    Default \code{TRUE}.}

  \item{...}{Additional arguments passed to
    \code{NormalizeData}, \code{FindVariableFeatures},
    \code{ScaleData}, or \code{SCTransform}.}
}
\value{
The input \code{Seurat} object with
normalised data, scaled data, and variable features
stored in the specified assay. If \code{use_sct = TRUE},
an additional \code{SCT} assay is created and set as
default.
}
\description{
Performs standard Seurat normalisation
(\code{LogNormalize}), identifies highly variable
features, and scales the data. Optionally runs
\code{SCTransform} for variance-stabilising
transformation.
}
\details{
When \code{use_sct = TRUE}, the function prints
downstream-analysis recommendations, e.g.:
\itemize{
  \item Differential expression --- use \code{RNA} assay,
    \code{data} slot.
  \item Visualisation / dimensionality reduction --- use
    \code{SCT} assay, \code{scale.data} slot.
  \item Cell-cell communication --- use \code{SCT} assay,
    \code{data} slot.
}
If the \pkg{glmGamPoi} package is installed, it is
automatically leveraged by \code{SCTransform} for speed.
}
\examples{
\dontrun{
  # Standard LogNormalize workflow
  sce_norm <- Easy_Normal(sce, use_sct = FALSE)

  # With SCTransform and mitochondrial regression
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
