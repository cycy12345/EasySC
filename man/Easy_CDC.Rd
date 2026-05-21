\name{Easy_CDC}
\alias{Easy_CDC}
\title{Ambient RNA Decontamination with scCDC and decontX}
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
  \item{object}{A \code{Seurat} object.}

  \item{cluster_col}{Column name in \code{object@meta.data}
    containing cluster labels (e.g.
    \code{"RNA_snn_res.0.1"}).}

  \item{assay}{Assay to use for extracting count data when
    \pkg{decontX} is invoked. Default \code{"RNA"}.}

  \item{decontX_threshold}{Numeric threshold on
    \code{decontX} contamination score. Cells with score
    \code{>=} this value are removed. If \code{NULL}, no
    filtering is performed. Default \code{0.2}.}

  \item{verbose}{Logical. Print progress messages.
    Default \code{TRUE}.}
}
\value{
A \code{Seurat} object with contamination
corrected. The \code{Corrected} assay from \pkg{scCDC}
is copied to the \code{RNA} assay slot. If
\pkg{decontX} is used, a \code{contamination} column
is added to metadata.
}
\description{
Detects and corrects ambient-RNA contamination using
\pkg{scCDC}. If the contamination ratio is low
(\code{<= 0.0003}), the function additionally runs
\pkg{decontX} for deeper decontamination and optionally
filters cells by a contamination threshold.
}
\details{
Workflow:
\enumerate{
  \item Detect contaminant genes with
    \code{scCDC::ContaminationDetection}.
  \item Quantify contamination ratio with
    \code{scCDC::ContaminationQuantification}.
  \item Correct counts with
    \code{scCDC::ContaminationCorrection}.
  \item If contamination ratio \code{> 0.0003}, use
    scCDC result directly.
  \item If \code{<= 0.0003}, run \pkg{decontX} on the
    corrected counts and optionally filter cells.
}

\strong{Important:} After running \code{Easy_CDC}, you
must re-run normalisation, PCA/Harmony, UMAP/t-SNE, and
clustering because the count matrix has changed.
}
\examples{
\dontrun{
  sce_CDC <- Easy_CDC(sce_c, cluster_col = "RNA_snn_res.0.1")
  # Re-run downstream steps after decontamination
  sce_CDC <- Easy_Normal(sce_CDC)
  sce_CDC <- Easy_harmony(sce_CDC, group.by = "Patient")
}
}
\seealso{
\code{\link[scCDC]{ContaminationDetection}},
\code{\link[scCDC]{ContaminationCorrection}},
\code{\link[decontX]{decontX}}
}
