\name{Easy_DbFinder}
\alias{Easy_DbFinder}
\title{Doublet Detection and Filtering with scDblFinder}
\usage{
Easy_DbFinder(object, dbr = 0.08, group_by = NULL)
}
\arguments{
  \item{object}{A \code{Seurat} object.}

  \item{dbr}{Expected doublet rate. Default is \code{0.08} (8\\%).
    Set to \code{NULL} to let \code{scDblFinder} estimate automatically.
    The rule of thumb for 10x data is ~1\\% per 1,000 captured cells
    (e.g., 4\\% for 4,000 cells).}

  \item{group_by}{Column name in \code{object@meta.data} indicating sample
    grouping. If provided, doublet detection is performed independently for
    each sample. If \code{NULL} (default), all cells are processed together.}
}
\value{
A named \code{list} with two elements:
\itemize{
  \item \code{sce_filter}: The input \code{Seurat} object after removing
    doublets (only \code{"singlet"} cells retained).
  \item \code{plot}: A \code{ggplot2} bar plot visualising singlet/doublet
    counts per sample (or overall when \code{group_by = NULL}).
}
}
\description{
Detects and removes doublets (multiplets) from a Seurat object using
\code{scDblFinder}. The function converts the Seurat object to
\code{SingleCellExperiment}, runs \code{scDblFinder}, and returns both
the filtered object and a bar plot showing singlet/doublet distribution.
}
\details{
\code{scDblFinder} is run with \code{BiocParallel::MulticoreParam(3)} for
moderate parallelisation. Progress messages are printed via \code{cli}.
}
\examples{
\dontrun{
  # Basic usage (all cells together)
  res <- Easy_DbFinder(sce, dbr = 0.08)
  res$plot
  sce_clean <- res$sce_filter

  # Per-sample doublet detection
  res <- Easy_DbFinder(sce, group_by = "Patient", dbr = 0.08)
  sce_clean <- res$sce_filter
}
}
\seealso{
\code{\link[scDblFinder]{scDblFinder}},
\code{\link[Seurat]{as.SingleCellExperiment}}
}
