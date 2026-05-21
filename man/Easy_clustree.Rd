\name{Easy_clustree}
\alias{Easy_clustree}
\title{Multi-Resolution Clustering and Cluster Tree Visualisation}
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
  \item{object}{A \code{Seurat} object. Must already
    contain a Harmony (or PCA) reduction.}

  \item{assay}{Assay used to name the SNN graph.
    One of \code{"RNA"} or \code{"SCT"}.
    Default \code{c("RNA", "SCT")} (the first match is used).}

  \item{resolutions}{Numeric vector of clustering
    resolutions to test. Default
    \code{c(0.01, 0.05, 0.1, 0.2, 0.3, 0.5, 0.8, 1)}.}

  \item{algorithm}{Clustering algorithm passed to
    \code{FindClusters}. Default \code{1} (Louvain).}

  \item{PC}{Number of principal components (Harmony dims)
    to use for \code{FindNeighbors}. Typically the output
    of \code{\link{Easy_PC}}.}

  \item{graph.name}{Name of the SNN graph. If \code{NULL},
    inferred from \code{assay} as \code{"{assay}_snn"}.}

  \item{prefix}{Prefix for cluster columns in metadata.
    If \code{NULL}, defaults to \code{"{graph.name}_res."}.}
}
\value{
A named \code{list}:
\itemize{
  \item \code{object}: the \code{Seurat} object with
    multiple \code{seurat_clusters} resolutions stored
    in metadata.
  \item \code{plot}: a \code{clustree} plot
    (\code{ggplot} object).
}
}
\description{
Builds a shared-nearest-neighbour (SNN) graph (if not
already present), runs \code{Seurat::FindClusters} at
multiple resolutions, and draws a cluster tree with
\pkg{clustree} to help choose the best resolution.
}
\details{
\code{FindNeighbors} is run using the \code{"harmony"}
reduction and dimensions \code{1:PC}. If the SNN graph
already exists, the function skips rebuilding it, but
the check is performed after \code{FindNeighbors} to
ensure compatibility.
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
