\name{Easy_harmony}
\alias{Easy_harmony}
\title{Batch Correction with Harmony}
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
  \item{object}{A \code{Seurat} object. Should already
    contain variable features (e.g. after
    \code{\link{Easy_Normal}}).}

  \item{group.by}{Column name in \code{object@meta.data}
    indicating the batch variable (e.g. \code{"Patient"},
    \code{"sample"}).}

  \item{features}{Vector of feature names to use for PCA.
    Defaults to \code{VariableFeatures(object)}.}

  \item{npcs}{Number of principal components to compute.
    Default \code{50}.}

  \item{verbose}{Logical. Print progress messages.
    Default \code{TRUE}.}

  \item{...}{Additional arguments passed to
    \code{harmony::RunHarmony}.}
}
\value{
The input \code{Seurat} object with PCA stored
in \code{reductions$pca} and Harmony-corrected
embeddings stored in \code{reductions$harmony}.
}
\description{
Runs PCA on variable features followed by Harmony batch
correction. This function is a convenience wrapper around
\code{Seurat::RunPCA} and \code{harmony::RunHarmony}.
}
\details{
The function first determines the optimal number of PCs
via \code{\link{Easy_PC}}, then runs Harmony on PCs
\code{1:optimal_PC}.
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
