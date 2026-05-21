\name{Easy_PC}
\alias{Easy_PC}
\title{Select Optimal Number of Principal Components}
\usage{
Easy_PC(seurat, reduction = "pca", cum = 90, var = 5)
}
\arguments{
  \item{seurat}{A \code{Seurat} object containing a PCA
    reduction.}

  \item{reduction}{Name of the reduction to inspect.
    Default \code{"pca"}.}

  \item{cum}{Cumulative variance threshold (\\%). The first
    PC where cumulative variance exceeds this value is
    candidate 1. Default \code{90}.}

  \item{var}{Individual variance threshold (\\%). Only PCs
    with contribution below this value are considered for
    \code{co1}. Default \code{5}.}
}
\value{
An integer: the recommended number of PCs to use.
A scatter plot is printed to the active graphics device.
}
\description{
Determines the optimal number of principal components
(PCs) to retain by combining two criteria:
cumulative variance > \code{cum}\\% and individual
variance contribution < \code{var}\\%; plus a minimum
drop-off criterion (difference between consecutive PCs
> 0.1\\%).
}
\details{
The two criteria are:
\enumerate{
  \item \code{co1}: first PC where cumulative percentage
    > \code{cum} AND individual percentage < \code{var}.
  \item \code{co2}: last PC where the drop from the
    previous PC > 0.1\\%.
}
The returned value is \code{min(co1, co2)}.
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
