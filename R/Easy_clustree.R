#' Multi-Resolution Clustering and Cluster Tree Visualisation
#'
#' @description
#' Builds a shared-nearest-neighbour (SNN) graph (if not
#' already present), runs \code{Seurat::FindClusters} at
#' multiple resolutions, and draws a cluster tree with
#' \pkg{clustree} to help choose the best resolution.
#'
#' @param object A \code{Seurat} object. Must already
#'   contain a Harmony (or PCA) reduction.
#' @param assay Assay used to name the SNN graph.
#'   One of \code{"RNA"} or \code{"SCT"}.
#'   Default \code{c("RNA", "SCT")} (the first match is used).
#' @param resolutions Numeric vector of clustering
#'   resolutions to test. Default
#'   \code{c(0.01, 0.05, 0.1, 0.2, 0.3, 0.5, 0.8, 1)}.
#' @param algorithm Clustering algorithm passed to
#'   \code{FindClusters}. Default \code{1} (Louvain).
#' @param PC Number of principal components (Harmony dims)
#'   to use for \code{FindNeighbors}. Typically the output
#'   of \code{\link{Easy_PC}}.
#' @param graph.name Name of the SNN graph. If \code{NULL},
#'   inferred from \code{assay} as \code{"{assay}_snn"}.
#' @param prefix Prefix for cluster columns in metadata.
#'   If \code{NULL}, defaults to \code{"{graph.name}_res."}.
#'
#' @return A named \code{list}:
#'   \itemize{
#'     \item \code{object}: the \code{Seurat} object with
#'       multiple \code{seurat_clusters} resolutions stored
#'       in metadata.
#'     \item \code{plot}: a \code{clustree} plot
#'       (\code{ggplot} object).
#'   }
#'
#' @details
#' \code{FindNeighbors} is run using the \code{"harmony"}
#' reduction and dimensions \code{1:PC}. If the SNN graph
#' already exists, the function skips rebuilding it, but
#' the check is performed after \code{FindNeighbors} to
#' ensure compatibility.
#'
#' @seealso \code{\link[Seurat]{FindNeighbors}},
#'   \code{\link[Seurat]{FindClusters}},
#'   \code{\link[clustree]{clustree}}
#'
#' @examples
#' \dontrun{
#'   sce <- Easy_harmony(sce, group.by = "Patient")
#'   PC <- Easy_PC(sce)
#'   cl <- Easy_clustree(sce, assay = "RNA", PC = PC)
#'   cl$plot
#'   sce <- cl$object
#' }
#'
#' @importFrom Seurat FindNeighbors FindClusters
#' @importFrom cli cli_progress_step cli_progress_done
#' @export
Easy_clustree <- function(object,
                          assay = c("RNA", "SCT"),
                          resolutions = c(0.01, 0.05, 0.1, 0.2, 0.3, 0.5, 0.8, 1),
                          algorithm = 1,
                          PC = 20,
                          graph.name = NULL,
                          prefix = NULL) {

  # Step 1: validate clustree dependency
  cli::cli_progress_step("Validating clustree dependency")
  if (!requireNamespace("clustree", quietly = TRUE)) {
    stop("Please install clustree: install.packages('clustree')")
  }

  # Step 2: determine graph name and prefix
  cli::cli_progress_step("Determining graph name and prefix")
  if (is.null(graph.name)) {
    assay <- match.arg(assay)
    graph.name <- paste0(assay, "_snn")
    if (is.null(prefix)) {
      prefix <- paste0(graph.name, "_res.")
    }
  } else {
    if (is.null(prefix)) {
      prefix <- paste0(graph.name, "_res.")
    }
  }

  # Step 2.5: build SNN graph
  cli::cli_progress_step("Building KNN / SNN graph")
  object <- Seurat::FindNeighbors(
    object,
    reduction = "harmony",
    dims = 1:PC,
    verbose = FALSE
  )

  # Step 3: verify graph exists
  cli::cli_progress_step("Checking SNN graph")
  if (!graph.name %in% names(object@graphs)) {
    stop(
      "SNN graph '", graph.name,
      "' not found. Please run FindNeighbors with the correct assay or graph.name."
    )
  }

  # Step 4: multi-resolution clustering
  cli::cli_progress_step("Multi-resolution clustering")
  for (res in resolutions) {
    message("  Clustering at resolution = ", res)
    object <- Seurat::FindClusters(
      object,
      resolution = res,
      verbose = TRUE
    )
  }

  # Step 5: draw cluster tree
  cli::cli_progress_step("Drawing cluster tree")
  p <- clustree::clustree(object@meta.data, prefix = prefix)

  # Step 6: finish
  cli::cli_progress_step("Returning results")
  cli::cli_progress_done()
  return(list(object = object, plot = p))
}
