#' Batch Correction with Harmony
#'
#' @description
#' Runs PCA on variable features followed by Harmony batch
#' correction. This function is a convenience wrapper around
#' \code{Seurat::RunPCA} and \code{harmony::RunHarmony}.
#'
#' @param object A \code{Seurat} object. Should already
#'   contain variable features (e.g. after
#'   \code{\link{Easy_Normal}}).
#' @param group.by Column name in \code{object@meta.data}
#'   indicating the batch variable (e.g. \code{"Patient"},
#'   \code{"sample"}).
#' @param features Vector of feature names to use for PCA.
#'   Defaults to \code{VariableFeatures(object)}.
#' @param npcs Number of principal components to compute.
#'   Default \code{50}.
#' @param verbose Logical. Print progress messages.
#'   Default \code{TRUE}.
#' @param ... Additional arguments passed to
#'   \code{harmony::RunHarmony}.
#'
#' @return The input \code{Seurat} object with PCA stored
#'   in \code{reductions$pca} and Harmony-corrected
#'   embeddings stored in \code{reductions$harmony}.
#'
#' @details
#' The function first determines the optimal number of PCs
#' via \code{\link{Easy_PC}}, then runs Harmony on PCs
#' \code{1:optimal_PC}.
#'
#' @seealso \code{\link[Seurat]{RunPCA}},
#'   \code{\link[harmony]{RunHarmony}},
#'   \code{\link{Easy_PC}}
#'
#' @examples
#' \dontrun{
#'   sce_norm <- Easy_Normal(sce)
#'   sce_norm <- Easy_harmony(sce_norm, group.by = "Patient")
#' }
#'
#' @importFrom Seurat RunPCA
#' @importFrom cli cli_progress_step cli_progress_done
#' @export
Easy_harmony <- function(object,
                         group.by,
                         features = Seurat::VariableFeatures(object),
                         npcs = 50,
                         verbose = TRUE,
                         ...) {

  # Step 1: validate dependencies and arguments
  cli::cli_progress_step("Validating dependencies and batch variable")
  if (!requireNamespace("harmony", quietly = TRUE)) {
    stop("Please install the harmony package: install.packages('harmony')")
  }
  if (!group.by %in% colnames(object[[]])) {
    stop("Column '", group.by, "' not found in object metadata.")
  }

  # Step 2: PCA
  cli::cli_progress_step("Running PCA")
  if (verbose) message("Running PCA using ", length(features), " variable features...")
  object <- Seurat::RunPCA(
    object,
    features = features,
    npcs = npcs,
    verbose = verbose
  )

  # Step 2.5: determine optimal PC number
  PC <- Easy_PC(object)

  # Step 3: Harmony batch correction
  cli::cli_progress_step("Running Harmony batch correction")
  if (verbose) message("Running Harmony correction, batch variable: ", group.by)
  object <- harmony::RunHarmony(
    object,
    dims = 1:PC,
    group.by.vars = group.by
  )

  # Step 4: finish
  cli::cli_progress_step("Dimensionality reduction complete")
  if (verbose) message("Done! Harmony embeddings stored in 'harmony' reduction.")

  cli::cli_progress_done()
  return(object)
}
