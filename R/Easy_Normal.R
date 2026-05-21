#' Normalisation, Variable Feature Selection, and Scaling
#'
#' @description
#' Performs standard Seurat normalisation
#' (\code{LogNormalize}), identifies highly variable
#' features, and scales the data. Optionally runs
#' \code{SCTransform} for variance-stabilising
#' transformation.
#'
#' @param object A \code{Seurat} object.
#' @param use_sct Logical. If \code{TRUE}, additionally
#'   runs \code{SCTransform} and stores the result in the
#'   \code{SCT} assay. Default is \code{FALSE}.
#' @param vars.to.regress Variables to regress out during
#'   \code{ScaleData} and/or \code{SCTransform} (e.g.
#'   \code{c("percent.mt", "nCount_RNA")}). Default is
#'   \code{NULL}.
#' @param scale.factor Scale factor for
#'   \code{NormalizeData}. Default \code{10000} (CP10K).
#' @param nfeatures Number of highly variable features to
#'   select. Default \code{2000}.
#' @param assay Assay to use. Default \code{"RNA"}.
#' @param verbose Logical. Print detailed messages.
#'   Default \code{TRUE}.
#' @param ... Additional arguments passed to
#'   \code{NormalizeData}, \code{FindVariableFeatures},
#'   \code{ScaleData}, or \code{SCTransform}.
#'
#' @return The input \code{Seurat} object with
#'   normalised data, scaled data, and variable features
#'   stored in the specified assay. If \code{use_sct = TRUE},
#'   an additional \code{SCT} assay is created and set as
#'   default.
#'
#' @details
#' When \code{use_sct = TRUE}, the function prints
#' downstream-analysis recommendations, e.g.:
#' \itemize{
#'   \item Differential expression — use \code{RNA} assay,
#'     \code{data} slot.
#'   \item Visualisation / dimensionality reduction — use
#'     \code{SCT} assay, \code{scale.data} slot.
#'   \item Cell-cell communication — use \code{SCT} assay,
#'     \code{data} slot.
#' }
#' If the \pkg{glmGamPoi} package is installed, it is
#' automatically leveraged by \code{SCTransform} for speed.
#'
#' @seealso \code{\link[Seurat]{NormalizeData}},
#'   \code{\link[Seurat]{FindVariableFeatures}},
#'   \code{\link[Seurat]{ScaleData}},
#'   \code{\link[Seurat]{SCTransform}}
#'
#' @examples
#' \dontrun{
#'   # Standard LogNormalize workflow
#'   sce_norm <- Easy_Normal(sce, use_sct = FALSE)
#'
#'   # With SCTransform and mitochondrial regression
#'   sce_norm <- Easy_Normal(sce, use_sct = TRUE,
#'                           vars.to.regress = "percent.mt")
#' }
#'
#' @importFrom Seurat NormalizeData FindVariableFeatures
#'   ScaleData DefaultAssay VariableFeatures
#' @importFrom cli cli_progress_step cli_progress_done
#' @export
Easy_Normal <- function(object,
                        use_sct = FALSE,
                        vars.to.regress = NULL,
                        scale.factor = 10000,
                        nfeatures = 2000,
                        assay = "RNA",
                        verbose = TRUE,
                        ...) {

  # Step 1: argument checks
  cli::cli_progress_step("Validating input arguments")
  if (!inherits(object, "Seurat")) {
    stop("object must be a Seurat object.")
  }
  if (!assay %in% names(object)) {
    stop("Specified assay '", assay, "' does not exist.")
  }

  # Step 2: normalisation
  cli::cli_progress_step("Running LogNormalize")
  if (verbose) message("Running Seurat standard normalisation (LogNormalize)...")
  object <- Seurat::NormalizeData(
    object,
    normalization.method = "LogNormalize",
    scale.factor = scale.factor,
    assay = assay,
    verbose = verbose,
    ...
  )

  # Step 3: highly variable features
  cli::cli_progress_step("Finding highly variable features")
  object <- Seurat::FindVariableFeatures(
    object,
    selection.method = "vst",
    nfeatures = nfeatures,
    assay = assay,
    verbose = verbose
  )

  # Step 4: scaling
  cli::cli_progress_step("Running ScaleData")
  object <- Seurat::ScaleData(
    object,
    vars.to.regress = vars.to.regress,
    assay = assay,
    verbose = verbose,
    ...
  )

  # Step 5: set default assay
  cli::cli_progress_step("Setting default assay")
  Seurat::DefaultAssay(object) <- assay

  if (verbose) {
    message("\n========== LogNormalize complete ==========")
    message("Normalised data stored in: ", assay, "$data")
    message("Scaled data stored in: ", assay, "$scale.data")
    message("Variable features: ", length(Seurat::VariableFeatures(object, assay = assay)))
    message("===========================================\n")
  }

  # Step 6: optional SCTransform
  if (use_sct) {
    cli::cli_progress_step("Running SCTransform")
    if (verbose) message("Additional SCTransform (stored in SCT assay)...")

    if (requireNamespace("glmGamPoi", quietly = TRUE) && verbose) {
      message("glmGamPoi detected; SCTransform will use it automatically.")
    }

    object <- Seurat::SCTransform(
      object,
      assay = assay,
      vars.to.regress = vars.to.regress,
      verbose = verbose,
      ...
    )

    if (verbose) {
      message("\n========== SCTransform complete ==========")
      message("Corrected data stored in: SCT$data (Pearson residuals)")
      message("SCT assay is now default; switch with DefaultAssay(object) <- 'RNA'")
      message("\n[Downstream recommendations]")
      message("Differential analysis: use RNA assay, data slot")
      message("  FindMarkers(object, assay = 'RNA', slot = 'data')")
      message("Visualisation / dim-red: use SCT assay, scale.data")
      message("  RunPCA(object, assay = 'SCT')")
      message("Cell communication: use SCT assay, data slot")
      message("  GetAssayData(object, assay = 'SCT', slot = 'data')")
      message("==========================================\n")
    }
  } else {
    if (verbose) {
      message("\n[Downstream recommendations]")
      message("Differential analysis: FindMarkers(object, assay = 'RNA', slot = 'data')")
      message("Visualisation / dim-red: based on RNA assay scale.data")
      message("  RunPCA(object, assay = 'RNA')")
      message("==========================================\n")
    }
  }

  cli::cli_progress_done()
  return(object)
}
