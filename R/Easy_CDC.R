#' Ambient RNA Decontamination with scCDC and decontX
#'
#' @description
#' Detects and corrects ambient-RNA contamination using
#' \pkg{scCDC}. If the contamination ratio is low
#' (\code{<= 0.0003}), the function additionally runs
#' \pkg{decontX} for deeper decontamination and optionally
#' filters cells by a contamination threshold.
#'
#' @param object A \code{Seurat} object.
#' @param cluster_col Column name in \code{object@meta.data}
#'   containing cluster labels (e.g.
#'   \code{"RNA_snn_res.0.1"}).
#' @param assay Assay to use for extracting count data when
#'   \pkg{decontX} is invoked. Default \code{"RNA"}.
#' @param decontX_threshold Numeric threshold on
#'   \code{decontX} contamination score. Cells with score
#'   \code{>=} this value are removed. If \code{NULL}, no
#'   filtering is performed. Default \code{0.2}.
#' @param verbose Logical. Print progress messages.
#'   Default \code{TRUE}.
#'
#' @return A \code{Seurat} object with contamination
#'   corrected. The \code{Corrected} assay from \pkg{scCDC}
#'   is copied to the \code{RNA} assay slot. If
#'   \pkg{decontX} is used, a \code{contamination} column
#'   is added to metadata.
#'
#' @details
#' Workflow:
#' \enumerate{
#'   \item Detect contaminant genes with
#'     \code{scCDC::ContaminationDetection}.
#'   \item Quantify contamination ratio with
#'     \code{scCDC::ContaminationQuantification}.
#'   \item Correct counts with
#'     \code{scCDC::ContaminationCorrection}.
#'   \item If contamination ratio \code{> 0.0003}, use
#'     scCDC result directly.
#'   \item If \code{<= 0.0003}, run \pkg{decontX} on the
#'     corrected counts and optionally filter cells.
#' }
#'
#' \strong{Important:} After running \code{Easy_CDC}, you
#' must re-run normalisation, PCA/Harmony, UMAP/t-SNE, and
#' clustering because the count matrix has changed.
#'
#' @seealso \code{\link[scCDC]{ContaminationDetection}},
#'   \code{\link[scCDC]{ContaminationCorrection}},
#'   \code{\link[decontX]{decontX}}
#'
#' @examples
#' \dontrun{
#'   sce_CDC <- Easy_CDC(sce_c, cluster_col = "RNA_snn_res.0.1")
#'   # Re-run downstream steps after decontamination
#'   sce_CDC <- Easy_Normal(sce_CDC)
#'   sce_CDC <- Easy_harmony(sce_CDC, group.by = "Patient")
#' }
#'
#' @importFrom Seurat GetAssayData
#' @importFrom cli cli_progress_step cli_progress_done
#' @export
Easy_CDC <- function(object, cluster_col, assay = "RNA",
                     decontX_threshold = 0.2, verbose = TRUE) {

  # Step 1: validate dependencies
  cli::cli_progress_step("Validating dependencies")
  if (!requireNamespace("scCDC", quietly = TRUE)) {
    stop("Please install the scCDC package.")
  }
  if (!requireNamespace("decontX", quietly = TRUE)) {
    stop("Please install the decontX package.")
  }

  # Step 2: detect contaminant genes
  cli::cli_progress_step("Detecting contaminant genes (scCDC)")
  Zc <- object[[cluster_col]]
  Seurat::Idents(object) <- Zc[[cluster_col]]
  GCGs <- scCDC::ContaminationDetection(object)

  # Step 3: quantify contamination ratio
  cli::cli_progress_step("Quantifying contamination ratio")
  contamination_ratio <- scCDC::ContaminationQuantification(object, rownames(GCGs))

  # Step 4: correct contamination
  cli::cli_progress_step("Correcting contamination (scCDC)")
  seuratobj_corrected <- scCDC::ContaminationCorrection(object, rownames(GCGs))
  seuratobj_corrected@assays$RNA <- seuratobj_corrected@assays$Corrected

  # Step 5: choose follow-up strategy based on contamination ratio
  if (contamination_ratio > 0.0003) {
    cli::cli_progress_step("Contamination ratio > 0.0003; using scCDC correction")
    if (verbose) {
      message(
        "Contamination ratio = ", contamination_ratio,
        " > 0.0003; using scCDC decontamination."
      )
    }
  } else {
    cli::cli_progress_step("Contamination ratio <= 0.0003; adding decontX")
    if (verbose) {
      message(
        "Contamination ratio = ", contamination_ratio,
        " <= 0.0003; adding decontX decontamination."
      )
    }

    counts <- Seurat::GetAssayData(
      seuratobj_corrected, assay = assay, layer = "counts"
    )
    res <- decontX::decontX(counts, z = Zc[[cluster_col]])
    seuratobj_corrected$contamination <- res$contamination

    if (verbose) {
      message(
        "Contamination range: ",
        round(min(res$contamination), 3), " ~ ",
        round(max(res$contamination), 3)
      )
    }

    # filter cells by threshold
    if (!is.null(decontX_threshold)) {
      keep <- seuratobj_corrected$contamination < decontX_threshold
      n_removed <- sum(!keep)
      if (n_removed > 0) {
        seuratobj_corrected <- seuratobj_corrected[, keep]
        if (verbose) {
          message(
            "Removed ", n_removed,
            " cells (contamination >= ", decontX_threshold, ")"
          )
        }
      } else {
        if (verbose) {
          message(
            "No cells removed (all contamination < ",
            decontX_threshold, ")"
          )
        }
      }
    }
  }

  # Step 6: finish
  cli::cli_progress_step("Decontamination complete")
  message(
    "Done. Please re-run normalisation, PCA/Harmony, ",
    "UMAP/t-SNE, and clustering."
  )

  cli::cli_progress_done()
  return(seuratobj_corrected)
}
