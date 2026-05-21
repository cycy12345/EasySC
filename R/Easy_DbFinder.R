#' Doublet Detection and Filtering with scDblFinder
#'
#' @description
#' Detects and removes doublets (multiplets) from a Seurat object using
#' \code{scDblFinder}. The function converts the Seurat object to
#' \code{SingleCellExperiment}, runs \code{scDblFinder}, and returns both
#' the filtered object and a bar plot showing singlet/doublet distribution.
#'
#' @param object A \code{Seurat} object.
#' @param dbr Expected doublet rate. Default is \code{0.08} (8\%).
#'   Set to \code{NULL} to let \code{scDblFinder} estimate automatically.
#'   The rule of thumb for 10x data is ~1\% per 1,000 captured cells
#'   (e.g., 4\% for 4,000 cells).
#' @param group_by Column name in \code{object@meta.data} indicating sample
#'   grouping. If provided, doublet detection is performed independently for
#'   each sample. If \code{NULL} (default), all cells are processed together.
#'
#' @return A named \code{list} with two elements:
#'   \itemize{
#'     \item \code{sce_filter}: The input \code{Seurat} object after removing
#'       doublets (only \code{"singlet"} cells retained).
#'     \item \code{plot}: A \code{ggplot2} bar plot visualising singlet/doublet
#'       counts per sample (or overall when \code{group_by = NULL}).
#'   }
#'
#' @details
#' \code{scDblFinder} is run with \code{BiocParallel::MulticoreParam(3)} for
#' moderate parallelisation. Progress messages are printed via \code{cli}.
#'
#' @seealso \code{\link[scDblFinder]{scDblFinder}},
#'   \code{\link[Seurat]{as.SingleCellExperiment}}
#'
#' @examples
#' \dontrun{
#'   # Basic usage (all cells together)
#'   res <- Easy_DbFinder(sce, dbr = 0.08)
#'   res$plot
#'   sce_clean <- res$sce_filter
#'
#'   # Per-sample doublet detection
#'   res <- Easy_DbFinder(sce, group_by = "Patient", dbr = 0.08)
#'   sce_clean <- res$sce_filter
#' }
#'
#' @importFrom Seurat as.SingleCellExperiment subset
#' @importFrom scDblFinder scDblFinder
#' @importFrom BiocParallel MulticoreParam
#' @importFrom ggplot2 ggplot aes geom_bar geom_text ggtitle xlab ylab labs
#'   theme element_text position_dodge scale_fill_manual
#' @importFrom cli cli_progress_step cli_progress_done
#' @export
Easy_DbFinder <- function(object, dbr = 0.08, group_by = NULL) {

  # Step 1
  cli::cli_progress_step("Converting to SingleCellExperiment")
  re <- Seurat::as.SingleCellExperiment(object)

  # Step 2: run scDblFinder (most time-consuming step)
  cli::cli_progress_step("Running scDblFinder (please wait)")
  re <- scDblFinder::scDblFinder(
    re,
    samples = group_by,
    BPPARAM = BiocParallel::MulticoreParam(3),
    dbr = dbr,
    verbose = TRUE
  )

  # Step 3: extract classification results
  cli::cli_progress_step("Extracting doublet classification")
  object$Doublets <- re$scDblFinder.class
  df <- object@meta.data %>% as.data.frame()

  # prepare count table for plotting
  if (!is.null(group_by) && group_by %in% colnames(df)) {
    tmp <- as.data.frame(table(df$Doublets, df[[group_by]]))
  } else {
    tmp <- as.data.frame(table(df$Doublets))
    tmp$Var2 <- "All"
  }

  # Step 4: generate plot
  cli::cli_progress_step("Generating visualisation")
  p <- ggplot(tmp, aes(x = Var2, y = Freq, fill = Var1)) +
    geom_bar(stat = "identity", position = position_dodge()) +
    geom_text(aes(label = Freq), size = 3.5,
              position = position_dodge(width = 0.9)) +
    ggtitle("Singlet/Doublet Distribution") +
    xlab(ifelse(is.null(group_by), "Sample", group_by)) +
    labs(fill = "") +
    ylab("Cell Count") +
    theme(legend.position = "bottom",
          axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))

  # Step 5: filter doublets
  cli::cli_progress_step("Filtering doublets")
  message("Cells before filtering: ", ncol(object))
  object <- subset(object, subset = Doublets == "singlet")
  message("Cells after filtering: ", ncol(object))

  cli::cli_progress_done()
  message("Done: sce_filter = filtered object, plot = doublet bar plot")

  list(sce_filter = object, plot = p)
}
