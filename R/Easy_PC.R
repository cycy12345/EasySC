#' Select Optimal Number of Principal Components
#'
#' @description
#' Determines the optimal number of principal components
#' (PCs) to retain by combining two criteria:
#' cumulative variance > \code{cum}\% and individual
#' variance contribution < \code{var}\%; plus a minimum
#' drop-off criterion (difference between consecutive PCs
#' > 0.1\%).
#'
#' @param seurat A \code{Seurat} object containing a PCA
#'   reduction.
#' @param reduction Name of the reduction to inspect.
#'   Default \code{"pca"}.
#' @param cum Cumulative variance threshold (\%). The first
#'   PC where cumulative variance exceeds this value is
#'   candidate 1. Default \code{90}.
#' @param var Individual variance threshold (\%). Only PCs
#'   with contribution below this value are considered for
#'   \code{co1}. Default \code{5}.
#'
#' @return An integer: the recommended number of PCs to use.
#'   A scatter plot is printed to the active graphics device.
#'
#' @details
#' The two criteria are:
#' \enumerate{
#'   \item \code{co1}: first PC where cumulative percentage
#'     > \code{cum} AND individual percentage < \code{var}.
#'   \item \code{co2}: last PC where the drop from the
#'     previous PC > 0.1\%.
#' }
#' The returned value is \code{min(co1, co2)}.
#'
#' @seealso \code{\link[Seurat]{RunPCA}}
#'
#' @examples
#' \dontrun{
#'   sce <- Seurat::RunPCA(sce, features = VariableFeatures(sce))
#'   npc <- Easy_PC(sce, reduction = "pca", cum = 90, var = 5)
#'   npc
#' }
#'
#' @importFrom cli cli_progress_step cli_progress_done
#' @importFrom ggplot2 ggplot aes geom_text geom_vline
#'   geom_hline scale_color_manual labs theme_bw theme
#'   element_text element_blank
#' @export
Easy_PC <- function(seurat, reduction = "pca", cum = 90, var = 5) {

  # Step 1: validate input
  cli::cli_progress_step("Validating input")
  if (missing(seurat)) {
    stop("A Seurat object must be provided.")
  }
  if (!reduction %in% names(seurat@reductions)) {
    stop("Specified reduction does not exist in the Seurat object.")
  }

  # Step 2: extract standard deviations
  cli::cli_progress_step("Extracting standard deviations")
  if (reduction == "pca") {
    stdevs <- seurat[[reduction]]@stdev
  } else {
    embeddings <- seurat[[reduction]]@cell.embeddings
    stdevs <- apply(embeddings, 2, sd)
  }
  if (is.null(stdevs) || length(stdevs) == 0) {
    stop("Specified reduction contains no valid standard deviations.")
  }

  # Step 3: compute contributions
  cli::cli_progress_step("Computing variance contributions")
  pct <- stdevs / sum(stdevs) * 100
  cumu <- cumsum(pct)

  # Step 4: determine optimal PC count
  cli::cli_progress_step("Determining optimal PC count")
  co1 <- which(cumu > cum & pct < var)[1]

  if (length(pct) > 1) {
    co2 <- sort(
      which((pct[1:(length(pct) - 1)] - pct[2:length(pct)]) > 0.1),
      decreasing = TRUE
    )[1] + 1
  } else {
    co2 <- NA
  }

  pcs <- min(co1, co2, na.rm = TRUE)

  # Step 5: plot
  cli::cli_progress_step("Generating PC selection plot")
  plot_df <- data.frame(pct = pct, cumu = cumu, rank = seq_along(pct))

  p <- ggplot(plot_df, aes(x = cumu, y = pct, label = rank, color = rank > pcs)) +
    geom_text(check_overlap = TRUE, size = 3.5, fontface = "bold") +
    geom_vline(xintercept = cum, color = "grey", linetype = "dashed") +
    geom_hline(
      yintercept = min(pct[pct > var], na.rm = TRUE),
      color = "grey",
      linetype = "dashed"
    ) +
    scale_color_manual(values = c("TRUE" = "red", "FALSE" = "skyblue")) +
    labs(
      title = "Optimal Principal Components Selection",
      x = "Cumulative Percentage",
      y = "Percentage of Variation"
    ) +
    theme_bw(base_rect_size = 1.5) +
    theme(
      plot.title = element_text(hjust = 0.5),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      legend.position = "none"
    )
  print(p)

  # Step 6: finish
  cli::cli_progress_step("Returning optimal PC count")
  message("Recommended number of PCs: ", pcs)

  cli::cli_progress_done()
  return(pcs)
}
