#' Quality Control for scRNA-seq Data
#'
#' @description
#' Computes mitochondrial (MT), ribosomal (RP), and haemoglobin (HB)
#' gene percentages; filters low-quality cells; and returns violin plots
#' for \code{nFeature_RNA}, \code{nCount_RNA}, and the three percentage
#' metrics.
#'
#' @param object A \code{Seurat} object. Must contain \code{nFeature_RNA}
#'   and \code{nCount_RNA} in its metadata.
#' @param species Species identifier. Either \code{"human"} (default) or
#'   \code{"mouse"}. Determines the gene-name patterns used for MT, RP,
#'   and HB detection.
#' @param MT Threshold for mitochondrial gene percentage (\%). Cells with
#'   \code{percent.mt >= MT} are removed. Default is \code{5}.
#' @param RP Threshold for ribosomal gene percentage (\%). Cells with
#'   \code{percent.rp >= RP} are removed. Default is \code{5}.
#' @param HB Threshold for haemoglobin gene percentage (\%). Cells with
#'   \code{percent.hb >= HB} are removed. Default is \code{1}.
#' @param Feature_high Upper bound for \code{nFeature_RNA}. Cells with
#'   \code{nFeature_RNA > Feature_high} are removed. Default \code{5000}.
#' @param Feature_low Lower bound for \code{nFeature_RNA}. Cells with
#'   \code{nFeature_RNA < Feature_low} are removed. Default \code{200}.
#' @param plot_group Column name in \code{object@meta.data} for grouping
#'   the violin plots (e.g. \code{"sample"}, \code{"Patient"}). If
#'   \code{NULL}, defaults to \code{"orig.ident"}.
#' @param color A vector of colours for \code{ggpubr::ggviolin}.
#'   Default is a 9-colour pastel palette.
#'
#' @return A named \code{list} with three elements:
#'   \itemize{
#'     \item \code{sce_filter}: The \code{Seurat} object after filtering.
#'     \item \code{plot_Count_Feature}: Violin plot of
#'       \code{nFeature_RNA} and \code{nCount_RNA}.
#'     \item \code{plot_percent}: Violin plot of \code{percent.mt},
#'       \code{percent.rp}, and \code{percent.hb}.
#'   }
#'
#' @details
#' Gene patterns:
#' \itemize{
#'   \item Human — MT: \code{^MT-}, RP: \code{^RP[sl]}, HB: \code{^HB[^(p)]}
#'   \item Mouse — MT: \code{^Mt-}, RP: \code{^Rp[sl]}, HB: \code{^Hb[^(p)]}
#' }
#'
#' @seealso \code{\link[Seurat]{PercentageFeatureSet}},
#'   \code{\link[Seurat]{FetchData}}, \code{\link[ggpubr]{ggviolin}}
#'
#' @examples
#' \dontrun{
#'   sce_qc <- Easy_QC(sce, species = "human", MT = 5, RP = 5, HB = 1,
#'                     Feature_high = 5000, Feature_low = 200,
#'                     plot_group = "Patient")
#'   sce_qc$plot_Count_Feature
#'   sce_qc$plot_percent
#'   sce_clean <- sce_qc$sce_filter
#' }
#'
#' @importFrom Seurat PercentageFeatureSet FetchData subset
#' @importFrom tidyr pivot_longer
#' @importFrom ggpubr ggviolin
#' @importFrom ggplot2 facet_wrap
#' @importFrom cli cli_progress_step cli_progress_done
#' @export
Easy_QC <- function(object, species = "human", MT = 5, RP = 5, HB = 1,
                    Feature_high = 5000, Feature_low = 200,
                    plot_group = "sample",
                    color = c("#fbb4ae", "#b3cde3", "#ccebc5", "#decbe4",
                              "#fed9a6", "#ffffcc", "#e5d8bd", "#fddaec",
                              "#f2f2f2")) {

  # Step 1: argument checks
  cli::cli_progress_step("Validating input arguments")
  if (missing(object)) {
    stop("You must provide a Seurat object.")
  }
  if (!(species %in% c("human", "mouse"))) {
    stop("Species must be either 'human' or 'mouse'.")
  }
  if (is.null(plot_group)) {
    plot_group <- "orig.ident"
  }

  # Step 2: compute QC metrics
  cli::cli_progress_step("Calculating MT / RP / HB percentages")
  if (species == "human") {
    object[["percent.mt"]] <- Seurat::PercentageFeatureSet(object, pattern = "^MT-")
    object[["percent.rp"]] <- Seurat::PercentageFeatureSet(object, pattern = "^RPS|^RPL")
    object[["percent.hb"]] <- Seurat::PercentageFeatureSet(object, pattern = "^HB[^(p)]")
  } else if (species == "mouse") {
    object[["percent.mt"]] <- Seurat::PercentageFeatureSet(object, pattern = "^Mt-")
    object[["percent.rp"]] <- Seurat::PercentageFeatureSet(object, pattern = "^Rps|^Rpl")
    object[["percent.hb"]] <- Seurat::PercentageFeatureSet(object, pattern = "^Hb[^(p)]")
  }

  # Step 3: filter cells
  cli::cli_progress_step("Filtering low-quality cells")
  message("Cells before filtering: ", ncol(object))
  object <- subset(
    x = object,
    subset = nFeature_RNA > Feature_low &
      nFeature_RNA < Feature_high &
      percent.mt < MT &
      percent.rp > RP &
      percent.hb < HB
  )
  message("Cells after filtering: ", ncol(object))

  # Step 4: extract and pivot data
  cli::cli_progress_step("Extracting QC data")
  df <- Seurat::FetchData(
    object,
    vars = c("nFeature_RNA", "nCount_RNA",
             "percent.mt", "percent.rp", "percent.hb",
             plot_group)
  )
  df_long <- tidyr::pivot_longer(
    df,
    cols = c(nFeature_RNA, nCount_RNA, percent.mt, percent.rp, percent.hb),
    names_to = "indicator",
    values_to = "value"
  )

  # Step 5: nFeature / nCount violin plot
  cli::cli_progress_step("Generating nFeature / nCount violin plot")
  p1 <- ggpubr::ggviolin(
    df_long[df_long$indicator %in% c("nFeature_RNA", "nCount_RNA"), ],
    x = plot_group, y = "value", fill = "indicator", palette = color,
    add = "boxplot", add.params = list(fill = "white", error.plot = "linerange")
  ) + facet_wrap(~indicator, scales = "free_y", nrow = 2)

  # Step 6: percentage metrics violin plot
  cli::cli_progress_step("Generating percentage metrics violin plot")
  p2 <- ggpubr::ggviolin(
    df_long[df_long$indicator %in% c("percent.mt", "percent.rp", "percent.hb"), ],
    x = plot_group, y = "value", fill = "indicator", palette = color,
    add = "boxplot", add.params = list(fill = "white", error.plot = "linerange")
  ) + facet_wrap(~indicator, scales = "free_y", nrow = 3)

  # Step 7: return results
  cli::cli_progress_step("Summarising results")
  res <- list(
    sce_filter = object,
    plot_Count_Feature = p1,
    plot_percent = p2
  )

  cli::cli_progress_done()
  message("QC done: sce_filter = filtered object; ",
          "plot_Count_Feature / plot_percent = QC plots")

  return(res)
}
