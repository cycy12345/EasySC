##使用scDblFinder去除双细胞
#object:Seurat对象
#dbr:预期的双细胞率（默认0.08，可填NULL）。默认情况下，假定每千个捕获的细胞中双细胞率为1%（即4000个细胞中为4%），这适用于10x数据集。将根据给定的比率对同型双细胞进行校正
#group_by:meta.data中的样本分类列名，若指定，将按样本分裂独立处理。若省略，将对所有细胞一起进行双细胞搜索。

##结果输出为一个list,其中sce_filter为过滤后Seurat对象，plot为双细胞图
##运行示例：res<- Easy_DbFinder(sce,group_by = "Patient",dbr=0.08)

library(cli)
library(Seurat)
library(scDblFinder)
library(BiocParallel)
library(ggplot2)

Easy_DbFinder <- function(object, dbr = 0.08, group_by = NULL) {
  
  # Step 1
  cli::cli_progress_step("转换为 SingleCellExperiment")
  re <- Seurat::as.SingleCellExperiment(object)
  
  # Step 2：最耗时，scDblFinder 内部 verbose=T 也会输出日志
  cli::cli_progress_step("运行 scDblFinder（请稍候）")
  re <- scDblFinder::scDblFinder(
    re,
    samples = group_by,
    BPPARAM = BiocParallel::MulticoreParam(3),
    dbr = dbr,
    verbose = TRUE
  )
  
  # Step 4
  cli::cli_progress_step("提取双细胞分类结果")
  object$Doublets <- re$scDblFinder.class
  df <- object@meta.data %>% as.data.frame()
  
  # 补全 tmp（原代码未定义）
  if (!is.null(group_by) && group_by %in% colnames(df)) {
    tmp <- as.data.frame(table(df$Doublets, df[[group_by]]))
  } else {
    tmp <- as.data.frame(table(df$Doublets))
    tmp$Var2 <- "All"
    tmp$Var1 <- tmp$Var1
  }
  
  # Step 5
  cli::cli_progress_step("生成可视化图表")
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
  
  # Step 6
  cli::cli_progress_step("过滤双细胞并输出结果")
  message("过滤前细胞数：", ncol(object))
  object <- subset(object, subset = Doublets == "singlet")
  message("过滤后细胞数：", ncol(object))
  
  # 结束
  cli::cli_progress_done()
  message("分析完成：sce_filter 为过滤后对象，plot 为双细胞图")
  
  list(sce_filter = object, plot = p)
}

species = "human", MT = 5, RP = 5, HB = 1,
                   Feature_high = 5000, Feature_low = 200


Easy_QC <- function(object, species = "human", MT = 5, RP = 5, HB = 1,
                    Feature_high = 5000, Feature_low = 200, plot_group = "sample",
                    color = c("#fbb4ae", "#b3cde3", "#ccebc5", "#decbe4",
                              "#fed9a6", "#ffffcc", "#e5d8bd", "#fddaec",
                              "#f2f2f2")) {
  
  # Step 1: 参数校验
  cli::cli_progress_step("验证输入参数")
  if (missing(object)) {
    stop("You must provide a Seurat object.")
  }
  if (!(species %in% c("human", "mouse"))) {
    stop("Species must be either 'human' or 'mouse'.")
  }
  if (is.null(plot_group)) {
    plot_group <- "orig.ident"
  }
  
  # Step 2: 计算质控指标
  cli::cli_progress_step("计算线粒体、核糖体和血红蛋白比例")
  if (species == "human") {
    object[["percent.mt"]] <- Seurat::PercentageFeatureSet(object, pattern = "^MT-")
    object[["percent.rp"]] <- Seurat::PercentageFeatureSet(object, pattern = "^RP[sl]")
    object[["percent.hb"]] <- Seurat::PercentageFeatureSet(object, pattern = "^HB[^(p)]")
  } else if (species == "mouse") {
    object[["percent.mt"]] <- Seurat::PercentageFeatureSet(object, pattern = "^Mt-")
    object[["percent.rp"]] <- Seurat::PercentageFeatureSet(object, pattern = "^Rp[sl]")
    object[["percent.hb"]] <- Seurat::PercentageFeatureSet(object, pattern = "^Hb[^(p)]")
  }
  
  # Step 3: 过滤细胞
  cli::cli_progress_step("按阈值过滤低质量细胞")
  message("过滤前细胞数：", ncol(object))
  object <- subset(x = object,
                   subset = nFeature_RNA > Feature_low &
                     nFeature_RNA < Feature_high &
                     percent.mt < MT &
                     percent.rp < RP &
                     percent.hb < HB)
  message("过滤后细胞数：", ncol(object))
  
  # Step 4: 提取数据
  cli::cli_progress_step("提取质控数据并转换长格式")
  df <- Seurat::FetchData(object, vars = c("nFeature_RNA", "nCount_RNA", 
                                           "percent.mt", "percent.rp", "percent.hb", 
                                           plot_group))
  df_long <- tidyr::pivot_longer(
    df, 
    cols = c(nFeature_RNA, nCount_RNA, percent.mt, percent.rp, percent.hb),
    names_to = "indicator",
    values_to = "value"
  )
  
  # Step 5: 绘制 nFeature / nCount
  cli::cli_progress_step("生成 nFeature / nCount 小提琴图")
  p1 <- ggpubr::ggviolin(
    df_long[df_long$indicator %in% c("nFeature_RNA", "nCount_RNA"), ],
    x = plot_group, y = "value", fill = "indicator", palette = color,
    add = "boxplot", add.params = list(fill = "white", error.plot = "linerange")
  ) + facet_wrap(~indicator, scales = "free_y", nrow = 2)
  
  # Step 6: 绘制百分比指标
  cli::cli_progress_step("生成百分比指标小提琴图")
  p2 <- ggpubr::ggviolin(
    df_long[df_long$indicator %in% c("percent.mt", "percent.rp", "percent.hb"), ],
    x = plot_group, y = "value", fill = "indicator", palette = color,
    add = "boxplot", add.params = list(fill = "white", error.plot = "linerange")
  ) + facet_wrap(~indicator, scales = "free_y", nrow = 3)
  
  # Step 7: 返回结果
  cli::cli_progress_step("汇总返回结果")
  res <- list(
    sce_filter = object,
    plot_Count_Feature = p1,
    plot_percent = p2
  )
  
  cli::cli_progress_done()
  message("质控完成：sce_filter 为过滤后对象，plot_Count_Feature / plot_percent 为质控图")
  
  return(res)
}


Easy_Normal <- function(object,
                        use_sct = FALSE,
                        vars.to.regress = NULL,
                        scale.factor = 10000,
                        nfeatures = 2000,
                        assay = "RNA",
                        verbose = TRUE,
                        ...) {
  
  # Step 1: 参数校验
  cli::cli_progress_step("验证输入参数")
  if (!inherits(object, "Seurat")) {
    stop("object 必须是 Seurat 对象")
  }
  if (!assay %in% names(object)) {
    stop("指定的 assay '", assay, "' 不存在")
  }
  
  # Step 2: 标准化
  cli::cli_progress_step("执行 LogNormalize 标准化")
  if (verbose) message("执行 Seurat 标准标准化 (LogNormalize)...")
  object <- NormalizeData(object,
                          normalization.method = "LogNormalize",
                          scale.factor = scale.factor,
                          assay = assay,
                          verbose = verbose,
                          ...)
  
  # Step 3: 高变基因
  cli::cli_progress_step("筛选高变基因")
  object <- FindVariableFeatures(object,
                                 selection.method = "vst",
                                 nfeatures = nfeatures,
                                 assay = assay,
                                 verbose = verbose)
  
  # Step 4: 缩放
  cli::cli_progress_step("执行 ScaleData")
  object <- ScaleData(object,
                      vars.to.regress = vars.to.regress,
                      assay = assay,
                      verbose = verbose,
                      ...)
  
  # Step 5: 设置默认 assay
  cli::cli_progress_step("设置默认 assay")
  DefaultAssay(object) <- assay
  
  if (verbose) {
    message("\n========== LogNormalize 完成 ==========")
    message("标准化数据存储在: ", assay, "$data")
    message("缩放数据存储在: ", assay, "$scale.data")
    message("高变基因数: ", length(VariableFeatures(object, assay = assay)))
    message("======================================\n")
  }
  
  # Step 6: 可选 SCTransform
  if (use_sct) {
    cli::cli_progress_step("执行 SCTransform")
    if (verbose) message("额外执行 SCTransform (结果保存在 SCT assay)...")
    
    if (requireNamespace("glmGamPoi", quietly = TRUE) && verbose) {
      message("检测到 glmGamPoi，将自动加速 SCTransform")
    }
    
    object <- SCTransform(object,
                          assay = assay,
                          vars.to.regress = vars.to.regress,
                          verbose = verbose,
                          ...)
    
    if (verbose) {
      message("\n========== SCTransform 完成 ==========")
      message("校正数据存储在: SCT$data (皮尔逊残差)")
      message("SCT assay 已设为默认 assay，可通过 DefaultAssay(object) <- 'RNA' 切换")
      message("\n【下游分析建议】")
      message("差异分析: 推荐使用 RNA assay 的 data slot")
      message("  FindMarkers(object, assay = 'RNA', slot = 'data')")
      message("可视化/降维: 推荐使用 SCT assay 的 scale.data")
      message("  RunPCA(object, assay = 'SCT')")
      message("细胞通讯: 使用 SCT assay 的 data slot")
      message("  GetAssayData(object, assay = 'SCT', slot = 'data')")
      message("====================================\n")
    }
  } else {
    if (verbose) {
      message("\n【下游分析建议】")
      message("差异分析: FindMarkers(object, assay = 'RNA', slot = 'data')")
      message("可视化/降维: 基于 RNA assay 的 scale.data")
      message("  RunPCA(object, assay = 'RNA')")
      message("====================================\n")
    }
  }
  
  cli::cli_progress_done()
  return(object)
}





Easy_harmony <- function(object,
                         group.by,
                         features = VariableFeatures(object),
                         npcs = 50,
                         verbose = TRUE,
                         ...) {
  
  # Step 1: 验证依赖与参数
  cli::cli_progress_step("验证依赖包与批次变量")
  if (!requireNamespace("harmony", quietly = TRUE)) {
    stop("请先安装 harmony 包: install.packages('harmony')")
  }
  if (!group.by %in% colnames(object[[]])) {
    stop("在 object 的 metadata 中找不到 '", group.by, "' 列")
  }
  
  # Step 2: PCA
  cli::cli_progress_step("运行 PCA")
  if (verbose) message("正在运行 PCA，使用 ", length(features), " 个高变基因...")
  object <- Seurat::RunPCA(object,
                           features = features,
                           npcs = npcs,
                           verbose = verbose)
  PC<- Easy_PC(object,)
  # Step 3: Harmony 批次校正
  cli::cli_progress_step("运行 Harmony 批次校正")
  if (verbose) message("正在运行 Harmony 批次校正，批次变量: ", group.by)
  object <- harmony::RunHarmony(object,dims =1:PC,
                                group.by.vars = group.by)
  
  # Step 4: 完成
  cli::cli_progress_step("完成降维整合")
  if (verbose) message("完成！Harmony 降维结果存储在 'harmony' 中。")
  
  cli::cli_progress_done()
  return(object)
}





Easy_PC <- function(seurat, reduction = "pca", cum = 90, var = 5) {
  
  # Step 1: 验证输入
  cli::cli_progress_step("验证输入参数")
  if (missing(seurat)) {
    stop("必须提供 Seurat 对象。")
  }
  if (!reduction %in% names(seurat@reductions)) {
    stop("指定的降维方法在 Seurat 对象中不存在。")
  }
  
  # Step 2: 提取标准差
  cli::cli_progress_step("提取降维标准差")
  if (reduction == "pca") {
    stdevs <- seurat[[reduction]]@stdev
  } else {
    embeddings <- seurat[[reduction]]@cell.embeddings
    stdevs <- apply(embeddings, 2, sd)
  }
  if (is.null(stdevs) || length(stdevs) == 0) {
    stop("指定的降维方法不包含有效的标准差。")
  }
  
  # Step 3: 计算贡献率
  cli::cli_progress_step("计算贡献率与累积贡献率")
  pct <- stdevs / sum(stdevs) * 100
  cumu <- cumsum(pct)
  
  # Step 4: 确定最优主成分数
  cli::cli_progress_step("确定最优主成分数")
  co1 <- which(cumu > cum & pct < var)[1]
  
  if (length(pct) > 1) {
    co2 <- sort(which((pct[1:(length(pct) - 1)] - pct[2:length(pct)]) > 0.1), decreasing = TRUE)[1] + 1
  } else {
    co2 <- NA
  }
  
  pcs <- min(co1, co2, na.rm = TRUE)
  
  # Step 5: 绘制可视化
  cli::cli_progress_step("生成 PC 选择图")
  plot_df <- data.frame(pct = pct, cumu = cumu, rank = 1:length(pct))
  
  p <- ggplot(plot_df, aes(x = cumu, y = pct, label = rank, color = rank > pcs)) +
    geom_text(check_overlap = TRUE, size = 3.5, fontface = "bold") +
    geom_vline(xintercept = cum, color = "grey", linetype = "dashed") +
    geom_hline(yintercept = min(pct[pct > var], na.rm = TRUE), color = "grey", linetype = "dashed") +
    scale_color_manual(values = c("TRUE" = "red", "FALSE" = "skyblue")) +
    labs(
      title = "Optimal Principal Components Selection",
      x = "Cumulative Percentage",
      y = "Percentage of Variation"
    ) +
    theme_test(base_rect_size = 1.5) +
    theme(
      plot.title = element_text(hjust = 0.5),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      legend.position = "none"
    )
  print(p)
  
  # Step 6: 完成
  cli::cli_progress_step("返回最优 PC 数")
  message("推荐使用主成分数: ", pcs)
  
  cli::cli_progress_done()
  return(pcs)
}








Easy_clustree <- function(object,
                          assay = c("RNA", "SCT"),
                          resolutions = c(0.01, 0.05, 0.1, 0.2, 0.3, 0.5, 0.8, 1),
                          algorithm = 1,
                          PC=20,
                          graph.name = NULL,
                          prefix = NULL) {
  
  # Step 1: 验证依赖包
  cli::cli_progress_step("验证 clustree 依赖包")
  if (!requireNamespace("clustree", quietly = TRUE)) {
    stop("请先安装 clustree 包：install.packages('clustree')")
  }
  
  # Step 2: 确定邻接图名称和列名前缀
  cli::cli_progress_step("确定邻接图名称与前缀")
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
  cli::cli_progress_step("构建KNN图")
  object<-FindNeighbors(
    object,
    reduction = "harmony",
    dims = 1:PC,
    verbose = FALSE
  )
  # Step 3: 检查邻接图是否存在
  cli::cli_progress_step("检查邻接图是否存在")
  if (!graph.name %in% names(object)) {
    stop("邻接图 '", graph.name, "' 不存在。请先运行 FindNeighbors 并指定正确的 assay 或 graph.name。")
  }
  
  # Step 4: 多分辨率聚类
  cli::cli_progress_step("多分辨率聚类")
  for (res in resolutions) {
    message("  正在聚类 resolution = ", res)
    object <- FindClusters(object,
                           # graph.name = graph.name,
                           resolution = res,
                           # algorithm = algorithm,
                           verbose = TRUE)
  }
  
  # Step 5: 绘制聚类树
  cli::cli_progress_step("绘制聚类树")
  p <- clustree::clustree(object@meta.data, prefix = prefix)
  
  # Step 6: 完成
  cli::cli_progress_step("返回结果")
  
  cli::cli_progress_done()
  return(list(object = object, plot = p))
}



Easy_CDC <- function(object, cluster_col, assay = "RNA", 
                     decontX_threshold = 0.2, verbose = TRUE) {
  
  # Step 1: 验证依赖包
  cli::cli_progress_step("验证依赖包")
  if (!requireNamespace("scCDC", quietly = TRUE)) {
    stop("请先安装 scCDC 包")
  }
  if (!requireNamespace("decontX", quietly = TRUE)) {
    stop("请先安装 decontX 包")
  }
  
  # Step 2: 检测污染基因
  cli::cli_progress_step("检测污染基因 (scCDC)")
  Zc <- object[[cluster_col]]
  Idents(object) <- Zc[[cluster_col]]
  GCGs <- scCDC::ContaminationDetection(object)
  
  # Step 3: 量化污染率
  cli::cli_progress_step("量化污染率")
  contamination_ratio <- scCDC::ContaminationQuantification(object, rownames(GCGs))
  
  # Step 4: 执行污染校正
  cli::cli_progress_step("执行污染校正 (scCDC)")
  seuratobj_corrected <- scCDC::ContaminationCorrection(object, rownames(GCGs))
  seuratobj_corrected@assays$RNA <- seuratobj_corrected@assays$Corrected
  
  # Step 5: 根据污染率选择后续策略
  if (contamination_ratio > 0.0003) {
    cli::cli_progress_step("污染率 > 0.0003，使用 scCDC 校正结果")
    if (verbose) message("污染率=", contamination_ratio, " 大于 0.0003，使用 scCDC 方法去除污染")
  } else {
    cli::cli_progress_step("污染率 ≤ 0.0003，联合 decontX 深度去污染")
    if (verbose) message("污染率=", contamination_ratio, " 小于 0.0003，增加 decontX 方法去除污染")
    
    counts <- Seurat::GetAssayData(seuratobj_corrected, assay = assay, layer = "counts")
    res <- decontX::decontX(counts, z = Zc[[cluster_col]])
    seuratobj_corrected$contamination <- res$contamination
    
    if (verbose) {
      message("污染比例范围：", round(min(res$contamination), 3), " ~ ", round(max(res$contamination), 3))
    }
    
    # 根据阈值过滤细胞
    if (!is.null(decontX_threshold)) {
      keep <- seuratobj_corrected$contamination < decontX_threshold
      n_removed <- sum(!keep)
      if (n_removed > 0) {
        seuratobj_corrected <- seuratobj_corrected[, keep]
        if (verbose) message("过滤掉 ", n_removed, " 个细胞（污染比例 >= ", decontX_threshold, "）")
      } else {
        if (verbose) message("没有细胞被过滤（所有细胞污染比例 < ", decontX_threshold, "）")
      }
    }
  }
  
  # Step 6: 完成
  cli::cli_progress_step("去污染完成")
  message("去污染已完成，请重新执行降维聚类（Easy_Normal、run_pca_harmony、RunUMAP、Easy_clustree）")
  
  cli::cli_progress_done()
  return(seuratobj_corrected)
}
