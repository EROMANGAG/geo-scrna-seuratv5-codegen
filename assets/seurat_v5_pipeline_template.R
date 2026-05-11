# Template only: Codex should tailor this script to the current GEO dataset.
# The generated script must replace all <PLACEHOLDER> values after inspecting GEO/SRA/article metadata.

# ---- 00. project metadata ----
project_id <- "<GSEID>"

meta <- c(
  project_GEOid = "<GSEID>",
  cancer_type   = "<TCGA_ABBREVIATION_OR_NA>",
  platform      = "<PLATFORM_OR_NA>",
  patient_id    = "<PATIENT_IDS_OR_NA>",
  sample_id     = "<SAMPLE_IDS_OR_NA>",
  sample_type   = "<T_N_L_OR_NA>"
)

sample_meta <- data.frame(
  sample_id = c("<GSM_OR_SAMPLE_ID>"),
  patient_id = c("NA"),
  sample_type = c("NA"),
  sra_run = c("NA"),
  platform = c("<PLATFORM_OR_NA>"),
  stringsAsFactors = FALSE
)

# ---- 01. packages ----
require_or_stop <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop(sprintf("Package '%s' is required. Install it before running this script.", pkg), call. = FALSE)
  }
}

required_packages <- c("Seurat", "qs", "Matrix", "data.table", "dplyr", "ggplot2", "DoubletFinder")
invisible(lapply(required_packages, require_or_stop))

library(Seurat)
library(qs)
library(Matrix)
library(data.table)
library(dplyr)
library(ggplot2)

# ---- 02. parameters ----
# 中文注释：这里必须写入用户明确提供的工作目录和随机种子
work_dir <- "<WORK_DIR>"
seed <- <SEED>
setwd(work_dir)
set.seed(seed)

input_dir <- file.path(work_dir, "data", "<GSEID>")
output_dir <- file.path(work_dir, "output", "<GSEID>")
fig_dir <- file.path(output_dir, "figures")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

min_features <- 200
max_features <- 6000
max_counts <- Inf
max_percent_mt <- 20
min_cells_per_gene <- 3
n_hvg <- 3000
dims_use <- NULL
npcs <- 40
vars_to_regress <- c("nFeature_RNA", "nCount_RNA", "percent.mt")
cluster_resolution <- 0.5
mt_pattern <- "^MT-"
contamination_method <- "auto" # auto, soupx, decontx, none
run_doubletfinder <- TRUE
doubletfinder_pN <- 0.25
doubletfinder_npcs <- 30
doubletfinder_resolution <- 0.1
axel_proxy <- "http://www.cirno999.cn:12306"
axel_connections <- 10

multiplet_rates_10x <- data.frame(
  Multiplet_rate = c(0.004, 0.008, 0.0160, 0.023, 0.031, 0.039, 0.046, 0.054, 0.061, 0.069, 0.076),
  Loaded_cells = c(800, 1600, 3200, 4800, 6400, 8000, 9600, 11200, 12800, 14400, 16000),
  Recovered_cells = c(500, 1000, 2000, 3000, 4000, 5000, 6000, 7000, 8000, 9000, 10000)
)

# ---- 03. helper functions ----
# 中文注释：当用户没有显式给出 dims_use 时，按拐点和累计贡献率自动确定主成分数
calculate_min_pc <- function(obj, reduction = "harmony") {
  stdv <- obj[[reduction]]@stdev
  percent_stdv <- (stdv / sum(stdv)) * 100
  cumulative <- cumsum(percent_stdv)
  co1 <- which(cumulative > 90 & percent_stdv < 5)[1]
  co2 <- sort(
    which((percent_stdv[1:(length(percent_stdv) - 1)] -
             percent_stdv[2:length(percent_stdv)]) > 0.1),
    decreasing = TRUE
  )[1] + 1
  min_pc <- min(co1, co2, na.rm = TRUE)
  if (!is.finite(min_pc)) {
    min_pc <- min(30, length(stdv))
  }
  max(2, min_pc)
}

resolve_dims_use <- function(obj, reduction = "harmony") {
  if (!is.null(dims_use)) {
    return(dims_use)
  }
  seq_len(calculate_min_pc(obj, reduction = reduction))
}

cleanup_vars <- function(var_names, envir = .GlobalEnv) {
  existing <- intersect(var_names, ls(envir = envir, all.names = TRUE))
  if (length(existing) > 0) {
    rm(list = existing, envir = envir)
  }
  invisible(gc())
}

# 中文注释：公共数据下载统一使用 axel，不使用 R 内置下载函数
download_with_axel <- function(url, dest_path) {
  dir.create(dirname(dest_path), recursive = TRUE, showWarnings = FALSE)
  status <- system2(
    "axel",
    args = c("-n", as.character(axel_connections), "-o", dest_path, url),
    env = c(
      sprintf("http_proxy=%s", axel_proxy),
      sprintf("https_proxy=%s", axel_proxy),
      sprintf("all_proxy=%s", axel_proxy)
    )
  )
  if (!identical(status, 0L) || !file.exists(dest_path)) {
    stop(sprintf("axel download failed for: %s", url), call. = FALSE)
  }
  dest_path
}

# 中文注释：标准预处理流程固定采用技能约定的参数
preprocess_standard <- function(obj, npcs_use = npcs) {
  obj <- NormalizeData(obj, normalization.method = "LogNormalize",
                       scale.factor = 10000, margin = 1, assay = "RNA")
  obj <- FindVariableFeatures(obj, nfeatures = 3000)
  obj <- ScaleData(obj, vars.to.regress = vars_to_regress)
  obj <- RunPCA(obj, features = VariableFeatures(obj), npcs = npcs_use)
  obj
}

estimate_multiplet_rate <- function(cell_number) {
  idx <- which(multiplet_rates_10x$Recovered_cells < cell_number)
  if (length(idx) == 0) {
    return(multiplet_rates_10x$Multiplet_rate[1])
  }
  multiplet_rates_10x$Multiplet_rate[max(idx)]
}

run_doubletfinder_one_sample <- function(o) {
  o <- NormalizeData(o, normalization.method = "LogNormalize",
                     scale.factor = 10000, margin = 1, assay = "RNA")
  o <- FindVariableFeatures(o, nfeatures = 3000)
  o <- ScaleData(o, vars.to.regress = vars_to_regress)
  o <- RunPCA(o, features = VariableFeatures(o), npcs = doubletfinder_npcs)

  min_pc <- calculate_min_pc(o, reduction = "pca")
  o <- RunUMAP(o, reduction = "pca", dims = 1:min_pc)
  o <- FindNeighbors(o, reduction = "pca", dims = 1:min_pc)
  o <- FindClusters(o, resolution = doubletfinder_resolution)

  sweep.res.hcc <- DoubletFinder::paramSweep(o, PCs = 1:min_pc, sct = FALSE)
  sweep.stats <- DoubletFinder::summarizeSweep(sweep.res.hcc, GT = FALSE)
  bcmvn <- DoubletFinder::find.pK(sweep.stats)
  pK_bcmvn <- as.numeric(as.vector(bcmvn$pK[which.max(bcmvn$BCmetric)]))

  homotypic.prop <- DoubletFinder::modelHomotypic(o$seurat_clusters)
  multiplet_rate <- estimate_multiplet_rate(nrow(o@meta.data))
  nExp.poi <- round(multiplet_rate * nrow(o@meta.data))
  nExp.poi.adj <- round(nExp.poi * (1 - homotypic.prop))

  DoubletFinder::doubletFinder(
    seu = o,
    PCs = 1:min_pc,
    pN = doubletfinder_pN,
    pK = pK_bcmvn,
    nExp = nExp.poi.adj,
    sct = FALSE
  )
}

remove_doublets_by_sample <- function(obj) {
  # 中文注释：双细胞按样本分别检测，可以减少不同样本混合带来的参数偏移
  sample_obj <- SplitObject(obj, split.by = "sample")

  h <- lapply(names(sample_obj), function(n) {
    run_doubletfinder_one_sample(sample_obj[[n]])
  })
  names(h) <- names(sample_obj)

  h.singlet <- lapply(names(sample_obj), function(n) {
    o <- h[[n]]
    colnames(o@meta.data)[grepl("DF.classifications", colnames(o@meta.data))] <- "doublet_finder"
    subset(o, subset = doublet_finder == "Singlet")
  })
  names(h.singlet) <- names(sample_obj)

  if (length(h.singlet) == 1) {
    singlet_merged <- h.singlet[[1]]
  } else {
    singlet_merged <- merge(h.singlet[[1]], y = h.singlet[-1], add.cell.ids = names(h.singlet))
  }
  singlet_merged <- JoinLayers(singlet_merged)

  # 中文注释：单细胞结果已经合并完成，可释放 DoubletFinder 中间对象
  cleanup_vars(c("h", "h.singlet", "sample_obj"), envir = environment())
  singlet_merged
}

# ---- 04. read counts ----
# 中文注释：这里必须替换为与 GEO 文件格式匹配的读取逻辑
# 中文注释：如果 GEO 补充文件需要先下载，请调用 download_with_axel(url, dest_path)
read_one_sample <- function(sample_id, path) {
  counts <- Read10X(data.dir = path)
  one_meta <- sample_meta[match(sample_id, sample_meta$sample_id), , drop = FALSE]
  obj <- CreateSeuratObject(
    counts = counts,
    project = project_id,
    min.cells = min_cells_per_gene,
    meta.data = one_meta
  )
  obj$sample_id <- sample_id
  obj$sample <- sample_id
  obj <- RenameCells(obj, add.cell.id = sample_id)
  obj
}

# object_list <- setNames(lapply(sample_paths, function(p) read_one_sample(...)), sample_meta$sample_id)
# 中文注释：如果读取阶段产生了额外的大矩阵或临时表，此处应及时删除
# cleanup_vars(c("raw_counts", "filtered_counts", "sample_paths"))

# ---- 05. contamination correction ----
run_contamination_correction <- function(obj) {
  # 中文注释：优先使用 SoupX；若只有 filtered matrix，可选 DecontX；都不可用时明确跳过
  obj
}

# 中文注释：污染校正阶段结束后，及时清理不再需要的原始临时对象
# cleanup_vars(c("raw_counts", "tod", "toc", "sce_tmp"))

# ---- 06. QC helpers ----
qc_filter_one <- function(obj) {
  obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = mt_pattern)
  subset(
    obj,
    subset = nFeature_RNA >= min_features &
      nFeature_RNA <= max_features &
      nCount_RNA <= max_counts &
      percent.mt <= max_percent_mt
  )
}

# ---- 07. merge and doublet removal ----
# qc_list <- lapply(object_list, qc_filter_one)
# obj_qc <- if (length(qc_list) == 1) qc_list[[1]] else merge(qc_list[[1]], y = qc_list[-1], add.cell.ids = names(qc_list))
# obj_qc <- JoinLayers(obj_qc)
# cleanup_vars(c("qc_list", "object_list"))

# 中文注释：按样本拆分后分别去除双细胞，再合并回单细胞对象
# obj_singlet_merged_raw <- if (run_doubletfinder) remove_doublets_by_sample(obj_qc) else obj_qc
# cleanup_vars(c("obj_qc"))

# 中文注释：去双细胞后的合并对象需要重新运行标准预处理流程
# obj_singlet_merged <- preprocess_standard(obj_singlet_merged_raw)

# ---- 08. integration templates ----
# 中文注释：merge-only 基线用于检查是否发生过度校正
obj_merged_no_correction <- obj_singlet_merged
dims_pca <- resolve_dims_use(obj_merged_no_correction, reduction = "pca")
obj_merged_no_correction <- RunUMAP(
  obj_merged_no_correction,
  reduction = "pca",
  dims = dims_pca,
  reduction.name = "umap.unintegrated"
)
obj_merged_no_correction <- FindNeighbors(obj_merged_no_correction, reduction = "pca", dims = dims_pca)
obj_merged_no_correction <- FindClusters(obj_merged_no_correction, resolution = cluster_resolution)

obj_rpca <- obj_singlet_merged_raw
obj_rpca[["RNA"]] <- split(obj_rpca[["RNA"]], f = obj_rpca$sample_id)
obj_rpca$sample <- obj_rpca$sample_id
obj_rpca <- preprocess_standard(obj_rpca)
obj_rpca <- IntegrateLayers(
  object = obj_rpca,
  method = RPCAIntegration,
  orig.reduction = "pca",
  new.reduction = "integrated.rpca",
  verbose = FALSE
)
obj_rpca <- JoinLayers(obj_rpca)
dims_rpca <- resolve_dims_use(obj_rpca, reduction = "integrated.rpca")
obj_rpca <- FindNeighbors(obj_rpca, reduction = "integrated.rpca", dims = dims_rpca)
obj_rpca <- RunUMAP(obj_rpca, reduction = "integrated.rpca", dims = dims_rpca, reduction.name = "umap.rpca")
obj_rpca <- FindClusters(obj_rpca, resolution = cluster_resolution)

obj_harmony <- NULL
if (!requireNamespace("harmony", quietly = TRUE)) {
  message("Package 'harmony' is not installed; skip integrated_harmony.")
} else {
  obj_harmony <- obj_singlet_merged_raw
  obj_harmony$sample <- obj_harmony$sample_id
  obj_harmony <- preprocess_standard(obj_harmony)
  obj_harmony <- RunHarmony(obj_harmony, "sample")
  dims_harmony <- resolve_dims_use(obj_harmony, reduction = "harmony")
  obj_harmony <- FindNeighbors(obj_harmony, reduction = "harmony", dims = dims_harmony)
  obj_harmony <- RunUMAP(obj_harmony, reduction = "harmony", dims = dims_harmony, reduction.name = "umap.harmony")
  obj_harmony <- FindClusters(obj_harmony, resolution = cluster_resolution)
}

# 中文注释：三个结果对象都生成后，去双流程中间对象和降维辅助变量可以释放
cleanup_vars(c("obj_singlet_merged_raw", "obj_singlet_merged", "dims_pca", "dims_rpca", "dims_harmony"))

# ---- 09. final object contract ----
seurat_results <- Filter(
  Negate(is.null),
  list(
    merged_no_correction = obj_merged_no_correction,
    integrated_rpca = obj_rpca,
    integrated_harmony = obj_harmony
  )
)

qsave(seurat_results, file.path(output_dir, paste0(project_id, "_seurat_results.qs")))

# 中文注释：最终结果已经收纳到 seurat_results 中，可移除重复绑定以降低后续内存压力
cleanup_vars(c("obj_merged_no_correction", "obj_rpca", "obj_harmony"))
