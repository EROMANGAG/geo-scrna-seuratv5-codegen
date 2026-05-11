---
name: geo-scrna-seuratv5-codegen
description: Generate a high-readability Seurat v5 R pipeline from a scRNA-seq GEO accession, GEO URL, or SRA Selector URL, including GEO/SRA/PubMed metadata extraction, a named meta R vector, RNA contamination correction, QC, normalization, dimensionality reduction, and at least two merge/integration Seurat objects saved with qs.
---

# GEO scRNA-seq → Seurat v5 code generator

Use this skill when the user gives a GEO accession, GEO page URL, SRA Selector URL, SRP/PRJNA/SRA accession, or article link and asks for R code to process a single-cell RNA-seq dataset into Seurat objects.

Do not use this skill for bulk RNA-seq, spatial transcriptomics, scATAC-only datasets, or non-R pipelines unless the user explicitly asks to adapt the workflow.

## Required user inputs before generation

Before generating the R script, require the user to provide both:

1. A working directory path that will be written into the script and applied with `setwd(work_dir)`.
2. A numeric seed value that will be written into the script and applied with `set.seed(seed)`.

If either value is missing, ask for it before writing the final code. Do not silently assume a default path or seed.

## Required final deliverable

Produce a complete, directly runnable R script tailored to the dataset. The generated script must:

1. Use **Seurat v5** as the main analysis framework.
2. Use **qs** for all saved R objects. Use `qsave()` / `qread()`. Do not use `saveRDS()` / `readRDS()` for generated intermediate or final objects.
3. Return and save a named list called `seurat_results` whose elements are Seurat objects.
4. Include at least two different merge/integration strategies, and save each corresponding Seurat object separately.
5. Include a readable parameter block at the top so the user can easily modify `work_dir`, `seed`, QC cutoffs, dimensions, resolution, contamination method, and integration methods. The script must call both `setwd(work_dir)` and `set.seed(seed)`.
6. Include a named R vector called `meta`, and when multiple samples exist also include `sample_meta` as a data frame.
7. Implement the full processing flow: RNA contamination correction → basic QC → DoubletFinder-based doublet removal → normalization → HVG selection → scaling → PCA → UMAP → neighbors → clustering.
8. Include verification code at the end: `stopifnot(inherits(..., "Seurat"))`, `print(seurat_results)`, and checks that `.qs` outputs exist.
9. Add concise Chinese comments at key analysis steps, helper functions, and any place where parameter choices or assumptions need explanation.
10. After completing each major section, remove large temporary variables that are no longer needed and run `gc()` to keep memory usage low.
11. If the script needs to download public supplementary files, do not use `download.file()` or other R built-in download helpers. Use `axel` with 10 connections and the fixed proxy `http://www.cirno999.cn:12306`.

## Required literature-derived integration logic

Use the user's uploaded benchmark article as the method-selection reference. Encode these principles in comments and in the method choices:

- Always generate an unintegrated / merge-only baseline object. It is required for detecting overcorrection and for inspecting whether biological differences are being erased.
- Use highly variable gene selection before integration by default. The benchmark reports that HVG selection generally improves integration performance.
- Treat scaling as useful but potentially overcorrecting: scaling can improve batch removal but may reduce biological signal conservation. Keep parameters explicit and avoid hiding this choice.
- For an R/Seurat-v5-only deliverable, prefer `RPCAIntegration` and Harmony via `RunHarmony(obj, "sample")` as the default two integration methods, with `merged_no_correction` as the baseline.
- For small/simple data with distinct biological signal, Harmony is acceptable and often convenient.
- For strong or nested batch effects, include RPCA and Harmony outputs and instruct the user to compare both against the merge-only baseline.
- Mention in comments that Scanorama, scVI/scANVI are strong choices in the benchmark, but do not make them default unless the user explicitly permits Python/reticulate or asks for non-Seurat-native methods.
- Do not integrate away biology: if sample type, tissue site, species, tumor/normal status, or disease group is the main biological contrast, do not set it as the only batch variable unless there are independent technical batch labels.

## Metadata retrieval requirements

Before writing code, inspect the GEO page, SRA Selector page, and the linked article if necessary. Use these canonical URL patterns:

- GEO page: `https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=<GSEID>`
- SRA Selector page: `https://www.ncbi.nlm.nih.gov/Traces/study/?acc=<GSEID>` or `https://www.ncbi.nlm.nih.gov/Traces/study/?acc=<SRP/PRJNA>`

Extract or infer:

- `project_GEOid`: GSE accession, e.g. `GSE123456`.
- `cancer_type`: TCGA tumor abbreviation, e.g. `PAAD`, `LUAD`, `BRCA`, `LIHC`. Use `NA` if not cancer or not inferable.
- `platform`: sequencing platform or library protocol, e.g. `10x Genomics Chromium 3'`, `Smart-seq2`, `Drop-seq`, `BD Rhapsody`, `snRNA-seq`, or the most specific available value.
- `patient_id`: patient/donor identifier. Use `NA` if absent.
- `sample_id`: GSM/sample identifier or author-provided sample ID. Use `NA` if absent.
- `sample_type`: `T` for tumor/primary tumor/metastasis tumor, `N` for normal/adjacent normal/non-tumor, `L` for lymph node/lymphatic tissue. Use `NA` if absent or ambiguous.

The generated R script must define:

```r
meta <- c(
  project_GEOid = "GSEXXXXXX",
  cancer_type   = "PAAD",
  platform      = "10x Genomics Chromium 3'",
  patient_id    = "P01;P02;P03",
  sample_id     = "GSMXXXXXXX;GSMXXXXXXY;GSMXXXXXXZ",
  sample_type   = "T;N;L"
)
```

For multi-sample projects, also define:

```r
sample_meta <- data.frame(
  sample_id = c("GSM..."),
  patient_id = c("P01"),
  sample_type = c("T"),
  sra_run = c("SRR..."),
  platform = c("10x Genomics Chromium 3'"),
  stringsAsFactors = FALSE
)
```

If a field cannot be found, set it exactly to `"NA"`, not an empty string.

## Data input handling

Tailor the code to the files that GEO/SRA provides. Do not output a generic reader if the GEO supplementary files clearly indicate the format.

## Download requirements

When the generated script needs to download public GEO/SRA supplementary files, enforce the following rules:

- Do not use `download.file()`, `curl::curl_download()`, `wget`, or other R-native download helpers.
- Use `axel` through `system2()` for all public file downloads.
- Fix the download concurrency at 10 connections.
- Fix the proxy to `http://www.cirno999.cn:12306`.
- Do not add logic that checks whether `axel` is installed. Assume the runtime environment already provides it.
- If `axel` returns a non-zero status or the target file is still missing after download, stop with a clear error message.

Use an explicit helper similar to:

```r
axel_proxy <- "http://www.cirno999.cn:12306"
axel_connections <- 10

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
```

Common cases:

1. 10x `filtered_feature_bc_matrix` folder: use `Read10X()`.
2. 10x `.h5`: use `Read10X_h5()`.
3. Matrix Market `.mtx.gz` + `barcodes.tsv.gz` + `features.tsv.gz` / `genes.tsv.gz`: use `ReadMtx()` or manually read Matrix Market if sample-specific prefixes require it.
4. Dense tabular count matrix `.csv.gz`, `.tsv.gz`, `.txt.gz`: use `data.table::fread()`, set gene symbols as row names, convert to sparse matrix with `Matrix::Matrix(..., sparse = TRUE)`.
5. Multiple sample matrices in one file: split columns by sample prefix if documented by GEO metadata.
6. FASTQ-only datasets: generate an upstream note and a separate alignment/counting placeholder. Do not pretend Seurat can directly read FASTQ. The downstream Seurat code should begin from Cell Ranger / STARsolo / kallisto-bustools / alevin-fry count matrices.

Cell barcodes must be made globally unique using sample IDs, e.g. `RenameCells(obj, add.cell.id = sample_id)` or `colnames(counts) <- paste(sample_id, colnames(counts), sep = "_")` before merging.

## R script structure to generate

Generate the R script in this order:

1. Header: dataset name, source URLs, article/PMID if used, and short assumptions.
2. Package block with `require_or_stop()` and optional package messages.
3. User parameter block, including `work_dir`, `seed`, `setwd(work_dir)`, and `set.seed(seed)`.
4. `meta` vector and `sample_meta` table.
5. Directory creation and qs path definitions.
6. Dataset-tailored count readers.
7. Seurat object construction per sample.
8. RNA contamination correction.
9. QC metrics and filtering.
10. DoubletFinder workflow on split samples.
11. Merge singlets and rerun the standard preprocessing workflow.
12. Merge/integration method 1: merge-only baseline.
13. Merge/integration method 2: Seurat v5 RPCA integration.
14. Merge/integration method 3: Harmony integration when available; otherwise explicit fallback message and skip.
15. Save each object with `qsave()`.
16. Build `seurat_results <- list(...)` and save it with `qsave()`.
17. Diagnostic plots saved to `output_dir/figures` with `ggsave()`.
18. Final validation checks.

## RNA contamination correction rules

Include this step before QC filtering.

- If 10x raw + filtered matrices are both available, prefer SoupX.
- If only filtered count matrices are available, include DecontX as optional fallback if `celda` and `SingleCellExperiment` are installed.
- If neither method is available or appropriate, do not fail silently. Set `contamination_method <- "none"` or print a clear message. The pipeline must still run.
- Save both pre- and post-contamination-correction objects or count matrices with qs when feasible.

Use comments to tell the user that contamination correction can alter low-expression markers and should be checked against canonical marker genes.

## QC defaults

Use explicit, editable defaults unless the GEO article gives specific values:

```r
min_features <- 200
max_features <- 6000
max_counts   <- Inf
max_percent_mt <- 20
min_cells_per_gene <- 3
n_hvg <- 3000
dims_use <- NULL
npcs <- 40
vars_to_regress <- c("nFeature_RNA", "nCount_RNA", "percent.mt")
cluster_resolution <- 0.5
```

For human datasets, mitochondrial genes generally match `^MT-`; for mouse, generally `^mt-`. Infer species if possible; otherwise expose `mt_pattern` in the parameter block.

## Standard preprocessing defaults

Unless the article gives a justified dataset-specific alternative, use the following preprocessing calls and parameters in generated code:

```r
obj <- NormalizeData(obj, normalization.method = "LogNormalize",
                     scale.factor = 10000, margin = 1, assay = "RNA")
obj <- FindVariableFeatures(obj, nfeatures = 3000)
obj <- ScaleData(obj, vars.to.regress = c("nFeature_RNA", "nCount_RNA", "percent.mt"))
obj <- RunPCA(obj, features = VariableFeatures(obj), npcs = 40)
```

Keep these defaults explicit in the script instead of hiding them behind vague wrappers.

## Dimension selection helper

Create a helper function that follows this logic. When `dims_use` is `NULL`, use this function to determine the effective dimensions for the relevant reduction:

```r
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
```

Also add a small resolver such as `resolve_dims_use()` so the generated script can do `seq_len(calculate_min_pc(...))` when `dims_use` is not provided.

## DoubletFinder rules

When the dataset has sample-level structure, generated scripts should remove doublets after QC and before the final merged/integration workflows. Use `sample` as the split field for this step.

- Split the post-QC object with `SplitObject(obj, split.by = "sample")`.
- Run DoubletFinder independently for each split sample, then extract only cells classified as `Singlet`.
- Merge the singlet-only objects back together, call `JoinLayers()` when appropriate, and rerun the standard preprocessing workflow on the merged singlet object before any final downstream Harmony/RPCA/baseline analyses.
- Keep this logic explicit in the script instead of hiding it in prose comments.

Use the following parameterization and workflow as the default:

```r
multiplet_rates_10x <- data.frame(
  Multiplet_rate = c(0.004, 0.008, 0.0160, 0.023, 0.031, 0.039, 0.046, 0.054, 0.061, 0.069, 0.076),
  Loaded_cells = c(800, 1600, 3200, 4800, 6400, 8000, 9600, 11200, 12800, 14400, 16000),
  Recovered_cells = c(500, 1000, 2000, 3000, 4000, 5000, 6000, 7000, 8000, 9000, 10000)
)

sample_obj <- SplitObject(obj, split.by = "sample")

h <- lapply(names(sample_obj), function(n) {
  o <- sample_obj[[n]]
  o <- NormalizeData(o, normalization.method = "LogNormalize",
                     scale.factor = 10000, margin = 1, assay = "RNA")
  o <- FindVariableFeatures(o, nfeatures = 3000)
  o <- ScaleData(o, vars.to.regress = c("nFeature_RNA", "nCount_RNA", "percent.mt"))
  o <- RunPCA(o, features = VariableFeatures(o), npcs = 30)
  min_pc <- calculate_min_pc(o, reduction = "pca")
  o <- RunUMAP(o, reduction = "pca", dims = 1:min_pc)
  o <- FindNeighbors(o, reduction = "pca", dims = 1:min_pc)
  o <- FindClusters(o, resolution = 0.1)
  sweep.res.hcc <- DoubletFinder::paramSweep(o, PCs = 1:min_pc, sct = FALSE)
  sweep.stats <- DoubletFinder::summarizeSweep(sweep.res.hcc, GT = FALSE)
  bcmvn <- DoubletFinder::find.pK(sweep.stats)
  pK_bcmvn <- as.numeric(as.vector(bcmvn$pK[which.max(bcmvn$BCmetric)]))
  homotypic.prop <- DoubletFinder::modelHomotypic(o$seurat_clusters)
  multiplet_rate <- multiplet_rates_10x %>%
    dplyr::filter(Recovered_cells < nrow(o@meta.data)) %>%
    dplyr::slice(which.max(Recovered_cells)) %>%
    dplyr::pull(Multiplet_rate)
  nExp.poi <- round(multiplet_rate * nrow(o@meta.data))
  nExp.poi.adj <- round(nExp.poi * (1 - homotypic.prop))
  DoubletFinder::doubletFinder(
    seu = o,
    PCs = 1:min_pc,
    pN = 0.25,
    pK = pK_bcmvn,
    nExp = nExp.poi.adj,
    sct = FALSE
  )
})

h.singlet <- lapply(names(sample_obj), function(n) {
  o <- h[[n]]
  colnames(o@meta.data)[grepl("DF.classifications", colnames(o@meta.data))] <- "doublet_finder"
  subset(o, subset = doublet_finder == "Singlet")
})

obj_singlet_merged_raw <- if (length(h.singlet) == 1) {
  h.singlet[[1]]
} else {
  merge(h.singlet[[1]], y = h.singlet[-1], add.cell.ids = names(h.singlet))
}
obj_singlet_merged_raw <- JoinLayers(obj_singlet_merged_raw)
cleanup_vars(c("h", "h.singlet", "sample_obj"), envir = environment())
obj_singlet_merged <- preprocess_standard(obj_singlet_merged_raw)
```

Additional requirements:

- Add `DoubletFinder` to the package checks whenever this workflow is emitted.
- Keep `pN = 0.25`, `sct = FALSE`, per-sample `RunPCA(..., npcs = 30)`, and `FindClusters(resolution = 0.1)` as the default DoubletFinder settings unless the source article justifies a change.
- Rename the detected classification column to a stable name such as `doublet_finder` before subsetting singlets.
- If only one sample exists, still run the same DoubletFinder logic on that single sample object.
- Use Chinese comments to explain why doublet removal is done per sample before the final merged analysis.

## Memory management requirements

Generated scripts should actively control memory usage instead of waiting until the end.

- Add a small helper such as `cleanup_vars()` that removes named objects if they exist and then calls `gc()`.
- After finishing a major section, clear temporary matrices, temporary Seurat objects, temporary dimension vectors, and other no-longer-needed intermediates.
- After per-sample objects have been merged into downstream objects and are no longer needed, clear `object_list` and any upstream raw-count containers when safe.
- After creating the final `seurat_results` object and saving it, clear redundant standalone bindings that are already captured inside `seurat_results` when that does not break the remaining validation code.
- Use concise Chinese comments to explain why a cleanup step is safe.

## Seurat v5 integration templates

Default baseline:

```r
obj_merged_no_correction <- obj_singlet_merged
dims_pca <- if (is.null(dims_use)) seq_len(calculate_min_pc(obj_merged_no_correction, reduction = "pca")) else dims_use
obj_merged_no_correction <- RunUMAP(obj_merged_no_correction, reduction = "pca", dims = dims_pca, reduction.name = "umap.unintegrated")
obj_merged_no_correction <- FindNeighbors(obj_merged_no_correction, reduction = "pca", dims = dims_pca)
obj_merged_no_correction <- FindClusters(obj_merged_no_correction, resolution = cluster_resolution)
```

Default Seurat v5 RPCA integration:

```r
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
dims_rpca <- if (is.null(dims_use)) seq_len(calculate_min_pc(obj_rpca, reduction = "integrated.rpca")) else dims_use
obj_rpca <- FindNeighbors(obj_rpca, reduction = "integrated.rpca", dims = dims_rpca)
obj_rpca <- RunUMAP(obj_rpca, reduction = "integrated.rpca", dims = dims_rpca, reduction.name = "umap.rpca")
obj_rpca <- FindClusters(obj_rpca, resolution = cluster_resolution)
```

Default Harmony integration:

```r
obj_harmony <- obj_singlet_merged_raw
obj_harmony$sample <- obj_harmony$sample_id
obj_harmony <- preprocess_standard(obj_harmony)
obj_harmony <- RunHarmony(obj_harmony, "sample")
dims_harmony <- if (is.null(dims_use)) seq_len(calculate_min_pc(obj_harmony, reduction = "harmony")) else dims_use
obj_harmony <- FindNeighbors(obj_harmony, reduction = "harmony", dims = dims_harmony)
obj_harmony <- RunUMAP(obj_harmony, reduction = "harmony", dims = dims_harmony, reduction.name = "umap.harmony")
obj_harmony <- FindClusters(obj_harmony, resolution = cluster_resolution)
```

If `RunHarmony()` is unavailable in the user's installed Seurat/Harmony setup, generate a clear message and skip the Harmony object rather than failing after the RPCA object has been created.

## Output object contract

The final generated script must end with a result list similar to:

```r
seurat_results <- list(
  merged_no_correction = obj_merged_no_correction,
  integrated_rpca = obj_rpca,
  integrated_harmony = obj_harmony
)

qsave(seurat_results, file.path(output_dir, paste0(project_id, "_seurat_results.qs")))
```

If one optional method was skipped, omit it from the list using `Filter(Negate(is.null), ...)` but keep at least two Seurat objects in the list.

## Readability requirements

- Use section dividers such as `# ---- 01. packages ----`.
- Use small helper functions with explicit names.
- Do not compress the workflow into one long pipe.
- Add comments explaining why each major step exists, and use concise Chinese comments where a human reader is likely to need clarification.
- Insert explicit memory-cleanup blocks between major sections rather than leaving all intermediates in memory.
- Keep all key thresholds in the parameter block; do not hard-code them in functions.
- Never overwrite output files without making it clear where they are saved.

## Safety and reproducibility requirements

- Use both `setwd(work_dir)` and `set.seed(seed)`.
- Use `rm()` plus `gc()` to clear large temporary objects after they are no longer needed.
- Print `sessionInfo()` at the end and save it with qs or text output.
- Do not auto-install packages unless the user explicitly asked for an installation script. Instead, stop with a clear package list.
- Do not download controlled-access data.
- Do not run shell commands that delete user files.
- For large datasets, include `future::plan()` only if `future` is installed, and expose `future_workers` / `future_max_size` as parameters.

## When uncertain

If metadata or file format is ambiguous after searching GEO/SRA/article sources, do not fabricate. Use `"NA"` for missing metadata and add a short `assumptions` comment in the R script explaining exactly what is uncertain and what the user should verify.
