---
name: geo-scrna-seuratv5-codegen
description: Generate a high-readability Seurat v5 R or R Markdown (Rmd) pipeline from a scRNA-seq GEO accession, GEO URL, or SRA Selector URL, including GEO/SRA/PubMed metadata extraction, named meta fields, RNA contamination correction, QC, DoubletFinder-based doublet removal, modular qread/qsave state handoff, selectable RPCA/Harmony/BBKNN integration outputs saved with qs, and an optional SingleR-based initial annotation module for res0.1 clustering results.
---

# GEO scRNA-seq to Seurat v5 code generator

Use this skill when the user gives a GEO accession, GEO page URL, SRA Selector URL, SRP/PRJNA/SRA accession, or article link and asks for Seurat v5 code to process a single-cell RNA-seq dataset into Seurat objects.

Use this skill for both `.R` and `.Rmd` outputs.

Do not use this skill for bulk RNA-seq, spatial transcriptomics, scATAC-only datasets, or non-R pipelines unless the user explicitly asks to adapt the workflow.

## Required user inputs before generation

Before generating code, require the user to provide:

1. A working directory path that will be written into the code and applied with `setwd(work_dir)`.
2. A numeric seed value that will be written into the code and applied with `set.seed(seed)`.
3. An output format choice: `r` or `rmd`.

If the user has not clearly chosen `r` or `rmd`, ask them to confirm the target format before generating code.

If the user explicitly asks for both, generate both templates.

Do not silently assume a default path, seed, or output format.

## Required final deliverable

Produce a complete, directly runnable Seurat v5 pipeline tailored to the dataset. The generated code must:

1. Use **Seurat v5** as the main analysis framework.
2. Use **qs** for all generated intermediate and final R objects. Use `qsave()` / `qread()`. Do not use `saveRDS()` / `readRDS()` for pipeline state handoff.
3. Return and save a named list called `seurat_results` whose elements are Seurat objects.
4. Include separate, user-selectable downstream modules for:
   - `merged_no_correction`
   - `integrated_rpca`
   - `integrated_harmony`
   - `integrated_bbknn`
   - `singleR_initial_annotation`
5. Make the downstream integration modules independent, so the user can run any one or multiple modules without them affecting each other.
6. Include a readable parameter block so the user can easily modify `work_dir`, `seed`, QC cutoffs, dimensions, resolution, contamination method, and integration settings.
7. Include a named R vector called `meta`, and when multiple samples exist also include `sample_meta` as a data frame.
8. Implement the full processing flow: RNA contamination correction -> basic QC -> DoubletFinder-based doublet removal -> normalization -> HVG selection -> scaling -> PCA -> UMAP -> neighbors -> clustering.
9. Include verification code at the end: `stopifnot(inherits(..., "Seurat"))`, `print(seurat_results)`, and checks that generated `.qs` outputs exist.
10. Add concise Chinese comments at key analysis steps, helper functions, and any place where parameter choices or assumptions need explanation.
11. Remove large temporary variables after each major module and run `gc()` to keep memory usage low.
12. If the script needs to download public supplementary files, do not use `download.file()` or other R built-in download helpers. Use `axel` with 10 connections and the fixed proxy `http://www.cirno999.cn:12306`.

## Output format and modular execution requirements

### R output requirements

For `.R` output:

- Use section dividers such as `# ---- 01. environment ----`.
- Treat the script as a stepwise modular pipeline rather than one long sequential run.
- The first environment/setup section should define paths, parameters, helper functions, package checks, and `state_paths`.
- The first raw-file-processing module may read directly from raw files and write the first `.qs` state.
- Every later module must:
  1. Load its required upstream state with `qread()`.
  2. Run only the logic for that module.
  3. Save its output state with `qsave()`.
  4. Clean up temporary objects with `cleanup_vars()` / `gc()`.
- Use clear short titles before each module so the user can quickly jump to the relevant section.

### Rmd output requirements

For `.Rmd` output:

- Start with a YAML header and a setup chunk.
- Use one dedicated code chunk for environment setup: paths, parameters, helper functions, package checks, and `state_paths`.
- Give every later analysis stage its own Markdown heading and its own code chunk.
- Except for the initial raw-file-processing chunk, every later chunk must contain complete `qread()` input and `qsave()` output logic.
- Include clear short Markdown titles before each chunk so the user can quickly move through the notebook.
- Include at least one optional Python chunk placeholder with `eval=FALSE` so the user can add Python code when needed.
- Make the RPCA, Harmony, and BBKNN chunks independent so the user can run any one or multiple integration chunks without cross-dependency.

### Shared modular state requirements

Both `.R` and `.Rmd` outputs must define a `state_paths` structure similar to:

```r
state_dir <- file.path(output_dir, "states")
annotation_dir <- file.path(output_dir, "annotations")
metrics_dir <- file.path(output_dir, "metrics")
state_paths <- list(
  object_list_raw = file.path(state_dir, paste0(project_id, "_01_object_list_raw.qs")),
  object_list_contam = file.path(state_dir, paste0(project_id, "_02_object_list_contam.qs")),
  obj_qc = file.path(state_dir, paste0(project_id, "_03_obj_qc.qs")),
  obj_singlet_merged_raw = file.path(state_dir, paste0(project_id, "_04_obj_singlet_merged_raw.qs")),
  obj_singlet_merged = file.path(state_dir, paste0(project_id, "_05_obj_singlet_merged.qs")),
  doubletfinder_cell_counts_csv = file.path(metrics_dir, paste0(project_id, "_doubletfinder_cell_counts.csv")),
  merged_no_correction = file.path(state_dir, paste0(project_id, "_06_merged_no_correction.qs")),
  integrated_rpca = file.path(state_dir, paste0(project_id, "_07_integrated_rpca.qs")),
  integrated_harmony = file.path(state_dir, paste0(project_id, "_08_integrated_harmony.qs")),
  integrated_bbknn = file.path(state_dir, paste0(project_id, "_09_integrated_bbknn.qs")),
  seurat_results = file.path(state_dir, paste0(project_id, "_10_seurat_results.qs")),
  annotation_dir = annotation_dir
)
```

## Required literature-derived integration logic

Use the user's uploaded benchmark article as the method-selection reference. Encode these principles in comments and in the method choices:

- Always generate an unintegrated / merge-only baseline object. It is required for detecting overcorrection and for inspecting whether biological differences are being erased.
- Use highly variable gene selection before integration by default. The benchmark reports that HVG selection generally improves integration performance.
- Treat scaling as useful but potentially overcorrecting: scaling can improve batch removal but may reduce biological signal conservation. Keep parameters explicit and avoid hiding this choice.
- Provide separate downstream modules for `RPCA`, `Harmony`, and `BBKNN` so the user can compare different correction strategies against the merge-only baseline.
- For small/simple data with distinct biological signal, Harmony is acceptable and often convenient.
- For strong or nested batch effects, instruct the user to compare RPCA, Harmony, and BBKNN against the merge-only baseline rather than trusting a single corrected view.
- Mention in comments that Scanorama, scVI/scANVI are strong choices in some benchmarks, but do not make them default unless the user explicitly permits Python/reticulate or asks for non-Seurat-native methods.
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

The generated code must define:

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

When the generated code needs to download public GEO/SRA supplementary files, enforce the following rules:

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
check_existing_downloads <- TRUE

download_with_axel <- function(url, dest_path, check_existing = check_existing_downloads) {
  dir.create(dirname(dest_path), recursive = TRUE, showWarnings = FALSE)
  if (isTRUE(check_existing) && file.exists(dest_path)) {
    return(dest_path)
  }
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

Expose `check_existing_downloads <- TRUE` in the parameter block. When it is `TRUE`, the generated code should skip re-downloading a public file if the target path already exists. When it is `FALSE`, allow the script to download again and overwrite the target file.

Common cases:

1. 10x `filtered_feature_bc_matrix` folder: use `Read10X()`.
2. 10x `.h5`: use `Read10X_h5()`.
3. Matrix Market `.mtx.gz` + `barcodes.tsv.gz` + `features.tsv.gz` / `genes.tsv.gz`: use `ReadMtx()` or manually read Matrix Market if sample-specific prefixes require it.
4. Dense tabular count matrix `.csv.gz`, `.tsv.gz`, `.txt.gz`: use `data.table::fread()`, set gene symbols as row names, convert to sparse matrix with `Matrix::Matrix(..., sparse = TRUE)`.
5. Multiple sample matrices in one file: split columns by sample prefix if documented by GEO metadata.
6. FASTQ-only datasets: generate an upstream note and a separate alignment/counting placeholder. Do not pretend Seurat can directly read FASTQ. The downstream Seurat code should begin from Cell Ranger / STARsolo / kallisto-bustools / alevin-fry count matrices.

Cell barcodes must be made globally unique using sample IDs, e.g. `RenameCells(obj, add.cell.id = sample_id)` or `colnames(counts) <- paste(sample_id, colnames(counts), sep = "_")` before merging.

## Pipeline structure to generate

Generate the pipeline in this order:

1. Header: dataset name, source URLs, article/PMID if used, and short assumptions.
2. Package block with `require_or_stop()`.
3. Environment block: `work_dir`, `seed`, `setwd(work_dir)`, `set.seed(seed)`, path creation, helper functions, and `state_paths`.
4. Raw input processing module.
5. RNA contamination correction module.
6. QC module.
7. DoubletFinder module.
8. Re-preprocessing module for merged singlets.
9. Merge-only baseline module.
10. RPCA module.
11. Harmony module.
12. BBKNN module.
13. SingleR initial annotation module.
14. Optional Python placeholder chunk for `.Rmd`.
15. Final results assembly module.
16. Diagnostic plot module.
17. Final validation module.

## RNA contamination correction rules

Include this step before QC filtering.

- If 10x raw + filtered matrices are both available, prefer SoupX.
- If only filtered count matrices are available, include DecontX as optional fallback if `celda` and `SingleCellExperiment` are installed.
- If neither method is available or appropriate, do not fail silently. Set `contamination_method <- "none"` or print a clear message. The pipeline must still run.
- Save both pre- and post-contamination-correction objects or count matrices with qs when feasible.
- For the DecontX branch, follow this structure and do not use `as.SingleCellExperiment()` there:

```r
corrected <- lapply(names(object_list), function(n) {
  obj <- object_list[[n]]
  sce <- SingleCellExperiment::SingleCellExperiment(
    assays = list(counts = LayerData(obj, assay = "RNA", layer = "counts"))
  )
  sce <- celda::decontX(sce)
  corrected_counts <- celda::decontXcounts(sce)
  obj[["RNA"]] <- SeuratObject::CreateAssay5Object(counts = corrected_counts)
  obj <- add_sample_metadata(obj, sample_meta[sample_meta$sample_id == n, , drop = FALSE][1, ])
  cleanup_vars(c("sce", "corrected_counts"))
  obj
})
```

- Keep the DecontX implementation at the `object_list` level, because metadata restoration depends on the sample name.
- Add a small helper such as `add_sample_metadata()` so each corrected object gets the expected sample metadata back after the RNA assay is replaced.

Use comments to tell the user that contamination correction can alter low-expression markers and should be checked against canonical marker genes.

## Seurat v5 compatibility rules

Generated code must target the validated package baseline below, derived from the user's reference `sessionInfo()`.

```r
package_version_requirements <- c(
  Seurat = "5.5.0",
  SeuratObject = "5.4.0",
  qs = "0.27.3",
  ggplot2 = "4.0.3",
  dplyr = "1.2.1",
  data.table = "1.17.8",
  Matrix = "1.7-5",
  DoubletFinder = "2.0.6",
  harmony = "2.0.2",
  SingleR = "2.8.0",
  SingleCellExperiment = "1.28.0"
)
```

At minimum, generated code should check the versions of the main packages it actually uses and stop early if the installed version is below the validated baseline.

- Add an explicit version-check helper near the package block and use it against the validated baseline.
- Do not generate `FetchData(..., slot = ...)`; always use `FetchData(..., layer = ..., assay = ...)` when expression data must be fetched.
- Do not generate `GetAssayData(..., slot = ...)` or `SetAssayData(..., slot = ...)`.
- Prefer `LayerData()` and `LayerData<-` for expression matrix access in Seurat v5 style code.
- Do not generate `CreateAssayObject()` fallback branches for legacy Seurat versions; use `SeuratObject::CreateAssay5Object()` directly.
- If a helper or code block was originally written for Seurat v4 or older, rewrite it to v5-compatible layer semantics before emitting the final script.
- Keep optional-module checks aligned with the same baseline; for example, `harmony`, `SingleR`, and `SingleCellExperiment` should also be checked before their modules run.

## QC defaults

Use explicit, editable defaults unless the GEO article gives specific values:

```r
min_features <- 200
max_features <- 6000
max_counts <- Inf
max_percent_mt <- 20
min_cells_per_gene <- 3
n_hvg <- 3000
dims_use <- NULL
npcs <- 40
vars_to_regress <- c("nFeature_RNA", "nCount_RNA", "percent.mt")
cluster_resolutions <- c(0.1, 0.5)
```

For human datasets, mitochondrial genes generally match `^MT-`; for mouse, generally `^mt-`. Infer species if possible; otherwise expose `mt_pattern` in the parameter block.

Keep `0.1` inside `cluster_resolutions`, because the SingleR initial annotation module depends on `res0.1` cluster labels.

## Standard preprocessing defaults

Unless the article gives a justified dataset-specific alternative, use the following preprocessing calls and parameters in generated code:

```r
obj <- NormalizeData(obj, normalization.method = "LogNormalize",
                     scale.factor = 10000, margin = 1, assay = "RNA")
obj <- FindVariableFeatures(obj, nfeatures = 3000)
obj <- ScaleData(obj, vars.to.regress = c("nFeature_RNA", "nCount_RNA", "percent.mt"))
obj <- RunPCA(obj, features = VariableFeatures(obj), npcs = 40)
```

Keep these defaults explicit in the generated code.

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

Also add a small resolver such as `resolve_dims_use()` so the generated code can do `seq_len(calculate_min_pc(...))` when `dims_use` is not provided.

## DoubletFinder rules

When the dataset has sample-level structure, generated code should remove doublets after QC and before the final merged/integration workflows. Use `sample` as the split field for this step.

- Split the post-QC object with `SplitObject(obj, split.by = "sample")`.
- Run DoubletFinder independently for each split sample, then extract only cells classified as `Singlet`.
- Merge the singlet-only objects back together, call `JoinLayers()` when appropriate, and rerun the standard preprocessing workflow on the merged singlet object before any final downstream analyses.
- Add `DoubletFinder` to the package checks whenever this workflow is emitted.
- Keep `pN = 0.25`, `sct = FALSE`, per-sample `RunPCA(..., npcs = 30)`, and `FindClusters(resolution = 0.1)` as the default DoubletFinder settings unless the source article justifies a change.
- Rename the detected classification column to a stable name such as `doublet_finder` before subsetting singlets.
- If only one sample exists, still run the same DoubletFinder logic on that single sample object.
- Use Chinese comments to explain why doublet removal is done per sample before the final merged analysis.
- Print the per-sample cell counts before DoubletFinder and after DoubletFinder.
- Save a stable per-sample summary table that includes `sample`, `cells_before_doubletfinder`, `cells_after_doubletfinder`, and `cells_removed_by_doubletfinder`.
- Save that summary table to a stable path such as `output/metrics/<project_id>_doubletfinder_cell_counts.csv`.
- Even if `run_doubletfinder <- FALSE`, still emit the same summary table with identical before/after counts so the output shape stays stable.

## SingleR rules

Add an independent SingleR initial annotation module for `res0.1` clustering results.

- Let the user choose the upstream result via a simple setting such as `singleR_annotation_input_key`.
- Add a nearby Chinese comment that explicitly lists the allowed values:
  - `merged_no_correction`
  - `integrated_rpca`
  - `integrated_harmony`
  - `integrated_bbknn`
- Save SingleR outputs into `output_dir/annotations`.
- Save both:
  1. a cluster-level annotation table, and
  2. an annotated Seurat object with per-cell propagated labels.
- Use `SingleR` cluster-level annotation on the `res0.1` clustering column rather than annotating every cell independently.
- For `merged_no_correction`, `integrated_rpca`, and `integrated_harmony`, use `RNA_snn_res.0.1` as the cluster column.
- For `integrated_bbknn`, use `bbknn_res.0.1` as the cluster column.
- If the chosen upstream object does not contain the required `res0.1` cluster column, stop with a clear error message instead of guessing.
- Use stable output naming derived from `project_id` and `singleR_annotation_input_key`, so different projects keep the same output structure.

Use a template similar to:

```r
# 中文注释：可选上游结果有 merged_no_correction / integrated_rpca / integrated_harmony / integrated_bbknn
singleR_annotation_input_key <- "integrated_harmony"
singleR_label_field <- "label.main"
singleR_reference <- NULL

get_singleR_cluster_column <- function(input_key) {
  switch(
    input_key,
    merged_no_correction = "RNA_snn_res.0.1",
    integrated_rpca = "RNA_snn_res.0.1",
    integrated_harmony = "RNA_snn_res.0.1",
    integrated_bbknn = "bbknn_res.0.1",
    stop(sprintf("Unsupported SingleR input key: %s", input_key), call. = FALSE)
  )
}

build_singleR_output_paths <- function(input_key) {
  list(
    table_csv = file.path(annotation_dir, paste0(project_id, "_", input_key, "_SingleR_res0.1.csv")),
    annotated_qs = file.path(annotation_dir, paste0(project_id, "_", input_key, "_SingleR_annotated.qs"))
  )
}

require_or_stop("SingleR")
obj_singleR <- qread(state_paths[[singleR_annotation_input_key]])
cluster_col <- get_singleR_cluster_column(singleR_annotation_input_key)
sce_singleR <- as.SingleCellExperiment(obj_singleR, assay = "RNA")
pred_singleR <- SingleR::SingleR(
  test = sce_singleR,
  ref = singleR_reference,
  labels = SummarizedExperiment::colData(singleR_reference)[[singleR_label_field]],
  clusters = obj_singleR[[cluster_col]][, 1],
  assay.type.test = 1
)
```

Then propagate the cluster-level labels back to each cell and save both the annotation table and the annotated object.

## BBKNN rules

Add an independent BBKNN integration module.

- Load the upstream singlet-merged raw object with `qread(state_paths$obj_singlet_merged_raw)`.
- Set `obj_bbknn$sample <- obj_bbknn$sample_id`.
- Rerun the standard preprocessing workflow inside the BBKNN module so the module is independent.
- Use `bbknnR::RunBBKNN()` as the default entry point.
- Keep BBKNN package loading module-local rather than forcing it into the global package block, so users who do not run the BBKNN block are not blocked.
- Save the BBKNN result to its own `.qs` file.

Use a template similar to:

```r
require_or_stop("bbknnR")
obj_bbknn <- qread(state_paths$obj_singlet_merged_raw)
obj_bbknn$sample <- obj_bbknn$sample_id
obj_bbknn <- preprocess_standard(obj_bbknn)
dims_bbknn <- if (is.null(dims_use)) calculate_min_pc(obj_bbknn, reduction = "pca") else max(dims_use)
obj_bbknn <- bbknnR::RunBBKNN(
  object = obj_bbknn,
  batch_key = "sample",
  reduction = "pca",
  n_pcs = dims_bbknn,
  graph_name = "bbknn",
  run_TSNE = FALSE,
  run_UMAP = TRUE,
  UMAP_name = "umap",
  seed = seed,
  verbose = FALSE
)
obj_bbknn <- FindClusters(obj_bbknn, graph.name = "bbknn", resolution = cluster_resolutions)
qsave(obj_bbknn, state_paths$integrated_bbknn)
```

## Output object contract

The final results-assembly module must collect whichever downstream result files the user has actually generated.

Do not require the user to run all downstream integration modules before assembling results.

Use a pattern similar to:

```r
result_files <- c(
  merged_no_correction = state_paths$merged_no_correction,
  integrated_rpca = state_paths$integrated_rpca,
  integrated_harmony = state_paths$integrated_harmony,
  integrated_bbknn = state_paths$integrated_bbknn
)

seurat_results <- Filter(
  Negate(is.null),
  lapply(names(result_files), function(nm) {
    if (file.exists(result_files[[nm]])) qread(result_files[[nm]]) else NULL
  })
)
names(seurat_results) <- names(result_files)
seurat_results <- Filter(Negate(is.null), seurat_results)
```

If no downstream result files exist, stop with a clear message.

## Memory management requirements

Generated code should actively control memory usage instead of waiting until the end.

- Add a small helper such as `cleanup_vars()` that removes named objects if they exist and then calls `gc()`.
- After finishing a major module, clear temporary matrices, temporary Seurat objects, temporary dimension vectors, and other no-longer-needed intermediates.
- After state has been handed off with `qsave()`, clean up the in-memory copy when safe.
- After creating the final `seurat_results` object and saving it, clear redundant standalone bindings that are already captured inside `seurat_results` when that does not break the remaining validation code.
- Use concise Chinese comments to explain why a cleanup step is safe.

## Output format stability requirements

The generated pipeline must keep a stable top-level format across different projects.

- Preserve the same module order unless a module is explicitly unavailable.
- Preserve the same section titles or chunk titles across projects.
- Preserve the same parameter names, helper function names, and `state_paths` keys across projects.
- Preserve the same output directory layout: `output/`, `output/figures/`, `output/states/`, `output/annotations/`, and `output/metrics/`.
- Preserve the same object names for core states: `object_list_raw`, `object_list_contam`, `obj_qc`, `obj_singlet_merged_raw`, `obj_singlet_merged`, `merged_no_correction`, `integrated_rpca`, `integrated_harmony`, `integrated_bbknn`, `seurat_results`.
- Preserve the same annotation-setting names: `singleR_annotation_input_key`, `singleR_label_field`, `singleR_reference`.
- Preserve the same shared UMAP reduction name, `umap`, across all downstream result objects so later plotting modules can run unchanged.
- Keep `DoubletFinder` and merged-singlet re-preprocessing as two separate modules; do not merge them back into one ad hoc block.
- Do not rename modules, reorder them, or switch between ad hoc naming schemes from one project to another.

## Readability requirements

- Use short, clear titles before every section or chunk.
- Use small helper functions with explicit names.
- Do not compress the workflow into one long pipe.
- Add comments explaining why each major step exists, and use concise Chinese comments where a human reader is likely to need clarification.
- Keep all key thresholds in the parameter block; do not hard-code them inside downstream modules.
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

If metadata or file format is ambiguous after searching GEO/SRA/article sources, do not fabricate. Use `"NA"` for missing metadata and add a short `assumptions` comment in the code explaining exactly what is uncertain and what the user should verify.
