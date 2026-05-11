# Integration decision rules for generated Seurat v5 scRNA-seq pipelines

Use these rules when generating the R script.

## Default output set

Always generate at least:

1. `merged_no_correction`: raw merge baseline, processed with standard PCA/UMAP/clustering.
2. `integrated_rpca`: Seurat v5 `RPCAIntegration`.
3. `integrated_harmony`: Harmony integration generated with `RunHarmony(obj, "sample")` when the user's installed packages support it.

The baseline is not optional. It is the visual and analytic reference for overcorrection.

## Method choice rationale

- HVG selection is default because the benchmark found it generally improves RNA integration performance.
- Scaling is explicit because it can improve batch removal while reducing biological variation conservation.
- RPCA is used as a Seurat-native integration method suitable for multi-sample R workflows.
- Harmony is included because it is usable and practical for simple or moderately complex batch structure, and should be written with an explicit `sample` metadata column for `RunHarmony(obj, "sample")`.
- Scanorama/scVI/scANVI may be mentioned in script comments as benchmark-favored options, but do not implement them unless the user authorizes Python/reticulate.

## Batch variable selection

Use `sample_id` as the canonical sample identifier. Also create `sample <- sample_id` in generated objects so Harmony can run with `RunHarmony(obj, "sample")`. Use `patient_id` only when patient is the technical batch and not the biological contrast. Do not use `sample_type`, disease group, tumor/normal status, species, or tissue site as a batch variable if that is the biological contrast to preserve.

## Comparison plots

For each result object, generate UMAP plots colored by:

- `sample_id`
- `patient_id` if available
- `sample_type` if available
- `nFeature_RNA`
- `percent.mt`

These plots help distinguish batch removal from biology removal.
