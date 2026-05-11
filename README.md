# geo-scrna-seuratv5-codegen

Codex skill for generating readable Seurat v5 `R` or `Rmd` pipelines from GEO/SRA scRNA-seq datasets.

## Install

Copy this folder to one of Codex's skill search locations, for example:

```bash
mkdir -p ~/.agents/skills
cp -r geo-scrna-seuratv5-codegen ~/.agents/skills/
```

For a repository-scoped workflow:

```bash
mkdir -p .agents/skills
cp -r geo-scrna-seuratv5-codegen .agents/skills/
```

Restart Codex if the skill does not appear.

## Example prompt

```text
Use $geo-scrna-seuratv5-codegen to generate a Seurat v5 R or Rmd pipeline for https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSEXXXXXX
```
