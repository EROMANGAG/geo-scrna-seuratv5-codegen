# GEO/SRA metadata extraction checklist

Before writing code, collect the following and preserve the source URLs in script comments.

## GEO page

Open `https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=<GSEID>` and collect:

- Title and summary
- Overall design
- Organism
- Platform(s), e.g. GPL IDs and sequencing technology text
- Supplementary file names and URLs
- GSM sample titles and characteristics
- Related BioProject/SRA/PubMed links

## SRA Selector

Open `https://www.ncbi.nlm.nih.gov/Traces/study/?acc=<GSEID>` and, when available, download or inspect RunInfo / RunTable fields:

- Run accession (`SRR...`)
- BioSample
- sample alias / title
- LibraryStrategy, LibrarySource, LibrarySelection
- Platform / Instrument
- sample attributes

## Article / PubMed

Use the linked article when GEO/SRA lacks patient/sample annotations. Extract:

- Cancer type and cohort description
- Patient/donor IDs
- Tissue/sample type: tumor, normal, lymph node, metastatic lesion
- Single-cell technology and count-generation method
- QC criteria reported by authors

## Coding rule

If metadata remains absent, write `"NA"`. Never infer patient IDs or sample types from column order alone.
