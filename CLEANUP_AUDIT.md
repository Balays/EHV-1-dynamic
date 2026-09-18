# Repository cleanup audit

This audit maps the repository to the analyses described in the EHV-1 dynamic-transcriptome manuscript and distinguishes reproducibility material from working files and generated output.

## Keep: directly relevant to the manuscript

### Core workflow

- `EHV-1.wf.LoRTIA.R`
- `_WF.part0.R`
- `_WF.part1.R`
- `project_config.txt`
- `functions/` after confirming that every helper is actually called

### CAGE integration and transcript annotation

- `CAGE.refgenes.Rmd`
- `Novel_TXs_from_CAGE_and_transfrags.Rmd`
- `Novel_Spliced_TXs_from_CAGE_and_transfrags.Rmd`
- `CAGEFighteR.SIG.tsv`
- final transcript GFF3/TSV files used in Supplementary Table S4
- compact CAGE and transfrag tables required as notebook inputs

### dRNA-Seq/NAGATA validation

- `NAGATA.R`
- `NAGATA_dRNA.gff3`
- `EHV_novel_NAGATA_merged.gff3`
- `EHV_dRNA_minion_20210129_pass_mapped_intron.gff3`

The complete NAGATA source tree should not be vendored here; the repository should instead document the exact upstream version and command.

### Kinetic analyses

- `TSS.only_VirusRead.norm_hclust.Rmd`
- `TES.only_VirusRead.norm_hclust.Rmd`
- `TSS.TES_VirusRead.norm_hclust.Rmd`
- `All_Methods_Combined.Rmd`
- clustering helpers that are demonstrably sourced by these notebooks

### Isoform, splicing and overlap analyses

- `TR.ratios.Rmd`
- `Transcript.Isoform.Ratios.refgenes.Rmd`
- `prime.plots.Rmd`
- `coverage_virus.R`
- `coverage_virus_plots.R`
- final compact tables used to construct manuscript figures

### CHX analysis

- `Gene.counts.R`
- the compact CHX count table used for Supplementary Table S3

### Reference and metadata

- `EHV-1_metadata.tsv`
- `NC_001491.2.fasta`
- one canonical EHV-1 gene annotation GFF3
- one canonical prior-transcript annotation GFF3/TSV
- PRV reference sequences/annotations only where required for the OriL/OriS comparison

## Remove: unequivocal repository clutter

- manuscript drafts and response-letter files under `Article/`
- knitted `*.html` files
- R Markdown cache and `*_files/figure-html/` directories
- raw or mapped sequencing data (`*.bam`, `*.bai`, FASTQ, POD5, FAST5)
- hundreds of automatically generated exploratory JPG files under `LoRTIA_virus/`
- copies of external software such as the `NAGATA/` source tree
- local IDE state and session files
- diagnostic plots not used in the paper

Raw sequencing data belong in ENA; large derived files that are not practical to regenerate should be deposited in a data archive rather than in the code repository.

## Consolidate before deleting

The following groups contain multiple versions or overlapping implementations and should be reduced only after tracing their consumers:

- `CAGEFighteR.SIG.tsv` and `CAGEFighteR.SIG.v1.tsv`
- the many `NC_001491.2.mRNA_corrected*` GFF2/GFF3/GTF files
- `TR.all.data.sp*.tsv`
- `final_tx_nov*.gff3`
- `customclust.R`, `dbscan.R`, `mclust.R`, `mcluster.R`, `cluster.R`, `hcluster.R`, and `hclust.wss.R`
- `All_Methods_Combined.Rmd` versus the three endpoint-specific clustering notebooks
- `TR.ratios.Rmd` versus `Transcript.Isoform.Ratios.refgenes.Rmd`

For each group, retain the version read by the final manuscript workflow and move superseded versions to a tagged archival release if they are scientifically useful.

## Main reproducibility blockers

1. `_WF.part0.R` loads helper functions and lookup tables through absolute local paths:
   - `C:/GitHub/Rlyeh/R`
   - `C:/GitHub/minitax/R`
2. The workflow does not provide a locked R package environment.
3. Execution order is implicit and distributed across notebooks.
4. Many notebooks contain disabled, exploratory, or duplicated code paths.
5. Inputs, derived data, final tables, manuscript figures and scratch output are mixed together.
6. Some repository filenames and notebook titles do not describe their current function.

## Recommended final structure

```text
EHV-1-dynamic/
├── README.md
├── LICENSE
├── renv.lock
├── analysis/
│   ├── 01_build_transfrags.R
│   ├── 02_integrate_cage.R
│   ├── 03_count_reference_transcripts.R
│   ├── 04_assemble_novel_transcripts.R
│   ├── 05_validate_drna_nagata.R
│   ├── 06_cluster_tss_tes_transcripts.R
│   ├── 07_isoform_dynamics.R
│   └── 08_chx_analysis.R
├── R/
│   └── project-specific helper functions
├── config/
│   └── project.yml
├── data/
│   ├── reference/
│   ├── metadata/
│   └── processed/          # compact, publication-facing inputs only
├── results/
│   └── manuscript_tables/  # small final TSV/GFF3 outputs
└── figures/
    └── manuscript/         # optional final figures only
```

## Recommended next technical pass

1. Identify every function imported from `Rlyeh` and `minitax` and vendor only those project-specific helpers under `R/`.
2. Replace absolute paths with paths relative to the repository root.
3. Add an `renv.lock` file and a short installation command.
4. Convert the final workflow into numbered scripts or a single workflow manager such as `targets`.
5. Move compact publication inputs to `data/processed/` and regenerate all plots into an ignored `outputs/` directory.
6. Remove large tracked files from Git history using `git filter-repo` after the cleaned branch is validated.

## Conservative cleanup policy

The cleanup branch removes only files that cannot affect the scientific results: manuscript drafts, knitted reports and standalone diagnostic images. Raw BAM removal and consolidation of competing annotation/script versions should follow after the retained workflow has been executed successfully from a fresh environment.
