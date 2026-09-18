# EHV-1 dynamic transcriptome analysis

This repository contains the analysis and figure-generation code for the manuscript **“Mapping the Temporal Transcriptomic Signature of a Viral Pathogen through CAGE and Nanopore Sequencing.”**

The study combines time-resolved Oxford Nanopore direct-cDNA sequencing (dcDNA-Seq), Illumina CAGE-Seq, previously generated native direct-RNA sequencing (dRNA-Seq), and cycloheximide-treated samples to refine the EHV-1 transcript annotation and characterize transcription start sites (TSSs), transcription end sites (TESs), full-length canonical transcripts, transcript isoforms, splicing, transcriptional overlaps, and kinetic expression classes.

## Data

Sequencing data are available from the European Nucleotide Archive under:

- `PRJEB52190`
- `PRJEB6233`

Large sequencing files and alignment files should be downloaded or regenerated from the archived data and are not intended to be stored in this code repository.

## Analysis map

| Manuscript analysis | Main repository files |
|---|---|
| dcDNA-Seq import, LoRTIA-derived transfrags, adapter status and coverage | `EHV-1.wf.LoRTIA.R`, `_WF.part0.R`, `_WF.part1.R` |
| CAGE integration and confidence assignment | `CAGE.refgenes.Rmd`, `CAGEFighteR.SIG.tsv` |
| Reference-transcript counting with iterative GFF-compare matching | `EHV-1.wf.LoRTIA.R`, `import.ref.TRs.R` |
| Novel transcript assembly from refined TSSs and TESs | `Novel_TXs_from_CAGE_and_transfrags.Rmd`, `Novel_Spliced_TXs_from_CAGE_and_transfrags.Rmd` |
| dRNA-Seq validation with NAGATA | `NAGATA.R`, `NAGATA_dRNA.gff3`, `EHV_novel_NAGATA_merged.gff3` |
| TSS, TES and canonical-transcript kinetic clustering | `TSS.only_VirusRead.norm_hclust.Rmd`, `TES.only_VirusRead.norm_hclust.Rmd`, `TSS.TES_VirusRead.norm_hclust.Rmd`, `All_Methods_Combined.Rmd` |
| Transcript isoform ratios and switching | `TR.ratios.Rmd`, `Transcript.Isoform.Ratios.refgenes.Rmd` |
| CHX-treated gene counting | `Gene.counts.R` |
| Genome-wide 5′/3′-end and coverage plots | `prime.plots.Rmd`, `coverage_virus.R`, `coverage_virus_plots.R` |

## Reference files

The principal EHV-1 reference files are the `NC_001491.2` FASTA and GFF3 files. Additional GFF3/TSV files contain the previously published transcript annotation, refined transcript models, CAGE support, and NAGATA validation used in the manuscript.

## Running the analysis

The scripts were developed as an interactive R workflow rather than as a single command-line pipeline. The intended order is:

1. Prepare mapped dcDNA-Seq BAM files and run LoRTIA.
2. Build transfrag-level counts and adapter-supported 5′/3′ ends with `EHV-1.wf.LoRTIA.R`.
3. Process CAGE-Seq data and assign CAGE clusters with `CAGE.refgenes.Rmd`.
4. Count previously annotated transcripts and assemble novel transcripts.
5. Validate selected TSSs and introns against dRNA-Seq/NAGATA output.
6. Generate TSS-, TES-, transcript-, isoform- and CHX-based analyses and figures.

`project_config.txt` currently points the workflows to `LoRTIA_virus/` as their output directory.

## Reproducibility status

The scientific analysis code is present, but the repository still requires consolidation before it can be considered fully portable:

- `_WF.part0.R` contains absolute paths to helper code in local `Rlyeh` and `minitax` repositories.
- The exact R package environment is not locked with `renv`.
- Several R Markdown notebooks overlap substantially and include exploratory or disabled code.
- Generated figures, knitted HTML, raw BAM files and intermediate outputs were historically committed alongside source code.

The `paper-reproducibility-cleanup` branch documents and applies a conservative cleanup. Ambiguous scientific files are retained until their dependencies are fully resolved.

## External tools

- LoRTIA: https://github.com/zsolt-balazs/LoRTIA
- NAGATA: use the software version cited in the manuscript
- minimap2
- STAR
- GFFCompare
- CAGEfightR
- pvclust

## License

MIT License. See `LICENSE`.
