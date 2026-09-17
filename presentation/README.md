# EHV-1 conference presentation

This directory contains the master conference slide deck for the EHV-1 dynamic transcriptome project.

- **Deck:** `EHV1_dynamic_transcriptome_conference_master.pptx`
- **Build source:** `build_ehv1_master.js`
- **Primary source:** Tombácz et al., *PLOS ONE* 20(4): e0320439 (2025), doi:10.1371/journal.pone.0320439
- **Analysis source:** this repository, including the final normalized transcript and isoform-ratio tables used for publication.

The deck is intentionally longer than a 5–10 minute talk. It is a master version with 15 main slides plus 2 backup slides, including presenter notes. A shorter conference version can be made by removing intermediate methods/results slides while preserving the central story: integrated long-read/CAGE annotation → temporal kinetic complexity → dynamic splicing and isoform switching.

The PowerPoint is generated with PptxGenJS by the GitHub Actions workflow in `.github/workflows/build-ehv1-presentation.yml`.
