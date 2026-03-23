# Chromatin Remodeling Drives Immune-Fibroblast Crosstalk in Heart Failure Pathogenesis

Computational replication of Alexanian et al., *Nature* 2024 — dry lab workflows for scRNA-seq, scATAC-seq, bulk RNA-seq, ChIP-seq, and CUT&RUN data from the heart failure chromatin study.

**Source:** Alexanian et al., Nature 2024 ([bioRxiv 2023.01.06.522937](https://www.biorxiv.org/content/10.1101/2023.01.06.522937))  
**Data:** [GEO GSE221699](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE221699) | [BioProject PRJNA915384](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA915384)

---

## Overview

The study identifies BRD4-dependent chromatin remodeling in cardiac myeloid cells that drives *Il1b* expression and immune–fibroblast crosstalk in heart failure. Key findings include:

- **BRD4** in Cx3cr1+ myeloid cells regulates super-enhancers at *Il1b*
- **IL1B** signals to cardiac fibroblasts, inducing *Meox1* and contractility
- **Anti-IL1B** treatment improves cardiac function and reduces fibroblast activation

---

## Quick Start

### Download data

| Script | Data | Disk space | Platform |
|--------|------|------------|----------|
| **download_data.py** | Bulk RNA-seq, ChIP-seq, CUT&RUN, sorted fibroblasts | ~150 GB | Any (ARM64, x86) |
| **download_cellranger_data.py** | scRNA-seq, snRNA-seq, scATAC-seq (10x) | ~500 GB | x86 only (Cell Ranger) |

```bash
# Main data (bulk, ChIP, CUT&RUN) — run on any machine
python download_data.py

# Single-cell data — run on x86 with Cell Ranger
python download_cellranger_data.py
```

**Requirements:** Python 3.9+, `curl`. Optional: `sra-tools` for `--use-sra` (faster for large downloads).

---

## Project Structure

```
chromatin-hf-crosstalk/
├── README.md                          # This file
├── STUDY_ANALYSIS_AND_REPLICATION_ROADMAP.md  # Full replication plan
├── download_data.py                   # Main download (bulk, ChIP, CUT&RUN)
└── download_cellranger_data.py        # Single-cell download (10x)
```

---

## Data Summary

| GEO sub-series | Type | Samples | Key conditions |
|----------------|------|---------|----------------|
| GSE221693 | Bulk RNA-seq | 16 | Unstim, IL1B, TGFB, TGFB+IL1B (iPSC-CF) |
| GSE221694 | ChIP-seq | 6 | H3K27Ac, p65/RELA (iPSC-CF) |
| GSE221695 | CUT&RUN | 4 | BRD4 in Cx3cr1+ Sham/TAC |
| GSE221696 | scATAC-seq | 10 | Whole heart + CD45+ nuclei |
| GSE221698 | scRNA-seq | ~25 | Sham/TAC/JQ1, CD45+, Brd4KO, IgG/anti-IL1B |

---

## Analysis Pipelines

| Phase | Tools | Target |
|-------|-------|--------|
| scRNA/snRNA | Cell Ranger → Seurat, Harmony | Figs 1, 2, 3, 5 |
| scATAC | Cell Ranger ATAC → ArchR | Figs 3, 4 |
| CUT&RUN | nf-core/cutandrun | Fig 4 |
| ChIP-seq | nf-core/chipseq | Fig 5 |
| Bulk RNA | Rsubread, edgeR, limma | Fig 5 |

See [STUDY_ANALYSIS_AND_REPLICATION_ROADMAP.md](STUDY_ANALYSIS_AND_REPLICATION_ROADMAP.md) for full pipeline details, QC parameters, and figure-to-analysis mapping.

---

## References

- Alexanian et al. *Nature* 2024 — This study
- Alexanian et al. *Nature* 595, 438–443 (2021) — MEOX1/fibroblast switch
- Kuppe et al. *Nature* 608, 766–777 (2022) — Human cardiac scATAC
- Link et al. *Cell* 173, 1796–1809 (2018) — RELA ChIP BMDM
