# Chromatin Remodeling Drives Immune-Fibroblast Crosstalk in Heart Failure Pathogenesis

Computational replication of Alexanian et al., *Nature* 2024 — dry lab workflows for scRNA-seq, scATAC-seq, bulk RNA-seq, ChIP-seq, and CUT&RUN data from the heart failure chromatin study.

**Source:** Alexanian et al., Nature 2024 ([bioRxiv 2023.01.06.522937](https://www.biorxiv.org/content/10.1101/2023.01.06.522937))  
**Data:** [GEO GSE221699](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE221699) | [BioProject PRJNA915384](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA915384)

*Presentation project for UC Berkeley Extension **MCELLBIX110-122 Immunology** (Dr. Ian Fernandopulle). Not required for the course — created to replicate the study.*

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
| **scripts/download_data.py** | Bulk RNA-seq, ChIP-seq, CUT&RUN, sorted fibroblasts | ~150 GB | Any (ARM64, x86) |
| **scripts/download_cellranger_data.py** | scRNA-seq, snRNA-seq, scATAC-seq (10x) | ~500 GB | x86 only (Cell Ranger) |

```bash
# Run from project root
cd chromatin-hf-crosstalk

# Main data (bulk, ChIP, CUT&RUN) — run on any machine
python scripts/download_data.py

# Single-cell data — run on x86 with Cell Ranger
python scripts/download_cellranger_data.py
python scripts/download_cellranger_refs.py   # refs → data/refs/

# Then run Cell Ranger (creates symlinks + cellranger count / cellranger-atac count)
python scripts/run_cellranger.py --ref-scrna data/refs/refdata-gex-GRCm39-2024-A \
                                --ref-atac data/refs/refdata-cellranger-arc-GRCm39-2024-A
```

**Requirements:** Python 3.9+, `curl`. Optional: `sra-tools` for `--use-sra` (faster for large downloads). Cell Ranger references: see [INSTALL.md](INSTALL.md#6-cell-ranger-x86-only).

### Setup environment

```bash
conda env create -f environment.yml
conda activate chromatin-hf
```

See [INSTALL.md](INSTALL.md) for full setup. ArchR requires `bioconductor-rhdf5` and `r-cairo` (now in environment.yml) plus a post-install R step for genome annotations.

**One environment for all machines.** Use the same `environment.yml` on ARM64 and x86. The only reason for multiple machines is that **Cell Ranger is x86-only** — if your main machine is ARM64 (e.g. Apple Silicon), you need an x86 machine for scRNA/scATAC processing. If you have a single x86 machine, you can run the full workflow there with this same environment.

---

## Project Structure

```
chromatin-hf-crosstalk/
├── README.md
├── environment.yml
├── INSTALL.md
├── STUDY_ANALYSIS_AND_REPLICATION_ROADMAP.md
├── scripts/
│   ├── download_data.py
│   ├── download_cellranger_data.py
│   ├── download_cellranger_refs.py   # Cell Ranger refs → data/refs/
│   ├── run_cellranger.py             # Cell Ranger wrapper (x86)
│   └── install_r_packages.R
└── notebooks/                  # Analysis notebooks
    ├── 01_data_overview.ipynb
    ├── 02_bulk_rnaseq_de.ipynb
    ├── 03_chipseq_meox1.ipynb
    ├── 04_cutrun_brd4.ipynb
    ├── 05_scrnaseq_seurat.ipynb
    ├── 06_scatac_archr.ipynb
    └── 07_integration_summary.ipynb
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

## Notebooks

Analysis workflows are in `notebooks/`. **Machine split:**

| Notebook | Run on | Content |
|----------|--------|---------|
| 01_data_overview | either | Data structure, sample metadata |
| 02_bulk_rnaseq_de | **ARM64** | edgeR/limma, MEOX1 (Fig 5) |
| 03_chipseq_meox1 | **ARM64** | ChIP at MEOX1 locus (Fig 5) |
| 04_cutrun_brd4 | **ARM64** | BRD4 at Il1b peaks (Fig 4) |
| 05_scrnaseq_seurat | **x86** | Seurat workflow (Figs 1–3, 5) |
| 06_scatac_archr | **x86** | ArchR workflow (Figs 3, 4) |
| 07_integration_summary | either | Key takeaways |

ARM64: `scripts/download_data.py` (bulk, ChIP, CUT&RUN). x86: `scripts/download_cellranger_data.py` + Cell Ranger (scRNA, scATAC).

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
