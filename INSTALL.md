# Installation Guide

Tools and libraries for the chromatin-HF showcase notebooks. R runs in Jupyter via IRkernel (no RStudio needed).

---

## 0. Conda environment (recommended)

**One environment for all machines.** Use the same `environment.yml` on ARM64 and x86. Cell Ranger (scRNA/scATAC) is x86-only, so if your main machine is ARM64 (e.g. Apple Silicon), you’ll use an x86 machine only for Cell Ranger steps. If you have a single x86 machine, run the full workflow there with this same environment.

```bash
cd chromatin-hf-crosstalk
conda env create -f environment.yml
conda activate chromatin-hf
```

Register the R kernel (if not auto-registered):

```bash
R -e "IRkernel::installspec()"
```

**Post-install — ArchR dependencies:** `bioconductor-rhdf5` and `r-cairo` are in environment.yml. If ArchR fails with rhdf5 or Cairo errors, install manually:

```bash
conda install -c bioconda bioconductor-rhdf5
conda install -c conda-forge r-cairo
```

**Post-install — ArchR and genome annotations (in R):** Run these in R (not bash). Choose **3: None** when prompted to update packages.

```r
BiocManager::install(c("BSgenome.Mmusculus.UCSC.mm10", "BSgenome.Hsapiens.UCSC.hg38", "org.Mm.eg.db", "org.Hs.eg.db"))
devtools::install_github("GreenleafLab/ArchR", ref = "master", repos = BiocManager::repositories())
library(ArchR)
ArchR::installExtraPackages()
```

---

## 1. System prerequisites (if not using conda)

### R (≥ 4.0)

```bash
# Ubuntu/Debian
sudo apt update
sudo apt install r-base r-base-dev

# Or use conda
conda install -c conda-forge r-base
```

### Python (≥ 3.9)

```bash
# Ubuntu/Debian
sudo apt install python3 python3-pip python3-venv

# Or use conda
conda create -n chromatin-hf python=3.11
conda activate chromatin-hf
```

### System libraries (for IRkernel, some R packages)

```bash
# Ubuntu/Debian
sudo apt install libzmq3-dev libcurl4-openssl-dev libssl-dev libxml2-dev libfontconfig1-dev libharfbuzz-dev libfribidi-dev libfreetype6-dev libpng-dev libtiff5-dev libjpeg-dev
```

---

## 2. Python packages

```bash
pip install jupyterlab jupyter
pip install pandas numpy matplotlib seaborn
pip install scanpy  # optional, for Python single-cell
```

Or with conda:

```bash
conda install -c conda-forge jupyterlab pandas numpy matplotlib seaborn
```

---

## 3. R kernel for Jupyter (IRkernel)

**1.** Start R (in terminal or Cursor):

```bash
R
```

**2.** In R:

```r
install.packages("IRkernel")
IRkernel::installspec()
```

**3.** Check kernels:

```bash
jupyter kernelspec list
```

You should see `ir` (or `R`) in the list. R will then be available as a kernel in Jupyter notebooks in Cursor.

---

## 4. R packages

### 4.1 CRAN packages

```r
install.packages(c(
  "Seurat",        # scRNA-seq
  "Harmony",       # batch correction
  "ggplot2",       # plotting
  "dplyr",         # data manipulation
  "tidyr",
  "tibble",
  "devtools",      # for GitHub installs
  "BiocManager"    # for Bioconductor
))
```

### 4.2 Bioconductor packages

```r
BiocManager::install(c(
  "edgeR",         # bulk RNA-seq DE
  "limma",         # bulk RNA-seq DE, voom
  "Rsubread",      # alignment (bulk RNA-seq)
  "GenomicRanges",
  "GenomicAlignments",
  "rtracklayer",
  "BSgenome.Mmusculus.UCSC.mm10",   # mouse ref (scATAC, etc.)
  "BSgenome.Hsapiens.UCSC.hg38",    # human ref (ChIP, bulk)
  "TxDb.Mmusculus.UCSC.mm10.knownGene",
  "org.Mm.eg.db",
  "org.Hs.eg.db"
))
```

### 4.3 ArchR (from GitHub)

ArchR is for scATAC-seq. **Prerequisites:** Install `bioconductor-rhdf5` and `r-cairo` via conda first (see Section 0), or ArchR will fail with rhdf5/Cairo errors.

Install after the above:

```r
devtools::install_github("GreenleafLab/ArchR", ref = "master", repos = BiocManager::repositories())
```

Then install ArchR extras (choose 3: None when prompted):

```r
library(ArchR)
ArchR::installExtraPackages()
```

**Note:** ArchR works on Linux/macOS, not Windows.

### 4.4 Optional R packages

```r
install.packages("enrichR")   # GO enrichment (or use web)
BiocManager::install("csaw")  # ChIP-seq differential binding
```

---

## 5. One-shot R install script

Run from project root: `Rscript scripts/install_r_packages.R`

```r
# CRAN
install.packages(c("Seurat", "Harmony", "ggplot2", "dplyr", "tidyr", "devtools", "BiocManager"))

# Bioconductor
BiocManager::install(c("edgeR", "limma", "Rsubread", "GenomicRanges", "rtracklayer",
  "BSgenome.Mmusculus.UCSC.mm10", "BSgenome.Hsapiens.UCSC.hg38",
  "org.Mm.eg.db", "org.Hs.eg.db"))

# IRkernel (for Jupyter)
install.packages("IRkernel")
IRkernel::installspec()

# ArchR (run separately if needed)
# devtools::install_github("GreenleafLab/ArchR", ref = "master", repos = BiocManager::repositories())
```

---

## 6. Cell Ranger (x86 only)

Needed for scRNA-seq and scATAC-seq pipelines. **Cell Ranger is x86-only** — that’s why some workflows require an x86 machine. Use the same `chromatin-hf` conda env there. Install from [10x Genomics support](https://support.10xgenomics.com/single-cell-gene-expression/software/downloads/latest).

### 6.1 Reference genomes (mouse)

Download mouse references (scRNA + scATAC) to `data/refs/`:

```bash
python scripts/download_cellranger_refs.py
```

This fetches:
- **scRNA:** refdata-gex-GRCm39-2024-A (~9.7 GB) for Cell Ranger 10
- **scATAC:** refdata-cellranger-arc-GRCm39-2024-A (~13 GB) for Cell Ranger ATAC 2.2

Options: `--scrna-only`, `--atac-only`, `-o /path/to/refs`

### 6.2 Download FASTQs

```bash
python scripts/download_cellranger_data.py
```

- **scRNA/snRNA:** Uses ENA by default (2 FASTQs per run).
- **scATAC:** ENA provides only 2 FASTQs; 10x needs R1, I2 (barcode), R2. The script auto-uses SRA (`fasterq-dump --split-files --include-technical`) for scATAC runs. Install `sra-tools` (`conda install -c bioconda sra-tools`).

If you previously downloaded scATAC via ENA, remove those run dirs and re-run to get the 3-file layout.

### 6.3 Run Cell Ranger wrapper

After `scripts/download_cellranger_data.py` and `scripts/download_cellranger_refs.py`:

```bash
python scripts/run_cellranger.py \
  --ref-scrna data/refs/refdata-gex-GRCm39-2024-A \
  --ref-atac data/refs/refdata-cellranger-arc-GRCm39-2024-A
```

The script creates 10x-compatible symlinks and runs `cellranger count` (22 scRNA runs) and `cellranger-atac count` (10 scATAC runs). Output goes to `output/cellranger/`. Use `--dry-run` to preview commands.

### 6.4 Other optional tools

| Tool | Purpose | Install |
|------|---------|---------|
| **HOMER** | Motif enrichment | `conda install -c bioconda homer` |

---

## 7. Quick checklist

**Conda path:**
- [ ] `conda env create -f environment.yml` + `conda activate chromatin-hf`
- [ ] `conda install -c bioconda bioconductor-rhdf5` (for ArchR)
- [ ] `conda install -c conda-forge r-cairo` (for ArchR plotting)
- [ ] R: genome annotations + ArchR (see Section 0)

**Manual path:**
- [ ] R ≥ 4.0, Python ≥ 3.9
- [ ] System libs: `libzmq3-dev libcurl4-openssl-dev libssl-dev`
- [ ] `install.packages("IRkernel")` + `IRkernel::installspec()`
- [ ] R: Seurat, Harmony, edgeR, limma, BiocManager packages
- [ ] R: ArchR (if doing scATAC; rhdf5 + Cairo required)

---

## 8. Troubleshooting

| Issue | Fix |
|-------|-----|
| **ArchR fails: "rhdf5 is not available"** | `conda install -c bioconda bioconductor-rhdf5` |
| **Cairo build fails (X11/X.h missing)** | `conda install -c conda-forge r-cairo` (use conda, not R) |
| **BSgenome download timeout** | `options(timeout=600)` in R, then retry; or install one at a time |
| **"Enter one or more numbers"** (update prompt) | Choose **3: None** to skip and continue |
| **R code in bash** | Run R first (`R`), then paste commands; or use `R -e 'command'` |

---

## 9. Verify

```bash
# Python + Jupyter
jupyter --version

# R kernel
jupyter kernelspec list
```

```r
# In R
library(Seurat)
library(ArchR)   # if installed
```
