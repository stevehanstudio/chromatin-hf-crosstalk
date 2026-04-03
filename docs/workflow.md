# End-to-end workflow

This page is a **map of inputs, scripts, tools, notebooks, and outputs** for the replication project. It complements the tables in [README.md](../README.md) and the detail in [STUDY_ANALYSIS_AND_REPLICATION_ROADMAP.md](../STUDY_ANALYSIS_AND_REPLICATION_ROADMAP.md).

**How to read it:** follow solid arrows for the default path. The scATAC notebook has two ways to obtain fragments (Cell Ranger from FASTQ vs. preprocessed GEO); only one is needed for a given run.

---

## 1. Executive pipeline (all modalities)

Shows how public records feed download scripts, where **Cell Ranger requires x86**, and how notebooks connect to the integration summary.

```mermaid
flowchart TB
  SRC[("GEO GSE221699 family\nPRJNA915384 / SRA")]

  subgraph pathA ["Path A — bulk & epigenomics (any machine)"]
    DA["scripts/python/download_data.py"]
    ASSET["data/ — alignments, peaks, counts,\nsorted-cell inputs (per series)"]
    DA --> ASSET
    ASSET --> N01["notebooks/01_data_overview.ipynb"]
    ASSET --> N02["notebooks/02_bulk_rnaseq_de.ipynb"]
    ASSET --> N03["notebooks/03_chipseq_meox1.ipynb"]
    ASSET --> N04["notebooks/04_cutrun_brd4.ipynb"]
  end

  subgraph pathB ["Path B — single-cell from FASTQ (x86)"]
    DC["scripts/python/download_cellranger_data.py"]
    RF["scripts/python/download_cellranger_refs.py\n→ data/refs/"]
    RC["scripts/python/run_cellranger.py\n+ cellranger count / cellranger-atac count"]
    FQ["FASTQ"]
    OUT["Per-run outs —\nmatrix / filtered matrix / fragments / peaks"]
    DC --> FQ
    RF --> RC
    FQ --> RC
    RC --> OUT
    OUT --> N05["notebooks/05_scrnaseq_seurat.ipynb"]
    OUT --> N06["notebooks/06_scatac_archr.ipynb"]
  end

  subgraph pathC ["Path C — scATAC fragments without Cell Ranger (optional)"]
    GEOFRAG["GSE221696 (or similar)\npreprocessed fragment files"]
    GEOFRAG --> N06
  end

  subgraph close ["Synthesis (either machine)"]
    N07["notebooks/07_integration_summary.ipynb"]
  end

  SRC --> DA
  SRC --> DC
  SRC --> GEOFRAG

  N01 --> N07
  N02 --> N07
  N03 --> N07
  N04 --> N07
  N05 --> N07
  N06 --> N07
```

**Notes**

- **Path A:** ChIP/CUT&RUN in this repo are often analyzed from files produced under `data/`; external pipelines (e.g. nf-core) may sit between download and the notebooks depending on your setup—see the roadmap doc.
- **Path B vs. Path C for `06`:** use **either** Cell Ranger fragments **or** GEO fragments, not both as redundant inputs for the same cells.

---

## 2. scATAC detail (ArchR notebook)

Zooms in on **06** inputs: raw 10x pipeline vs. fragment files, and typical ArchR artifacts (regenerated when you re-run the notebook).

```mermaid
flowchart TB
  FQ2["Path B — FASTQ + refs"]
  CRATAC["cellranger-atac count"]
  ATACOUT["outs: fragments.tsv.gz, singlecell.csv, peaks"]
  FR["Path C — GEO fragments.tsv.gz (alternative to Path B)"]
  MERGE["ArchR reads fragments → notebooks/06_scatac_archr.ipynb"]
  ARROW["*.arrow per sample (on disk)"]
  PROJ["ArchRProject"]
  FIGS["UMAP, embeddings, QC plots"]
  N07["07_integration_summary.ipynb"]

  FQ2 --> CRATAC --> ATACOUT --> MERGE
  FR --> MERGE
  MERGE --> ARROW --> PROJ --> FIGS --> N07
```

**Notes**

- **`.arrow` files** are ArchR’s on-disk Arrow format for fragments and derived matrices; they are large and reproducible from fragments—typically gitignored.
- **Spaces in paths:** ArchR is picky about spaces in `outputDirectory`; the notebook uses a no-space path where needed.

---

## Machine split (quick reference)

| Where | Typical work |
| ----- | ------------ |
| **ARM64 or x86** | Path A notebooks (01–04), 07 |
| **x86 only** | Path B downloads, Cell Ranger, 05–06 from FASTQ |

Same Conda environment can be used on both architectures; only Cell Ranger is x86-bound.
