# Chromatin Remodeling Drives Immune-Fibroblast Crosstalk in Heart Failure Pathogenesis

## Dry Lab Replication Plan

**Source:** Alexanian et al., Nature 2024 (bioRxiv 2023.01.06.522937)  
**Data:** GEO **GSE221699** | BioProject **PRJNA915384** | SRA Run Selector: [PRJNA915384](https://www.ncbi.nlm.nih.gov/Traces/study/?acc=PRJNA915384)

This document focuses on **dry lab (computational) replication**. Wet lab experiments are summarized to provide context for how each dataset was generated and how it relates to the analyses.

---

## Part 1: Wet Lab Context (Data Provenance)

Understanding what wet lab experiments produced each dataset is essential for correct interpretation and downstream analysis.

### 1.1 In Vivo Mouse Experiments → sc/snRNA-seq, scATAC-seq

| Experiment | What was done | Output data |
|------------|---------------|-------------|
| **Sham / TAC / TAC+JQ1** | C57BL/6J mice: sham surgery or TAC (aortic constriction). TAC mice received vehicle or JQ1 (50 mg/kg IP daily) for 30 days. Hearts perfused (Langendorff), non-cardiomyocytes isolated, 10X scRNA-seq. | Non-CM scRNA-seq (Sham, TAC, TAC_JQ1) → Fig 1 |
| **TAC / TAC-Brd4KO** | Cx3cr1-CreERT2;Brd4flox/flox vs Brd4flox/flox. Tamoxifen to delete Brd4 in Cx3cr1+ cells. CD45+ cells FACS-sorted from heart → 10X scRNA-seq. | CD45+ scRNA-seq (Sham, TAC, TAC_Brd4KO) → Fig 2 |
| **Whole heart nuclei** | Same genotypes, TAC vs TAC-Brd4KO. Nuclei isolated (PAN-INTACT), snRNA-seq and scATAC-seq. | snRNA-seq + scATAC-seq (all cardiac cells) → Fig 3 |
| **CD45+ nuclei** | CD45+ cells FACS-sorted, nuclei → scATAC-seq. | CD45+ scATAC-seq (Sham, TAC, TAC_Brd4KO) → Fig 4 |
| **TAC + anti-IL1B** | TAC mice treated with 500 µg IgG or anti-IL1B IP every 3 days. Non-CM scRNA-seq at day 30. | Non-CM scRNA-seq (IgG vs anti-IL1B) → Fig 5 |

### 1.2 Epigenomics → CUT&RUN, ChIP-seq

| Experiment | What was done | Output data |
|------------|---------------|-------------|
| **BRD4 CUT&RUN** | Brd4flag/flag mice (Sham or TAC). Cx3cr1+ cells FACS-sorted, anti-FLAG CUT&RUN. | BRD4 occupancy in Cx3cr1+ cells → Fig 4 |
| **ChIP-seq (iPSC-CF)** | Human iPSC-derived cardiac fibroblasts, unstimulated. H3K27Ac and p65/RELA ChIP-seq. | Active enhancers, RELA binding at MEOX1 locus → Fig 5 |

### 1.3 Bulk RNA-seq (iPSC-CF)

| Experiment | What was done | Output data |
|------------|---------------|-------------|
| **iPSC-CF bulk RNA-seq** | iPSC cardiac fibroblasts treated: Unstim, IL1B, TGFB, TGFB+IL1B. 4 replicates each. | DE analysis → MEOX1 as TGFB+IL1B–responsive gene → Fig 5 |

### 1.4 External Public Data (Used by Paper)

| Dataset | Source | Use |
|---------|--------|-----|
| **Human cardiac scATAC-seq** | Kuppe et al. Nature 2022 (GEO) | IL1B locus accessibility in control vs MI patients → Fig 4N |
| **RELA ChIP-seq (BMDM)** | Link et al. Cell 2018 | p65/RELA binding at Il1b enhancer peaks → Fig 4I |

### 1.5 Wet Lab Validation (No Sequencing — Informs Dry Lab Interpretation)

These assays do not produce data in GEO but validate findings from the computational pipelines:

| Assay | Purpose | Dry lab finding validated |
|-------|---------|---------------------------|
| **CRISPR Il1b Peak 5/6 deletion** (RAW 264.7) | Deletion of peaks identified by scATAC + CUT&RUN ablates LPS-induced Il1b | Il1b enhancer peaks 5 & 6 are functional |
| **Luciferase (Il1b Peak 5/6, MEOX1 Peak 9/10)** | BRD4 + p65/RELA co-regulation of cloned enhancers | Enhancers identified by ChIP/ATAC are responsive |
| **IL1B neutralization in TAC mice** | Anti-IL1B improves cardiac function, reduces fibroblast activation | scRNA shows fibroblast cluster shift (IgG vs anti-IL1B) |
| **BMDM conditioned media → iPSC-CF contraction** | Paracrine IL1B from macrophages drives fibroblast contractility | Supports IL1B as the myeloid→fibroblast signal |

---

## Part 2: GEO Data Structure (GSE221699)

The super-series **GSE221699** contains 57 samples across 6 sub-series:

| Sub-series | Data type | Samples | Key conditions |
|------------|-----------|---------|----------------|
| **GSE221693** | Bulk RNA-seq | 16 | Unstim, IL1B, TGFB, TGFB+IL1B (iPSC-CF) |
| **GSE221694** | ChIP-seq | 6 | H3K27Ac, p65/RELA, input (iPSC-CF) |
| **GSE221695** | CUT&RUN | 4 | Sham/TAC Cx3cr1+ FLAG, IgG control |
| **GSE221696** | scATAC-seq | 10 | Whole heart (TAC WT/KO); CD45+ (Sham, TAC, TAC_Brd4KO) |
| **GSE221698** | scRNA-seq | ~25 | Sham/TAC/TAC_JQ1; CD45+; TAC WT/KO; TAC IgG/anti-IL1B |
| **GSE247261** | scRNA-seq | (alt) | Additional scRNA-seq |

**Platforms:** Illumina NextSeq 500 (mouse, human), NovaSeq 6000 (mouse).

---

## Part 3: Dry Lab Replication Plan (Detailed)

### Phase 1: Environment & Data Acquisition

#### 1.1 Software Stack

| Tool | Version | Purpose |
|------|---------|---------|
| Cell Ranger | 3.1.x | sc/snRNA-seq alignment, count matrices |
| Cell Ranger ATAC | 2.0.x | scATAC-seq alignment, fragments |
| Seurat | 4.0.1 | scRNA-seq analysis |
| ArchR | 1.0.1 | scATAC-seq analysis |
| R | ≥4.0 | |
| Harmony | (R package) | Batch correction |
| nf-core/cutandrun | 2.4.2 | CUT&RUN processing |
| nf-core/chipseq | 1.2.2 | ChIP-seq processing |
| HOMER | (findMotifsGenome) | Motif enrichment |
| Enrichr | (web or R) | GO analysis |
| Rsubread, edgeR, limma | (Bioconductor) | Bulk RNA-seq |

#### 1.2 Download Strategy

| Script | Data | Disk space | Run on |
|--------|------|------------|--------|
| **scripts/download_data.py** | Bulk RNA-seq, ChIP-seq, CUT&RUN, sorted fibroblasts | **~150 GB** | Any machine (ARM64, x86) |
| **scripts/download_cellranger_data.py** | scRNA-seq, snRNA-seq, scATAC-seq (10x) | ~500 GB | x86 only (Cell Ranger) |

```bash
# Main download (bulk, ChIP, CUT&RUN) — run on any machine
python scripts/download_data.py

# Cell Ranger data (scRNA, scATAC) — run on x86 where Cell Ranger is installed
python scripts/download_cellranger_data.py
```

**Note:** For scRNA/scATAC, you typically need raw FASTQs to reproduce Cell Ranger → Seurat/ArchR workflows. GEO may host count matrices; check each sub-series.

---

### Phase 2: scRNA-seq Analysis (Seurat)

**Target figures:** 1, 2, 3, 5 (scRNA parts), S1, S2, S3, S8

#### 2.1 Pipeline (Per Paper Methods)

1. **Cell Ranger**
   - `cellranger mkfastq` (if BCL) or use pre-demultiplexed FASTQs
   - `cellranger count` → mm10, `--include-introns` for snRNA compatibility
   - `cellranger aggr` → normalize to least-sequenced sample

2. **Seurat**
   - `Read10X` / `CreateSeuratObject`
   - Metadata: `gem.group` = condition + replicate
   - QC filters:
     - `nFeature_RNA`: 2000–7500
     - `nCount_RNA`: < 80,000
     - `percent.mt`: < 15%
   - `SCTransform` (or standard normalization)
   - `RunPCA` → Harmony (`split.by = "gem.group"`)
   - `RunUMAP`, `FindNeighbors`, `FindClusters`
   - `FindAllMarkers` (Wilcoxon), `FindMarkers` for pairwise DE

3. **Replicate per dataset**
   - **Fig 1:** Sham, TAC, TAC_JQ1 non-CM → myeloid sub-clustering, stress-correlated genes
   - **Fig 2:** CD45+ Sham, TAC, TAC_Brd4KO → monocyte/macrophage clusters, Il1b DE
   - **Fig 3:** snRNA TAC vs TAC_Brd4KO → fibroblast clusters, DE
   - **Fig 5:** TAC_IgG vs TAC_antiIL1B → fibroblast clusters, Postn/Meox1

---

### Phase 3: Correlation Analysis (Cardiac Function ↔ Gene Expression)

**Target:** Fig 1D–G (stress-correlated genes in myeloid cells)

#### 3.1 Method (From Paper)

- For each cell type (≥250 cells): fit **Poisson GLM**: `reads ~ ejection_fraction`
- Score = correlation coefficient of EF with gene expression
- Normalize scores by IQR; compute variance → p-value
- Threshold: `|score| > 5` (p < 1e−6) for strong correlation

#### 3.2 Implementation Notes

- Need **ejection fraction per condition** (Sham, TAC, TAC_JQ1). Paper reports means; use those or extract from figure/supplement.
- Downsample to 1,000 cells per cell type per condition when fitting.
- Expect ~22 genes in myeloid with score < −5 (anticorrelated with function).

---

### Phase 4: scATAC-seq Analysis (ArchR)

**Target figures:** 3I–O, 4A–N, S4, S5, S6

#### 4.1 Pipeline

1. **Cell Ranger ATAC**
   - `cellranger-atac count` → mm10

2. **ArchR**
   - Create project, add fragments
   - QC: `minTSS = 15–16`, `minFrags = 1584–3163`, `maxFrags = 1e6`
   - `addDoubletScores`, `filterDoublets`
   - `addIterativeLSI` (genome-wide tiling)
   - `addHarmony` (batch correction)
   - `addClusters`, `addUMAP`
   - Gene scores: `getMarkerFeatures(useMatrix = "GeneScoreMatrix")`
   - Peak calling: `addGroupCoverages`, `addReproduciblePeakSet` (MACS2)
   - DARs: `getMarkerFeatures` (useGroups vs bgdGroups); FDR<0.1, Log2FC>1
   - Motifs: `peakAnnoEnrichment`; HOMER `findMotifsGenome.pl -size given -mask`

3. **Super-enhancer discovery (custom)**
   - Modified ROSE: stitch distal peaks (12.5 kb gap), exclude gene bodies
   - Score = sum TSS-normalized accessibility
   - Tangent at elbow of rank vs score → 749 super-enhancers
   - Wilcoxon (p < 1e−5): TAC vs Sham, TAC vs TAC_Brd4KO
   - Identify regions that open with TAC and close with Brd4 KO

---

### Phase 5: CUT&RUN Processing

**Target:** Fig 4H–I, S6

#### 5.1 Pipeline

```bash
nextflow run nf-core/cutandrun \
  --genome GRCm38 \
  --input design_matrix.csv \
  -profile singularity
```

- Input: design CSV with sample IDs, paths to FASTQ
- Output: peaks, bigWig coverage
- Compare Sham vs TAC BRD4 signal at Il1b peaks 1–7
- Use `csaw` (R) for 500 bp bin-level differential binding if needed

---

### Phase 6: ChIP-seq Processing

**Target:** Fig 5C (H3K27Ac, p65/RELA at MEOX1)

#### 6.1 Pipeline

```bash
nextflow run nf-core/chipseq \
  --max_memory 80.GB --single_end --narrow_peak \
  --skipBiotypeQC --skipTrimming \
  --genome GRCh38 \
  --input design_matrix.csv \
  -profile singularity
```

- Genome: **hg38** (human iPSC-CF)
- Extract MEOX1 locus coverage; confirm Peak 9/10 syntenic region

---

### Phase 7: Bulk RNA-seq (iPSC-CF)

**Target:** Fig 5B (MEOX1 in TGFB+IL1B quadrant)

#### 7.1 Pipeline

- Align: Rsubread to **hg19** (paper uses GRCm37/hg19)
- Count: featureCounts (Ensembl hg19)
- Normalize: edgeR `calcNormFactors`
- DE: limma-voom; filter CPM as in paper
- Identify genes: LogFC > 2 for TGFB response AND increased in TGFB+IL1B vs TGFB alone
- MEOX1 should fall in top-right quadrant

---

### Phase 8: Human scATAC (External Data)

**Target:** Fig 4L–N

- Download Kuppe et al. 2022 human cardiac scATAC (GEO)
- Process with ArchR (or use preprocessed)
- Subset myeloid cells → cluster 3 (CX3CR1/IL1B high)
- Plot IL1B locus coverage: control vs MI patients
- Lift mouse Il1b locus to human (mm10 chr2:129.3–129.4 Mb → hg38 chr2:112.8–112.84 Mb)

---

### Phase 9: Locus-Specific Tracks & Co-accessibility

**Target:** Fig 3O, 4I, 4N

#### 9.1 TSS-Normalized Coverage

- Per-cell TSS-normalized accessibility (from ArchR bins)
- Mean per condition at Postn, Meox1, Il1b loci
- Confidence intervals from cell-to-cell variance

#### 9.2 Co-accessibility

- Jaccard similarity: binarized accessibility at peak vs promoter
- Highlight Postn Peak 11, Meox1 Peak 9/10, Il1b Peaks 5–6

---

### Phase 10: Gene Ontology & Motif Enrichment

| Analysis | Tool | Input |
|----------|------|-------|
| GO (DE genes) | Enrichr | Gene lists from FindMarkers |
| TF motifs (DARs) | HOMER findMotifsGenome | Peak BED, `-size given -mask` |
| TF motifs (ArchR) | peakAnnoEnrichment | ArchR peak set |

---

## Part 4: Figure-to-Analysis Mapping

| Figure | Primary data | Key analyses |
|--------|--------------|--------------|
| Fig 1 | scRNA Sham/TAC/TAC_JQ1 | UMAP, correlation, myeloid clusters, Il1b/Cx3cr1 |
| Fig 2 | CD45+ scRNA | Monocyte/macrophage clusters, Il1b DE, cluster 4 |
| Fig 3 | snRNA + scATAC TAC vs Brd4KO | Fibroblast DE, chromatin accessibility, motifs |
| Fig 4 | CD45+ scATAC, CUT&RUN | Super-enhancers, Il1b peaks, BRD4/RELA |
| Fig 5 | Bulk RNA, ChIP, scRNA IgG/anti-IL1B | MEOX1 DE, ChIP at MEOX1, fibroblast resolution |

---

## Part 5: Key Parameters Quick Reference

| Parameter | Value |
|-----------|-------|
| scRNA QC nFeature | 2000–7500 |
| scRNA QC mito | < 15% |
| scRNA QC nCount | < 80,000 |
| scATAC minTSS | 15–16 |
| scATAC minFrags | 1584–3163 |
| DAR FDR | < 0.1 |
| DAR Log2FC | > 1 |
| Correlation threshold | \|score\| > 5 (p < 1e−6) |
| Super-enhancer Wilcoxon | p < 1e−5 |
| Bulk DE (MEOX1) | LogFC > 2 |

---

## Part 6: Expected Outcomes & Validation

- **Fig 1:** ~41,626 cells, 20 clusters; myeloid cluster 1 enriched for Cx3cr1/Il1b; 22 genes score < −5
- **Fig 2:** Cluster 4 = monocyte-derived, Il1b high; 13 genes TAC↑ and Brd4KO↓
- **Fig 3:** 310 fibroblast DE genes TAC vs Brd4KO; Postn/Meox1 accessibility ↓ in Brd4KO
- **Fig 4:** 749 super-enhancers; Il1b peaks 5 & 6 BRD4-bound; CRISPR validates peaks 5 & 6
- **Fig 5:** MEOX1 LogFC > 2 in TGFB+IL1B; fibroblast clusters 4 & 6 depleted in anti-IL1B

---

## Part 7: References

- Alexanian et al. *Nature* 595, 438–443 (2021) — MEOX1/fibroblast switch, Peak 9/10
- Alexanian et al. *Nature* 2024 — This study
- GEO GSE221699; BioProject PRJNA915384
- Kuppe et al. *Nature* 608, 766–777 (2022) — Human cardiac scATAC
- Link et al. *Cell* 173, 1796–1809 (2018) — RELA ChIP BMDM

