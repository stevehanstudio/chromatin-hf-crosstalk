#!/usr/bin/env Rscript
# Install R packages for chromatin-HF showcase notebooks
# Run from project root: Rscript scripts/install_r_packages.R

message("Installing CRAN packages...")
install.packages(
  c("Seurat", "Harmony", "ggplot2", "dplyr", "tidyr", "tibble", "devtools", "BiocManager"),
  repos = "https://cloud.r-project.org/"
)

message("Installing Bioconductor packages...")
BiocManager::install(
  c("edgeR", "limma", "Rsubread", "GenomicRanges", "rtracklayer",
    "BSgenome.Mmusculus.UCSC.mm10", "BSgenome.Hsapiens.UCSC.hg38",
    "org.Mm.eg.db", "org.Hs.eg.db"),
  update = TRUE,
  ask = FALSE
)

message("Installing IRkernel for Jupyter...")
install.packages("IRkernel", repos = "https://cloud.r-project.org/")
IRkernel::installspec()

message("Done. To install ArchR (scATAC-seq), run in R:")
message('  devtools::install_github("GreenleafLab/ArchR", ref = "master", repos = BiocManager::repositories())')
message('  library(ArchR); ArchR::installExtraPackages()')
