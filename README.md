Aging Gene Expression Analysis with R
Overview

This project explores aging-associated differential gene expression and functional enrichment using human transcriptomic data from the GEO database. 
It applies traditional bioinformatics pipelines and demonstrates reproducible analysis in R, aligning with the research interests of longevity and healthspan.

Dataset and Tools

- Dataset: GSE11882 – Human tissues from young vs. old individuals.
- Tools: R, Bioconductor packages (limma, EnhancedVolcano, clusterProfiler, org.Hs.eg.db)
- Goals:
  - Identify differentially expressed genes (DEGs) in aging
  - Perform GO enrichment analysis to reveal biological processes involved
  - Visualize findings for scientific communication

Repository Structure

aging-gene-expression-analysis/
├── data/                   # Processed expression and metadata
├── results/                # DEGs, GO enrichment, plots
├── scripts/                # R scripts for reproducible analysis
├── README.md               # Project overview

Reproducibility Instructions

1. Install dependencies:
```r
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(c("GEOquery", "limma", "EnhancedVolcano", "clusterProfiler", "org.Hs.eg.db"))
```

2. Run scripts in order:
```r
source("scripts/01_data_download.R")
source("scripts/02_DEG_analysis.R")
source("scripts/03_GO_enrichment.R")
```

Outputs

- DEGs.csv: Differentially expressed genes
- GO_enrichment.csv: Enriched biological processes (GO terms)
- Volcano plot and GO barplot for publication-ready figures

Key Findings

- Significant transcriptomic changes associated with age
- Enrichment in biological processes such as immune regulation, oxidative stress response, and cellular senescence

Purpose

This project was prepared as part of my application to the Institute for Biostatistics and Informatics in Medicine and Ageing Research (IBIMA) at 
University Medicine Rostock to demonstrate technical competence and interest in longevity science.

Author

Maruf Hasan
Email: mhasanmaruf@gmail.com
Interested in molecular aging, bioinformatics, and translational health research

License
MIT License – use freely with citation.
