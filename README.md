🧬 Aging Gene Expression Analysis with R

<h3>Overvie<h3>

This project investigates aging-associated changes in gene expression using human transcriptomic data. Applying a reproducible R-based bioinformatics pipeline, it identifies differentially expressed genes (DEGs) and performs gene ontology (GO) enrichment analysis to uncover biological processes related to aging. The analysis aligns with current research interests in longevity, healthspan, and intervention analytics.

📊 Dataset and Tools

Dataset: GSE11882 – Gene expression profiles from young vs. old human tissues (GEO database)
Tools & Libraries:
R (v4.3 or later)
Bioconductor packages:
GEOquery
limma
EnhancedVolcano
clusterProfiler
org.Hs.eg.db
🎯 Goals

Identify differentially expressed genes (DEGs) associated with aging.
Conduct GO enrichment analysis to highlight impacted biological pathways.
Create visualizations for effective scientific communication and publication-ready figures.
📁 Repository Structure

aging-gene-expression-analysis/
├── data/             # Processed expression data and sample metadata
├── results/          # DEG lists, GO enrichment results, plots
├── scripts/          # Modular R scripts for each analysis step
├── README.md         # Project overview and instructions
🔁 Reproducibility Instructions

Install dependencies in R:
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install(c(
  "GEOquery", "limma", "EnhancedVolcano",
  "clusterProfiler", "org.Hs.eg.db"
))
Run scripts in sequence:
source("scripts/01_data_download.R")     # Download and process GEO data
source("scripts/02_DEG_analysis.R")      # Perform differential gene expression analysis
source("scripts/03_GO_enrichment.R")     # Run functional enrichment analysis
📂 Outputs

DEGs.csv: Table of differentially expressed genes
GO_enrichment.csv: Enriched GO biological processes
volcano_plot.png: Volcano plot of DEGs
GO_barplot.png: Barplot of top GO terms
🔬 Key Findings

Aging is associated with significant transcriptomic alterations across tissues.
Enrichment in pathways related to:
Immune regulation
Oxidative stress response
Cellular senescence
Results support known hallmarks of aging and provide insights for potential interventions.
🧠 Purpose

This project was developed in the context of my application to the
Institute for Biostatistics and Informatics in Medicine and Ageing Research (IBIMA)
at University Medicine Rostock, to demonstrate both technical competence and genuine interest in aging and longevity research using bioinformatics and AI.

👤 Author

Maruf Hasan
📧 mhasanmaruf@gmail.com
🔬 Interests: Molecular aging | Bioinformatics | Translational health research

📜 License

MIT License — feel free to use, modify, and cite appropriately.
