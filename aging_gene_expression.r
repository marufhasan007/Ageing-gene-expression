# Aging Gene Expression Analysis Using R

## 📁 Project Structure

```
aging-gene-expression-analysis/
├── data/                   # Processed data files
├── results/                # Output: DEGs, plots, enrichment results
├── scripts/                # R scripts for each step
├── README.md               # Project documentation
```

---

## 🔧 Requirements

- R ≥ 4.0
- Bioconductor packages: `GEOquery`, `Biobase`, `limma`, `EnhancedVolcano`, `clusterProfiler`, `org.Hs.eg.db`

---

## 📥 01_data_download.R

```r
# 01_data_download.R
# Download and prepare expression data from GEO

if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(c("GEOquery", "Biobase"))

library(GEOquery)

# Download dataset GSE11882
gse <- getGEO("GSE11882", GSEMatrix = TRUE)[[1]]

# Extract expression data and phenotype
expr_data <- exprs(gse)
pheno_data <- pData(gse)

# Save data
dir.create("data", showWarnings = FALSE)
save(expr_data, pheno_data, file = "data/processed_data.RData")
```

---

## 📊 02_DEG_analysis.R

```r
# 02_DEG_analysis.R
# Perform differential expression analysis

BiocManager::install(c("limma", "EnhancedVolcano"))
library(limma)
library(EnhancedVolcano)

load("data/processed_data.RData")

# Create design matrix
group <- factor(ifelse(grepl("young", pheno_data$title, ignore.case = TRUE), "Young", "Old"))
design <- model.matrix(~0 + group)
colnames(design) <- levels(group)

# Fit model
fit <- lmFit(expr_data, design)
contrast_matrix <- makeContrasts(Old-Young, levels = design)
fit2 <- contrasts.fit(fit, contrast_matrix)
fit2 <- eBayes(fit2)

# Get results
results <- topTable(fit2, adjust = "fdr", number = Inf)
dir.create("results", showWarnings = FALSE)
write.csv(results, "results/DEGs.csv")

# Volcano plot
EnhancedVolcano(results,
                lab = rownames(results),
                x = 'logFC',
                y = 'P.Value',
                pCutoff = 0.05,
                FCcutoff = 1.0,
                title = 'Differential Expression: Old vs. Young')
```

---

## 🧠 03_GO_enrichment.R

```r
# 03_GO_enrichment.R
# Perform GO enrichment on significant DEGs

BiocManager::install(c("clusterProfiler", "org.Hs.eg.db"))
library(clusterProfiler)
library(org.Hs.eg.db)

results <- read.csv("results/DEGs.csv", row.names = 1)
sig_genes <- rownames(results[results$adj.P.Val < 0.05 & abs(results$logFC) > 1, ])

# Map to Entrez IDs
entrez_ids <- mapIds(org.Hs.eg.db, keys = sig_genes, column = "ENTREZID",
                     keytype = "SYMBOL", multiVals = "first")

# Remove NAs
entrez_ids <- na.omit(entrez_ids)

# GO enrichment
ego <- enrichGO(gene = entrez_ids,
                OrgDb = org.Hs.eg.db,
                keyType = "ENTREZID",
                ont = "BP",
                pAdjustMethod = "fdr",
                pvalueCutoff = 0.05)

# Save results
write.csv(as.data.frame(ego), "results/GO_enrichment.csv")
barplot(ego, showCategory = 20)
```

---

## 📊 Example Outputs

- `results/DEGs.csv`: List of all DEGs with statistics.
- `results/GO_enrichment.csv`: Top enriched GO terms.
- Volcano plot of DEGs.
- Barplot of top GO categories.

---

## 🔍 Dataset Info

- **GEO Accession:** [GSE11882](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE11882)
- **Study:** Gene expression in young and old human tissues to understand aging.

---

## 👨‍🔬 Author

**Maruf Hasan**  
GitHub: [https://github.com/marufhasan007]

---

## 🧠 Purpose

This project was built to demonstrate bioinformatics proficiency and domain interest in aging/longevity research, particularly in the context of a postdoctoral application to the Institute for Biostatistics and Informatics in Medicine and Ageing Research (IBIMA), University Medicine Rostock.

---

## 📜 License

MIT License – Feel free to use, modify, and cite this project.
