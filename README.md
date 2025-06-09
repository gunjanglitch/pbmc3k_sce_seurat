# 🔬 PBMC3k Single-Cell RNA-Seq Analysis Using Seurat

This project presents a full downstream analysis pipeline of the **PBMC3k** single-cell RNA-seq dataset using the Seurat package in R. The dataset contains ~3,000 peripheral blood mononuclear cells from 10X Genomics.

---

## 📌 Objectives
- Perform quality control and filtering
- Normalize data and identify variable genes
- Run PCA, UMAP for dimensionality reduction
- Cluster cells and identify marker genes
- Annotate clusters using SingleR
- Perform functional enrichment analysis (GO)
- Visualize cell type proportions and pathways

---

## 🧰 Tools & Packages
- `Seurat`
- `TENxPBMCData`, `SingleR`, `celldex`
- `clusterProfiler`, `org.Hs.eg.db`
- `biomaRt`, `ggplot2`, `dplyr`

---

## 📁 Repository Structure

| Folder       | Description                           |
|--------------|---------------------------------------|
| `data/`      | Contains output files and results     |
| `figures/`   | Plots: UMAP, heatmap, QC, GO analysis |
| `scripts/`   | R scripts used for full pipeline      |

---

## 📚 Dataset Info

- **Source**: `TENxPBMCData::pbmc3k` (Bioconductor)
- ~3,000 PBMCs processed using 10X Genomics Chromium

---

## 🧪 Author

**Gunjan Sarode**   
🌐 [LinkedIn](https://www.linkedin.com/in/gunjan-sarode/) | 📫 gunjansarode.bioinfo@gmail.com

---

## ✅ License

This project is for educational and research use.
