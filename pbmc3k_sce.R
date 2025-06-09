#installing bioconductor packages
library(BiocManager)
BiocManager::install(c("TENxPBMCData", "SingleCellExperiment"))
BiocManager::install(c("SingleR", "celldex"))
#install seurat
install.packages("Seurat")

#load packages
library(TENxPBMCData)
library(SingleCellExperiment)
library(Seurat)
library(patchwork)
library(Matrix)
library(scater)
library(biomaRt)
library(SingleR)
library(celldex)
library(clusterProfiler)
library(org.Hs.eg.db)

#load pbmc 3k dataset (singlecellexperiment format)
sce = TENxPBMCData("pbmc3k")

#extract counts assay as matirx
counts_delayed = assay(sce, "counts")

#convert delayedarray to sparse matrix
counts_mat = as(counts_delayed, "dgCMatrix")

#convert SCE to seurat object
pbmc = as.Seurat(sce, counts = counts_mat, data = NULL)
# Create Seurat object from sparse matrix counts
pbmc = CreateSeuratObject(counts = counts_mat)

# 1. QC: add mitochondrial % and visualize
#get gene symbols
gene_symbols = rowData(sce)$Symbol_TENx
#find mitochondrial genes by symbol
mt_gene = rownames(pbmc)[which(gene_symbols %in% grep("^MT-", gene_symbols, value = TRUE))]
#calculate percent.mt manually
pbmc[["percent.mt"]] = PercentageFeatureSet(pbmc, features = mt_gene)
#visualize
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"))


#filter cells
pbmc = subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)

#normalize
pbmc = NormalizeData(pbmc)

#find variable genes
pbmc = FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
VariableFeatures(pbmc)

#scaling the data
all.genes = rownames(pbmc)
pbmc = ScaleData(pbmc, features = all.genes)

#PCA
pbmc = RunPCA(pbmc, features = VariableFeatures(object = pbmc))

#determine PCs to use
ElbowPlot(pbmc, ndims = 20) + 
  scale_x_continuous(breaks = 1:20)

#clustering
pbmc = FindNeighbors(pbmc, dims = 1:7)
pbmc = FindClusters(pbmc, resolution = 0.5)

#umap
pbmc = RunUMAP(pbmc, dims = 1:7)

#visualization
DimPlot(pbmc, reduction = "umap", label = TRUE, pt.size = 0.7)

#find marker genes for each cluster
markers = FindAllMarkers(pbmc, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
head(markers)

#example markers to plot
names(gene_symbols) = rownames(sce)
rownames(pbmc) = gene_symbols[rownames(pbmc)]
FeaturePlot(pbmc, features = c("CD3D", "MS4A1", "NKG7", "CD14"))
VlnPlot(pbmc, features = c("CD3D", "MS4A1"))

write.csv(markers, "pbmc3k_cluster_markers.csv", row.names = FALSE)

#changine ensemble id to gene symbols
getwd()
setwd("C:/Users/Asus/Desktop/resume_project/pbmc3k_sce")
ensembl = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
gene_data = read.csv("pbmc3k_cluster_markers.csv", stringsAsFactors = FALSE)
colnames(gene_data)
ensembl_ids = gene_data$gene
result <- getBM(attributes = c('ensembl_gene_id', 'hgnc_symbol'),
                filters = 'ensembl_gene_id',
                values = ensembl_ids,
                mart = ensembl)
colnames(result)
merged_data = merge(gene_data, result, by.x = 'gene', by.y = 'ensembl_gene_id', all.x = TRUE)
write.csv(merged_data, 'pbmc3k_cluster_marker_with_symbols.csv', row.names = FALSE)


#annotation
ref = celldex::HumanPrimaryCellAtlasData()
expr= GetAssayData(pbmc, assay = "RNA", slot = "data")
pred = SingleR(test = expr, ref = ref, labels =  ref$label.main)
pbmc$SingleR_label = pred$pruned.labels
DimPlot(pbmc, group.by = "SingleR_label", label = TRUE, repel = TRUE)

gene_markers = read.csv("pbmc3k_cluster_marker_with_symbols.csv", stringsAsFactors = FALSE)
table(pbmc$seurat_clusters, pbmc$SingleR_label)

cluster_annotation = as.data.frame(table(pbmc$seurat_clusters, pbmc$SingleR_label))
library(dplyr)
cluster_map <- cluster_annotation %>%
  group_by(Var1) %>%
  slice_max(order_by = Freq, n = 1) %>%
  select(cluster = Var1, Predicted_Cell_Type = Var2)

cluster_map$cluster = as.numeric(as.character(cluster_map$cluster))
final_markers = merge(gene_markers, cluster_map, by = "cluster", all.x = TRUE)
write.csv(final_markers, "Pbmc3k_marker_with_predicted_cell_types.csv", row.names = FALSE)

# annot_data =  pbmc@meta.data
# head(colnames(annot_data))
# colnames(annot_data)[which(names(annot_data)== "SingleR_label")] = "Predicted_Cell_Type"
# annot_data$Cell_Barcode = rownames(annot_data)
# write.csv(annot_data[, c("Cell_Barcode", "Predicted_Cell_Type")],
#           "pbmc3k_SingleR_celltype_annotation.csv", row.names = FALSE)


#functional enrichment analysis
gene_path = read.csv("pbmc3k_cluster_marker_with_symbols.csv", stringsAsFactors = FALSE)
colnames(gene_path)
marker_genes = gene_path$hgnc_symbol 
entrez_id = bitr(marker_genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)

ego = enrichGO(gene = entrez_id$ENTREZID,
               OrgDb = org.Hs.eg.db,
               ont = "BP",
               pAdjustMethod = "BH",
               pvalueCutoff = 0.05,
               qvalueCutoff = 0.05,
               readable = TRUE)
ego_df = as.data.frame(ego)
write.csv(ego_df, "pathway_involved_in_pbmc3k_marker_genes.csv", row.names = FALSE)
head(ego)

barplot(ego, showCategory = 10)


#cell type proportion
pbmc@meta.data %>%
  count(SingleR_label) %>%
  ggplot(aes(x = reorder(SingleR_label, -n), y = n, fill = SingleR_label)) +
  geom_bar(stat = "identity") +
  labs(title = "Cell Type Proportions", x = "Cell Type", y = "Number of Cells") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_fill_brewer(palette = "Paired")

saveRDS(pbmc, file = "pbmc3k_sce_seurat.rds")
