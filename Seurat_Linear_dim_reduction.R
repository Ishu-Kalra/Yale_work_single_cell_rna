
## Perform linear dimensional reduction

# --> refered to Seurat v2: Next we perform PCA on the scaled data. By default, the genes in `object@var.genes` are used as
# input, but can be defined using pc.genes. We have typically found that running dimensionality reduction on highly variable genes
# can improve performance. However, with UMI data - particularly after regressing out technical variables, we often see that PCA returns 
# similar (albeit slower) results when run on much larger subsets of genes, including the whole transcriptome.

# --> refered to Seurat v3 (latest): high variable features are accessed through the function HVFInfo(object). Despite RunPCA
# has a features argument where to specify the features to compute PCA on, I've been modifying its values and the output PCA graph
# has always the same dimensions, indicating that the provided genes in the features argument are not exactly the ones used to compute
# PCA. Wether the function gets the HVG directly or does not take them into account, I don't know.

pbmc <- RunPCA(object = pbmc,  npcs = 30, verbose = FALSE)

# --> refered to Seurat v2: Seurat provides several useful ways of visualizing both cells and genes that define the PCA,
# including `PrintPCA`, `VizPCA`,  `PCAPlot`, and `PCHeatmap`

#--> refered to Seurat v3 (latest):
#   Seurat v3 provides functions for visualizing:
#   - PCA
# - PCA plot coloured by a quantitative feature
# - Scatter plot across single cells
# - Scatter plot across individual features
# - Variable Feature Plot
# - Violin and Ridge plots
# - Heatmaps

# Examine and visualize PCA results a few different ways
DimPlot(object = pbmc, reduction = "pca")

# Dimensional reduction plot, with cells colored by a quantitative feature
FeaturePlot(object = pbmc, features = "MS4A1")


# Scatter plot across single cells, replaces GenePlot
FeatureScatter(object = pbmc, feature1 = "MS4A1", feature2 = "PC_1")
FeatureScatter(object = pbmc, feature1 = "MS4A1", feature2 = "CD3D")

# Scatter plot across individual features, repleaces CellPlot
CellScatter(object = pbmc, cell1 = "AGTCTACTAGGGTG", cell2 = "CACAGATGGTTTCT")

VariableFeaturePlot(object = pbmc)

# Violin and Ridge plots
VlnPlot(object = pbmc, features = c("LYZ", "CCL5", "IL32"))
RidgePlot(object = pbmc, feature = c("LYZ", "CCL5", "IL32"))

# In particular `DimHeatmap` allows for easy exploration of the primary sources of heterogeneity in a dataset,
# and can be useful when trying to decide which PCs to include for further downstream analyses. Both cells and 
# genes are ordered according to their PCA scores. Setting cells.use to a number plots the ‘extreme’ cells on both
# ends of the spectrum, which dramatically speeds plotting for large datasets. Though clearly a supervised analysis,
# we find this to be a valuable tool for exploring correlated gene sets.

# Heatmaps
DimHeatmap(object = pbmc, reduction = "pca", cells = 200, balanced = TRUE)

# ProjectPCA function is no loger available in Seurat 3.0.