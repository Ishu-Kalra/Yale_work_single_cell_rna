
## Finding differentially expressed genes (cluster biomarkers)

# Seurat can help you find markers that define clusters via differential expression. By default, it identifes positive and
# negative markers of a single cluster (specified in `ident.1`), compared to all other cells. `FindAllMarkers` automates this
# process for all clusters, but you can also test groups of clusters vs. each other, or against all cells.

# The `min.pct` argument requires a gene to be detected at a minimum percentage in either of the two groups of cells, and the 
# thresh.test argument requires a gene to be differentially expressed (on average) by some amount between the two groups. You can
# set both of these to 0, but with a dramatic increase in time - since this will test a large number of genes that are unlikely 
# to be highly discriminatory. As another option to speed up these computations, `max.cells.per.ident` can be set. This will downsample 
# each identity class to have no more cells than whatever this is set to. While there is generally going to be a loss in power, the
# speed increases can be significiant and the most highly differentially expressed genes will likely still rise to the top.


# find all markers of cluster 1
cluster1.markers <- FindMarkers(object = pbmc, ident.1 = 1, min.pct = 0.25)
print(x = head(x = cluster1.markers, n = 5))

# find all markers distinguishing cluster 5 from clusters 0 and 3
cluster5.markers <- FindMarkers(object = pbmc, ident.1 = 2, ident.2 = c(0, 3), min.pct = 0.25)
print(x = head(x = cluster5.markers, n = 5))


# find markers for every cluster compared to all remaining cells, report
# only the positive ones
pbmc.markers <- FindAllMarkers(object = pbmc, only.pos = TRUE, min.pct = 0.25, thresh.use = 0.25)
pbmc.markers %>% group_by(cluster) %>% top_n(2, avg_logFC)

# Seurat has several tests for differential expression which can be set with the test.use parameter
# (see our [DE vignette](http://satijalab.org/seurat/de_vignette.html) for details). For example, the ROC test returns 
# the ‘classification power’ for any individual marker (ranging from 0 - random, to 1 - perfect).


cluster1.markers <- FindMarkers(object = pbmc, ident.1 = 0, thresh.use = 0.25, test.use = "roc", only.pos = TRUE)


# We include several tools for visualizing marker expression. 
# • `VlnPlot` (shows expression probability distributions across clusters), 
# • and `FeaturePlot` (visualizes gene expression on a tSNE or PCA plot) are our most commonly used visualizations. 
# We also suggest exploring:
#   • `RidgePlot`, 
# • `CellPlot`, and 
# • `DotPlot` as additional methods to view your dataset.

VlnPlot(object = pbmc, features =c("NKG7", "PF4"))


FeaturePlot(object = pbmc, features = c("MS4A1", "GNLY", "CD3E", "CD14", "FCER1A", "FCGR3A", "LYZ", "PPBP", "CD8A"), cols = c("grey", "blue"), reduction = "tsne")

# `DoHeatmap` generates an expression heatmap for given cells and genes. In this case, we are plotting the top 20 markers
# (or all markers if less than 20) for each cluster.

top10 <- pbmc.markers %>% group_by(cluster) %>% top_n(10, avg_logFC)
# setting slim.col.label to TRUE will print just the cluster IDS instead of
# every cell name
DoHeatmap(object = pbmc, features = top10$gene, label = TRUE)

## Assigning cell type identity to clusters

# Fortunately in the case of this dataset, we can use canonical markers to easily match the unbiased clustering to known cell types.
