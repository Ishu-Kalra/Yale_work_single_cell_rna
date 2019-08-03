
pbmc <- StashIdent(object = pbmc, save.name = "ClusterNames_0.6")

# Note that if you set save.snn=T above, you don't need to recalculate the
# SNN, and can simply put: pbmc <- FindClusters(pbmc,resolution = 0.8)
pbmc <- FindClusters(object = pbmc, reduction.type = "pca", dims.use = 1:10, resolution = 0.8, print.output = FALSE)

# Demonstration of how to plot two tSNE plots side by side, and how to color
# points based on different criteria
plot1 <- DimPlot(object = pbmc, reduction = "tsne", do.return = TRUE, no.legend = TRUE, do.label = TRUE)
plot2 <- DimPlot(object = pbmc, reduction = "tsne", do.return = TRUE, group.by = "ClusterNames_0.6", no.legend = TRUE, do.label = TRUE)
plot_grid(plot1, plot2)

# Find discriminating markers
tcell.markers <- FindMarkers(object = pbmc, ident.1 = 0, ident.2 = 1)

# Most of the markers tend to be expressed in C1 (i.e. S100A4). However, we
# can see that CCR7 is upregulated in C0, strongly indicating that we can
# differentiate memory from naive CD4 cells.  cols.use demarcates the color
# palette from low to high expression
FeaturePlot(object = pbmc, features = c("S100A4", "CCR7"), cols = c("green", "blue"))

# The memory/naive split is bit weak, and we would probably benefit from looking at more cells to see if this becomes
# more convincing. In the meantime, we can restore our old cluster identities for downstream processing.

pbmc <- SetIdent(object = pbmc, value = "ClusterNames_0.6")
saveRDS(pbmc, file = "data/pbmc3k_final.rds")
