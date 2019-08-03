## Run UMAP

# To visualize the two conditions side-by-side, we can use the split.by argument to show each condition colored by cluster.
pbmc <- RunUMAP(pbmc, reduction = "pca", dims = 1:20)
DimPlot(pbmc, reduction = "umap", split.by = "seurat_clusters")

# You can save the object at this point so that it can easily be loaded back in without having to rerun the computationally intensive steps
# performed above, or easily shared with collaborators.

saveRDS(pbmc, file = "data/pbmc_tutorial.rds")