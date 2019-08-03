

## Run Non-linear dimensional reduction (tSNE)

# Seurat continues to use tSNE as a powerful tool to visualize and explore these datasets. While we no longer advise 
# clustering directly on tSNE components, cells within the graph-based clusters determined above should co-localize on
# the tSNE plot. This is because the tSNE aims to place cells with similar local neighborhoods in high-dimensional space
# together in low-dimensional space. As input to the tSNE, we suggest using the same PCs as input to the clustering analysis,
# although computing the tSNE based on scaled gene expression is also supported using the genes.use argument.

pbmc <- RunTSNE(object = pbmc, dims.use = 1:10, do.fast = TRUE)
# note that you can set do.label=T to help label individual clusters
DimPlot(object = pbmc, reduction = "tsne")
