
## Determine statistically significant principal components


# To overcome the extensive technical noise in any single gene for scRNA-seq data, Seurat clusters cells based on their PCA scores,
# with each PC essentially representing a ‘metagene’ that combines information across a correlated gene set. Determining how many PCs
# to include downstream is therefore an important step.

# In [Macosko et al](http://www.cell.com/abstract/S0092-8674(15)00549-8), we implemented a resampling test inspired by the jackStraw
# procedure. We randomly permute a subset of the data (1% by default) and rerun PCA, constructing a ‘null distribution’ of gene scores,
# and repeat this procedure. We identify ‘significant’ PCs as those who have a strong enrichment of low p-value genes.

# NOTE: This process can take a long time for big datasets, comment out for
# expediency.  More approximate techniques such as those implemented in
# PCElbowPlot() can be used to reduce computation time
pbmc <- JackStraw(object = pbmc, reduction = "pca", dims = 20, num.replicate = 100,  prop.freq = 0.1, verbose = FALSE)

# The `JackStrawPlot` function provides a visualization tool for comparing the distribution of p-values for each PC with a uniform
# distribution (dashed line). ‘Significant’ PCs will show a strong enrichment of genes with low p-values (solid curve above the dashed line).
# In this case it appears that PCs 1-10 are significant.

pbmc <- ScoreJackStraw(object = pbmc, dims = 1:20, reduction = "pca")

JackStrawPlot(object = pbmc, dims = 1:20, reduction = "pca")

# A more ad hoc method for determining which PCs to use is to look at a plot of the standard deviations of the principle components
# and draw your cutoff where there is a clear elbow in the graph. This can be done with `ElbowPlot`. In this example, it looks like
# the elbow would fall around PC 5.

ElbowPlot(object = pbmc)

# PC selection – identifying the true dimensionality of a dataset – is an important step for Seurat, but can be challenging/uncertain
# for the user. We therefore suggest these three approaches to consider. The first is more supervised, exploring PCs to determine relevant
# sources of heterogeneity, and could be used in conjunction with GSEA for example. The second implements a statistical test based on a
# random null model, but is time-consuming for large datasets, and may not return a clear PC cutoff. The third is a heuristic that is 
# commonly used, and can be calculated instantly. In this example, all three approaches yielded similar results, but we might have been
# justified in choosing anything between PC 7-10 as a cutoff. We followed the jackStraw here, admittedly buoyed by seeing the PCHeatmap
# returning interpretable signals (including canonical dendritic cell markers) throughout these PCs. Though the results are only subtly
# affected by small shifts in this cutoff (you can test below), we strongly suggest always explore the PCs they choose to include downstream.
