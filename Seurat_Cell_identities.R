
current.cluster.ids <- c(0, 1, 2, 3, 4, 5, 6, 7)
new.cluster.ids <- c("CD4 T cells", "CD14+ Monocytes", "B cells", "CD8 T cells", "FCGR3A+ Monocytes", "NK cells", "Dendritic cells", "Megakaryocytes")
pbmc@active.ident <- plyr::mapvalues(x = pbmc@active.ident, from = current.cluster.ids, to = new.cluster.ids)
DimPlot(object = pbmc, reduction = "tsne", do.label = TRUE, pt.size = 0.5)
## Further subdivisions within cell types

# If you perturb some of our parameter choices above (for example, setting `resolution=0.8` or changing the number of PCs),
# you might see the CD4 T cells subdivide into two groups. You can explore this subdivision to find markers separating the two
# T cell subsets. However, before reclustering (which will overwrite `object@ident`), we can stash our renamed identities to be
# easily recovered later.
# First lets stash our identities for later
