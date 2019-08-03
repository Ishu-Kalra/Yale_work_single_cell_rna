## Detection of variable genes across the single cells

# Seurat calculates highly variable genes and focuses on these for downstream analysis. `FindVariableGenes`
# calculates the average expression and dispersion for each gene, places these genes into bins, and then calculates
# a z-score for dispersion within each bin. This helps control for the relationship between variability and average expression.
# This function is unchanged from (Macosko et al.), but new methods for variable gene expression identification are coming soon.
# We suggest that users set these parameters to mark visual outliers on the dispersion plot, but the exact parameter settings may
# vary based on the data type, heterogeneity in the sample, and normalization strategy. The parameters here identify ~2,000 variable genes,
# and represent typical parameter settings for UMI data that is normalized to a total of 1e4 molecules.
pbmc <- FindVariableFeatures(object = pbmc, mean.function = ExpMean, dispersion.function = LogVMR, x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5, nfeatures = 2000)


# To view the output of the FindVariableFeatures output we use this function. The genes appear not to be stored in the object,
# but can be accessed this way.
head(x = HVFInfo(object = pbmc))
