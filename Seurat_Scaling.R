
## Scaling the data and removing unwanted sources of variation

# Your single cell dataset likely contains ‘uninteresting’ sources of variation. This could include not only technical noise, but
# batch effects, or even biological sources of variation (cell cycle stage). As suggested in [Buettner et al, NBT, 2015]
# (https://www.nature.com/articles/nbt.3102), regressing these signals out of the analysis can improve downstream dimensionality reduction
# and clustering. To mitigate the effect of these signals, Seurat constructs linear models to predict gene expression based on user-defined
# variables. The scaled z-scored residuals of these models are stored in the scale.data slot, and are used for dimensionality reduction and
# clustering.

# We can regress out cell-cell variation in gene expression driven by batch (if applicable), cell alignment rate (as provided by
# Drop-seq tools for Drop-seq data), the number of detected molecules, and mitochondrial gene expression. For cycling cells, we can
# also learn a ‘cell-cycle’ score (see example [here](http://satijalab.org/seurat/cell_cycle_vignette.html)) and regress this out as
# well. In this simple example here for post-mitotic blood cells, we regress on the number of detected molecules per cell as well as
# the percentage mitochondrial gene content.

# Seurat v2.0 implements this regression as part of the data scaling process. Therefore, the `RegressOut` function has been deprecated,
# and replaced with the vars.to.regress argument in `ScaleData`.

pbmc <- ScaleData(object = pbmc, vars.to.regress = c("nCounts_RNA", "percent.mito"))
