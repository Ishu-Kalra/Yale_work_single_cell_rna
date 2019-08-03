## Normalizing the data

# After removing unwanted cells from the dataset, the next step is to normalize the data. By default, we employ a
# global-scaling normalization method “LogNormalize” that normalizes the gene expression measurements for each cell
# by the total expression, multiplies this by a scale factor (10,000 by default), and log-transforms the result.

pbmc <- NormalizeData(object = pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
