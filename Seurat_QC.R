## QC and selecting cells for further analysis

# While the `CreateSeuratObject` imposes a basic minimum gene-cutoff, you may want to filter out cells at this stage based
# on technical or biological parameters. Seurat allows you to easily explore QC metrics and filter cells based on any user-defined
# criteria. In the example below, we visualize gene and molecule counts, plot their relationship, and exclude cells with a clear outlier
# number of genes detected as potential multiplets. Of course this is not a guaranteed method to exclude cell doublets, but we include
# this as an example of filtering user-defined outlier cells. We also filter cells based on the percentage of mitochondrial genes present.

# The number of genes and UMIs (nGene and nUMI) are automatically calculated
# for every object by Seurat.  For non-UMI data, nUMI represents the sum of
# the non-normalized values within a cell We calculate the percentage of
# mitochondrial genes here and store it in percent.mito using AddMetaData.
# We use object@raw.data since this represents non-transformed and
# non-log-normalized counts The % of UMI mapping to MT-genes is a common
# scRNA-seq QC metric.
mito.genes <- grep(pattern = "^MT-", x = rownames(pbmc@assays[["RNA"]]), value = TRUE)

percent.mito <- Matrix::colSums(pbmc@assays[["RNA"]][mito.genes, ])/Matrix::colSums(pbmc@assays[["RNA"]])

# AddMetaData adds columns to object@meta.data, and is a great place to
# stash QC stats

#Seurat v2 function, but shows compatibility in Seurat v3
pbmc <- AddMetaData(object = pbmc, metadata = percent.mito, col.name = "percent.mito") 
#in case the above function does not work simply do:
pbmc$percent.mito <- percent.mito

VlnPlot(object = pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3)

# GenePlot is typically used to visualize gene-gene relationships, but can
# be used for anything calculated by the object, i.e. columns in
# object@meta.data, PC scores etc.  Since there is a rare subset of cells
# with an outlier level of high mitochondrial percentage and also low UMI
# content, we filter these as well
par(mfrow = c(1, 2))
FeatureScatter(object = pbmc, feature1 = "nCount_RNA", feature2 = "percent.mito")
FeatureScatter(object = pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

# We filter out cells that have unique gene counts (nFeature_RNA) over 2,500 or less than
# 200 Note that > and < are used to define a'gate'.  
#-Inf and Inf should be used if you don't want a lower or upper threshold.
pbmc <- subset(x = pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mito >  -Inf & percent.mito < 0.05 )
