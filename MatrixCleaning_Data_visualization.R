
library(SingleCellExperiment)
library(scater)
options(stringsAsFactors = FALSE)
umi <- readRDS("tung/umi.rds")
umi.qc <- umi[rowData(umi)$use, colData(umi)$use]
endog_genes <- !rowData(umi.qc)$is_feature_control

tmp <- runPCA(
  umi[endog_genes, ],
  exprs_values = "counts"
)
plotPCA(
  tmp,
  colour_by = "batch",
  size_by = "total_features_by_counts",
  shape_by = "individual"
)

# ###WARNING Highly recommended to not use logcounts_raw for downstream analysis
# However, note that just a log-transformation is not enough to account for different technical factors 
# between the cells (e.g. sequencing depth). Therefore, please do not use logcounts_raw for 
# your downstream analysis, instead as a minimum suitable data use the logcounts slot of the 
# SingleCellExperiment object, which not just log-transformed, but also normalised by library size 
# (e.g. CPM normalisation). In the course we use logcounts_raw only for demonstration purposes!


tmp <- runPCA(
  umi[endog_genes, ],
  exprs_values = "logcounts_raw"
)

plotPCA(
  tmp,
  colour_by = "batch",
  size_by = "total_features_by_counts",
  shape_by = "individual"
)

# Clearly log-transformation is benefitial for our data - it reduces the variance on the first principal 
# component and already separates some biological effects. Moreover, it makes the distribution of the 
# expression values more normal. In the following analysis and chapters we will be using log-transformed 
# raw counts by default.

tmp <- runPCA(
  umi.qc[endog_genes, ],
  exprs_values = "logcounts_raw"
)

plotPCA(
  tmp,
  colour_by = "batch",
  size_by = "total_features_by_counts",
  shape_by = "individual"
)

##By default plotPCA takes top 500 genes. We can set the number of genes by "ntop" args as demonstrated below
## based on ntop value, the plotting time of PCA varies
tmp <- runPCA(
  umi.qc[endog_genes, ],
  ntop = 10000,
  exprs_values = "logcounts_raw"
)

plotPCA(
  tmp,
  colour_by = "batch",
  size_by = "total_features_by_counts",
  shape_by = "individual"
)

install.packages("Rtsne")
library("Rtsne")

####WARNING tsne is stochastic so if you want same plot everywhere, set seed #####
# Furthermore tSNE requires you to provide a value of perplexity which reflects the 
# number of neighboursused to build the nearest-neighbour network; a high value 
# creates a dense network which clumps cells together while a low value makes the
# network more sparse allowing groups of cells to separate from each other. scater
# uses a default perplexity of the total number of cells divided by five (rounded down).

###HIGHLY RECOMMEND #####
# You can read more about the pitfalls of using tSNE here -> http://distill.pub/2016/misread-tsne/

tmp <- runTSNE(
  umi[endog_genes, ],
  exprs_values = "logcounts_raw",
  perplexity = 130
)

plotTSNE(
  tmp,
  colour_by = "batch",
  size_by = "total_features_by_counts",
  shape_by = "individual"
)

tmp <- runTSNE(
  umi.qc[endog_genes, ],
  exprs_values = "logcounts_raw",
  perplexity = 130
)
plotTSNE(
  tmp,
  colour_by = "batch",
  size_by = "total_features_by_counts",
  shape_by = "individual"
)