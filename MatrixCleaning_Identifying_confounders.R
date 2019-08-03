
library(scater)
options(stringsAsFactors = FALSE)
reads <- readRDS("tung/reads.rds")
reads.qc <- reads[rowData(reads)$use, colData(reads)$use]
endog_genes <- !rowData(reads.qc)$is_feature_control


tmp <- runPCA(
  reads[endog_genes, ],
  exprs_values = "counts"
)
plotPCA(
  tmp,
  colour_by = "batch",
  size_by = "total_features_by_counts",
  shape_by = "individual"
)

tmp <- runPCA(
  reads[endog_genes, ],
  exprs_values = "logcounts_raw"
)
plotPCA(
  tmp,
  colour_by = "batch",
  size_by = "total_features_by_counts",
  shape_by = "individual"
)

tmp <- runPCA(
  reads.qc[endog_genes, ],
  exprs_values = "logcounts_raw"
)
plotPCA(
  tmp,
  colour_by = "batch",
  size_by = "total_features_by_counts",
  shape_by = "individual"
)

tmp <- runTSNE(
  reads[endog_genes, ],
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
  reads.qc[endog_genes, ],
  exprs_values = "logcounts_raw",
  perplexity = 130
)
plotTSNE(
  tmp,
  colour_by = "batch",
  size_by = "total_features_by_counts",
  shape_by = "individual"
)