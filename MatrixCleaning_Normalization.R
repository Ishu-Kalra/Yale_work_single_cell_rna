##Normalization  UMI
##This function was developed by Hemberg lab itself. 
install.packages("devtools")
devtools::install_github("hemberg-lab/scRNA.seq.funcs")
install.packages("scRNA.seq.funcs")
library(scRNA.seq.funcs)
library(scater)
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("scran")
a
no

library(scran)
options(stringsAsFactors = FALSE)
umi <- readRDS("tung/umi.rds")
umi.qc <- umi[rowData(umi)$use, colData(umi)$use]
endog_genes <- !rowData(umi.qc)$is_feature_control

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

logcounts(umi.qc) <- log2(calculateCPM(umi.qc, use_size_factors = FALSE) + 1)
plotPCA(
  umi.qc[endog_genes, ],
  colour_by = "batch",
  size_by = "total_features_by_counts",
  shape_by = "individual"
)

plotRLE(
  umi.qc[endog_genes, ], 
  exprs_values = "logcounts_raw",
  colour_by = "batch"
)

plotRLE(
  umi.qc[endog_genes, ], 
  exprs_values = "logcounts",
  colour_by = "batch"
)

##SCRAN

qclust <- quickCluster(umi.qc, min.size = 30)
umi.qc <- computeSumFactors(umi.qc, sizes = 15, clusters = qclust)
umi.qc <- normalize(umi.qc)

plotPCA(
  umi.qc[endog_genes, ],
  colour_by = "batch",
  size_by = "total_features_by_counts",
  shape_by = "individual"
)

plotRLE(
  umi.qc[endog_genes, ], 
  exprs_values = "logcounts_raw",
  colour_by = "batch"
)

plotRLE(
  umi.qc[endog_genes, ], 
  exprs_values = "logcounts",
  colour_by = "batch"
)

summary(sizeFactors(umi.qc))

##Downsampling of UMI Normalization
Down_Sample_Matrix <- 
  function (expr_mat) {
    min_lib_size <- min(colSums(expr_mat))
    down_sample <- function(x) {
      prob <- min_lib_size/sum(x)
      return(unlist(lapply(x, function(y) {
        rbinom(1, y, prob)
      })))
    }
    down_sampled_mat <- apply(expr_mat, 2, down_sample)
    return(down_sampled_mat)
  }

logcounts(umi.qc) <- log2(Down_Sample_Matrix(counts(umi.qc)) + 1)

plotPCA(
  umi.qc[endog_genes, ],
  colour_by = "batch",
  size_by = "total_features_by_counts",
  shape_by = "individual"
)

plotRLE(
  umi.qc[endog_genes, ], 
  exprs_values = "logcounts_raw",
  colour_by = "batch"
)

plotRLE(
  umi.qc[endog_genes, ], 
  exprs_values = "logcounts",
  colour_by = "batch"
)

# Some methods combine library size and fragment/gene length normalization such as:
#
# RPKM - Reads Per Kilobase Million (for single-end sequencing)
# FPKM - Fragments Per Kilobase Million (same as RPKM but for paired-end sequencing, makes sure that
# paired ends mapped to the same fragment are not counted twice)
# TPM - Transcripts Per Kilobase Million (same as RPKM, but the order of normalizations is reversed
# - length first and sequencing depth second)
# These methods are not applicable to our dataset since the end of the transcript which contains the
# UMI was preferentially sequenced. Furthermore in general these should only be calculated using appropriate
# quantification software from aligned BAM files not from read counts since often only a portion of the
# entire gene/transcript is sequenced, not the entire length. If in doubt check for a relationship between
# gene/transcript length and expression level.
#  
# However, here we show how these normalisations can be calculated using scater. First, we need to
# find the effective transcript length in Kilobases. However, our dataset containes only gene IDs,
# therefore we will be using the gene lengths instead of transcripts. scater uses the biomaRt package,
# which allows one to annotate genes by other attributes:
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("grimbough/biomaRt")
a
no

library(biomaRt)
install.packages("xml2")
##Press enter twice as it will retrive the data from net. Takes time so be patient
####DEBUG -> Added ensembleRedirect to avoid returning an invalid result when a character string of length 1 is expected###

#####  BUG  #####
umi.qc <- getBMFeatureAnnos(
  umi.qc,
  filters = "ensembl_gene_id", 
  attributes = c(
    "ensembl_gene_id",
    "hgnc_symbol",
    "chromosome_name",
    "start_position",
    "end_position"
  ), 
  biomart = "ENSEMBL_MART_ENSEMBL", 
  dataset = "hsapiens_gene_ensembl",
  host = "www.ensembl.org",
  #ensemblRedirect = FALSE,
)


# If you have mouse data, change the arguments based on this example:
# getBMFeatureAnnos(
#     object,
#     filters = "ensembl_transcript_id",
#     attributes = c(
#         "ensembl_transcript_id",
#         "ensembl_gene_id", 
#         "mgi_symbol",
#         "chromosome_name",
#         "transcript_biotype",
#         "transcript_start",
#         "transcript_end",
#         "transcript_count"
#     ),
#     biomart = "ENSEMBL_MART_ENSEMBL",
#     dataset = "mmusculus_gene_ensembl",
#     host = "www.ensembl.org"
# )

umi.qc.ann <- umi.qc[!is.na(rowData(umi.qc)$ensembl_gene_id), ]

eff_length <- abs(rowData(umi.qc.ann)$end_position - rowData(umi.qc.ann)$start_position) / 1000
plot(eff_length, rowMeans(counts(umi.qc.ann)))

tpm(umi.qc.ann) <- log2(calculateTPM(umi.qc.ann, eff_length) + 1)

tmp <- runPCA(
  umi.qc.ann,
  exprs_values = "tpm",
)
plotPCA(
  tmp,
  colour_by = "batch",
  size_by = "total_features_by_counts",
  shape_by = "individual"
)

tpm(umi.qc.ann) <- log2(calculateFPKM(umi.qc.ann, eff_length) + 1)

tmp <- runPCA(
  umi.qc.ann,
  exprs_values = "tpm",
)
plotPCA(
  tmp,
  colour_by = "batch",
  size_by = "total_features_by_counts",
  shape_by = "individual"
)

# Note The PCA looks for differences between cells. Gene length is the same across cells
# for each gene thus FPKM is almost identical to the CPM plot (it is just rotated) since
# it performs CPM first then normalizes gene length. Whereas, TPM is different because it
# weights genes by their length before performing CPM.

###############################################
####TODO Normalization for reads#####
###############################################

# Technical confounders (aka batch effects) can arise from difference in reagents, isolation
# methods, the lab/experimenter who performed the experiment, even which day/time the experiment
# was performed. Accounting for technical confounders, and batch effects particularly, is a
# large topic that also involves principles of experimental design. Here we address approaches
# that can be taken to account for confounders when the experimental design is appropriate.
# 
# Fundamentally, accounting for technical confounders involves identifying and, ideally, removing
# sources of variation in the expression data that are not related to (i.e. are confounding) the
# biological signal of interest. Various approaches exist, some of which use spike-in or housekeeping
# genes, and some of which use endogenous genes.
