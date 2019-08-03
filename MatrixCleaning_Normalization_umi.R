library(scater, quietly = TRUE)
options(stringsAsFactors = FALSE)
umi <- readRDS("tung/umi.rds")
umi.qc <- umi[rowData(umi)$use, colData(umi)$use]
endog_genes <- !rowData(umi.qc)$is_feature_control

# scater allows one to identify principal components that correlate with experimental and QC variables of interest (it ranks principle components by R^2 from a linear model regressing PC value against the variable of interest).

logcounts(umi.qc) <- assay(umi.qc, "logcounts_raw")
plotExplanatoryPCs(
  umi.qc[endog_genes, ],
  variables = "total_features_by_counts"
)
logcounts(umi.qc) <- NULL

# Indeed, we can see that PC1 can be almost completely explained by the number of detected genes.
# In fact, it was also visible on the PCA plot above. This is a well-known issue in scRNA-seq and
# was described here. http://biorxiv.org/content/early/2015/12/27/025528

plotExplanatoryVariables(
  umi.qc[endog_genes, ],
  exprs_values = "logcounts_raw",
  variables = c(
    "total_features_by_counts",
    "total_counts",
    "batch",
    "individual",
    "pct_counts_ERCC",
    "pct_counts_MT"
  )
)

# This analysis indicates that the number of detected genes (again) and also the 
# sequencing depth (number of counts) have substantial explanatory power for many
# genes, so these variables are good candidates for conditioning out in a normalisation
# step, or including in downstream statistical models. Expression of ERCCs also appears
# to be an important explanatory variable and one notable feature of the above plot
# is that batch explains more than individual.

# In addition to correcting for batch, there are other factors that one may want to compensate for.
# As with batch correction, these adjustments require extrinsic information. One popular method is
# scLVM which allows you to identify and subtract the effect from processes such as cell-cycle or apoptosis.
# scLVM -> https://github.com/PMBio/scLVM
# In addition, protocols may differ in terms of their coverage of each transcript, their bias based
# on the average content of A/T nucleotides, or their ability to capture short transcripts. Ideally,
# we would like to compensate for all of these differences and biases.

###############################################
####TODO Identifying confounders for reads#####
###############################################

# Library sizes vary because scRNA-seq data is often sequenced on highly multiplexed platforms the total
# reads which are derived from each cell may differ substantially. Some quantification methods
# (eg. Cufflinks, RSEM) incorporated library size when determining gene expression estimates thus do
# not require this normalization.
# 
# Cufflinks -> http://cole-trapnell-lab.github.io/cufflinks/
# RSEM -> http://deweylab.github.io/RSEM/
#
# However, if another quantification method was used then library size must be corrected for by multiplying
# or dividing each column of the expression matrix by a “normalization factor” which is an estimate of
# the library size relative to the other cells. Many methods to correct for library size have been
# developped for bulk RNA-seq and can be equally applied to scRNA-seq (eg. UQ, SF, CPM, RPKM, FPKM, TPM).


##CPM -> Counts Per Million

calc_cpm <-
  function (expr_mat, spikes = NULL) 
  {
    norm_factor <- colSums(expr_mat[-spikes, ])
    return(t(t(expr_mat)/norm_factor)) * 10^6
  }

# One potential drawback of CPM is if your sample contains genes that are both very highly
# expressed and differentially expressed across the cells. In this case, the total molecules
# in the cell may depend of whether such genes are on/off in the cell and normalizing by 
# total molecules may hide the differential expression of those genes and/or falsely create
# differential expression for the remaining genes.

##RLE -> Relative Log Expression

calc_sf <-
  function (expr_mat, spikes = NULL) 
  {
    geomeans <- exp(rowMeans(log(expr_mat[-spikes, ])))
    SF <- function(cnts) {
      median((cnts/geomeans)[(is.finite(geomeans) & geomeans > 
                                0)])
    }
    norm_factor <- apply(expr_mat[-spikes, ], 2, SF)
    return(t(t(expr_mat)/norm_factor))
  }

# The size factor (SF) was proposed and popularized by DESeq (Anders and Huber 2010).
# First the geometric mean of each gene across all cells is calculated. The size factor
# for each cell is the median across genes of the ratio of the expression to the gene’s
# geometric mean. A drawback to this method is that since it uses the geometric mean only
# genes with non-zero expression values across all cells can be used in its calculation,
# making it unadvisable for large low-depth scRNASeq experiments. edgeR & scater call this
# method RLE for “relative log expression”.

##UQ -> Upper Quantile

calc_uq <-
  function (expr_mat, spikes = NULL) 
  {
    UQ <- function(x) {
      quantile(x[x > 0], 0.75)
    }
    uq <- unlist(apply(expr_mat[-spikes, ], 2, UQ))
    norm_factor <- uq/median(uq)
    return(t(t(expr_mat)/norm_factor))
  }

# The upperquartile (UQ) was proposed by (Bullard et al. 2010). Here each column is divided by
# the 75% quantile of the counts for each library. Often the calculated quantile is scaled by
# the median across cells to keep the absolute level of expression relatively consistent. A 
# drawback to this method is that for low-depth scRNASeq experiments the large number of undetected
# genes may result in the 75% quantile being zero (or close to it). This limitation can be overcome
# by generalizing the idea and using a higher quantile (eg. the 99% quantile is the default in
# scater) or by excluding zeros prior to calculating the 75% quantile.

##DOWNSAMPLING

Down_Sample_Matrix <-
  function (expr_mat) 
  {
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

# A final way to correct for library size is to downsample the expression matrix so that
# each cell has approximately the same total number of molecules. The benefit of this method
# is that zero values will be introduced by the down sampling thus eliminating any biases
# due to differing numbers of detected genes. However, the major drawback is that the process
# is not deterministic so each time the downsampling is run the resulting expression matrix
# is slightly different. Thus, often analyses must be run on multiple downsamplings to ensure
# results are robust.

##Effectiveness with Relative Log Expression

calc_cell_RLE <-
  function (expr_mat, spikes = NULL) 
  {
    RLE_gene <- function(x) {
      if (median(unlist(x)) > 0) {
        log((x + 1)/(median(unlist(x)) + 1))/log(2)
      }
      else {
        rep(NA, times = length(x))
      }
    }
    if (!is.null(spikes)) {
      RLE_matrix <- t(apply(expr_mat[-spikes, ], 1, RLE_gene))
    }
    else {
      RLE_matrix <- t(apply(expr_mat, 1, RLE_gene))
    }
    cell_RLE <- apply(RLE_matrix, 2, median, na.rm = T)
    return(cell_RLE)
  }

# to compare the efficiency of different normalization methods we will use visual
# inspection of PCA plots and calculation of cell-wise relative log expression via
# scater’s plotRLE() function. Namely, cells with many (few) reads have higher (lower)
# than median expression for most genes resulting in a positive (negative) RLE across
# the cell, whereas normalized cells have an RLE close to zero.

######WARNING######
# Note The RLE, TMM, and UQ size-factor methods were developed for bulk RNA-seq data and,
# depending on the experimental context, may not be appropriate for single-cell RNA-seq data,
# as their underlying assumptions may be problematically violated.