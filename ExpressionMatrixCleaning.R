###WARNING Update R to 3.6 (latest as of now) to avoid any errors###

###SOURCE: https://scrnaseq-course.cog.sanger.ac.uk/website/cleaning-the-expression-matrix.html
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("scater")
a
no

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("DelayedArray")
a
no

#####WARNING: RStudio > 3.6 required for scater ##########
library("DelayedArray")
library("SingleCellExperiment")
library("scater")
options(stringsAsFactors = FALSE)

molecules <- read.table("tung/molecules.txt", sep = "\t")
anno <- read.table("tung/annotation.txt", sep = "\t", header = TRUE)

head(molecules[ , 1:3])
head(anno)

#//The data consists of 3 individuals and 3 replicates and therefore has 9 batches in total.
#//
#//We standardize the analysis by using both SingleCellExperiment (SCE) and scater packages. First, create the SCE object:

umi <- SingleCellExperiment(
  assays = list(counts = as.matrix(molecules)),
  colData = anno
)

keep_feature <- rowSums(counts(umi) > 0) > 0
umi <- umi[keep_feature, ]

##isSpike for ERCC
isSpike(umi, "ERCC") <- grepl("^ERCC-", rownames(umi))

##isSpike for mitochondiral genes as provided by authors
isSpike(umi, "MT") <- rownames(umi) %in%  c("ENSG00000198899", "ENSG00000198727", "ENSG00000198888","ENSG00000198886", "ENSG00000212907", "ENSG00000198786","ENSG00000198695", "ENSG00000198712", "ENSG00000198804","ENSG00000198763", "ENSG00000228253", "ENSG00000198938","ENSG00000198840")

##Quality Metrics
umi <- calculateQCMetrics(
  umi,
  feature_controls = list(
    ERCC = isSpike(umi, "ERCC"), 
    MT = isSpike(umi, "MT")
  )
)

#Next we consider the total number of RNA molecules detected per sample (if we were using read counts rather than UMI counts this would be the total number of reads). Wells with few reads/molecules are likely to have been broken or failed to capture a cell, and should thus be removed.

hist(
  umi$total_counts,
  breaks = 100
)
abline(v = 25000, col = "red")

filter_by_total_counts <- (umi$total_counts > 25000)
table(filter_by_total_counts)


#In addition to ensuring sufficient sequencing depth for each sample, we also want to make sure that the reads are distributed across the transcriptome. Thus, we count the total number of unique genes detected in each sample.

hist(
  umi$total_features_by_counts,
  breaks = 100
)
abline(v = 7000, col = "red")


#From the plot we conclude that most cells have between 7,000-10,000 detected genes, which is normal for high-depth scRNA-seq. However, this varies by experimental protocol and sequencing depth. For example, droplet-based methods or samples with lower sequencing-depth typically detect fewer genes per cell. The most notable feature in the above plot is the “heavy tail” on the left hand side of the distribution. If detection rates were equal across the cells then the distribution should be approximately normal. Thus we remove those cells in the tail of the distribution (fewer than 7,000 detected genes).


filter_by_expr_features <- (umi$total_features_by_counts > 7000)
table(filter_by_expr_features)

#Another measure of cell quality is the ratio between ERCC spike-in RNAs and endogenous RNAs. This ratio can be used to estimate the total amount of RNA in the captured cells. Cells with a high level of spike-in RNAs had low starting amounts of RNA, likely due to the cell being dead or stressed which may result in the RNA being degraded.

plotColData(
  umi,
  x = "total_features_by_counts",
  y = "pct_counts_MT",
  colour = "batch"
)

plotColData(
  umi,
  x = "total_features_by_counts",
  y = "pct_counts_ERCC",
  colour = "batch"
)

filter_by_ERCC <- umi$batch != "NA19098.r2"
table(filter_by_ERCC)
filter_by_MT <- umi$pct_counts_MT < 10
table(filter_by_MT)

#Using those values which we thresholded. This is manual thresholding
umi$use <- (
  # sufficient features (genes)
  filter_by_expr_features &
    # sufficient molecules counted
    filter_by_total_counts &
    # sufficient endogenous RNA
    filter_by_ERCC &
    # remove cells with unusual number of reads in MT genes
    filter_by_MT
)
table(umi$use)

install.packages("mvoutlier")
install.packages("mvtnorm")

library("mvtnorm")
no
library("mvoutlier")

install.packages("rjags")
library("rjags")
####WARNING: Iif it throws error that "rjags" is not installed. If on mac: run "brew install jags" and checkout######
##https://gist.github.com/casallas/8411082

#Another option available in scater is to conduct PCA on a set of QC metrics and then use automatic outlier detection to identify potentially problematic cells.

#By default, the following metrics are used for PCA-based outlier detection:

#pct_counts_top_100_features
#total_features
#pct_counts_feature_controls
#n_detected_feature_controls
#log10_counts_endogenous_features
#log10_counts_feature_controls
#scater first creates a matrix where the rows represent cells and the columns represent the different QC metrics. Then, outlier cells can also be identified by using the mvoutlier package on the QC metrics for all cells. This will identify cells that have substantially different QC metrics from the others, possibly corresponding to low-quality cells. We can visualize any outliers using a principal components plot as shown below:

umi <- runPCA(
  umi, 
  use_coldata = TRUE, 
  detect_outliers = TRUE
)
reducedDimNames(umi)

table(umi$outlier)

#Then, we can use a PCA plot to see a 2D representation of the cells ordered by their quality metrics.

plotReducedDim(
  umi,
  use_dimred = "PCA_coldata",
  size_by = "total_features_by_counts", 
  shape_by = "use", 
  colour_by = "outlier"
)


if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("limma")
a
no

library("limma")
auto <- colnames(umi)[umi$outlier]
man <- colnames(umi)[!umi$use]
venn.diag <- vennCounts(
  cbind(colnames(umi) %in% auto,
        colnames(umi) %in% man)
)

##Comapring manual and automatic filtering.

vennDiagram(
  venn.diag,
  names = c("Automatic", "Manual"),
  circle.col = c("blue", "green")
)

###Takes time. May not show
#In addition to removing cells with poor quality, it is usually a good idea to exclude genes where we suspect that technical artefacts may have skewed the results. Moreover, inspection of the gene expression profiles may provide insights about how the experimental procedures could be improved.

#It is often instructive to consider the number of reads consumed by the top 50 expressed genes.

plotHighestExprs(umi, exprs_values = "counts")

#It is typically a good idea to remove genes whose expression level is considered “undetectable”. We define a gene as detectable if at least two cells contain more than 1 transcript from the gene. If we were considering read counts rather than UMI counts a reasonable threshold is to require at least five reads in at least two cells. However, in both cases the threshold strongly depends on the sequencing depth. It is important to keep in mind that genes must be filtered after cell filtering since some genes may only be detected in poor quality cells (note colData(umi)$use filter applied to the umi dataset).

keep_feature <- nexprs(
  umi[,colData(umi)$use], 
  byrow = TRUE, 
  detection_limit = 1
) >= 2
rowData(umi)$use <- keep_feature

table(keep_feature)
dim(umi[rowData(umi)$use, colData(umi)$use])

assay(umi, "logcounts_raw") <- log2(counts(umi) + 1)
reducedDim(umi) <- NULL

saveRDS(umi, file = "tung/umi.rds")

library(SingleCellExperiment)
library(scater)
options(stringsAsFactors = FALSE)

reads <- read.table("tung/reads.txt", sep = "\t")
anno <- read.table("tung/annotation.txt", sep = "\t", header = TRUE)
head(reads[ , 1:3])
head(anno)

reads <- SingleCellExperiment(
  assays = list(counts = as.matrix(reads)), 
  colData = anno
)
keep_feature <- rowSums(counts(reads) > 0) > 0
reads <- reads[keep_feature, ]
isSpike(reads, "ERCC") <- grepl("^ERCC-", rownames(reads))
isSpike(reads, "MT") <- rownames(reads) %in% 
  c("ENSG00000198899", "ENSG00000198727", "ENSG00000198888",
    "ENSG00000198886", "ENSG00000212907", "ENSG00000198786",
    "ENSG00000198695", "ENSG00000198712", "ENSG00000198804",
    "ENSG00000198763", "ENSG00000228253", "ENSG00000198938",
    "ENSG00000198840")
reads <- calculateQCMetrics(
  reads,
  feature_controls = list(
    ERCC = isSpike(reads, "ERCC"), 
    MT = isSpike(reads, "MT")
  )
)

hist(
  reads$total_counts,
  breaks = 100
)
abline(v = 1.3e6, col = "red")

filter_by_total_counts <- (reads$total_counts > 1.3e6)
table(filter_by_total_counts)

hist(
  reads$total_features_by_counts,
  breaks = 100
)
abline(v = 7000, col = "red")

filter_by_expr_features <- (reads$total_features_by_counts > 7000)
table(filter_by_expr_features)

plotColData(
  reads,
  x = "total_features_by_counts",
  y = "pct_counts_MT",
  colour = "batch"
)

plotColData(
  reads,
  x = "total_features_by_counts",
  y = "pct_counts_ERCC",
  colour = "batch"
)

filter_by_ERCC <- 
  reads$batch != "NA19098.r2" & reads$pct_counts_ERCC < 25
table(filter_by_ERCC)

filter_by_MT <- reads$pct_counts_MT < 30
table(filter_by_MT)

reads$use <- (
  # sufficient features (genes)
  filter_by_expr_features &
    # sufficient molecules counted
    filter_by_total_counts &
    # sufficient endogenous RNA
    filter_by_ERCC &
    # remove cells with unusual number of reads in MT genes
    filter_by_MT
)
table(reads$use)

reads <- runPCA(
  reads,
  use_coldata = TRUE, 
  detect_outliers = TRUE
)
reducedDimNames(reads)

table(reads$outlier)

plotReducedDim(
  reads,
  use_dimred = "PCA_coldata",
  size_by = "total_features_by_counts", 
  shape_by = "use", 
  colour_by = "outlier"
)

library(limma)

auto <- colnames(reads)[reads$outlier]
man <- colnames(reads)[!reads$use]
venn.diag <- vennCounts(
  cbind(colnames(reads) %in% auto,
        colnames(reads) %in% man)
)
vennDiagram(
  venn.diag,
  names = c("Automatic", "Manual"),
  circle.col = c("blue", "green")
)

###Takes time. May not show
plotHighestExprs(reads, exprs_values = "counts")

keep_feature <- nexprs(
  reads[,colData(reads)$use], 
  byrow = TRUE, 
  detection_limit = 1
) >= 2
rowData(reads)$use <- keep_feature
table(keep_feature)

dim(reads[rowData(reads)$use, colData(reads)$use])
assay(reads, "logcounts_raw") <- log2(counts(reads) + 1)
reducedDim(reads) <- NULL

saveRDS(reads, file = "tung/reads.rds")

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

####################################################################################
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

##Counfounders

library(scRNA.seq.funcs)

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("RUVSeq")
a
no

library(RUVSeq)
library(scater)
library(SingleCellExperiment)
library(scran)
library(devtools)
install_github('theislab/kBET')
library(kBET)
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("sva")
a
no

library(devtools)
install_github("immunogenomics/harmony")
###If any error then checkout -> https://github.com/immunogenomics/harmony/issues/10
library(harmony)

library(sva) # Combat
library(edgeR)
options(stringsAsFactors = FALSE)
umi <- readRDS("tung/umi.rds")
umi.qc <- umi[rowData(umi)$use, colData(umi)$use]
endog_genes <- !rowData(umi.qc)$is_feature_control
erccs <- rowData(umi.qc)$is_feature_control

qclust <- quickCluster(umi.qc, min.size = 30)
umi.qc <- computeSumFactors(umi.qc, sizes = 15, clusters = qclust)
umi.qc <- normalize(umi.qc)

##Checkout https://scrnaseq-course.cog.sanger.ac.uk/website/cleaning-the-expression-matrix.html#dealing-with-confounders
## For Theory

###RUVg
ruvg <- RUVg(counts(umi.qc), erccs, k = 1)
assay(umi.qc, "ruvg1") <- log2(
  t(t(ruvg$normalizedCounts) / colSums(ruvg$normalizedCounts) * 1e6) + 1
)
ruvg <- RUVg(counts(umi.qc), erccs, k = 10)
assay(umi.qc, "ruvg10") <- log2(
  t(t(ruvg$normalizedCounts) / colSums(ruvg$normalizedCounts) * 1e6) + 1
)

#####RUVs

scIdx <- matrix(-1, ncol = max(table(umi.qc$individual)), nrow = 3)
tmp <- which(umi.qc$individual == "NA19098")
scIdx[1, 1:length(tmp)] <- tmp
tmp <- which(umi.qc$individual == "NA19101")
scIdx[2, 1:length(tmp)] <- tmp
tmp <- which(umi.qc$individual == "NA19239")
scIdx[3, 1:length(tmp)] <- tmp
cIdx <- rownames(umi.qc)
ruvs <- RUVs(counts(umi.qc), cIdx, k = 1, scIdx = scIdx, isLog = FALSE)
assay(umi.qc, "ruvs1") <- log2(
  t(t(ruvs$normalizedCounts) / colSums(ruvs$normalizedCounts) * 1e6) + 1
)
ruvs <- RUVs(counts(umi.qc), cIdx, k = 10, scIdx = scIdx, isLog = FALSE)
assay(umi.qc, "ruvs10") <- log2(
  t(t(ruvs$normalizedCounts) / colSums(ruvs$normalizedCounts) * 1e6) + 1
)

####Combat

###Warning -> data containing multiple experimental replicates rather than a balanced design ->
# so using mod1 to preserve biological variability will result in an error.

combat_data <- logcounts(umi.qc)
mod_data <- as.data.frame(t(combat_data))
# Basic batch removal
mod0 = model.matrix(~ 1, data = mod_data) 
# Preserve biological variability
mod1 = model.matrix(~ umi.qc$individual, data = mod_data) 
# adjust for total genes detected
mod2 = model.matrix(~ umi.qc$total_features_by_counts, data = mod_data)

##Covariate is Basic
assay(umi.qc, "combat") <- ComBat(
  dat = t(mod_data), 
  batch = factor(umi.qc$batch), 
  mod = mod0,
  par.prior = TRUE,
  prior.plots = FALSE
)

##Covariate is total features by counts. Accordingly covariates can be changed
assay(umi.qc, "combat_tf") <- ComBat(
  dat = t(mod_data), 
  batch = factor(umi.qc$batch), 
  mod = mod2,
  par.prior = TRUE,
  prior.plots = FALSE
)

# mnnCorrect (Haghverdi et al. 2017) assumes that each batch shares at least one biological condition
# with each other batch. Thus it works well for a variety of balanced experimental designs. However,
# the Tung data contains multiple replicates for each invidividual rather than balanced batches, thus
# we will normalized each individual separately. Note that this will remove batch effects between batches
# within the same individual but not the batch effects between batches in different individuals, due to
# the confounded experimental design.

##For any other type of data not similar to Tung, we don't need to normalize it separately.
# This is given below as balanced design


do_mnn <- function(data.qc) {
  batch1 <- logcounts(data.qc[, data.qc$replicate == "r1"])
  batch2 <- logcounts(data.qc[, data.qc$replicate == "r2"])
  batch3 <- logcounts(data.qc[, data.qc$replicate == "r3"])
  
  if (ncol(batch2) > 0) {
    x = mnnCorrect(
      batch1, batch2, batch3,  
      k = 20,
      sigma = 0.1,
      cos.norm.in = TRUE,
      svd.dim = 2
    )
    res1 <- x$corrected[[1]]
    res2 <- x$corrected[[2]]
    res3 <- x$corrected[[3]]
    dimnames(res1) <- dimnames(batch1)
    dimnames(res2) <- dimnames(batch2)
    dimnames(res3) <- dimnames(batch3)
    return(cbind(res1, res2, res3))
  } else {
    x = mnnCorrect(
      batch1, batch3,  
      k = 20,
      sigma = 0.1,
      cos.norm.in = TRUE,
      svd.dim = 2
    )
    res1 <- x$corrected[[1]]
    res3 <- x$corrected[[2]]
    dimnames(res1) <- dimnames(batch1)
    dimnames(res3) <- dimnames(batch3)
    return(cbind(res1, res3))
  }
}

indi1 <- do_mnn(umi.qc[, umi.qc$individual == "NA19098"])

indi2 <- do_mnn(umi.qc[, umi.qc$individual == "NA19101"])

indi3 <- do_mnn(umi.qc[, umi.qc$individual == "NA19239"])

assay(umi.qc, "mnn") <- cbind(indi1, indi2, indi3)

# For a balanced design: 
#assay(umi.qc, "mnn") <- mnnCorrect(
#    list(B1 = logcounts(batch1), B2 = logcounts(batch2), B3 = logcounts(batch3)),  
#    k = 20,
#    sigma = 0.1,
#    cos.norm = TRUE,
#    svd.dim = 2
#)

###GLM
# A general linear model is a simpler version of Combat. It can correct for batches while
# preserving biological effects if you have a balanced design. In a confounded/replicate 
# design biological effects will not be fit/preserved. Similar to mnnCorrect we could remove
# batch effects from each individual separately in order to preserve biological (and technical)
# variance between individuals. For demonstation purposes we will naively correct all cofounded batch effects:

glm_fun <- function(g, batch, indi) {
  model <- glm(g ~ batch + indi)
  model$coef[1] <- 0 # replace intercept with 0 to preserve reference batch.
  return(model$coef)
}
effects <- apply(
  logcounts(umi.qc), 
  1, 
  glm_fun, 
  batch = umi.qc$batch, 
  indi = umi.qc$individual
)
corrected <- logcounts(umi.qc) - t(effects[as.numeric(factor(umi.qc$batch)), ])
assay(umi.qc, "glm") <- corrected


##Performing GLM correction for each slot separately

glm_fun1 <- function(g, batch) {
  model <- glm(g ~ batch)
  model$coef[1] <- 0 # replace intercept with 0 to preserve reference batch.
  return(model$coef)
}
do_glm <- function(data.qc) {
  effects <- apply(
    logcounts(data.qc), 
    1, 
    glm_fun1, 
    batch = data.qc$batch
  )
  corrected <- logcounts(data.qc) - t(effects[as.numeric(factor(data.qc$batch)), ])
  return(corrected)
}
indi1 <- do_glm(umi.qc[, umi.qc$individual == "NA19098"])
indi2 <- do_glm(umi.qc[, umi.qc$individual == "NA19101"])
indi3 <- do_glm(umi.qc[, umi.qc$individual == "NA19239"])
assay(umi.qc, "glm_indi") <- cbind(indi1, indi2, indi3)


#### Harmony ####
# Harmony [Korsunsky2018fast] is a newer batch correction method, which is designed to operate on PC space.
# The algorithm proceeds to iteratively cluster the cells, with the objective function formulated to promote
# cells from multiple datasets within each cluster. Once a clustering is obtained, the positions of the
# centroids of each dataset are obtained on a per-cluster basis and the coordinates are corrected. This
# procedure is iterated until convergence. Harmony comes with a theta parameter that controls the degree 
# of batch correction (higher values lead to more dataset integration), and can account for multiple
# experimental and biological factors on input.
# 
# Seeing how the end result of Harmony is an altered dimensional reduction space created on the basis
# of PCA, we plot the obtained manifold here and exclude it from the rest of the follow-ups in the section.

umi.qc.endog = umi.qc[endog_genes,]
umi.qc.endog = runPCA(umi.qc.endog, exprs_values = 'logcounts', ncomponents = 20)
pca <- as.matrix(umi.qc.endog@reducedDims@listData[["PCA"]])
harmony_emb <- HarmonyMatrix(pca, umi.qc.endog$batch, theta=2, do_pca=FALSE)

umi.qc.endog@reducedDims@listData[['harmony']] <- harmony_emb
plotReducedDim(
  umi.qc.endog,
  use_dimred = 'harmony',
  colour_by = "batch",
  size_by = "total_features_by_counts",
  shape_by = "individual"
)

## How to evaluate and compare confounder removal strategies

# A key question when considering the different methods for removing confounders is how to 
# quantitatively determine which one is the most effective. The main reason why comparisons
# are challenging is because it is often difficult to know what corresponds to technical
# counfounders and what is interesting biological variability. Here, we consider three 
# different metrics which are all reasonable based on our knowledge of the experimental design.
# Depending on the biological question that you wish to address, it is important to choose a
# metric that allows you to evaluate the confounders that are likely to be the biggest concern
# for the given situation.

####   Effectiveness 1
# We evaluate the effectiveness of the normalization by inspecting the PCA plot where colour
# corresponds the technical replicates and shape corresponds to different biological samples
# (individuals). Separation of biological samples and interspersed batches indicates that 
# technical variation has been removed. We always use log2-cpm normalized data to match the
# assumptions of PCA.


## Takes time. Be patient
for(n in assayNames(umi.qc)) {
  tmp <- runPCA(
    umi.qc[endog_genes, ],
    exprs_values = n
  )
  print(
    plotPCA(
      tmp,
      colour_by = "batch",
      size_by = "total_features_by_counts",
      shape_by = "individual"
    ) +
      ggtitle(n)
  )
}

#### Effectiveness 2

#We can also examine the effectiveness of correction using the relative log expression (RLE)
# across cells to confirm technical noise has been removed from the dataset. Note RLE only
# evaluates whether the number of genes higher and lower than average are equal for each cell
# - i.e. systemic technical effects. Random technical noise between batches may not be detected by RLE.

###### BUG ############
res <- list()
for(n in assayNames(umi.qc)) {
  res[[n]] <- suppressWarnings(calc_cell_RLE(assay(umi.qc, n), erccs))
}
par(mar=c(6,4,1,1))
boxplot(res, las=2)



#We can repeat the analysis from Chapter 12 to check whether batch effects have been removed.

#### BUG ###
library(scater)
for(n in assayNames(umi.qc)) {
  print(
    plotQC(
      umi.qc[endog_genes, ],
      type = "expl",
      exprs_values = n,
      variables = c(
        "total_features",
        "total_counts",
        "batch",
        "individual",
        "pct_counts_ERCC",
        "pct_counts_MT"
      )
    ) +
      ggtitle(n)
  )
}

###############################

##Effectiveness 3

# Another method to check the efficacy of batch-effect correction is to consider the intermingling of points
# from different batches in local subsamples of the data. If there are no batch-effects then proportion of
# cells from each batch in any local region should be equal to the global proportion of cells in each batch. 
# 
# `kBET` [@Buttner2017-ds] takes `kNN` networks around random cells and tests the number of cells from each
# batch against a binomial distribution. The rejection rate of these tests indicates the severity of 
# batch-effects still present in the data (high rejection rate = strong batch effects). `kBET` assumes each
# batch contains the same complement of biological groups, thus it can only be applied to the entire dataset
# if a perfectly balanced design has been used. However, `kBET` can also be applied to replicate-data if it
# is applied to each biological group separately. In the case of the Tung data, we will apply `kBET` to each
# individual independently to check for residual batch effects. However, this method will not identify residual
# batch-effects which are confounded with biological conditions. In addition, `kBET` does not determine if
# biological signal has been preserved. 

compare_kBET_results <- function(sce){
  indiv <- unique(sce$individual)
  norms <- assayNames(sce) # Get all normalizations
  results <- list()
  for (i in indiv){ 
    for (j in norms){
      tmp <- kBET(
        df = t(assay(sce[,sce$individual== i], j)), 
        batch = sce$batch[sce$individual==i], 
        heuristic = TRUE, 
        verbose = FALSE, 
        addTest = FALSE, 
        plot = FALSE)
      results[[i]][[j]] <- tmp$summary$kBET.observed[1]
    }
  }
  return(as.data.frame(results))
}

eff_debatching <- compare_kBET_results(umi.qc)

require("reshape2")
require("RColorBrewer")
# Plot results
dod <- melt(as.matrix(eff_debatching),  value.name = "kBET")
colnames(dod)[1:2] <- c("Normalisation", "Individual")

colorset <- c('gray', brewer.pal(n = 9, "RdYlBu"))

ggplot(dod, aes(Normalisation, Individual, fill=kBET)) +  
  geom_tile() +
  scale_fill_gradient2(
    na.value = "gray",
    low = colorset[2],
    mid=colorset[6],
    high = colorset[10],
    midpoint = 0.5, limit = c(0,1)) +
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) + 
  theme(
    axis.text.x = element_text(
      angle = 45, 
      vjust = 1, 
      size = 12, 
      hjust = 1
    )
  ) + 
  ggtitle("Effect of batch regression methods per individual")

## Dealing with confounders (reads)

library(scRNA.seq.funcs)
library(RUVSeq)
library(scater)
library(SingleCellExperiment)
library(scran)
library(kBET)
library(sva) # Combat
library(harmony)
library(edgeR)
options(stringsAsFactors = FALSE)
reads <- readRDS("tung/reads.rds")
reads.qc <- reads[rowData(reads)$use, colData(reads)$use]
endog_genes <- !rowData(reads.qc)$is_feature_control
erccs <- rowData(reads.qc)$is_feature_control

qclust <- quickCluster(reads.qc, min.size = 30)
reads.qc <- computeSumFactors(reads.qc, sizes = 15, clusters = qclust)
reads.qc <- normalize(reads.qc)

ruvg <- RUVg(counts(reads.qc), erccs, k = 1)
assay(reads.qc, "ruvg1") <- log2(
  t(t(ruvg$normalizedCounts) / colSums(ruvg$normalizedCounts) * 1e6) + 1
)
ruvg <- RUVg(counts(reads.qc), erccs, k = 10)
assay(reads.qc, "ruvg10") <- log2(
  t(t(ruvg$normalizedCounts) / colSums(ruvg$normalizedCounts) * 1e6) + 1
)
scIdx <- matrix(-1, ncol = max(table(reads.qc$individual)), nrow = 3)
tmp <- which(reads.qc$individual == "NA19098")
scIdx[1, 1:length(tmp)] <- tmp
tmp <- which(reads.qc$individual == "NA19101")
scIdx[2, 1:length(tmp)] <- tmp
tmp <- which(reads.qc$individual == "NA19239")
scIdx[3, 1:length(tmp)] <- tmp
cIdx <- rownames(reads.qc)
ruvs <- RUVs(counts(reads.qc), cIdx, k = 1, scIdx = scIdx, isLog = FALSE)
assay(reads.qc, "ruvs1") <- log2(
  t(t(ruvs$normalizedCounts) / colSums(ruvs$normalizedCounts) * 1e6) + 1
)
ruvs <- RUVs(counts(reads.qc), cIdx, k = 10, scIdx = scIdx, isLog = FALSE)
assay(reads.qc, "ruvs10") <- log2(
  t(t(ruvs$normalizedCounts) / colSums(ruvs$normalizedCounts) * 1e6) + 1
)
combat_data <- logcounts(reads.qc)
mod_data <- as.data.frame(t(combat_data))
# Basic batch removal
mod0 = model.matrix(~ 1, data = mod_data) 
# Preserve biological variability
mod1 = model.matrix(~ reads.qc$individual, data = mod_data) 
# adjust for total genes detected
mod2 = model.matrix(~ reads.qc$total_features_by_counts, data = mod_data)
assay(reads.qc, "combat") <- ComBat(
  dat = t(mod_data), 
  batch = factor(reads.qc$batch), 
  mod = mod0,
  par.prior = TRUE,
  prior.plots = FALSE
)


