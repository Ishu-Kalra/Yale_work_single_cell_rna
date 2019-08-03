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

molecules <- read.table("/Users/ishukalra/OneDrive - Indian Institute of Technology Guwahati/YaleWork/data/tung/molecules.txt", sep = "\t")
anno <- read.table("/Users/ishukalra/OneDrive - Indian Institute of Technology Guwahati/YaleWork/data/tung/annotation.txt", sep = "\t", header = TRUE)

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
