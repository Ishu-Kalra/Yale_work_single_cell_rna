if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("pcaMethods")
a
no
library(pcaMethods)

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("SC3")
a
no
library(SC3)

library(scater)
library(SingleCellExperiment)
install.packages("pheatmap")
library(pheatmap)
library(mclust)

deng <- readRDS("/Users/ishukalra/YaleWork/Pipeline_crude/deng/deng-reads.rds")
deng

kidney <- readRDS("/Users/ishukalra/YaleWork/Pipeline_crude/kidney_droplet.rds")
table(colData(deng)$cell_type2)
plotPCA(deng)
  
#Estimating the number of SC3 clusters
deng <- sc3_estimate_k(deng)
metadata(deng)$sc3$k_estimation
plotPCA(deng, colour_by = "X1")

deng <- sc3(deng, ks = 10, biology = TRUE, n_cores = 1)

#SC3 result consists of several different outputs (please look in (Kiselev et al. 2017) and
#SC3 vignette for more details). 
#Thing to look at -> https://scrnaseq-course.cog.sanger.ac.uk/website/biological-analysis.html#ref-Kiselev2016-bq
#Another thing to look at -> http://bioconductor.org/packages/release/bioc/vignettes/SC3/inst/doc/my-vignette.html

sc3_plot_consensus(deng, k = 10, show_pdata = "cell_type2") #consensus matrix
sc3_plot_silhouette(deng, k = 10) #Silhoutte plot
sc3_plot_expression(deng, k = 10, show_pdata = "cell_type2") #heat map
sc3_plot_markers(deng, k = 10, show_pdata = "cell_type2")  #Identified marker genes
plotPCA(deng, colour_by = "sc3_10_clusters") #PCA plot with highlighted SC3 clusters

#Compare the results of SC3 clustering with the original publication cell type labels:
##Note To compare two sets of clustering labels we can use adjusted Rand index. The 
# index is a measure of the similarity between two data clusterings. Values of the 
# adjusted Rand index lie in [0, 1] interval, where 1means that two clusterings are
# identical and 0 means the level of similarity expected by chance.
adjustedRandIndex(colData(deng)$cell_type2, colData(deng)$sc3_10_clusters)

##Opens in web browser
sc3_interactive(deng)
##Press ESC in RStudio to continue

##tSNE + kMeans <- Stochastic so it is better to run multiple times to get ccurate results
deng <- runTSNE(deng, rand_seed = 1)
plotTSNE(deng)
# Note that all points on the plot above are black. This is different from what we saw before,
# when the cells were coloured based on the annotation. Here we do not have any annotation
# and all cells come from the same batch, therefore all dots are black.

# k = 8 is a hyperparamaeter. You can adjust that
colData(deng)$tSNE_kmeans <- as.character(kmeans(deng@reducedDims$TSNE, centers = 8)$clust)
plotTSNE(deng, colour_by = "tSNE_kmeans")

##Comparing the result with the original publication
colData(deng)$tSNE_kmeans <- as.character(kmeans(deng@reducedDims$TSNE, centers = 10)$clust)
adjustedRandIndex(colData(deng)$cell_type2, colData(deng)$tSNE_kmeans)

##SINCERA
# SINCERA is based on hierarchical clustering. One important thing to keep in mind is that it
# performs a gene-level z-score transformation before doing clustering:

# use the same gene filter as in SC3
input <- logcounts(deng[rowData(deng)$sc3_gene_filter, ])

# perform gene-by-gene per-sample z-score transformation
dat <- apply(input, 1, function(y) scRNA.seq.funcs::z.transform.helper(y))
# hierarchical clustering
dd <- as.dist((1 - cor(t(dat), method = "pearson"))/2)
hc <- hclust(dd, method = "average")

# If the number of cluster is not known SINCERA can identify k as the minimum height of the
# hierarchical tree that generates no more than a specified number of singleton clusters
# (clusters containing only 1 cell)

num.singleton <- 0
kk <- 1
for (i in 2:dim(dat)[2]) {
  clusters <- cutree(hc, k = i)
  clustersizes <- as.data.frame(table(clusters))
  singleton.clusters <- which(clustersizes$Freq < 2)
  if (length(singleton.clusters) <= num.singleton) {
    kk <- i
  } else {
    break;
  }
}
cat(kk)

##Result of SINCERA on a heatmap
pheatmap(
  t(dat),
  cluster_cols = hc,
  cutree_cols = kk,
  kmeans_k = 100,
  show_rownames = FALSE
)

##Comparing the result with the original publication. For analysis from the scratch, we would not be using this function
colData(deng)$SINCERA <- as.character(cutree(hc, k = kk))
adjustedRandIndex(colData(deng)$cell_type2, colData(deng)$SINCERA)

##Feature Selection

library(scRNA.seq.funcs)
library(matrixStats)
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("M3Drop")
a
no
library(M3Drop)
library(RColorBrewer)
library(SingleCellExperiment)

# Single-cell RNASeq is capable of measuring the expression of many thousands of genes in every cell.
# However, in most situations only a portion of those will show a response to the biological condition
# of interest, e.g. differences in cell-type, drivers of differentiation, respond to an environmental
# stimulus. Most genes detected in a scRNASeq experiment will only be detected at different levels due
# to technical noise. One consequence of this is that technical noise and batch effects can obscure the
# biological signal of interest.
# 
# Thus, it is often advantageous to perform feature selection to remove those genes which only exhibit
# technical noise from downstream analysis. Not only does this generally increase the signal:noise ratio
# in the data; it also reduces the computational complexity of analyses, by reducing the total amount of
# data to be processed.
# 
# For scRNASeq data, we will be focusing on unsupervised methods of feature selection which don’t
# require any a priori information, such as cell-type labels or biological group, since they are not
# available, or may be unreliable, for many experiments. In contrast, differential expression
# can be considered a form of supervised feature selection since it uses the known biological label
# of each sample to identify features (i.e. genes) which are expressed at different levels across groups.

deng <- readRDS("data/deng/deng-reads.rds")
celltype_labs <- colData(deng)$cell_type2  ##Filtered cell type labels
cell_colors <- brewer.pal(max(3,length(unique(celltype_labs))), "Set3")


# Feature selection is performed after QC, however this data (deng dataset) has already been QCed so we can skip that
# step here. M3Drop contain two different feature selection methods “M3DropFeatureSelection” which is
# based on a Michaelis-Menten curve and is designed for full-transcript single-cell RNA-seq data (such
# as Smartseq2) and “NBumiFeatureSelectionCombinedDrop” which is based on a negative binomial model and
# is designed for UMI count data. We will demonstrate both on the Deng Smartseq2 data.

# M3Drop feature selection is runs direction on a normalized (but not log-transformed) expression matrix.
# This can be extracted from our SingleCellExperiment object using the command below.

expr_matrix <- M3Drop::M3DropConvertData(deng)
# This function is compatible with most single-cell RNA-seq analysis packages including: scater,
# SingleCellExperiment, monocle, and Seurat. It can also convert an existing expression matrix to
# the correct form (removing undetected genes & normalizing/delogging) if you specify whether the
# matrix is raw counts, or log transformed.
?M3Drop::M3DropConvertData

##Confirming that the conversion function has removed undected genes
nrow(counts(deng)) - nrow(expr_matrix)
summary( rowSums(counts(deng))[! rownames(counts(deng)) %in% rownames(expr_matrix) ] )

##Identifying genes as a null model
# There are two main approaches to unsupervised feature selection. The first is to identify genes which
# behave differently from a null model describing just the technical noise expected in the dataset.
# 
# If the dataset contains spike-in RNAs they can be used to directly model technical noise. However,
# measurements of spike-ins may not experience the same technical noise as endogenous transcripts
# (Svensson et al., 2017). In addition, scRNASeq experiments often contain only a small number of
# spike-ins which reduces our confidence in fitted model parameters.

## Highly variable genes
# The first method proposed to identify features in scRNASeq datasets was to identify highly variable
# genes (HVG). HVG assumes that if genes have large differences in expression across cells some of those
# differences are due to biological difference between the cells rather than technical noise. However,
# because of the nature of count data, there is a positive relationship between the mean expression of
# a gene and the variance in the read counts across cells. This relationship must be corrected for to
# properly identify HVGs.

#Plotting the means and variance in gene expression on a log scale
plot(
  rowMeans(expr_matrix), 
  rowVars(expr_matrix), 
  log="xy", 
  pch=16,
  xlab="Mean Expression", 
  ylab="Variance", 
  main=""
)

# A popular method to correct for the relationship between variance and mean expression was proposed
# by Brennecke et al.. To use the Brennecke method, we first normalize for library size then calculate
# the mean and the square coefficient of variation (variation divided by the squared mean expression).
# A quadratic curve is fit to the relationship between these two variables for the ERCC spike-in, and then
# a chi-square test is used to find genes significantly above the curve. This method is included in the
# M3Drop package as the Brennecke_getVariableGenes(counts, spikes) function. However, this dataset does
# not contain spike-ins so we will use the entire dataset to estimate the technical noise.
# 
# In the figure below the red curve is the fitted technical noise model and the dashed line is the 95% CI.
# Pink dots are the genes with significant biological variability after multiple-testing correction.

Brennecke_HVG <- BrenneckeGetVariableGenes(
  expr_matrix,
  fdr = 0.01,
  minBiolDisp = 0.5
)

# This function returns a matrix of significant genes as well as their estimated effect size
# (difference between observed and expected coefficient of variation), and their significance
# as raw p.values and FDR corrected q.values. For now we will just keep the names of the significant HVG genes.

HVG_genes <- Brennecke_HVG$Gene

length(HVG_genes)  ##How many genes were significant using Brennecke Method

##High dropout genes
# An alternative to finding HVGs is to identify genes with unexpectedly high numbers of zeros
# . The frequency of zeros, known as the “dropout rate”, is very closely related to expression
# level in scRNASeq data. Zeros are the dominant feature of single-cell RNASeq data, typically
# accounting for over half of the entries in the final expression matrix. These zeros predominantly
# result from the failure of mRNAs failing to be reversed transcribed (Andrews and Hemberg, 2016).
# Reverse transcription is an enzyme reaction thus can be modelled using the Michaelis-Menten equation:
#   
# P_dropout = 1 - S/(K + S)
# where  
# S is the mRNA concentration in the cell (we will estimate this as average expression) and  
# K is the Michaelis-Menten constant.
# Because the Michaelis-Menten equation is a convex non-linear function, genes which are differentially
# expression across two or more populations of cells in our dataset will be shifted up/right of the
# Michaelis-Menten model. MM Plot is demonstrated by the code below

K <- 49
S_sim <- 10^seq(from = -3, to = 4, by = 0.05) # range of expression values
MM <- 1 - S_sim / (K + S_sim)
plot(
  S_sim,
  MM,       ##Add log = "x" here to visualize on a log scale
  type = "l", 
  lwd = 3, 
  xlab = "Expression", 
  ylab = "Dropout Rate", 
  xlim = c(1,1000)
)
S1 <- 10 # Mean expression in population 1
P1 <- 1 - S1 / (K + S1) # Dropouts for cells in condition 1
S2 <- 750 # Mean expression in population 2
P2 <- 1 - S2 / (K + S2) # Dropouts for cells in condition 2
points(
  c(S1, S2),
  c(P1, P2), 
  pch = 16, 
  col = "grey85", 
  cex = 3
)
mix <- 0.5 # proportion of cells in condition 1
points(
  S1 * mix + S2 * (1 - mix), 
  P1 * mix + P2 * (1 - mix), 
  pch = 16, 
  col = "grey35", 
  cex = 3
)

##MM plot for different exp levels or for mixture
plot(
  S_sim, 
  MM, 
  type = "l", 
  lwd = 3, 
  xlab = "Expression", 
  ylab = "Dropout Rate", 
  xlim = c(1, 1000), 
  log = "x"
)
S1 <- 100
P1 <- 1 - S1 / (K + S1) # Expression & dropouts for cells in condition 1
S2 <- 1000
P2 <- 1 - S2 / (K + S2) # Expression & dropouts for cells in condition 2
points(
  c(S1, S2),
  c(P1, P2), 
  pch = 16, 
  col = "grey85", 
  cex = 3
)
mix <- 0.75 # proportion of cells in condition 1. Mixture
points(
  S1 * mix + S2 * (1 - mix), 
  P1 * mix + P2 * (1 - mix), 
  pch = 16, 
  col = "grey35", 
  cex = 3
)

# We use M3Drop to identify significant outliers to the right of the MM curve. We also apply
# 1% FDR multiple testing correction:

M3Drop_genes <- M3DropFeatureSelection(
  expr_matrix,
  mt_method = "fdr",
  mt_threshold = 0.01
)

M3Drop_genes <- M3Drop_genes$Gene

# An alternative method is contained in the M3Drop package that is tailored specifically for
# UMI-tagged data which generally contains many zeros resulting from low sequencing coverage
# in addition to those resulting from insufficient reverse-transcription. This model is the
# Depth-Adjusted Negative Binomial (DANB). This method describes each expression observation
# as a negative binomial model with a mean related to both the mean expression of the respective
# gene and the sequencing depth of the respective cell, and a variance related to the mean-expression
# of the gene.
# 
# This method is designed to model the raw counts in a dataset directly, and we can extract the appropriate
# matrix using the “NBumiConvertData” function similar to M3Drop. However, we have an extra step for fitting
# the model since that is the slowest step of the method and we are currently working on additional methods
# that can use this model information for other things (such as normalization, co-expression testing, highly
# variable gene detection).
# 
# This method includes a binomial test of the significance of each feature, but since the Deng data is not UM
# I counts the model does not fit the noise sufficiently and far too many genes will be called as significant.
# Thus we will take the top 1500 by effect size.

deng_int <- NBumiConvertData(deng)

DANB_fit <- NBumiFitModel(deng_int) # DANB is fit to the raw count matrix
# Perform DANB feature selection
DropFS <- NBumiFeatureSelectionCombinedDrop(DANB_fit, method="fdr", qval.thresh=0.01, suppress.plot=FALSE)

DANB_genes <- DropFS[1:1500,]$Gene

##Number of genes significant by DANB genes
print(nrow((DropFS)))

## Correlated Expression
# A completely different approach to feature selection is to use gene-gene correlations. This method is
# based on the idea that multiple genes will be differentially expressed between different cell-types or
# cell-states. Genes which are expressed in the same cell-population will be positively correlated with each
# other where as genes expressed in different cell-populations will be negatively correated with each other.
# Thus important genes can be identified by the magnitude of their correlation with other genes.
# 
# The limitation of this method is that it assumes technical noise is random and independent for each cell,
# thus shouldn’t produce gene-gene correlations, but this assumption is violated by batch effects which are
# generally systematic between different experimental batches and will produce gene-gene correlations. As a
# result it is more appropriate to take the top few thousand genes as ranked by gene-gene correlation than
# consider the significance of the correlations.

cor_feat <- M3Drop::corFS(expr_matrix)
Cor_genes <- names(cor_feat)[1:1500]

# Lastly, another common method for feature selection in scRNASeq data is to use PCA loadings. Genes with
# high PCA loadings are likely to be highly variable and correlated with many other variable genes, thus
# may be relevant to the underlying biology. However, as with gene-gene correlations PCA loadings tend to
# be susceptible to detecting systematic variation due to batch effects; thus it is recommended to plot the PCA
# results to determine those components corresponding to the biological variation rather than batch effects.

# PCA is typically performed on log-transformed expression data
pca <- prcomp(log(expr_matrix + 1) / log(2))

# plot projection
plot(
  pca$rotation[,1], 
  pca$rotation[,2], 
  pch = 16, 
  col = cell_colors[as.factor(celltype_labs)]
) 
# calculate loadings for components 1 and 2
score <- rowSums(abs(pca$x[,c(1,2)])) 
names(score) <- rownames(expr_matrix)
score <- score[order(-score)]
PCA_genes <- names(score[1:1500])

#Considerinng the top 5 principal componennts to find out which is the most biologically relevant.

#Exercise Answer
plot(
  pca$rotation[,2], 
  pca$rotation[,3], 
  pch = 16, 
  col = cell_colors[as.factor(celltype_labs)]
)
plot(
  pca$rotation[,3], 
  pca$rotation[,4], 
  pch = 16, 
  col = cell_colors[as.factor(celltype_labs)]
)
# calculate loadings for components 1 and 2
score <- rowSums(abs(pca$x[,c(2, 3, 4)]))
names(score) <- rownames(expr_matrix)
score <- score[order(-score)]
PCA_genes2 = names(score[1:1500])
### Comparing Methods

# We can check whether the identified features really do represent genes differentially expressed between
# cell-types in this dataset.

M3DropExpressionHeatmap(
  M3Drop_genes,
  expr_matrix,
  cell_labels = celltype_labs
)


# We can also consider how consistent each feature selection method is with the others using the Jaccard Index:
J <- sum(M3Drop_genes %in% HVG_genes)/length(unique(c(M3Drop_genes, HVG_genes)))


# __Exercise 7__
# 
# Plot the expression of the features for each of the other methods. Which appear to be differentially expressed? How consistent are the different methods for this dataset?

#Answer
M3DropExpressionHeatmap(
  HVG_genes,
  expr_matrix,
  cell_labels = celltype_labs
)

M3DropExpressionHeatmap(
  Cor_genes,
  expr_matrix,
  cell_labels = celltype_labs
)


M3DropExpressionHeatmap(
  PCA_genes,
  expr_matrix,
  cell_labels = celltype_labs
)

M3DropExpressionHeatmap(
  PCA_genes2,
  expr_matrix,
  cell_labels = celltype_labs
)


list_of_features <- list(
  M3Drop_genes, 
  HVG_genes, 
  Cor_genes, 
  PCA_genes, 
  PCA_genes2
)

Out <- matrix(
  0, 
  ncol = length(list_of_features), 
  nrow = length(list_of_features)
)

for(i in 1:length(list_of_features) ) {
  for(j in 1:length(list_of_features) ) {
    Out[i,j] <- sum(list_of_features[[i]] %in% list_of_features[[j]])/
      length(unique(c(list_of_features[[i]], list_of_features[[j]])))
  }
}
colnames(Out) <- rownames(Out) <- c("M3Drop", "HVG", "Cor", "PCA", "PCA2")


## Pseudotime analysis

library(SingleCellExperiment)
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("TSCAN")
a
no
library(TSCAN)
library(M3Drop)
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("monocle")
a
no
library(monocle)
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("destiny")
a
no
library(destiny)
library("devtools")
install_github("jw156605/SLICER")
library(SLICER)
devtools::install_github("kieranrcampbell/ouija")
library(ouija)
library(scater)
library(ggplot2)
library(ggthemes)
library(ggbeeswarm)
install.packages("corrplot")
library(corrplot)


# In many situations, one is studying a process where cells change
# continuously. This includes, for example, many differentiation processes
# taking place during development: following a stimulus, cells
# will change from one cell-type to another. Ideally, we would like to
# monitor the expression levels of an individual cell over
# time. Unfortunately, such monitoring is not possible with scRNA-seq
# since the cell is lysed (destroyed) when the RNA is extracted.
# 
# Instead, we must sample at multiple time-points and obtain snapshots
# of the gene expression profiles. Since some of the cells will proceed
# faster along the differentiation than others, each snapshot may
# contain cells at varying points along the developmental
# progression. We use statistical methods to order the cells along one
# or more trajectories which represent the underlying developmental
# trajectories, this ordering is referred to as "pseudotime".
# 
# In this chapter we will consider five different tools: Monocle, TSCAN,
# destiny, SLICER and ouija for ordering cells according to their pseudotime
# development. To illustrate the methods we will be using a dataset on
# mouse embryonic development [@Deng2014-mx]. The dataset consists of
# 268 cells from 10 different time-points of early mouse development. In this case, there is no need for pseudotime alignment since the cell labels provide information about the development trajectory. Thus, the labels allow us to establish a ground truth so that we can evaluate and compare the different methods.
# 
# A recent review by Cannoodt et al provides a detailed summary of the
# various computational methods for trajectory inference from
# single-cell transcriptomics [@Cannoodt2016-uj]. They discuss several
# tools, but unfortunately for our purposes many of these tools do not
# have complete or well-maintained implementations, and/or are not
# implemented in R.

# Cannoodt et al cover:
#   
#   * SCUBA - Matlab implementation
# * Wanderlust - Matlab (and requires registration to even download)
# * Wishbone - Python
# * SLICER - R, but package only available on Github
# * SCOUP - C++ command line tool
# * Waterfall - R, but one R script in supplement
# * Mpath - R pkg, but available as tar.gz on Github; function
# documentation but no vignette/workflow
# * Monocle - Bioconductor package
# * TSCAN - Bioconductor package
# 
# Unfortunately only two tools discussed (Monocle and TSCAN) meet the
# gold standard of open-source software hosted in a reputable repository.
# 
# The following figures from the paper summarise some of the features of
# the various tools.


### First look at Deng data

# Let us take a first look at the Deng data, without yet applying sophisticated pseudotime methods.
# As the plot below shows, simple PCA does a very good job of displaying the structure in these data.
# It is only once we reach the blast cell types ("earlyblast", "midblast", "lateblast") that PCA struggles
# to separate the distinct cell types.

deng_SCE <- readRDS("deng/deng-reads.rds")
deng_SCE$cell_type2 <- factor(
  deng_SCE$cell_type2,
  levels = c("zy", "early2cell", "mid2cell", "late2cell",
             "4cell", "8cell", "16cell", "earlyblast",
             "midblast", "lateblast")
)
cellLabels <- deng_SCE$cell_type2
deng <- counts(deng_SCE)
colnames(deng) <- cellLabels
deng_SCE <- runPCA(deng_SCE)
plotPCA(deng_SCE, colour_by = "cell_type2")


# PCA, here, provides a useful baseline for assessing different pseudotime methods. For a very
# naive pseudotime we can just take the co-ordinates of the first principal component.

deng_SCE$PC1 <- reducedDim(deng_SCE, "PCA")[,1]
ggplot(as.data.frame(colData(deng_SCE)), aes(x = PC1, y = cell_type2, 
                                             colour = cell_type2)) +
  geom_quasirandom(groupOnX = FALSE) +
  scale_color_tableau() + theme_classic() +
  xlab("First principal component") + ylab("Timepoint") +
  ggtitle("Cells ordered by first principal component")


# As the plot above shows, PC1 struggles to correctly order cells early and late in the developmental
# timecourse, but overall does a relatively good job of ordering cells by developmental time.   
# 
# Can bespoke pseudotime methods do better than naive application of PCA?
#   
  
  ### TSCAN
#   
#   TSCAN combines clustering with pseudotime analysis. First it clusters the cells using `mclust`,
#   which is based on a mixture of normal distributions. Then it builds a minimum spanning tree to 
#   connect the clusters. The branch of this tree that connects the largest number of clusters is the
#   main branch which is used to determine pseudotime.
# 
# First we will try to use all genes to order the cells.

procdeng <- TSCAN::preprocess(deng)
colnames(procdeng) <- 1:ncol(deng)
dengclust <- TSCAN::exprmclust(procdeng, clusternum = 10)
TSCAN::plotmclust(dengclust)
dengorderTSCAN <- TSCAN::TSCANorder(dengclust, orderonly = FALSE)
pseudotime_order_tscan <- as.character(dengorderTSCAN$sample_name)
deng_SCE$pseudotime_order_tscan <- NA
deng_SCE$pseudotime_order_tscan[as.numeric(dengorderTSCAN$sample_name)] <- 
  dengorderTSCAN$Pseudotime

# 
# Frustratingly, TSCAN only provides pseudotime values for 221 of 268 cells, silently returning
# missing values for non-assigned cells.
# 
# Again, we examine which timepoints have been assigned to each state:
  
#{r tscan-vs-truth}
cellLabels[dengclust$clusterid == 10]
ggplot(as.data.frame(colData(deng_SCE)), 
       aes(x = pseudotime_order_tscan, 
           y = cell_type2, colour = cell_type2)) +
  geom_quasirandom(groupOnX = FALSE) +
  scale_color_tableau() + theme_classic() +
  xlab("TSCAN pseudotime") + ylab("Timepoint") +
  ggtitle("Cells ordered by TSCAN pseudotime")



# TSCAN gets the development trajectory the "wrong way around", in the sense that later pseudotime 
# values correspond to early timepoints and vice versa. This is not inherently a problem (it is easy
# enough to reverse the ordering to get the intuitive interpretation of pseudotime), but overall it 
# would be a stretch to suggest that TSCAN performs better than PCA on this dataset. (As it is a 
# PCA-based method, perhaps this is not entirely surprising.)


# __Exercise 1__ Compare results for different numbers of clusters (`clusternum`).

### monocle

# Monocle skips the clustering stage of TSCAN and directly builds a
# minimum spanning tree on a reduced dimension representation of the
# cells to connect all cells. Monocle then identifies the longest path
# in this tree to determine pseudotime. If the data contains diverging
# trajectories (i.e. one cell type differentiates into two different
# cell-types), monocle can identify these. Each of the resulting forked paths is
# defined as a separate cell state.
# 
# Unfortunately, Monocle does not work when all the genes are used, so
# we must carry out feature selection. First, we use M3Drop:
# m3d-select-genes}
m3dGenes <- as.character(
  M3DropFeatureSelection(deng)$Gene
)
d <- deng[which(rownames(deng) %in% m3dGenes), ]
d <- d[!duplicated(rownames(d)), ]


# Now run monocle:
# monocle-all-genes, message=FALSE, warning=FALSE}
colnames(d) <- 1:ncol(d)
geneNames <- rownames(d)
rownames(d) <- 1:nrow(d)
pd <- data.frame(timepoint = cellLabels)
pd <- new("AnnotatedDataFrame", data=pd)
fd <- data.frame(gene_short_name = geneNames)
fd <- new("AnnotatedDataFrame", data=fd)

dCellData <- newCellDataSet(d, phenoData = pd, featureData = fd, expressionFamily = tobit())
dCellData <- setOrderingFilter(dCellData, which(geneNames %in% m3dGenes))
dCellData <- estimateSizeFactors(dCellData)
dCellDataSet <- reduceDimension(dCellData, pseudo_expr = 1)
dCellDataSet <- orderCells(dCellDataSet, reverse = FALSE)
plot_cell_trajectory(dCellDataSet)
# Store the ordering
pseudotime_monocle <-
  data.frame(
    Timepoint = phenoData(dCellDataSet)$timepoint,
    pseudotime = phenoData(dCellDataSet)$Pseudotime,
    State = phenoData(dCellDataSet)$State
  )
rownames(pseudotime_monocle) <- 1:ncol(d)
pseudotime_order_monocle <-
  rownames(pseudotime_monocle[order(pseudotime_monocle$pseudotime), ])


# We can again compare the inferred pseudotime to the known sampling timepoints.
#monocle-vs-truth}
deng_SCE$pseudotime_monocle <- pseudotime_monocle$pseudotime
ggplot(as.data.frame(colData(deng_SCE)), 
       aes(x = pseudotime_monocle, 
           y = cell_type2, colour = cell_type2)) +
  geom_quasirandom(groupOnX = FALSE) +
  scale_color_tableau() + theme_classic() +
  xlab("monocle pseudotime") + ylab("Timepoint") +
  ggtitle("Cells ordered by monocle pseudotime")

# Monocle - at least with its default settings - performs poorly on these data. The "late2cell" group
# is completely separated from the "zy", "early2cell" and "mid2cell" cells (though these are correctly
# ordered), and there is no separation at all of "4cell", "8cell", "16cell" or any blast cell groups.


### Diffusion maps

# [Diffusion maps](https://en.wikipedia.org/wiki/Diffusion_map) were introduced by [Ronald Coifman and
# Stephane Lafon](http://www.sciencedirect.com/science/article/pii/S1063520306000546), and the underlying
# idea is to assume that the data are samples from a diffusion process. The method infers the low-dimensional
# manifold by estimating the eigenvalues and eigenvectors for the diffusion operator related to the data.

# [Angerer et al](https://academic.oup.com/bioinformatics/article/32/8/1241/1744143) have applied the diffusion
# maps concept to the analysis of single-cell RNA-seq data to create an R package called [destiny]
# (http://bioconductor.org/packages/destiny).

# We will take the ranko prder of cells in the first diffusion map component as "diffusion map pseudotime" here.

#destiny-deng}
deng <- logcounts(deng_SCE)
colnames(deng) <- cellLabels
dm <- DiffusionMap(t(deng))

tmp <- data.frame(DC1 = eigenvectors(dm)[,1],
                  DC2 = eigenvectors(dm)[,2],
                  Timepoint = deng_SCE$cell_type2)
ggplot(tmp, aes(x = DC1, y = DC2, colour = Timepoint)) +
  geom_point() + scale_color_tableau() + 
  xlab("Diffusion component 1") + 
  ylab("Diffusion component 2") +
  theme_classic()

deng_SCE$pseudotime_diffusionmap <- rank(eigenvectors(dm)[,1])
ggplot(as.data.frame(colData(deng_SCE)), 
       aes(x = pseudotime_diffusionmap, 
           y = cell_type2, colour = cell_type2)) +
  geom_quasirandom(groupOnX = FALSE) +
  scale_color_tableau() + theme_classic() +
  xlab("Diffusion map pseudotime (first diffusion map component)") +
  ylab("Timepoint") +
  ggtitle("Cells ordered by diffusion map pseudotime")

# Like the other methods, using the first diffusion map component from destiny as pseudotime does a good job
# at ordering the early time-points (if we take high values as "earlier" in developement), but it is unable
# to distinguish the later ones.

# __Exercise 2__ Do you get a better resolution between the later time points by considering additional
# eigenvectors?
#   
#   __Exercise 3__ How does the ordering change if you only use the genes identified by M3Drop?

  ### SLICER
  
# The SLICER method is an algorithm for constructing trajectories that
# describe gene expression changes during a sequential biological
# process, just as Monocle and TSCAN are. SLICER is designed to capture
# highly nonlinear gene expression changes, automatically select genes
# related to the process, and detect multiple branch and loop features
# in the trajectory [@Welch2016-jr]. The SLICER R package is available
# from its [GitHub repository](https://github.com/jw156605/SLICER) and
# can be installed from there using the `devtools` package.
# 
# We use the `select_genes` function in SLICER to automatically select
# the genes to use in builing the cell trajectory. The function uses
# "neighbourhood variance" to identify genes that vary smoothly, rather
# than fluctuating randomly, across the set of cells. Following this, we
# determine which value of "k" (number of nearest neighbours) yields an embedding that
# most resembles a trajectory. Then we estimate the [locally linear
# embedding](https://en.wikipedia.org/wiki/Nonlinear_dimensionality_reduction) of the cells.

# slicer-analyis, message=FALSE, warning=FALSE
install.packages("lle")
library("lle")
slicer_genes <- select_genes(t(deng))
k <- select_k(t(deng[slicer_genes,]), kmin = 30, kmax=60)
slicer_traj_lle <- lle(t(deng[slicer_genes,]), m = 2, k)$Y
reducedDim(deng_SCE, "LLE") <- slicer_traj_lle
plotReducedDim(deng_SCE, use_dimred = "LLE", colour_by = "cell_type2") +
  xlab("LLE component 1") + ylab("LLE component 2") +
  ggtitle("Locally linear embedding of cells from SLICER")


# With the locally linear embedding computed we can construct a
# k-nearest neighbour graph that is fully connected. This plot displays
# a (yellow) circle for each cell, with the cell ID number overlaid in
# blue. Here we show the graph computed using 10 nearest
# neighbours. Here, SLICER appears to detect one major trajectory with
# one branch.

# slicer-build-graph
slicer_traj_graph <- conn_knn_graph(slicer_traj_lle, 10)
plot(slicer_traj_graph, main = "Fully connected kNN graph from SLICER")

# From this graph we can identify "extreme" cells that are candidates
# for start/end cells in the trajectory.

# slicer
ends <- find_extreme_cells(slicer_traj_graph, slicer_traj_lle)
start <- ends[1]


# Having defined a start cell we can order the cells in the estimated pseudotime.

pseudotime_order_slicer <- cell_order(slicer_traj_graph, start)
branches <- assign_branches(slicer_traj_graph, start)

pseudotime_slicer <-
  data.frame(
    Timepoint = cellLabels,
    pseudotime = NA,
    State = branches
  )
pseudotime_slicer$pseudotime[pseudotime_order_slicer] <-
  1:length(pseudotime_order_slicer)
deng_SCE$pseudotime_slicer <- pseudotime_slicer$pseudotime


# We can again compare the inferred pseudotime to the known sampling
# timepoints. SLICER does not provide a pseudotime value per se, just an
# ordering of cells.

# slicer-vs-truth}
ggplot(as.data.frame(colData(deng_SCE)), 
       aes(x = pseudotime_slicer, 
           y = cell_type2, colour = cell_type2)) +
  geom_quasirandom(groupOnX = FALSE) +
  scale_color_tableau() + theme_classic() +
  xlab("SLICER pseudotime (cell ordering)") +
  ylab("Timepoint") +
  theme_classic()

# Like the previous method, SLICER here provides a good ordering for the
# early time points. It places "16cell" cells before "8cell" cells, but provides better ordering
# for blast cells than many of the earlier methods.

# __Exercise 4__ How do the results change for different k? (e.g. k = 5) What about changing the number
# of nearest neighbours in the call to `conn_knn_graph`?
  
#   __Exercise 5__ How does the ordering change if you use a different set
# of genes from those chosen by SLICER (e.g. the genes identified by M3Drop)?
  
  ### Ouija
  
  # Ouija (http://kieranrcampbell.github.io/ouija/) takes a different approach from the pseudotime estimation
  # methods we have looked at so far. Earlier methods have all been "unsupervised", which is to say that apart
  # from perhaps selecting informative genes we do not supply the method with any prior information about how
  # we expect certain genes or the trajectory as a whole to behave. 
# 
# Ouija, in contrast, is a probabilistic framework that allows for interpretable learning of single-cell
# pseudotimes using only small panels of marker genes. This method:
#   
# * infers pseudotimes from a small number of marker genes letting you understand why the pseudotimes have
# been learned in terms of those genes;
# * provides parameter estimates (with uncertainty) for interpretable gene regulation behaviour (such as the
# peak time or the upregulation time); 
# * has a Bayesian hypothesis test to find genes regulated before others along the trajectory; 
# * identifies metastable states, ie discrete cell types along the continuous trajectory.
# 
# We will supply the following marker genes to Ouija (with timepoints where they are expected to be highly
# expressed):
#   
# * Early timepoints: Dazl, Rnf17, Sycp3, Nanog, Pou5f1, Fgf8, Egfr, Bmp5, Bmp15
# * Mid timepoints: Zscan4b, Foxa1, Prdm14, Sox21
# * Late timepoints: Creb3, Gpx4, Krt8, Elf5, Eomes, Cdx2, Tdgf1, Gdf3

# With Ouija we can model genes as either exhibiting monotonic up or down regulation (known as switch-like
# behaviour), or transient behaviour where the gene briefly peaks. By default, Ouija assumes all genes exhibit
# switch-like behaviour (the authors assure us not to worry if we get it wrong - the noise model means incorrectly
# specifying a transient gene as switch-like has minimal effect).

# Here we can "cheat" a little and check that our selected marker genes do actually identify different
# timepoints of the differentiation process.

# ouija-response-type, fig.height=11}
ouija_markers_down <- c("Dazl", "Rnf17", "Sycp3", "Fgf8", 
                        "Egfr", "Bmp5", "Bmp15", "Pou5f1")
ouija_markers_up <- c("Creb3", "Gpx4", "Krt8", "Elf5", "Cdx2", 
                      "Tdgf1", "Gdf3", "Eomes")
ouija_markers_transient <- c("Zscan4b", "Foxa1", "Prdm14", "Sox21")
ouija_markers <- c(ouija_markers_down, ouija_markers_up, 
                   ouija_markers_transient)
plotExpression(deng_SCE, ouija_markers, x = "cell_type2", colour_by = "cell_type2") +
  theme(axis.text.x = element_text(angle = 60, hjust = 1))


# In order to fit the pseudotimes wesimply call `ouija`, passing in the expected response types. Note that if
# no response types are provided then they are all assumed to be switch-like by default, which we will do here.
# The input to Ouija can be a cell-by-gene matrix of non-negative expression values, or an ExpressionSet object,
# or, happily, by selecting the `logcounts` values from a SingleCellExperiment object.

# We can apply prior information about whether genes are up- or down-regulated across the differentiation
# process, and also provide prior information about when the switch in expression or a peak in expression
# is likely to occur. 

# We can fit the Ouija model using either:
  
# * Hamiltonian Monte Carlo (HMC) - full MCMC inference where gradient information of the log-posterior is
# used to “guide” the random walk through the parameter space, or
# * Automatic Differentiation Variational Bayes (ADVI or simply VI) - approximate inference where the KL
# divergence to an approximate distribution is minimised.

# In general, HMC will provide more accurate inference with approximately correct posterior variance for all
# parameters. However, VB is orders of magnitude quicker than HMC and while it may underestimate posterior
# variance, the Ouija authors suggest that anecdotally it often performs as well as HMC for discovering posterior
# pseudotimes.

# To help the Ouija model, we provide it with prior information about the strength of switches for up- 
# and down-regulated genes. By setting switch strength to -10 for down-regulated genes and 10 for 
# up-regulated genes with a prior strength standard deviation of 0.5 we are telling the model that we
# are confident about the expected behaviour of these genes across the differentiation process.


# ouija-fit, warning=FALSE, message=FALSE, result='hide'}
options(mc.cores = parallel::detectCores())
response_type <- c(rep("switch", length(ouija_markers_down) + 
                         length(ouija_markers_up)), 
                   rep("transient", length(ouija_markers_transient)))
switch_strengths <- c(rep(-10, length(ouija_markers_down)),
                      rep(10, length(ouija_markers_up)))
switch_strength_sd <- c(rep(0.5, length(ouija_markers_down)),
                        rep(0.5, length(ouija_markers_up)))
garbage <- capture.output(
  oui_vb <- ouija(deng_SCE[ouija_markers,],
                  single_cell_experiment_assay = "logcounts", 
                  response_type = response_type,
                  switch_strengths = switch_strengths,
                  switch_strength_sd = switch_strength_sd,
                  inference_type = "vb")
)

print(oui_vb)


# We can plot the gene expression over pseudotime along with the maximum a posteriori (MAP)
# estimates of the mean function (the sigmoid or Gaussian transient function) using the 
# plot_expression function. 

plot_expression(oui_vb)


# We can also visualise when in the trajectory gene regulation behaviour occurs, either in the
# form of the switch time or the peak time (for switch-like or transient genes) using the
# plot_switch_times and plot_transient_times functions:
  
# ouija-plot-switch-times}
plot_switch_times(oui_vb)
plot_peak_times(oui_vb)


# Identify metastable states using consistency matrices.

# ouija-consistency}
cmo <- consistency_matrix(oui_vb)
plot_consistency(oui_vb)
cell_classifications <- cluster_consistency(cmo)


# ouija-pseudotime}
map_pst <- map_pseudotime(oui_vb)
ouija_pseudotime <- data.frame(map_pst, cell_classifications)

ggplot(ouija_pseudotime, aes(x = map_pst, y = cell_classifications)) +
  geom_point() +
  xlab("MAP pseudotime") +
  ylab("Cell classification")

deng_SCE$pseudotime_ouija <- ouija_pseudotime$map_pst
deng_SCE$ouija_cell_class <- ouija_pseudotime$cell_classifications

ggplot(as.data.frame(colData(deng_SCE)), 
       aes(x = pseudotime_ouija, 
           y = cell_type2, colour = cell_type2)) +
  geom_quasirandom(groupOnX = FALSE) +
  scale_color_tableau() + theme_classic() +
  xlab("Ouija pseudotime") +
  ylab("Timepoint") +
  theme_classic()


 
# Ouija does quite well in the ordering of the cells here, although it can be sensitive to the choice 
# of marker genes and prior information supplied. How do the results change if you select different marker
# genes or change the priors?

# Ouija identifies four metastable states here, which we might annotate as  "zygote/2cell", "4/8/16 cell",
# "blast1" and "blast2".

# ouija-states}
ggplot(as.data.frame(colData(deng_SCE)), 
       aes(x = as.factor(ouija_cell_class), 
           y = pseudotime_ouija, colour = cell_type2)) +
  geom_boxplot() + 
  coord_flip() +
  scale_color_tableau() + theme_classic() +
  xlab("Ouija cell classification") +
  ylab("Ouija pseudotime") +
  theme_classic()


# A common analysis is to work out the regulation orderings of genes. For example, is gene A 
# upregulated before gene B? Does gene C peak before the downregulation of gene D? Ouija answers these
# questions in terms of a Bayesian hypothesis test of whether the difference in regulation timing (either
# switch time or peak time) is significantly different to 0. This is collated using the gene_regulation function.

# ouija-regulation}
gene_regs <- gene_regulation(oui_vb)
head(gene_regs)



# What conclusions can you draw from the gene regulation output from Ouija?
  
# If you have time, you might try the HMC inference method and see if that changes the Ouija results
# in any way.


### Comparison of the methods

# How do the trajectories inferred by TSCAN, Monocle, Diffusion Map, SLICER and Ouija compare?
  
# TSCAN and Diffusion Map methods get the trajectory the "wrong way round", so we'll adjust that for
# these comparisons.

# compare-results, fig.width=10}
df_pseudotime <- as.data.frame(
    colData(deng_SCE)[, grep("pseudotime", colnames(colData(deng_SCE)))]
)
colnames(df_pseudotime) <- gsub("pseudotime_", "", 
                                colnames(df_pseudotime))
df_pseudotime$PC1 <- deng_SCE$PC1
df_pseudotime$order_tscan <- -df_pseudotime$order_tscan
df_pseudotime$diffusionmap <- -df_pseudotime$diffusionmap

corrplot.mixed(cor(df_pseudotime, use = "na.or.complete"), 
               order = "hclust", tl.col = "black",
               main = "Correlation matrix for pseudotime results",
               mar = c(0, 0, 3.1, 0))



# We see here that Ouija, TSCAN and SLICER all give trajectories that are similar and strongly
# correlated with PC1. Diffusion Map is less strongly correlated with these methods, and Monocle
# gives very different results.


### Expression of genes through time

# Each package also enables the visualization of expression through pseudotime. Following individual genes
# is very helpful for identifying genes that play an important role in the differentiation process. We
# illustrate the procedure using the `Rhoa` gene.

# We have added the pseudotime values computed with all methods here to
# the `colData` slot of an `SCE` object. Having done that, the full
# plotting capabilities of the `scater` package can be used to
# investigate relationships between gene expression, cell populations
# and pseudotime. This is particularly useful for the packages such as
# SLICER that do not provide plotting functions.


# __Principal components__
plotExpression(deng_SCE, "Rhoa", x = "PC1", 
               colour_by = "cell_type2", show_violin = FALSE,
               show_smooth = TRUE)


# __TSCAN__
plotExpression(deng_SCE, "Rhoa", x = "pseudotime_order_tscan", 
               colour_by = "cell_type2", show_violin = FALSE,
               show_smooth = TRUE)

#__Monocle__
plotExpression(deng_SCE, "Rhoa", x = "pseudotime_monocle", 
               colour_by = "cell_type2", show_violin = FALSE,
               show_smooth = TRUE)


#__Diffusion Map__
plotExpression(deng_SCE, "Rhoa", x = "pseudotime_diffusionmap", 
               colour_by = "cell_type2", show_violin = FALSE,
               show_smooth = TRUE)



#__SLICER__
plotExpression(deng_SCE, "Rhoa", x = "pseudotime_slicer", 
               colour_by = "cell_type2", show_violin = FALSE,
               show_smooth = TRUE)


#__Ouija__
plotExpression(deng_SCE, "Rhoa", x = "pseudotime_ouija", 
               colour_by = "cell_type2", show_violin = FALSE,
               show_smooth = TRUE)


# How many of these methods outperform the naive approach of using the first principal component to
# crepresent pseudotime for these data?

# We can Repeat the exercise using a subset of the genes, e.g. the set of highly variable genes that can be obtained using `Brennecke_getVariableGenes()`

  ## Imputation
install.packages("devtools")
library(devtools)
install_github("Vivianstats/scImpute")
library(scImpute)
library(SC3)
library(scater)
library(SingleCellExperiment)
library(mclust)
install.packages("DrImpute")
library(DrImpute)
install.packages("Rmagic")
library(Rmagic)



# As discussed previously, one of the main challenges when analyzing scRNA-seq 
# data is the presence of zeros, or dropouts. The dropouts are assumed to have 
# arisen for three possible reasons:
#   
# * The gene was not expressed in the cell and hence there are no transcripts to sequence
# * The gene was expressed, but for some reason the transcripts were lost somewhere prior to sequencing
# * The gene was expressed and transcripts were captured and turned into cDNA, but the sequencing depth was
# not sufficient to produce any reads.
# 
# Thus, dropouts could be result of experimental shortcomings, and if this is 
# the case then we would like to provide computational corrections. One possible
# solution is to impute the dropouts in the expression matrix. To be able to 
# impute gene expression values, one must have an underlying model. However, 
# since we do not know which dropout events are technical artefacts and which 
# correspond to the transcript being truly absent, imputation is a difficult 
# challenge and prone to creating [false-positive results in downstream analysis]
# (https://f1000research.com/articles/7-1740/v2).
# 
# There are many different imputation methods available we will consider three 
# fast, published methods:
#   [MAGIC](https://github.com/pkathail/magic) [@Van_Dijk2017-bh], 
# [DrImpute](https://github.com/gongx030/DrImpute) and 
# [scImpute](https://github.com/Vivianstats/scImpute) [@Li2017-tz]. 
# 
# DrImpute and scImpute both use a model to determine which zeros are technical
# and impute only those values. Both use clustering to identify a group of cells
# that are assumed to have homogenous expression. DrImpute imputes all values that
# are not consistently zero in all cells of a cluster. Whereas, scImpute uses a
# zero-inflated normal distribution fit to log-normalized expression values and
# imputed all inflated zeros. 

### scImpute

# To test `scImpute`, we use the default parameters and we apply it to the Deng dataset that we have worked
# with before. scImpute takes a .csv or .txt file as an input:
  
deng <- readRDS("data/deng/deng-reads.rds")
write.csv(counts(deng), "deng.csv")
scimpute(
  count_path = "deng.csv",
  infile = "csv",
  outfile = "txt", 
  out_dir = "./",
  Kcluster = 10,
  ncores = 2
)


#Now we can compare the results with original data by considering a PCA plot

res <- read.table("scimpute_count.txt")
colnames(res) <- NULL
res <- SingleCellExperiment(
  assays = list(logcounts = log2(as.matrix(res) + 1)), 
  colData = colData(deng)
)
rowData(res)$feature_symbol <- rowData(deng)$feature_symbol
plotPCA(
  res, 
  colour_by = "cell_type2"
)


# Compare this result to the original data in Chapter \@ref(clust-methods). What are the most
# significant differences?
  
# We can examine the expression of specific genes to directly see the effect of
# imputation on the expression distribution.
plotExpression(res, c("Sox2", "Eomes", "Zscan4d", "Fgf4"))
plotExpression(deng, c("Sox2", "Eomes", "Zscan4d", "Fgf4"))



# To evaluate the impact of the imputation, we use `SC3` to cluster the imputed matrix
res <- sc3_estimate_k(res)
metadata(res)$sc3$k_estimation
res <- sc3(res, ks = 10, n_cores = 1, gene_filter = FALSE)
adjustedRandIndex(colData(deng)$cell_type2, colData(res)$sc3_10_clusters)
plotPCA(
  res, 
  colour_by = "sc3_10_clusters"
)


# __Exercise:__ Based on the PCA and the clustering results, do you think that imputation using `scImpute`
# is a good idea for the Deng dataset?
  
  ### DrImpute
  
# We can do the same for DrImpute. DrImpute runs on a log-normalized expression matrix directly in R, 
# we generate this matrix using scater, then run DrImpute. Unlike scImpute, DrImpute considers the consensus
# imputation across a range of ks using two differ correlation distances:
  
deng <- normalize(deng)
res <- DrImpute(deng@assays[["logcounts"]], ks=8:12)
colnames(res) <- colnames(deng)
rownames(res) <- rownames(deng)
res <- SingleCellExperiment(
  assays = list(logcounts = as.matrix(res)), 
  colData = colData(deng)
)
rowData(res)$feature_symbol <- rowData(deng)$feature_symbol
plotPCA(
  res, 
  colour_by = "cell_type2"
)
plotExpression(res, c("Sox2", "Eomes", "Zscan4d", "Fgf4"))

# __Exercise:__ Check the sc3 clustering of the DrImpute matrix, do you think that imputation using
# `DrImpute` is a good idea for the Deng dataset?
  
  
  ### MAGIC
  

# MAGIC is a python package but the authors have provided an R package wrapper,
# so it can be run seemlessly from R.
# 
# Unlike scImpute and DrImpute, MAGIC smoothes the entire dataset. It imputes 
# zeros but also smoothes non-zero values to further exaggerate any structure
# within the dataset. Since it is based on a diffusion process, it specifically
# enhances trajectory-like structure in a dataset, in contrast to scImpute and 
# DrImpute which assume a cluster-like structure to the underlying data.

##@@@@@@@@@@@@ BUG | DOESN'T WORK@@@@@@@@@@@@@###
# res <- magic(t(deng@assays[["logcounts"]]), genes="all_genes", k=10, t="auto")
# Yes  ##Use pip3 install magic-impute if error occurs
# 
# res <- t(as.matrix(res))
# rownames(res) <- rownames(deng)
# colnames(res) <- colnames(deng)
# res <- SingleCellExperiment(
#   assays = list(logcounts = res), 
#   colData = colData(deng)
# )
# rowData(res)$feature_symbol <- rownames(res)
# plotPCA(
#   res, 
#   colour_by = "cell_type2"
# )
##@@@@@@@@@@@@DOESN'T WORK@@@@@@@@@@@@@###

# Compare this result to the original data in Chapter \@ref(clust-methods). What are the most significant differences?
  
# To evaluate the impact of the imputation, we use `SC3` to cluster the imputed matrix
res <- sc3_estimate_k(res)
metadata(res)$sc3$k_estimation
res <- sc3(res, ks = 10, n_cores = 1, gene_filter = FALSE)
adjustedRandIndex(colData(deng)$cell_type2, colData(res)$sc3_10_clusters)
plotPCA(
  res, 
  colour_by = "sc3_10_clusters"
)

### Bulk RNA-seq

# One of the most common types of analyses when working with bulk RNA-seq
# data is to identify differentially expressed genes. By comparing the
# genes that change between two conditions, e.g. mutant and wild-type or
# stimulated and unstimulated, it is possible to characterize the
# molecular mechanisms underlying the change.
# 
# Several different methods,
# e.g. [DESeq2](https://bioconductor.org/packages/DESeq2) and
# [edgeR](https://bioconductor.org/packages/release/bioc/html/edgeR.html),
# have been developed for bulk RNA-seq. Moreover, there are also
# extensive
# [datasets](http://genomebiology.biomedcentral.com/articles/10.1186/gb-2013-14-9-r95)
# available where the RNA-seq data has been validated using
# RT-qPCR. These data can be used to benchmark DE finding algorithms and the available evidence suggests
# that the algorithms are performing quite well.

### Single cell RNA-seq
# In contrast to bulk RNA-seq, in scRNA-seq we usually do not have a defined
# set of experimental conditions. Instead, as was shown in a previous chapter
# (\@ref(clust-methods)) we can identify the cell groups by using an unsupervised
# clustering approach. Once the groups have been identified one can find differentially
# expressed genes either by comparing the differences in variance between the groups (like the
# Kruskal-Wallis test implemented in SC3), or by comparing gene expression between clusters in a pairwise
# manner. In the following chapter we will mainly consider tools developed for pairwise comparisons.

### Differences in Distribution

# Unlike bulk RNA-seq, we generally have a large number of samples (i.e. cells) for each group we are
# comparing in single-cell experiments. Thus we can take advantage of the whole distribution of expression
# values in each group to identify differences between groups rather than only comparing estimates of 
# mean-expression as is standard for bulk RNASeq.

# There are two main approaches to comparing distributions. Firstly, we can use existing statistical
# models/distributions and fit the same type of model to the expression in each group then test for
# differences in the parameters for each model, or test whether the model fits better if a particular
# paramter is allowed to be different according to group. For instance in (dealing-with-confounders) we used
# edgeR to test whether allowing mean expression to be different in different batches significantly improved 
# the fit of a negative binomial model of the data.

# Alternatively, we can use a non-parametric test which does not assume that expression values follow
# any particular distribution, e.g. the [Kolmogorov-Smirnov test (KS-test)]
# (https://en.wikipedia.org/wiki/Kolmogorov%E2%80%93Smirnov_test). Non-parametric tests generally convert
# observed expression values to ranks and test whether the distribution of ranks for one group are signficantly
# different from the distribution of ranks for the other group. However, some non-parametric methods fail
# in the presence of a large number of tied values, such as the case for dropouts (zeros) in single-cell
# RNA-seq expression data. Moreover, if the conditions for a parametric test hold, then it will typically
# be more powerful than a non-parametric test.

### Models of single-cell RNASeq data

# The most common model of RNASeq data is the negative binomial model:
hist(
  rnbinom(
    1000, 
    mu = 10, 
    size = 100), 
  col = "grey50", 
  xlab = "Read Counts", 
  main = "Negative Binomial"
)

# Mean:
#   $\mu = mu$
#   
#   Variance:
#   $\sigma^2 = mu + mu^2/size$
#   
#   It is parameterized by the mean expression (mu) and the dispersion (size), which is inversely
# related to the variance. The negative binomial model fits bulk RNA-seq data very well and it is used
# for most statistical methods designed for such data. In addition, it has been show to fit the distribution
# of molecule counts obtained from data tagged by unique molecular identifiers (UMIs) quite well
# ([Grun et al. 2014](http://www.nature.com/nmeth/journal/v11/n6/full/nmeth.2930.html),
#   [Islam et al. 2011](http://genome.cshlp.org/content/21/7/1160)).

# However, a raw negative binomial model does not fit full-length transcript data as well due
# to the high dropout rates relative to the non-zero read counts. For this type of data a variety
# of zero-inflated negative binomial models have been proposed 
# (e.g. [MAST](https://bioconductor.org/packages/release/bioc/html/MAST.html), 
#   [SCDE](https://bioconductor.org/packages/release/bioc/html/scde.html)).

# inflation-plot, fig.cap="Zero-inflated Negative Binomial distribution"}
d <- 0.5;
counts <- rnbinom(
  1000, 
  mu = 10, 
  size = 100
)
counts[runif(1000) < d] <- 0
hist(
  counts, 
  col = "grey50", 
  xlab = "Read Counts", 
  main = "Zero-inflated NB"
)
# Mean:
#   $\mu = mu \cdot (1 - d)$
#   
#   Variance:
#   $\sigma^2 = \mu \cdot (1-d) \cdot (1 + d \cdot \mu + \mu / size)$
#   
#   These models introduce a new parameter $d$, for the dropout rate, to the negative binomial model.
# As we saw in Chapter 19, the dropout rate of a gene is strongly correlated with the mean expression
# of the gene. Different zero-inflated negative binomial models use different relationships between mu
# and d and some may fit $\mu$ and $d$ to the expression of each gene independently.

# Finally, several methods use a Poisson-Beta distribution which is based on a mechanistic model of
# transcriptional bursting. There is strong experimental support for this model ([Kim and Marioni, 2013]
# (https://genomebiology.biomedcentral.com/articles/10.1186/gb-2013-14-1-r7)) and it provides a good fit
# to scRNA-seq data but it is less easy to use than the negative-binomial models and much less existing
# methods upon which to build than the negative binomial model.

# pois-beta-plot, fit.cap="Poisson-Beta distribution"}
a <- 0.1
b <- 0.1
g <- 100
lambdas <- rbeta(1000, a, b)
counts <- sapply(g*lambdas, function(l) {rpois(1, lambda = l)})
hist(
  counts, 
  col = "grey50", 
  xlab = "Read Counts", 
  main = "Poisson-Beta"
)

# Mean:
#   $\mu = g \cdot a / (a + b)$
#   
#   Variance:
#   $\sigma^2 = g^2 \cdot a \cdot b/((a + b + 1) \cdot (a + b)^2)$
#   
#   This model uses three parameters: $a$ the rate of activation of transcription; $b$ the
# rate of inhibition of transcription; and $g$ the rate of transcript production while transcription
# is active at the locus. Differential expression methods may test each of the parameters for differences
# across groups or only one (often $g$).
#
# All of these models may be further expanded to explicitly account for other sources of gene
# expression differences such as batch-effect or library depth depending on the particular DE algorithm.

### Bulk RNA-seq

# One of the most common types of analyses when working with bulk RNA-seq
# data is to identify differentially expressed genes. By comparing the
# genes that change between two conditions, e.g. mutant and wild-type or
# stimulated and unstimulated, it is possible to characterize the
# molecular mechanisms underlying the change.
# 
# Several different methods,
# e.g. [DESeq2](https://bioconductor.org/packages/DESeq2) and
# [edgeR](https://bioconductor.org/packages/release/bioc/html/edgeR.html),
# have been developed for bulk RNA-seq. Moreover, there are also
# extensive
# [datasets](http://genomebiology.biomedcentral.com/articles/10.1186/gb-2013-14-9-r95)
# available where the RNA-seq data has been validated using
# RT-qPCR. These data can be used to benchmark DE finding algorithms and the available evidence
# suggests that the algorithms are performing quite well.

### Single cell RNA-seq

# In contrast to bulk RNA-seq, in scRNA-seq we usually do not have a defined
# set of experimental conditions. Instead, as was shown in a previous chapter
# (\@ref(clust-methods)) we can identify the cell groups by using an unsupervised
# clustering approach. Once the groups have been identified one can find differentially
# expressed genes either by comparing the differences in variance between the groups 
# (like the Kruskal-Wallis test implemented in SC3), or by comparing gene expression between
# clusters in a pairwise manner. In the following chapter we will mainly consider tools developed
# for pairwise comparisons.

### Differences in Distribution

# Unlike bulk RNA-seq, we generally have a large number of samples (i.e. cells) for each group we are
# comparing in single-cell experiments. Thus we can take advantage of the whole distribution of expression
# values in each group to identify differences between groups rather than only comparing estimates of
# mean-expression as is standard for bulk RNASeq.

# There are two main approaches to comparing distributions. Firstly, we can use existing statistical
# models/distributions and fit the same type of model to the expression in each group then test for 
# differences in the parameters for each model, or test whether the model fits better if a particular
# paramter is allowed to be different according to group. For instance in (dealing-with-confounders) we used
# edgeR to test whether allowing mean expression to be different in different batches significantly improved 
# the fit of a negative binomial model of the data.

# Alternatively, we can use a non-parametric test which does not assume that expression values follow any
# particular distribution, e.g. the [Kolmogorov-Smirnov test (KS-test)]
# (https://en.wikipedia.org/wiki/Kolmogorov%E2%80%93Smirnov_test). Non-parametric tests generally convert
# observed expression values to ranks and test whether the distribution of ranks for one group are
# signficantly different from the distribution of ranks for the other group. However, some non-parametric
# methods fail in the presence of a large number of tied values, such as the case for dropouts (zeros)
# in single-cell RNA-seq expression data. Moreover, if the conditions for a parametric test hold, then
# it will typically be more powerful than a non-parametric test.

### Models of single-cell RNASeq data

# The most common model of RNASeq data is the negative binomial model:
  
  
# nb-plot, fig.cap="Negative Binomial distribution of read counts for a single gene across 1000 cells"}
hist(
  rnbinom(
    1000, 
    mu = 10, 
    size = 100), 
  col = "grey50", 
  xlab = "Read Counts", 
  main = "Negative Binomial"
)

# Mean:
#   $\mu = mu$
#   
#   Variance:
#   $\sigma^2 = mu + mu^2/size$
  
#   It is parameterized by the mean expression (mu) and the dispersion (size), which is inversely related
# to the variance. The negative binomial model fits bulk RNA-seq data very well and it is used for most
# statistical methods designed for such data. In addition, it has been show to fit the distribution of
# molecule counts obtained from data tagged by unique molecular identifiers (UMIs) quite well 
# ([Grun et al. 2014](http://www.nature.com/nmeth/journal/v11/n6/full/nmeth.2930.html), 
#   [Islam et al. 2011](http://genome.cshlp.org/content/21/7/1160)).

# However, a raw negative binomial model does not fit full-length transcript data as well due to the high
# dropout rates relative to the non-zero read counts. For this type of data a variety of zero-inflated
# negative binomial models have been proposed
# (e.g. [MAST](https://bioconductor.org/packages/release/bioc/html/MAST.html),
#   [SCDE](https://bioconductor.org/packages/release/bioc/html/scde.html)).

# zero-inflation-plot, fig.cap="Zero-inflated Negative Binomial distribution"}
d <- 0.5;
counts <- rnbinom(
  1000, 
  mu = 10, 
  size = 100
)
counts[runif(1000) < d] <- 0
hist(
  counts, 
  col = "grey50", 
  xlab = "Read Counts", 
  main = "Zero-inflated NB"
)

# Mean:
#   $\mu = mu \cdot (1 - d)$
#   
#   Variance:
#   $\sigma^2 = \mu \cdot (1-d) \cdot (1 + d \cdot \mu + \mu / size)$
#   
#   These models introduce a new parameter $d$, for the dropout rate, to the negative binomial model. As we saw
# in Chapter 19, the dropout rate of a gene is strongly correlated with the mean expression of the gene.
# Different zero-inflated negative binomial models use different relationships between mu and d and some 
# may fit $\mu$ and $d$ to the expression of each gene independently.
# 
# Finally, several methods use a Poisson-Beta distribution which is based on a mechanistic model of
# transcriptional bursting. There is strong experimental support for this model 
# ([Kim and Marioni, 2013](https://genomebiology.biomedcentral.com/articles/10.1186/gb-2013-14-1-r7)) and
# it provides a good fit to scRNA-seq data but it is less easy to use than the negative-binomial models
# and much less existing methods upon which to build than the negative binomial model.

# rpois-beta-plot, fit.cap="Poisson-Beta distribution"}
a <- 0.1
b <- 0.1
g <- 100
lambdas <- rbeta(1000, a, b)
counts <- sapply(g*lambdas, function(l) {rpois(1, lambda = l)})
hist(
  counts, 
  col = "grey50", 
  xlab = "Read Counts", 
  main = "Poisson-Beta"
)
# 
# Mean:
#   $\mu = g \cdot a / (a + b)$
#   
#   Variance:
#   $\sigma^2 = g^2 \cdot a \cdot b/((a + b + 1) \cdot (a + b)^2)$
#   
#   This model uses three parameters: $a$ the rate of activation of transcription; $b$ the rate of
# inhibition of transcription; and $g$ the rate of transcript production while transcription is active
# at the locus. Differential expression methods may test each of the parameters for differences across
# groups or only one (often $g$).
# 
# All of these models may be further expanded to explicitly account for other sources of gene expression
# differences such as batch-effect or library depth depending on the particular DE algorithm.

#### DE in a real dataset

library(scRNA.seq.funcs)
library(edgeR)
library(monocle)
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("MAST")
a
no
library(MAST)
library(ROCR)
### Introduction

# To test different single-cell differential expression methods we will be using the Blischak dataset
# For this experiment bulk RNA-seq data for each cell-line was generated in addition to single-cell data.
# We will use the differentially expressed genes identified using standard methods on the respective bulk
# data as the ground truth for evaluating the accuracy of each single-cell method. To save time we have
# pre-computed these for you. You can run the commands below to load these data.

DE <- read.table("data/tung/TPs.txt")
notDE <- read.table("data/tung/TNs.txt")
GroundTruth <- list(
  DE = as.character(unlist(DE)), 
  notDE = as.character(unlist(notDE))
)
 
# This ground truth has been produce for the comparison of individual NA19101 to NA19239.
# Now load the respective single-cell data:
molecules <- read.table("data/tung/molecules.txt", sep = "\t")
anno <- read.table("data/tung/annotation.txt", sep = "\t", header = TRUE)
keep <- anno[,1] == "NA19101" | anno[,1] == "NA19239"
data <- molecules[,keep]
group <- anno[keep,1]
batch <- anno[keep,4]
# remove genes that aren't expressed in at least 6 cells
gkeep <- rowSums(data > 0) > 5;
counts <- data[gkeep,]
# Library size normalization
lib_size = colSums(counts)
norm <- t(t(counts)/lib_size * median(lib_size)) 
# Variant of CPM for datasets with library sizes of fewer than 1 mil molecules
# Now we will compare various single-cell DE methods. Note that we will only be running methods
# which are available as R-packages and run relatively quickly.

### Kolmogorov-Smirnov test

# The types of test that are easiest to work with are non-parametric
# ones. The most commonly used non-parametric test is the
# [Kolmogorov-Smirnov test](https://en.wikipedia.org/wiki/Kolmogorov%E2%80%93Smirnov_test) (KS-test) 
# and we can use it to compare the distributions for each gene in the two individuals.
# 
# The KS-test quantifies the distance between the empirical cummulative distributions of the expression
# of each gene in each of the two populations. It is sensitive to changes in mean experession and changes 
# in variability. However it assumes data is continuous and may perform poorly when data contains a large
# number of identical values (eg. zeros). Another issue with the KS-test is that it can be very sensitive
# for large sample sizes and thus it may end up as significant even though the magnitude of the difference
# is very small.

# "Illustration of the two-sample Kolmogorov–Smirnov statistic. Red and blue lines each correspond to an
# empirical distribution function, and the black arrow is the two-sample KS statistic.
# (taken from [here](https://en.wikipedia.org/wiki/Kolmogorov%E2%80%93Smirnov_test))"}
# Now run the test:
  
pVals <- apply(
  norm, 1, function(x) {
    ks.test(
      x[group == "NA19101"], 
      x[group == "NA19239"]
    )$p.value
  }
)
# multiple testing correction
pVals <- p.adjust(pVals, method = "fdr")

# This code "applies" the function to each row (specified by 1) of the expression matrix,
# data. In the function we are returning just the p.value from the ks.test output. We can
# now consider how many of the ground truth positive and negative DE genes are detected by the KS-test:
  
  #### Evaluating Accuracy
  
sigDE <- names(pVals)[pVals < 0.05]
length(sigDE) 
# Number of KS-DE genes
sum(GroundTruth$DE %in% sigDE) 
# Number of KS-DE genes that are true DE genes
sum(GroundTruth$notDE %in% sigDE)
# Number of KS-DE genes that are truly not-DE

# As you can see many more of our ground truth negative genes were identified as DE by the KS-test
# (false positives) than ground truth positive genes (true positives), however this may be due to
# the larger number of notDE genes thus we typically normalize these counts as the True positive rate
# (TPR), TP/(TP + FN), and False positive rate (FPR), FP/(FP+TP).
tp <- sum(GroundTruth$DE %in% sigDE)
fp <- sum(GroundTruth$notDE %in% sigDE)
tn <- sum(GroundTruth$notDE %in% names(pVals)[pVals >= 0.05])
fn <- sum(GroundTruth$DE %in% names(pVals)[pVals >= 0.05])
tpr <- tp/(tp + fn)
fpr <- fp/(fp + tn)
cat(c(tpr, fpr))

# Now we can see the TPR is much higher than the FPR indicating the KS test is identifying DE genes.

# So far we've only evaluated the performance at a single significance threshold. Often it is
# informative to vary the threshold and evaluate performance across a range of values. This is then
# plotted as a receiver-operating-characteristic curve (ROC) and a general accuracy statistic can be
# calculated as the area under this curve (AUC). We will use the ROCR package to facilitate this plotting.

# ROC curve for KS-test."}
# Only consider genes for which we know the ground truth
pVals <- pVals[names(pVals) %in% GroundTruth$DE | 
               names(pVals) %in% GroundTruth$notDE] 
truth <- rep(1, times = length(pVals));
truth[names(pVals) %in% GroundTruth$DE] = 0;
pred <- ROCR::prediction(pVals, truth)
perf <- ROCR::performance(pred, "tpr", "fpr")
ROCR::plot(perf)
aucObj <- ROCR::performance(pred, "auc")
aucObj@y.values[[1]] # AUC

# Finally to facilitate the comparisons of other DE methods let's put this code into a function
# so we don't need to repeat it:

DE_Quality_AUC <- function(pVals) {
    pVals <- pVals[names(pVals) %in% GroundTruth$DE | 
                   names(pVals) %in% GroundTruth$notDE]
    truth <- rep(1, times = length(pVals));
    truth[names(pVals) %in% GroundTruth$DE] = 0;
    pred <- ROCR::prediction(pVals, truth)
    perf <- ROCR::performance(pred, "tpr", "fpr")
    ROCR::plot(perf)
    aucObj <- ROCR::performance(pred, "auc")
    return(aucObj@y.values[[1]])
}

### Wilcox/Mann-Whitney-U Test
# 
# The Wilcox-rank-sum test is another non-parametric test, but tests specifically if values in
# one group are greater/less than the values in the other group. Thus it is often considered a
# test for difference in median expression between two groups; whereas the KS-test is sensitive 
# to any change in distribution of expression values.

# ROC curve for Wilcox test
pVals <- apply(
    norm, 1, function(x) {
        wilcox.test(
            x[group == "NA19101"], 
            x[group == "NA19239"]
        )$p.value
    }
)
# multiple testing correction
pVals <- p.adjust(pVals, method = "fdr")

DE_Quality_AUC <- function(pVals) {
  pVals <- pVals[names(pVals) %in% GroundTruth$DE | 
                   names(pVals) %in% GroundTruth$notDE]
  truth <- rep(1, times = length(pVals));
  truth[names(pVals) %in% GroundTruth$DE] = 0;
  pred <- ROCR::prediction(pVals, truth)
  perf <- ROCR::performance(pred, "tpr", "fpr")
  ROCR::plot(perf)
  aucObj <- ROCR::performance(pred, "auc")
  return(aucObj@y.values[[1]])
}
DE_Quality_AUC(pVals)


### edgeR

# We've already used edgeR for differential expression in Chapter dealing-with-confounders). edgeR is based
# on a negative binomial model of gene expression and uses a generalized linear model (GLM) framework,
# the enables us to include other factors such as batch to the model.

# ROC curve for edgeR.", message=FALSE}
dge <- DGEList(
  counts = counts, 
  norm.factors = rep(1, length(counts[1,])), 
  group = group
)
group_edgeR <- factor(group)
design <- model.matrix(~ group_edgeR)
dge <- estimateDisp(dge, design = design, trend.method = "none")
fit <- glmFit(dge, design)
res <- glmLRT(fit)
pVals <- res$table[,4]
names(pVals) <- rownames(res$table)

pVals <- p.adjust(pVals, method = "fdr")

DE_Quality_AUC(pVals)

### Monocle

# [Monocle](https://bioconductor.org/packages/release/bioc/html/monocle.html) can use several different models
# for DE. For count data it recommends the Negative Binomial model (negbinomial.size). For normalized data
# it recommends log-transforming it then using a normal distribution (gaussianff). Similar to edgeR this
# method uses a GLM framework so in theory can account for batches, however in practice the model 
# fails for this dataset if batches are included.

# ROC curve for Monocle.", message=FALSE, warning=FALSE}
pd <- data.frame(group = group, batch = batch)
rownames(pd) <- colnames(counts)
pd <- new("AnnotatedDataFrame", data = pd)

Obj <- newCellDataSet(
  as.matrix(counts), 
  phenoData = pd, 
  expressionFamily = negbinomial.size()
)
Obj <- estimateSizeFactors(Obj)
Obj <- estimateDispersions(Obj)

## Intensive command
res <- differentialGeneTest(Obj, fullModelFormulaStr = "~group")

pVals <- res[,3]
names(pVals) <- rownames(res)
pVals <- p.adjust(pVals, method = "fdr")
DE_Quality_AUC(pVals)

# Comparing the results using the negative binomial model on counts and those from using the 
# normal/gaussian model (`gaussianff()`) on log-transformed normalized counts.

# ROC curve for Monocle-gaussian.", message=FALSE, echo=FALSE, warning=FALSE}
pd <- data.frame(group = group, batch = batch)
rownames(pd) <- colnames(norm)
pd <- new("AnnotatedDataFrame", data = pd)

Obj_log <- newCellDataSet(
  as.matrix(log(norm + 1) / log(2)), 
  phenoData = pd, 
  expressionFamily = gaussianff()
)
Obj_log <- estimateSizeFactors(Obj_log)
# Obj_log <- estimateDispersions(Obj_log)

## Intensive Command 
res <- differentialGeneTest(Obj_log, fullModelFormulaStr = "~group")

pVals <- res[,3]
names(pVals) <- rownames(res)
pVals <- p.adjust(pVals, method = "fdr")
DE_Quality_AUC(pVals)

### MAST

# [MAST](https://bioconductor.org/packages/release/bioc/html/MAST.html) is based on a zero-inflated
# negative binomial model. It tests for differential expression using a hurdle model to combine tests
# of discrete (0 vs not zero) and continuous (non-zero values) aspects of gene expression. Again this
# uses a linear modelling framework to enable complex models to be considered.

# ROC curve for MAST.", message=FALSE}
log_counts <- log(counts + 1) / log(2)
fData <- data.frame(names = rownames(log_counts))
rownames(fData) <- rownames(log_counts);
cData <- data.frame(cond = group)
rownames(cData) <- colnames(log_counts)

obj <- FromMatrix(as.matrix(log_counts), cData, fData)
colData(obj)$cngeneson <- scale(colSums(assay(obj) > 0))
cond <- factor(colData(obj)$cond)

# Model expression as function of condition & number of detected genes
zlmCond <- zlm.SingleCellAssay(~ cond + cngeneson, obj) 

summaryCond <- summary(zlmCond, doLRT = "condNA19101")
summaryDt <- summaryCond$datatable

summaryDt <- as.data.frame(summaryDt)
pVals <- unlist(summaryDt[summaryDt$component == "H",4]) # H = hurdle model
names(pVals) <- unlist(summaryDt[summaryDt$component == "H",1])
pVals <- p.adjust(pVals, method = "fdr")
DE_Quality_AUC(pVals)

### Slow Methods (>1h to run) 

# These methods are too slow to run today but we encourage you to try them out on your own:
  
  ### BPSC
  
# [BPSC](https://academic.oup.com/bioinformatics/article/32/14/2128/2288270/Beta-Poisson-model-for-single-cell-RNA-seq-data) uses the Poisson-Beta model
# of single-cell gene expression, which we discussed in the previous chapter, and combines it with generalized
# linear models which we've already encountered when using edgeR. BPSC performs comparisons of one or more
# groups to a reference group ("control") and can include other factors such as batches in the model.

#####@@@@@@@@@@@@@@@@@@Commenting them out so that I don't accidently run them @@@@@@@@@@@@@@@@@@####
install.packages("devtools")
library("devtools")
install_github("nghiavtr/BPSC")
library("BPSC")
library(BPSC)
# bpsc_data <- norm[,batch=="NA19101.r1" | batch=="NA19239.r1"]
# bpsc_group = group[batch=="NA19101.r1" | batch=="NA19239.r1"]
# 
# control_cells <- which(bpsc_group == "NA19101")
# design <- model.matrix(~bpsc_group)
# coef=2 # group label
# res=BPglm(data=bpsc_data, controlIds=control_cells, design=design, coef=coef, 
#                 estIntPar=FALSE, useParallel = FALSE)
# pVals = res$PVAL
# pVals <- p.adjust(pVals, method = "fdr")
# DE_Quality_AUC(pVals)


### SCDE
# [SCDE](http://hms-dbmi.github.io/scde/) is the first single-cell specific DE method. It
# fits a zero-inflated negative binomial model to expression data using Bayesian statistics.
# The usage below tests for differences in mean expression of individual genes across groups
# but recent versions include methods to test for differences in mean expression or dispersion
# of groups of genes, usually representing a pathway.
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("scde")
a
no
library(scde)

cnts <- apply(
    counts,
    2,
    function(x) {
        storage.mode(x) <- 'integer'
        return(x)
    }
)
names(group) <- 1:length(group)
colnames(cnts) <- 1:length(group)
########@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# Not working -> Try this for a pluasible soln. https://github.com/hms-dbmi/scde/issues/48
########@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# o.ifm <- scde::scde.error.models(
#     counts = cnts,
#     groups = group,
#     n.cores = 1,
#     threshold.segmentation = TRUE,
#     save.crossfit.plots = FALSE,
#     save.model.plots = FALSE,
#     verbose = 0,
#     min.size.entries = 2
# )
# priors <- scde::scde.expression.prior(
#     models = o.ifm,
#     counts = cnts,
#     length.out = 400,
#     show.plot = FALSE
# )
# resSCDE <- scde::scde.expression.difference(
#     o.ifm,
#     cnts,
#     priors,
#     groups = group,
#     n.randomizations = 100,
#     n.cores = 1,
#     verbose = 0
# )
# # Convert Z-scores into 2-tailed p-values
# pVals <- pnorm(abs(resSCDE$cZ), lower.tail = FALSE) * 2
# DE_Quality_AUC(pVals)



  






