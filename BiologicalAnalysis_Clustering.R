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

deng <- readRDS("data/deng/deng-reads.rds")
deng

table(colData(deng)$cell_type2)
plotPCA(deng, colour_by = "cell_type2")

#Estimating the number of SC3 clusters
deng <- sc3_estimate_k(deng)
metadata(deng)$sc3$k_estimation
plotPCA(deng, colour_by = "cell_type1")

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
