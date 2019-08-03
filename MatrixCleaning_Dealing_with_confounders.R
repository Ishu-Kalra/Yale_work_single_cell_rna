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
