
# Seurat {#seurat-chapter}

# [Seurat](http://satijalab.org/seurat/) was originally developed as a clustering tool for scRNA-seq data, however in the last
# few years the focus of the package has become less specific and at the moment `Seurat` is a popular R package that can perform 
# QC, analysis, and exploration of scRNA-seq data, i.e. many of the tasks covered in this course.

# __Note__ We recommend using `Seurat` for datasets with more than $5000$ cells. For smaller dataset a good alternative will be `SC3`.

# __Note__ In this chapter we use an exact copy of [this tutorial](https://satijalab.org/seurat/pbmc3k_tutorial.html).

## Setup the Seurat Object

# We will be analyzing the a dataset of Peripheral Blood Mononuclear Cells (PBMC) freely available from 10X Genomics. There are
# 2,700 single cells that were sequenced on the Illumina NextSeq 500. The raw data can be found [here]
# (https://s3-us-west-2.amazonaws.com/10x.files/samples/cell/pbmc3k/pbmc3k_filtered_gene_bc_matrices.tar.gz).

# We start by reading in the data. All features in Seurat have been configured to work with sparse matrices which results in
# significant memory and speed savings for Drop-seq/inDrop/10x data.

library(Seurat)
library(dplyr)
library(cowplot)

# Load the PBMC dataset
pbmc.data <- Read10X(data.dir = "data/pbmc3k_filtered_gene_bc_matrices/hg19/")

# Examine the memory savings between regular and sparse matrices
dense.size <- object.size(x = as.matrix(x = pbmc.data))
dense.size

sparse.size <- object.size(x = pbmc.data)
sparse.size

dense.size/sparse.size

# Initialize the Seurat object with the raw (non-normalized data).  Keep all
# genes expressed in >= 3 cells (~0.1% of the data). Keep all cells with at
# least 200 detected genes
pbmc <- CreateSeuratObject(counts = pbmc.data, min.cells = 3, min.features  = 200, project = "10X_PBMC", assay = "RNA")

## Standard pre-processing workflow

# The steps below encompass the standard pre-processing workflow for scRNA-seq data in Seurat. These represent the creation
# of a Seurat object, the selection and filtration of cells based on QC metrics, data normalization and scaling, and the detection
# of highly variable genes. 