dat = read.delim("FACS/FACS/Kidney-counts.csv", sep = ",", header = TRUE)
dat[1:5, 1:5]
dim(dat)

#first columns is gene names. We move them to become row names
rownames(dat) <- dat[, 1]
dat <- dat[,-1]

dat[1:5, 1:5]

#Removing spike-ins of ERCC
rownames(dat)[grep("^ERCC-", rownames(dat))]

cellIDs <- colnames(dat)
cell_info <- strsplit(cellIDs, "\\.")
Well <- lapply(cell_info, function(x){x[1]})
Well <- unlist(Well)
Plate <- unlist(lapply(cell_info, function(x){x[2]}))
Mouse <- unlist(lapply(cell_info, function(x){x[3]}))

summary(factor(Mouse))

table(Mouse, Plate)

ann <- read.table("FACS/FACS_annotations.csv", sep=",", header=TRUE)
ann <- ann[match(cellIDs, ann[,1]),]
celltype <- ann[,3]
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("scater")
a
no

library("SingleCellExperiment")
library("scater")
cell_anns <- data.frame(mouse = Mouse, well=Well, type=celltype)
rownames(cell_anns) <- colnames(dat)
sceset <- SingleCellExperiment(assays = list(counts = as.matrix(dat)), colData=cell_anns)
isSpike(sceset, "ERCC") <- grepl("ERCC-", rownames(sceset))

###########################################################################################################################################
###### 10x Genomics Data 
###########################################################################################################################################
library("Matrix")
cellbarcodes <- read.table("droplet/droplet/Kidney-10X_P4_5/barcodes.tsv")
genenames <- read.table("droplet/droplet/Kidney-10X_P4_5/genes.tsv")
molecules <- readMM("droplet/droplet/Kidney-10X_P4_5/matrix.mtx")

rownames(molecules) <- genenames[,1]
colnames(molecules) <- paste("10X_P4_5", cellbarcodes[,1], sep="_")

meta <- read.delim("droplet/droplet_metadata.csv", sep=",", header = TRUE)
head(meta)
meta[meta$channel == "10X_P4_5",]
mouseID <- "3_8_M"

ann <- read.delim("droplet/droplet_annotation.csv", sep=",", header=TRUE)
head(ann)

ann[,1] <- paste(ann[,1], "-1", sep="")
ann_subset <- ann[match(colnames(molecules), ann[,1]),]
celltype <- ann_subset[,3]

cell_anns <- data.frame(mouse = rep(mouseID, times=ncol(molecules)), type=celltype)
rownames(cell_anns) <- colnames(molecules);

all_cell_anns <- as.data.frame(rbind(cell_anns))
all_cell_anns$batch <- rep(c("10X_P4_5"), times = c(nrow(cell_anns)))

all_molecules <- as.matrix(molecules)
sceset <- SingleCellExperiment(
  assays = list(counts = as.matrix(all_molecules)),
  colData = all_cell_anns
)

saveRDS(sceset, "kidney_droplet.rds")
