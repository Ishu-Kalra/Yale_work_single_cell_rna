library("Matrix")
library("SingleCellExperiment")

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
rownames(cell_anns) <- colnames(molecules)


molecules1 <- molecules
cell_anns1 <- cell_anns

cellbarcodes <- read.table("droplet/droplet/Kidney-10X_P4_6/barcodes.tsv")
genenames <- read.table("droplet/droplet/Kidney-10X_P4_6/genes.tsv")
molecules <- Matrix::readMM("droplet/droplet/Kidney-10X_P4_6/matrix.mtx")
rownames(molecules) <- genenames[,1]
colnames(molecules) <- paste("10X_P4_6", cellbarcodes[,1], sep="_")
mouseID <- "3_9_M"
ann_subset <- ann[match(colnames(molecules), ann[,1]),]
celltype <- ann_subset[,3]
cell_anns <- data.frame(mouse = rep(mouseID, times=ncol(molecules)), type=celltype)
rownames(cell_anns) <- colnames(molecules)

molecules2 <- molecules
cell_anns2 <- cell_anns

cellbarcodes <- read.table("droplet/droplet/Kidney-10X_P7_5/barcodes.tsv")
genenames <- read.table("droplet/droplet/Kidney-10X_P7_5/genes.tsv")
molecules <- Matrix::readMM("droplet/droplet/Kidney-10X_P7_5/matrix.mtx")
rownames(molecules) <- genenames[,1]
colnames(molecules) <- paste("10X_P7_5", cellbarcodes[,1], sep="_")
mouseID <- "3_57_F"
ann_subset <- ann[match(colnames(molecules), ann[,1]),]
celltype <- ann_subset[,3]
cell_anns <- data.frame(mouse = rep(mouseID, times=ncol(molecules)), type=celltype)
rownames(cell_anns) <- colnames(molecules)

molecules3 <- molecules
cell_anns3 <- cell_anns

all_cell_anns <- as.data.frame(rbind(cell_anns1, cell_anns2, cell_anns3))
all_cell_anns$batch <- rep(c("10X_P4_5", "10X_P4_6","10X_P7_5"), times = c(nrow(cell_anns1), nrow(cell_anns2), nrow(cell_anns3)))

identical(rownames(molecules1), rownames(molecules2))
identical(rownames(molecules1), rownames(molecules3))

sum(colnames(molecules1) %in% colnames(molecules2))
sum(colnames(molecules1) %in% colnames(molecules3))
sum(colnames(molecules2) %in% colnames(molecules3))

all_molecules <- cbind(molecules1, molecules2, molecules3)
all_cell_anns <- as.data.frame(rbind(cell_anns1, cell_anns2, cell_anns3))
all_cell_anns$batch <- rep(c("10X_P4_5", "10X_P4_6","10X_P7_5"), times = c(nrow(cell_anns1), nrow(cell_anns2), nrow(cell_anns3)))

dim(all_molecules)[2]

all_molecules <- as.matrix(all_molecules)
sceset <- SingleCellExperiment(
  assays = list(counts = as.matrix(all_molecules)),
  colData = all_cell_anns
)

saveRDS(sceset, "kidney_droplet.rds")


