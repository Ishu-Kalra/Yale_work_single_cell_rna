install.packages("Matrix")
library(Matrix)

PATH = "hg19"

#Just the list of names
molecules_list <- list.files(path=PATH, pattern="*.mtx", full.names=TRUE, recursive=TRUE)
cellbarcodes_list <- list.files(path=PATH, pattern="barcodes.tsv", full.names=TRUE, recursive=TRUE)
genenames_list <- list.files(path=PATH, pattern="genes.tsv", full.names=TRUE, recursive=TRUE)
filenames_list <- list.files(path=PATH, pattern="*", full.names=FALSE, recursive=FALSE)

#Sorting for the downstram processing to be easier
# molecules_list <- sort(molecules_list)
# genenames_list <- sort(genenames_list)
# cellbarcodes_list <- sort(cellbarcodes_list)
# filenames_list <- sort(filenames_list)
n = length(molecules_list)

channel <- vector("list", n)
#length of each list should be same. If not then files are missing

#empty tables being initialized for each of the files
cellbarcodes_tables <- vector("list", n)
genenames_tables <- vector("list", n)
molecules_tables <- vector("list", n)

##10X data has no spike ins so that part is skipped

library(stringr)

i = 1

while(i <= n) {
  molecules_tables[[i]] <-readMM(molecules_list[i])
  genenames_tables[[i]] <-read.table(genenames_list[i])
  cellbarcodes_tables[[i]] <-read.table(cellbarcodes_list[i])
  
  molecule <- molecules_tables[i]
  genename <- genenames_tables[[i]][,1]
  file <- filenames_list[[i]] #eg. Kidney-10X_P4_5
  cellbarcode <- cellbarcodes_tables[[i]][,1]
  
  #file <- (str_match(string = file, pattern = '-(.*)')[2]) #eg 10X_P4_5 #In 10x genomics data, it will give NA so commenting it out
  channel[[i]] <- file
  print(file)
  
  rownames(molecules_tables[[i]])<-genename
  colnames(molecules_tables[[i]])<-paste(file, cellbarcode, sep = "_") #1<file>_AAAGTAGAGATGCCAG-1
  
  i <- i + 1
}

col_no = ncol(molecules_tables[[1]])
anns <- data.frame(matrix(nrow = col_no, ncol = col_no))
rownames(anns) <- colnames(molecules_tables[[1]])

sceset <- SingleCellExperiment(
  assays = list(counts = as.matrix(molecules_tables[[1]])),
  colData = anns
)

name <- paste(file, ".rds", sep="")

saveRDS(sceset, name)




