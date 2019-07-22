install.packages("Matrix")
library(Matrix)

PATH = "droplet/droplet"

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
  file <- filenames_list[i] #eg. Kidney-10X_P4_5
  cellbarcode <- cellbarcodes_tables[[i]][,1]
  
  file <- (str_match(string = file, pattern = '-(.*)')[2]) #eg 10X_P4_5
  channel[[i]] <- file
  print(file)
  
  rownames(molecules_tables[[i]])<-genename
  colnames(molecules_tables[[i]])<-paste(file, cellbarcode, sep = "_") #10X_P4_5_AAAGTAGAGATGCCAG-1
  
  i <- i + 1
}


meta <- read.delim("/Users/ishukalra/YaleWork/Pipeline_crude/Droplet/droplet_metadata.csv", sep=",", header = TRUE)
# meta <- meta[, order(3)]
head(meta)

# meta[meta$channel == "10X_P4_5",]
# cellbarcodes <- read.table("droplet/droplet/Kidney-10X_P4_5/barcodes.tsv")
# genenames <- read.table("droplet/droplet/Kidney-10X_P4_5/genes.tsv")
# molecules <- readMM("droplet/droplet/Kidney-10X_P4_5/matrix.mtx")
# rownames(molecules) <- genenames[,1]
# colnames(molecules) <- paste("10X_P4_5", cellbarcodes[,1], sep="_")
# meta[meta$channel == "10X_P4_5",]
# mouseID<-"3_8_M"

# ann <- read.delim("droplet/droplet_annotation.csv", sep=",", header=TRUE)
# head(ann)
# ann[,1] <- paste(ann[,1], "-1", sep="")
# ann_subset <- ann[match(colnames(molecules), ann[,1]),]
# celltype <- ann_subset[,3]
# print(celltype)
# print(ncol(molecules))
# cell_anns <- data.frame(mouse = rep(mouseID, times=ncol(molecules)), type=celltype)
# rownames(cell_anns) <- colnames(molecules);



ann_tables <- vector("list", n)
celltype_tables <- vector("list", n)
cell_anns_tables <- vector("list", n)

ann <- read.delim("/Users/ishukalra/YaleWork/Pipeline_crude/Droplet/droplet_annotation.csv", sep=",", header=TRUE)
head(ann)
# tissue_tables <- vector("list", n)
ann[,1] <- paste(ann[,1], "-1", sep="") #Making it consistent. Can change from one 10x dataset to another 

i <- 1
while (i <= n) {
  if (ncol(molecules_tables[[i]]) < 10000) {  ##So that ma does not get overloaded. Can remove it once we need to run the script in the FARNAM
    ann_tables[[i]] <- ann[match(colnames(molecules_tables[[i]]), ann[,1]), ]  ##Annotation of i-th element
    celltype_tables[[i]] <- ann_tables[[i]][,3]  ##Cell type of i-th element
    mouseID <- meta[meta$channel == channel[[i]],][2]   ##Mouse ID of i-th element
    # tissue <- meta[meta$channel == channel[[i]],][3]    ##Tissue type of i-th element
    #print(tissue[[1]][[1]])
    #print(paste(i, tissue, sep = "    ")) #Getting the tissue
    # tissue_tables[[i]] <- tissue 
    cell_anns_tables[[i]] <- data.frame(mouse = rep(mouseID, times = ncol(molecules_tables[[i]])), type = celltype_tables[[i]])
    rownames(cell_anns_tables[[i]]) <- colnames(molecules_tables[[i]]);
  }
  i <- i + 1
}

tissue <- meta[, 3]
  # if (ncol(molecules_tables[[i]]) < 10000) {
  # ann_tables[[i]] <- ann[match(colnames(molecules_tables[[i]]), ann[,1]), ]  #If col of molecules is annotated
  # celltype_tables[[i]] <- ann_tables[[i]][,3]
  # mouseID <- meta[meta$channel == channel[[i]],][2] #Comparing channel and getting that particular mouse ID
  # print(paste("Instance", i, sep = " "))
  # print(mouseID)
  # print(ncol(molecules_tables[[i]]))
  # #print(ncol(celltype_tables[[i]]))
  # if (ncol(celltype))
  # cell_anns_tables[[i]] <- data.frame(mouse = rep(mouseID, times = ncol(molecules_tables[[i]])), type = celltype_tables[[i]])
  # #rownames(cell_anns_tables[[i]]) <- colnames(molecules_tables[[i]]);
  # }
  # 
  # i <- i + 1
initial_row = molecules_tables[[1]]
##Checking if identical rows or not
correct_rows <- function(molecules_list) {
  init <- molecules_list[[1]]
  n <- length(molecules_list)
  i <- 1
  flag <- TRUE
  while (i <= n) {
    if (!identical(rownames(init), rownames(molecules_list[[i]]))) {
      flag <- FALSE
      break
    }
    i <- i + 1
  }
  return (flag)
}

##Checking for col
correct_cols <- function(molecules_list) {
  init <- molecules_list[[1]]
  n <- length(molecules_list)
  i <- 2
  flag <- TRUE
  while (i <= n) {
    if ((sum(colnames(init) %in% colnames(molecules_list[[i]]))) != 0) {
      flag <- FALSE
      break
    }
    i <- i + 1
  }
  return (flag)
}

install.packages("dplyr")
library("dplyr")

i <- 1
j <- 1
processed_tissue <- vector("list")
while (i <= n) {
  initial_tissue <- NULL
  molecules_tissue <- vector("list")
  cell_anns_tissue <- vector("list")
  rows_tissue <- vector()
  channels_tissue <- vector("list")
  print(tissue[[i]])
  #print(!is.null(tissue_tables[[i]][[1]][[1]]))
  if (!is.null(tissue[[i]])) {
    if (!(tissue[[i]] %in% processed_tissue)) {
      initial_tissue <- tissue[[i]]
      j <- i
      k <- 1
      while (j <= n) {
        if (!is.null(cell_anns_tables[[i]])) {
          if (initial_tissue == tissue[[j]]) {
          molecules_tissue[[k]] <- molecules_tables[[j]]
          cell_anns_tissue[[k]] <- cell_anns_tables[[j]]
          rows_tissue[[k]] <- nrow(molecules_tables[[j]])
          channels_tissue[[k]] <- channel[[j]]
          k <- k + 1
          }
        }
        j <- j + 1
      }
    }
  }
  
  molecules_tissue <- unlist(molecules_tissue)
  cell_anns_tissue <- unlist(cell_anns_tissue)
  rows_tissue <- unlist(rows_tissue)
  channels_tissue <- unlist(channels_tissue)
  #print((molecules_tissue))
  
  
  init <- molecules_tissue[[1]]
  ans <- init
  print(molecules_tissue)
  l <- 2
  while (l <= length(molecules_tissue)) {
    ans <- cbind(ans, molecules_tissue[[l]])
    l <- l + 1
  }
  
  
  print(class(ans))
  all_molecules_tissue <- ans
  #print(all_molecules_tissue)
  all_cell_anns_tissue <- as.data.frame(rbind((cell_anns_tissue)))
  all_cell_anns_tissue$batch <-  rep(channels_tissue, times = rows_tissue)
  all_molecules_tissue <- as.matrix(all_molecules_tissue)
  sceset <- SingleCellExperiment(
     assays = list(counts = as.matrix(all_molecules_tissue)),
     colData = all_cell_anns_tissue
   )
   name <- paste(tissue[[i]], "droplet.rds", sep = "_")
   saveRDS(sceset, name)
  i <- i + 1
}






# processed_tissue <- vector("list")
# while (i < n) {
#   initial_tissue = tissue_tables[[i]]
#   j <- i
#   molecules_same_tissue_list <- vector("list")
#   cellanns_same_tissue_list <- vector("list")
#   channels_same_tissue_list <- vector("list")
#   rows_same_tissue_list <- vector("list")
#   k <- 1
#   while (j <= n) {
#     if (!is.null(tissue_tables[[j]]) && tissue_tables[[j]] == initial_tissue) {
#     #if (!is.null(cell_anns_tables[[j]])) {
#     molecules_same_tissue_list[[k]] <- molecules_tables[[j]]
#     cellanns_same_tissue_list[[k]] <- cell_anns_tables[[j]]
#     channels_same_tissue_list[[k]] <- channel[[j]]
#     rows_same_tissue_list[[k]] <- nrow(cell_anns_tables[[j]])
#     k <- k + 1
#     print(j)
#     #}
#     }
#     else {
#       break
#     }
#     j <- j + 1
#   }
#   if (correct_cols(molecules_same_tissue_list) && correct_rows((molecules_same_tissue_list))) {
#   all_molecules_same_type <- as.matrix(cbind(molecules_same_tissue_list))
#   #all_cell_anns_same_type <- as.data.frame(rbind(cellanns_same_tissue_list))
#   #all_cell_anns_same_type$batch <- rep(c(channels_same_tissue_list), times = c(rows_same_tissue_list))
#   #sceset <- SingleCellExperiment(
#   #   assays = list(counts = as.matrix(all_molecules)),
#   #   colData = all_cell_anns
#   # )
#   # 
#   # name = paste(initial_tissue[[1]][[1]], "droplet.rds", sep = "_")
#   # print(name)
#   # saveRDS(sceset, name)
#   }
#   i <- j
# }

#######TO-DO#######
# i = 2
# while (i <= n) {
#   identical(rownames(initial_row), rownames(molecules_tables[[i]]))
#   i <- i + 2
# }
# 
# #######TO-DO#######
# rows <- vector("list", n)
# i = 1
# while (i <= n) {
#   rows[[i]] <- nrow(cell_anns_tables[[i]])
#   i <- i + 1
# }
# 
# install.packages("data.table")
# library(data.table)
# 
# ##Assuming no batches here
# all_molecules <- do.call(cbind, molecules_tables)
# 
# #For vector memory exhausted, see this 
# #https://stackoverflow.com/questions/51295402/r-on-macos-error-vector-memory-exhausted-limit-reached
# all_cell_anns <- as.data.frame(rbindlist(cell_anns_tables, fill = TRUE))
# 
# all_cell_anns$batch <- rep(channel, times = rows)




