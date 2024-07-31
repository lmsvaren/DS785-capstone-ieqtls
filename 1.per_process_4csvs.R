# This script is to process the downloaded human retina files from https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE237204
# This script is run using SLURM due to high memory demands
# These csv files contain the raw gene counts
# Selecting peripheral cells

# Load libraries
library(Seurat)
library(dplyr)
library(Matrix)

# Cell names for peripheral cells are selected from the meta file (none are in matrix 1)
meta2 <- read.csv("periph_merged_meta_071024.csv") %>% filter(matrix == "matrix 2")
meta3 <- read.csv("periph_merged_meta_071024.csv") %>% filter(matrix == "matrix 3")
meta4 <- read.csv("periph_merged_meta_071024.csv") %>% filter(matrix == "matrix 4")

# Select the peripheral cell names
cells2 <- meta2$meta_name
cells3 <- meta3$meta_name
cells4 <- meta4$meta_name

# Count matrix csv files loaded
csv2 <- read.csv("GSE237204_Human_count_mat_2.csv", row.names=1)
csv3 <- read.csv("GSE237204_Human_count_mat_3.csv", row.names=1)
csv4 <- read.csv("GSE237204_Human_count_mat_4.csv", row.names=1)

# Peripheral cell names are subset from the csv files
csv2 <- csv2 %>% select(cells2)
csv3 <- csv3 %>% select(cells3)
csv4 <- csv4 %>% select(cells4)

# Join all the macula cells into one file
csv_all <- cbind(csv2, csv3, csv4)

# Push dataframe into matrix format
mtx = as.matrix(csv_all)

# Push matrix into sparse matrix format
sparse = as(mtx, "sparseMatrix") 

# Save the macula counts as a sparse matrix
writeMM(sparse, "per_counts_071024.mtx")

