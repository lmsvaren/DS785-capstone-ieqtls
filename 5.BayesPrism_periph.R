# This script is adapted from the BayesPrism tutorial on GitHub
# https://github.com/Danko-Lab/BayesPrism/tree/main
# Purspose: Run BayesPrism

##### WHAT TO CHANGE #####
# Any of the inputs
# sink() file name
# bp.res, theta.cv, cellTypeFractions file names
# Set GEP/count.matrix in 2 places

##### PATHS #####
# Variables
input.type = "count.matrix"       #count.matrix or GEP
gene.group = c("Rb","Mrp","other_Rb","chrM","MALAT1","chrY", "chrX")
gene.type = c("protein_coding", "lincRNA")      # c("protein_coding", "pseudogene", "lincRNA")

# Files
dir = "/data/Segre_Lab/users/lsvaren/"
file = "periph_big.class_071424"

rdata.file = paste0(dir, file,".rdata")
frac.file = paste0(dir, "cellTypeFractions_", file, ".csv")
cv.file = paste0(dir, "thetaCV_", file, ".csv")

# Start writing to file
sink(file="periph_big.class_071424.txt")

# Load BayesPrism
library(BayesPrism)


##### INPUTS #####
## Bulk data ##
# Rows = genes; Columns = sampleIDs
bk.dat <- read.csv("bulk_periph_subset_counts.tsv", sep="\t", row.names=1)

## Single cell reference data ##
# Rows = bulk cell IDs; Columns = gene names or IDs
sc.dat <- read.csv("downsample5k_periph_071124", row.names=1)

## Cell type labels ##
# Should be the same length as nrow(sc.dat)
meta = read.csv("meta_downsample5k_periph_071124.csv", row.names=1)
cell.type.labels <- meta$big.class


## Cell state labels ##
# Same length as nrow(sc.dat)
cell.state.labels <- cell.type.labels


# DEFAULTS KEPT
##### FILTER #####
# Filter genes based on gene.group and number cells expressed
sc.dat.filtered = cleanup.genes(input = sc.dat,
                                input.type = input.type,
                                species = "hs",
                                gene.group = gene.group,
                                exp.cells = 50) 


# Subset based on gene type
sc.dat.filtered.pc <- select.gene.type(sc.dat.filtered, 
                                       gene.type = gene.type) 


##### SELECT MARKER GENES #####
# Pairwise t test for cell states from different cell types
diff.exp.stat <- get.exp.stat(sc.dat=sc.dat[,colSums(sc.dat>0)>3],
                              cell.type.labels=cell.type.labels,
                              cell.state.labels=cell.state.labels,
                              pseudo.count=0.1,
                              cell.count.cutoff=3, 
                              n.cores=1) 


# Subset count matrix over signature genes
sc.dat.filtered.pc.sig <- select.marker (sc.dat=sc.dat.filtered.pc,
                                         stat=diff.exp.stat,
                                         pval.max=0.01, 
                                         lfc.min=0.1) 


##### MAKE PRISM OBJECT #####
myPrism <- new.prism(reference = sc.dat.filtered.pc.sig, 
                     mixture = bk.dat,
                     input.type = input.type, 
                     cell.type.labels = cell.type.labels, 
                     cell.state.labels = cell.state.labels, 
                     key = NULL, 
                     outlier.cut = 0.01, 
                     outlier.fraction = 0.1) 



##### RUN PRISM ####
bp.res <- run.prism(prism = myPrism, n.cores = 50)

#save
save(bp.res, file=rdata.file)

##### EXTRACT FINAL FRACTIONS AND CV'S #####
theta <- get.fraction(bp = bp.res,
                      which.theta = "final",
                      state.or.type = "type")

theta.cv <- bp.res@posterior.theta_f@theta.cv

#Save
write.csv(theta, frac.file)
write.csv(theta.cv, cv.file)

##### End file writing #####
sink()

