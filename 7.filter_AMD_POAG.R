library(dplyr)

###### PERIPHERAL ######
periph <- read.csv("periph_0.05.cis_qtl_top_assoc.csv") 

# Find duplicated genes
dup_genes <- periph[duplicated(periph$V5),]$V5
per_duplicated <- periph %>% filter(V5 %in% dup_genes) %>% select(V5, cell_type, pval_adj_bh, b_gi, b_gi_se)

# Get tables
table(periph$cell_type)
table(periph$V9)

# Create final table, selecting only 5 columns
final_per <- periph %>% select(V5, cell_type, pval_adj_bh, b_gi, b_gi_se)
colnames(final_per) <- c("Gene", "Cell Type", "p-value adjusted", "beta g*i", "beta g*i SE")

# Save
write.csv(final_per, "periph_results.csv", row.names=F)


###### MACULA ######
macula <- read.csv("macula_0.1.cis_qtl_top_assoc.csv")

# Find duplicated genes
dup_genes <- macula[duplicated(macula$V5),]$V5
per_duplicated <- macula %>% filter(V5 %in% dup_genes) %>% select(V5, cell_type, pval_adj_bh, b_gi, b_gi_se)

# Get tables
table(macula$cell_type)
table(macula$V9)

# Create final table, selecting only 5 columns
final_mac <- macula %>% select(V5, cell_type, pval_adj_bh, b_gi, b_gi_se)
colnames(final_mac) <- c("Gene", "Cell Type", "p-value adjusted", "beta g*i", "beta g*i SE")

# Save
write.csv(final_mac, "macula_results.csv", row.names=F)


###### LOOK FOR AMD AND POAG OVERLAP ######
poag_genes <- c("RERE", "TRAPPC3", "RSPO1", "GLIS1", "ELOCP18", "RPE65", "LINC01364", "PKN2-AS1", "CDC7", "TGFBR3",
               "GPR88", "LINC01349", "COL11A1", "MOV10", "DDR2", "TMCO1", "MYOC", "LYPLAL1-DT", "TRIB2", "BABAM2",
               "PRKCE", "PNPT1", "EFEMP1", "LOC105374754", "SPRED2", "MIR4778", "ANTXR1", "ZNF638", "ACOXL", "FMNL2",
               "MIR4776-1", "THRB", "RARB", "ARHGEF3", "CADM2", "ALCAM", "LINC01214", "TSC22D2", "MECOM", "FNDC3B",
               "DGKG-LOC253573", "LPP", "AFAP1", "SCFD2", "FAM13A", "LOC105377364", "PITX2", "STOX2", "ANKH",
               "LOC102467147", "C5orf67", "FOXF2", "FOXCUT", "HLA-G", "SRSF3", "CLIC5", "TFAP2B", "PKHD1", "GJA1",
               "HSF2", "SLC2A12", "PDE7B", "TMEM181", "LINC00574", "LOC102724511", "THSD7A", "CREB5", "SEPT7",
               "LOC101928618", "POU6F2", "SEMA3C", "LOC105369146", "PCLO", "SEMA3E", "RELN", "CAV1", "CTTNBP2",
               "LSM8", "CALD1", "PRKAG2", "ANGPT2", "MCPH1", "GTF2E2", "ANGPT1", "FBXO32", "CDKN2B-AS1", "ABCA1",
               "SVEP1", "LMX1B", "ABO", "CELF2", "BICC1", "LRMDA", "CYP26A1", "MYOF", "PLCE1", "MIR4483", "PLEKHS1",
               "LHPP", "PLEKHA7", "RAPSN", "ZNRD2", "SSSCA1", "ME3", "YAP1", "CADM1", "TLCD5", "ARHGEF12", "TMEM136",
               "ETS1", "ADAMTS8", "PTHLH", "CCDC91", "TMTC2", "SLC6A15", "RIC8B", "ATXN2", "LINC00424", "LINC00540",
               "KLF5", "LMO7", "LINC00396", "COL4A1", "LOC105370504", "MIR5580", "SIX6", "LTBP2", "AREL1", "TTLL5",
               "SYNE3", "SNHG10", "TCF12", "LOC107984782", "RORA", "VPS13C", "SMAD6", "LOXL1", "SLCO3A1", "SV2B",
               "SALL1", "LINC01571", "LOC101927580", "LINC02141", "LOC105371299", "CDH11", "ADAMTS18", "NUDT7",
               "SMG6", "GAS7", "MAPT", "NPEPPS", "BCAS3", "CASC20", "LOC339568", "LINC01370", "EYA2", "GABPA", "APP",
               "LOC101928435", "PSMG1", "TXNRD2", "CHEK2", "TRIOBP", "MXRA5", "PRKX", "GPM6B", "NDP", "EFHC2",
               "TDGF1P3", "CHRDL1")

amd_genes <- c("CFH", "ADAMTS9-AS2", "COL8A1", "FILIP1L", "CFI", "C2", "CFB", "VEGFA", "COL10A1", "IER3", "DDR1",
               "TNFRSF10A", "TGFBR1", "ARMS2", "HTRA1", "B3GALTL", "RAD51B", "LIPC", "CETP", "C3", "APOE", "TIMP3",
               "SLC16A8")

# Check in periphery
POAG_per <- periph[periph$V5 %in% poag_genes, ] %>% select(V5, cell_type, pval_adj_bh, b_gi, b_gi_se)
AMD_per <- periph[periph$V5 %in% amd_genes, ] # NONE

# Check in macula
POAG_mac <- macula[macula$V5 %in% poag_genes, ] %>% select(V5, cell_type, pval_adj_bh, b_gi, b_gi_se) #NONE
AMD_per <- macula[macula$V5 %in% amd_genes, ] #NONE

# Select columns and save significant results
colnames(POAG_per) <- c("Gene", "Cell Type", "p-value adjusted", "beta g*i", "beta g*i SE")
write.csv(POAG_per, "per_POAG_loci.csv")



