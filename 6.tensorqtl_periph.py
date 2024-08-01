import pandas as pd
import torch
import tensorqtl
from tensorqtl import pgen, cis, trans, post

# Code adapted from tutorial at https://github.com/broadinstitute/tensorqtl
# Run code in command line python3.9
# Use virtual machine with GPU
# Change cell_type and prefix until iterated through all cell types

# Peripheral covariate and phenotype files
plink_prefix_path = 'files/final_403_SHAPEIT_phased_IMPUTE2_imputed_filtered_variants_no_X'
expression_bed = 'files/peripheral.rnaseqc_tmm.expression.bed.gz'
covariates_file = 'files/peripheral_eqtl_15_peers_10_genot_PCs_other_covariates.txt'

# File with name to ID information
merged = pd.read_csv('files/ids.csv', index_col=0)

# load phenotypes and covariates: 21324 rows x 403 cols; 403 rows x 28 cols
phenotype_df, phenotype_pos_df = tensorqtl.read_phenotype_bed(expression_bed)
covariates_df = pd.read_csv(covariates_file, sep='\t', index_col=0).T

# PLINK reader for genotypes 7057211 rows X 403 cols
pgr = pgen.PgenReader(plink_prefix_path)
genotype_df = pgr.load_genotypes()
variant_df = pgr.variant_df

# Make columns same using 
geno_cols = pd.DataFrame(genotype_df.columns.tolist(), columns=['cols'])
maps = geno_cols.merge(merged, left_on='cols', right_on='genotype_id')
genotype_df.columns = maps.X.tolist()

# Change variant_df['chrom'] to match phenotype_pos_df['chr'] format
variant_df['chrom'] = variant_df['chrom'].apply(lambda x: f'chr{x}')

# Sort indexes to make sure they match in the two dataframes
covariates_df = covariates_df.sort_index(axis=0)
phenotype_df = phenotype_df.sort_index(axis=1)



######### RUN CODE BELOW FOR EACH CELL TYPE #########
# ['RGC', 'AC', 'Rods', 'Cones', 'BC-OFF', 'MG', 'HC', 'BC-ON', 'Astro']
cell_type = 'RGC'
prefix = 'rgc_periph_0.05'
maf = 0.05

# Load interaction dataframe
interaction_df = pd.read_csv('files/cellTypeFractions_periph_big.class_071424.csv', index_col=0)

# load interaction_df
interaction_df = interaction_df[[cell_type]]

# Remove X at the beginning of each index name
interaction_df.index = interaction_df.index.str.lstrip('X')

# run file
cis.map_nominal(genotype_df, variant_df, phenotype_df, phenotype_pos_df, prefix,
                covariates_df=covariates_df,
                interaction_df=interaction_df, maf_threshold_interaction=maf, #default is 0.05
                run_eigenmt=True, output_dir='.', write_top=True, write_stats=True)