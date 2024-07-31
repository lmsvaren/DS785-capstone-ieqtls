# This notebook takes the output from 1.mac_process_4csv.R and 1.per_process_4csv.R and saves in anndata format
# The saved anndata object is unfiltered and is filtered in the next step

import scanpy as sc
import pandas as pd
import anndata as ad
import numpy as np
import scipy
from scipy.sparse import csr_matrix

# Load sparse matrices and transpose
mac_counts = scipy.io.mmread("mac_counts_071024.mtx").T
per_counts = scipy.io.mmread("per_counts_071024.mtx").T

# Load meta files
mac_meta = pd.read_csv("macula_barcodes_no.D2QGV_043024.csv", index_col=0)
per_meta = pd.read_csv("periph_merged_meta_071024.csv", index_col=0)

# Load gene names
genes = pd.read_csv("genes.csv")
genes_list = genes.x.tolist()

# Get cell lists
cells_list_mac = mac_meta.index.tolist()
cells_list_per = per_meta.index.tolist()

############################################################
# Make anndata for macula
mac_counts = csr_matrix(mac_counts)
anndata = ad.AnnData(mac_counts)

# Add layers to the anndata object
anndata.obs_names = cells_list_mac
anndata.var_names = genes_list
anndata.obs = mac_meta

# Save as h5ad
anndata.write_h5ad('macula_counts_unfilter_071024.h5ad')

############################################################
# Make anndata for periph
per_counts = csr_matrix(per_counts)
anndata = ad.AnnData(per_counts)

# Add layers to the anndata object
anndata.obs_names = cells_list_per
anndata.var_names = genes_list
anndata.obs = per_meta

# Save as h5ad
anndata.write_h5ad('periph_counts_unfilter_071124.h5ad')