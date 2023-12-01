import scanpy as sc
import numpy as np
import pandas as pd
from scipy.stats import median_abs_deviation as mad

import anndata2ri
import logging

import rpy2.rinterface_lib.callbacks as rcb
import rpy2.robjects as ro

import scanorama
from matplotlib.pyplot import rc_context

#https://appsilon.com/use-r-and-python-together/


####################### PREPROCESSING MY DATA
def mad_outlier(adata, metric, nmads):
    M = adata.obs[metric]
    return (M < np.median(M) - nmads * mad(M)) | (M > np.median(M) + nmads * mad(M))
    

def pp(sample_id):
    adata = sc.read_10x_mtx(sample_id + '/outs/filtered_feature_bc_matrix')
    adata.obs['sample_id'] = sample_id
    
    
    #calculate QC metrics
    adata.var["mt"] = adata.var_names.str.startswith("mt-")
    sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"],
                               inplace=True, percent_top=[20], log1p=True)
    
    
    #filter outliers
    bool_vector = mad_outlier(adata, 'log1p_total_counts', 5) +\
        mad_outlier(adata, 'log1p_n_genes_by_counts', 5) +\
        mad_outlier(adata, 'pct_counts_in_top_20_genes', 5) +\
        mad_outlier(adata, 'pct_counts_mt', 3)
    
    adata = adata[~bool_vector]
    
    return adata


########################### SOUP IN r
rcb.logger.setLevel(logging.ERROR)
ro.pandas2ri.activate()
anndata2ri.activate()
%load_ext rpy2.ipython


def get_soupx_group(adata):
    adata_pp = adata.copy()
    sc.pp.normalize_per_cell(adata_pp)
    sc.pp.log1p(adata_pp)
    sc.pp.pca(adata_pp)
    sc.pp.neighbors(adata_pp)
    sc.tl.leiden(adata_pp, key_added="soupx_groups")
    return adata_pp.obs['soupx_groups']



def get_soupx_group(adata):
    adata_pp = adata.copy()
    sc.pp.normalize_per_cell(adata_pp)
    sc.pp.log1p(adata_pp)
    sc.pp.pca(adata_pp)
    sc.pp.neighbors(adata_pp)
    sc.tl.leiden(adata_pp, key_added="soupx_groups")
    return adata_pp.obs['soupx_groups']
    
    

def prepare_broth(adata):
    # Make into individual components to pass to R
    cells = adata.obs_names
    genes = adata.var_names
    data = adata.X.T
    
    #get raw data
    sample_id = adata.obs.iloc[0]['sample_id']
    raw = sc.read_10x_mtx(sample_id + '/outs/raw_feature_bc_matrix/').X.T
    
    #get leiden clusters
    soupx_groups = get_soupx_group(adata)

    return data, raw, genes, cells, soupx_groups
    
################################## r CODE    

%%R
library(SoupX)

make_soup <- function(data, raw, genes, cells, soupx_groups){
    # specify row and column names of data
    rownames(data) = genes
    colnames(data) = cells
    # ensure correct sparse format for table of counts and table of droplets
    data <- as(data, "sparseMatrix")
    raw <- as(raw, "sparseMatrix")

    # Generate SoupChannel Object for SoupX 
    sc = SoupChannel(raw, data, calcSoupProfile = FALSE)

    # Add extra meta data to the SoupChannel object
    soupProf = data.frame(row.names = rownames(data), est = rowSums(data)/sum(data), counts = rowSums(data))
    sc = setSoupProfile(sc, soupProf)
    # Set cluster information in SoupChannel
    sc = setClusters(sc, soupx_groups)

    # Estimate contamination fraction
    sc  = autoEstCont(sc, doPlot=FALSE)
    # Infer corrected table of counts and round to integer
    out = adjustCounts(sc, roundToInt = TRUE)
    
    return(out)
}


################################# COOK SOUP
def cook_soup(adata):
    data, raw, genes, cells, soupx_groups = prepare_broth(adata)

    # Execute the R code and get the corrected counts
    %R -i data -i raw -i genes -i cells -i soupx_groups -o out out = make_soup(data, raw, genes, cells, soupx_groups)


    adata.layers["raw_counts"] = adata.X
    adata.layers["soupX_counts"] = out.T
    adata.X = adata.layers["soupX_counts"]
    
    return adata

################################# DATA INTEGRATION
def norm_for_integration(adata):
    sc.pp.filter_genes(adata, min_cells = 1, inplace = True)
    sc.pp.normalize_total(adata, inplace=True)
    sc.pp.log1p(adata)
    sc.pp.highly_variable_genes(adata, flavor="seurat", n_top_genes=4000, inplace=True, subset = True)
    return adata




################################### MAIN


sample_ids = ['Lung1', 'Lung2', 'Lung3', 'Lung4', 'Lung6']
adatas = [pp(adata) for adata in sample_ids]
adatas = [cook_soup(adata) for adata in adatas]
adatas[0].layers['soupX_counts'].sum()/adatas[0].layers['raw_counts'].sum()

norm_adatas = [norm_for_integration(adata.copy()) for adata in adatas]

scanorama.integrate_scanpy(norm_adatas)
scanorama_int = [adata.obsm['X_scanorama'] for adata in norm_adatas] #extracting from integrated objects

adata = sc.concat(adatas, index_unique='_')
adata.obsm["X_scanorama"] = np.concatenate(scanorama_int) #adding scanormama embeddings on to raw data

sc.pp.neighbors(adata, use_rep = "X_scanorama")
sc.tl.umap(adata)


before = adata.layers['raw_counts'].copy()
after = adata.layers['soupX_counts'].copy()

before.data = np.where(before.data > 0, 1, 0) #binarized before soupX
after.data = np.where(after.data > 0, 1, 0) #binarized soupX

diff = before - after #if went from 1 to 0 the diff will be 1, all other 0 or -1

changed = (diff == 1).astype(int) #binary matrix, will be 1 if was taken out, otherwise 0




