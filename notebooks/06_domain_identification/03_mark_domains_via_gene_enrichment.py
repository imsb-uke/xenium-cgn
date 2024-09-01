import numpy as np
import pandas as pd
import json
from scipy.signal import find_peaks
from tqdm import tqdm
import anndata as ad

import matplotlib.pyplot as plt
import seaborn as sns

import scanpy as sc
import squidpy as sq

import src.utilities as utils
from src.slide_utilities import *

import warnings
warnings.filterwarnings('ignore')


adata_all = sc.read_h5ad("data/adata_nichepca_clustered_per_sample_tmp.h5ad")
adata_all = ad.AnnData(
    X=adata_all.X, 
    obs=adata_all.obs,
    var=adata_all.var,
    uns={'spatial': adata_all.uns['spatial']},
    obsm={'spatial': adata_all.obsm['spatial']}
)
adata_all

df_res = pd.read_csv("data/best_resolution.csv", index_col=0)
sample_set = adata_all.obs['sample'].unique()
sample_set

adata_all.obs['domains'] = 'na'

for sample in tqdm(sample_set):

    adata = adata_all[adata_all.obs['sample'] == sample].copy()
    
    # add cluster markers
    resolution_best = df_res.loc[sample, 'resolution']
    if resolution_best == 'none':
        continue
    resolution_best = f"per_sample_leiden_res_{resolution_best}"

    adata = add_markers(adata, 
                        markers = KidneyCellMarkers, 
                        groupby = resolution_best,
                        colormap = colormap,
                        marker_subset = marker_subset, 
                        keyadded = 'final_clusters',
                        verbos=False,
                        return_adata = True
                       )
    
    adata_all.obs.loc[adata_all.obs['sample'] == sample ,'domains'] = adata.obs['final_clusters_cat']
    
    
adata_all.write("data/adata_nichepca_per_sample_with_domain_1.h5ad")