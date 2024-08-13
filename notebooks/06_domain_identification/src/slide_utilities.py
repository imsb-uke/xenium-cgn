import numpy as np
import pandas as pd
import scanpy as sc

# Markers ----------------------------
## Only for "Renal Corpuscle"
# kidney_atlas = pd.read_csv("annotations/kidney_atlas_map.csv", index_col='group')
# renal_corpuscle = kidney_atlas.loc['renal corpuscle', 'names'].tolist()

KidneyCellMarkers = {
    # 'Renal Corpuscle'   : renal_corpuscle,
    'Glom.'             : ['NPHS1', 'NPHS2', 'SYNPO', 'CDKN1C', 'WT1', 'FOXC2', 'MAFB', 'EFNB2', 'FOXL1', 
                          'CD2AP', 'PLCE1', 'MYH9'],
    'Glom. Endothelium' : ['PLAT', 'EMCN', 'TSPAN7', 'MAPT', 'KDR', 'SMAD6', 'EHD3', 'FLT1', 'KDR', 'BMX'],
    'Other Vasculature' : ['NRP1', 'CDH5', 'ELN', 'SMAD6', 'LPL', 'FBLN2', 'EDN1', 'FBLN5', 'KLF4', 'CAS6'],
    'Mesangium'         : ['SERPINE2', 'DES', 'TAGLN'],
    'Prox. Tubule'      : ['SLC34A1', 'LRP2', 'HXYD2', 'HRSP12', 'ACSM1', 'ACSM2', 'ATP11A', 'CPT1A', 'NOTCH2', 'VCAM1',
                           'SLC13A3', 'SLC1A1', 'SLC5A2', 'SLC5A12', 'SLC6A19', 'ADRA1A'],
    'DistaL Con. Tubule': ['PVALB', 'SLC12A3', 'CALB1', 'SLC8A1','KLK1' 'WNK1', 'FXYD2', 'TRPM7'],
    # "Collecting Duct' : ['SCNN1B', 'ALC25A4', 'AQP3', 'FOXI', 'ATP4A', 'HMX2', 'SPINK8', 'APELA'],
    }

ImmuneMarkers = {
    'Leukocytes'      : ['PTPRC'],
    'Monocytes'       : ['CD14', 'LYZ', 'IL1B', 'CSF2RA'],
    'Macrophages'     : ['C1QC', 'MRC1'],
    'M1Macrophages'   : ['CD86', 'NOS2', 'IL6'],
    'M2Macrophages'   : ['CD163', 'ARG1'],
    'DendriticCells'  : ['ITGAX'],
    'Tcells'          : ['CD3E', 'CD3D', 'CD3G'],
    'CytotoxicTcells' : ['CD8A', 'GZMB', 'PRF1', 'GZMK', 'ZNF683', 'GZMA', 'GZMH', 'GNLY'],
    'HelperTcells'    : ['CD4', 'CD40LG'],
    'Bcells'          : ['CD79A', 'CD19', 'MS4A1'],
    'Plasmacells'     : ['TENT5C', 'MZB1', 'FCRL5', 'CD38', 'JCHAIN', 'CD27'],
    }

marker_subset = ['Glom.', 'Prox. Tubule', 'DistaL Con. Tubule']
colormap = {'Glom.': 'red', 'Prox. Tubule' : 'green', 'DistaL Con. Tubule': 'blue', 'Other' : 'lightgray'}

# marker_subset = ['Glom.']
# colormap = {'Glom.': 'red', 'Other' : 'lightgray'}

# Functions to cut samples in slides ----------------------------
def bound_check(x, x0, x1):
    if x0 == None:
        result = (x < x1)
    elif x1 == None:
        result = (x > x0)
    else:
        result = (x > x0) and (x < x1)
    return result

# xy_check = lambda x, y : (bound_check(x, x0, x1) and bound_check(y, y0, y1))
def xy_check(x, y, x0, x1, y0, y1):
    return bound_check(x, x0, x1) and bound_check(y, y0, y1)

# Functions to add marker for clusteres -------------------------------
def remove_extra_genes(markers, adata):
    for key, val in markers.items():
        markers[key] = [v for v in val if v in adata.var['gene_ids']]
    return(markers)


def add_markers(adata, markers, groupby, colormap, marker_subset=None, 
                keyadded='clusters', verbos=True, return_adata=True, cluster2marker=None):

    if marker_subset == None:
        marker_subset = markers.keys()
    
    # reduce markers to their the avaulabe ones
    markers = remove_extra_genes(markers, adata)

    # get dotplot as a dataframe
    dp = sc.pl.dotplot(adata, markers, groupby=groupby, use_raw=False, return_fig=True)
    df = dp.dot_color_df

    # cluster2marker
    if cluster2marker == None:
        cluster2marker = dict()
        for set in marker_subset:
            genes = markers[set]
            df_set = df[genes]
            df_set.mean()
            cluster = df_set.mean(axis=1).idxmax()
            max = df_set.mean(axis=1).max()
            top2 = df_set.mean(axis=1).sort_values(ascending=False)[0: 2]
            cluster2marker[cluster] = {'set' : set, 'max' : max, 'diff' : top2[0] / (top2[1] + 1)}

    if verbos:
        print("Cluster number to marker dict:")
        print(cluster2marker)


    # Add selected_clusters column ()
    adata.obs[keyadded] = [cluster2marker[str(i)]['set'] if (str(i) in cluster2marker.keys() and cluster2marker[str(i)]['max'] > 0.01) else 'Other' for i in adata.obs[groupby]]

    adata.obs[keyadded + '_cat'] = pd.Categorical(adata.obs[keyadded], 
                                                        categories=list(colormap.keys()), 
                                                        ordered=True)
    adata.uns[keyadded + '_cat_colors'] = list(colormap.values())

    return adata if return_adata else cluster2marker