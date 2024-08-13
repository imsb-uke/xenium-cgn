import numpy as np
from scipy.spatial import KDTree
import torch
import scanpy as sc
from joblib import Parallel, delayed
from sklearn.metrics import silhouette_score, davies_bouldin_score, calinski_harabasz_score, adjusted_rand_score


# Cell centers to edgeindex
def centers2edgeindex(centers, thershold):
    
    """ This function coverts centers to an edge index matrix (as a pytorch tensor)
    given a threshold"""
    
    # Craete a KDTree object
    kdt = KDTree(data=centers)
    
    # Find the edges with distances smaller than 'thershold' and build undirected'edge_index'
    edge_index = kdt.query_pairs(thershold)
    edge_index = list(edge_index)
    edge_index = torch.tensor(edge_index)
    edge_index = edge_index.T
    edge_index = torch.cat([edge_index, edge_index[(1,0),:]], dim=1)
    
    # Calculate the pairwise distance and obtain the 'edge_weight' for 'edge_index'
    distance_mat = kdt.sparse_distance_matrix(kdt, thershold)
    edge_weight = [distance_mat[row[0], row[1]] for row in edge_index.T]
    edge_weight = torch.tensor(edge_weight)

    return edge_index, edge_weight


# Perform multiple paralell leiden clustering
def multiple_lediden(adata, resolutions, key_added='cluster', n_jobs=10, verbos=True, random_state=0):
    """Perform Lediden clustering with different resolutions in paralell and
    add result as colums to adata.obs 'in_place' """

    # for r in resolutions:
    def loop(r, adata):
        if verbos:
            print(f"Resolution = {r} Started!")
        sc.tl.leiden(adata, resolution = r, key_added = key_added + "_res_" + str(r), random_state = random_state)
        if verbos:
            print(f"Resolution = {r} Done!")
        return adata.obs[key_added + "_res_" + str(r)]
    
    Clusters = Parallel(n_jobs=n_jobs)(delayed(loop)(r, adata) for r in resolutions)

    for clusters in Clusters:
        adata.obs[clusters.name] = clusters
    
    return adata, Clusters


# Calculate clustering metrics
def clustering_metrics(adata, obs_key, obsm_key, obsp_key):

    """
    Calculate clusterigng metrics.
    adata: adata with: 1) clusterings in obs, 2) embedding in obsp, 3) and spatial connectivity
    calculated in obsm (see below).
    obs_key: the key in adata.obsm, where the clusterig result is saved.
    obsm_key: the key in adata.obsm, where the embedding space is saved.
    obsp_key: the key in adata.obsp, where the spatial connectivity is saved. spatial connectivity
        can be obtained by ```sc.pp.neighbors(adata, use_rep=obsm_key)``` and is used only for the
        morans_i method. 
    """
    
    clusters = adata.obs[obs_key]
    n_cluster = len(clusters.unique())
    
    try:
        sil_score = silhouette_score(adata.obsm[obsm_key], clusters, sample_size=1000)
    except:
        sil_score = 0
    
    try:
        db_score = davies_bouldin_score(adata.obsm[obsm_key], clusters)
    except:
        db_score = 0
    
    try:
        ch_score = calinski_harabasz_score(adata.obsm[obsm_key], clusters)
    except:
        ch_score = 0

    try:
        a = np.array([clusters, clusters])    
        mi_score = sc.metrics.morans_i(adata.obsp[obsp_key], a).mean()
    except:
        mi_score = 0

    return n_cluster, sil_score, db_score, ch_score, mi_score


# Apply moving average on ARI
def smooth(x, N=10):
    x_smoothed = x.copy()
    for i in range(N):
        for j in range(len(x) - 1):
            x_smoothed[j] = x_smoothed[j+1] = (x_smoothed[j] + x_smoothed[j+1]) / 2
    return x_smoothed


# windowed ARI index
def windowed_ari(Clusters, window_len_half = 3, do_smooth = True, N=10):

    ARI = []
    window_half = np.linspace(0, 1, window_len_half)[-1::-1]
    window_half = np.sqrt(window_half + .0000001)
    print(f"half of the window is: {window_half}")
    
    for i in range(len(Clusters)):
        ari = 0
        n = 0
        for l in range(1, window_len_half):
            try:
                ari += adjusted_rand_score(Clusters[i], Clusters[i+l]) * window_half[l-1]
                n += 1
            except:
                pass
            try:
                ari += adjusted_rand_score(Clusters[i], Clusters[i-l]) * window_half[l-1]
                n += 1
            except:
                pass
        ari = ari / n
        ARI.append(ari)
    if do_smooth:
        ARI = smooth(ARI, N=N)
    return ARI


# ---
