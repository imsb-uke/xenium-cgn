import numpy as np
import scipy
import torch
import torch_geometric as pyg
from torch_geometric.nn import MessagePassing
from torch_geometric.data import Data
import scanpy as sc
from scipy.sparse import csr_matrix
from scipy.spatial import Delaunay
from tqdm import tqdm

from src.data import toarray, check_for_raw_counts


def print_graph_stats(adata=None, edge_index=None, num_nodes=None, verbose=True):
    if adata is not None:
        edge_index = torch.from_numpy(adata.uns["graph"]["edge_index"])
        num_nodes = adata.shape[0]
    else:
        assert edge_index is not None, "Either adata or edge_index must be provided."

    if not verbose:
        return None

    if num_nodes is None:
        num_nodes = edge_index.max().item() + 1

    num_edges = edge_index.shape[1]
    in_degrees = pyg.utils.degree(edge_index[1], num_nodes=num_nodes, dtype=torch.float)
    out_degrees = pyg.utils.degree(
        edge_index[0], num_nodes=num_nodes, dtype=torch.float
    )
    # check for self loops
    # n_self_loops = (edge_index[0] == edge_index[1]).sum().item()

    print("----------- Graph Stats -----------")
    print(f"Number of nodes: {num_nodes}")
    print(f"Number of edges: {num_edges}")
    print(f"Average in-degree: {in_degrees.mean().item()}")
    print(f"Average out-degree: {out_degrees.mean().item()}")
    print(f"Contains self-loops: {pyg.utils.contains_self_loops(edge_index)}")
    # print(f"Is undirected: {pyg.utils.is_undirected(edge_index)}")


def draw_graph(adata):
    raise NotImplementedError


def save_graph_to_adata(adata, edge_index, edge_weight):
    adata.uns["graph"] = {}
    adata.uns["graph"]["edge_index"] = toarray(edge_index)
    adata.uns["graph"]["edge_weight"] = toarray(edge_weight)


def knn_graph(
    adata,
    knn,
    obsm_key="spatial",
    undirected=True,
    remove_self_loops=False,
    p=2,
    verbose=True,
):
    if obsm_key is not None:
        coords = adata.obsm[obsm_key]
    else:
        coords = toarray(adata.X)

    kdtree = scipy.spatial.KDTree(coords)
    distances, indices = kdtree.query(coords, k=knn + 1, p=p)
    edge_index = torch.cat(
        [
            torch.tensor(indices.flatten())[None, :],  # source
            torch.arange(0, coords.shape[0]).repeat_interleave(knn + 1)[
                None, :
            ],  # target
        ],
        axis=0,
    )
    edge_weight = torch.tensor(distances.flatten()).unsqueeze(-1)

    if undirected:
        edge_index, edge_weight = pyg.utils.to_undirected(edge_index, edge_weight)

    if remove_self_loops:
        edge_index, edge_weight = pyg.utils.remove_self_loops(edge_index, edge_weight)

    save_graph_to_adata(adata, edge_index, edge_weight)

    print_graph_stats(adata=adata, verbose=verbose)


def distance_graph(
    adata, radius=50, obsm_key="spatial", remove_self_loops=False, p=2, verbose=True
):
    if obsm_key is not None:
        coords = adata.obsm[obsm_key]
    else:
        coords = toarray(adata.X)
    kdtree = scipy.spatial.KDTree(coords)
    dist_mat = kdtree.sparse_distance_matrix(kdtree, radius, p=p)
    dist_mat = scipy.sparse.csr_matrix(dist_mat)

    edge_index, edge_weight = pyg.utils.from_scipy_sparse_matrix(dist_mat)
    edge_weight = edge_weight.unsqueeze(-1)

    if remove_self_loops:
        edge_index, edge_weight = pyg.utils.remove_self_loops(edge_index, edge_weight)

    save_graph_to_adata(adata, edge_index, edge_weight)

    print_graph_stats(adata=adata, verbose=verbose)


def knn_distance_graph(adata):
    # combines distance graph with knn graph, alyways selecting at least the k neighbors
    raise NotImplementedError


def delaunay_graph(
    adata,
    obsm_key="spatial",
    add_self_loops=False,
    remove_long_links_=True,
    verbose=True,
):
    coords = adata.obsm[obsm_key]
    N = coords.shape[0]

    # creates a delaunay triangulation graph
    tri = Delaunay(coords)
    indptr, indices = tri.vertex_neighbor_vertices
    adj_mat = csr_matrix(
        (np.ones_like(indices, dtype=np.float64), indices, indptr), shape=(N, N)
    )
    edge_index, _ = pyg.utils.from_scipy_sparse_matrix(adj_mat)
    # calc euclidan distance
    edge_weight = torch.tensor(
        np.linalg.norm(coords[edge_index[0]] - coords[edge_index[1]], axis=1)
    )
    edge_weight = edge_weight.unsqueeze(-1)

    if add_self_loops:
        edge_index, edge_weight = pyg.utils.add_self_loops(
            edge_index, edge_weight, fill_value=0
        )

    save_graph_to_adata(adata, edge_index, edge_weight)

    if remove_long_links_:
        remove_long_links(adata)

    print_graph_stats(adata=adata, verbose=verbose)


def remove_long_links(adata, dist_percentile=99.0, copy=False):
    edge_index = adata.uns["graph"]["edge_index"]
    edge_weight = adata.uns["graph"]["edge_weight"]

    if copy:
        edge_index, edge_weight = edge_index.copy(), edge_weight.copy()
    threshold = np.percentile(
        np.array(edge_weight[edge_weight != 0]).squeeze(), dist_percentile
    )
    mask = (edge_weight <= threshold).squeeze()
    edge_weight = edge_weight[mask]
    edge_index = edge_index[:, mask]

    if copy:
        return edge_index, edge_weight
    else:
        adata.uns["graph"]["edge_index"] = edge_index
        adata.uns["graph"]["edge_weight"] = edge_weight


def to_squidpy(adata):
    N = adata.shape[0]
    edge_index = adata.uns["graph"]["edge_index"]
    edge_weight = adata.uns["graph"]["edge_weight"]

    row_indices, col_indices = edge_index
    distances = edge_weight.squeeze()

    adj_mat = csr_matrix(
        (np.ones_like(distances), (row_indices, col_indices)), shape=(N, N)
    )
    dist_mat = csr_matrix((distances, (row_indices, col_indices)), shape=(N, N))
    adata.obsp["spatial_connectivities"] = adj_mat
    adata.obsp["spatial_distances"] = dist_mat


def from_squidpy(adata):
    pass


class GraphAggregation(MessagePassing):
    def __init__(self, aggr="mean"):
        super(GraphAggregation, self).__init__(aggr=aggr)

    def forward(self, x, edge_index, **kwargs):
        return self.propagate(edge_index, x=x)

    def message(self, x_j):
        return x_j

    def update(self, aggr_out):
        return aggr_out


def aggregate(adata, obsm_key=None, suffix="_agg", n_layers=1, out_key=None, **kwargs):
    aggr_fn = GraphAggregation(**kwargs)
    if obsm_key is None:
        X = toarray(adata.X)
    else:
        X = adata.obsm[obsm_key]

    X = torch.tensor(X).float()
    edge_index = torch.tensor(adata.uns["graph"]["edge_index"])
    edge_weight = torch.tensor(adata.uns["graph"]["edge_weight"])

    for _ in range(n_layers):
        X = aggr_fn(X, edge_index, edge_weight=edge_weight)

    if obsm_key is None:
        adata.X = toarray(X)
    elif out_key is not None:
        adata.obsm[out_key] = toarray(X)
    else:
        adata.obsm[obsm_key + suffix] = toarray(X)


def run_domain_pca(adata, knn, radius=None, n_comps=30, preprocess=True, **kwargs):

    if preprocess:
        check_for_raw_counts(adata)
        sc.pp.normalize_total(adata)
        sc.pp.log1p(adata)

    if knn is not None:
        knn_graph(adata, knn, **kwargs)
    elif radius is not None:
        distance_graph(adata, radius, **kwargs)
    else:
        raise ValueError("Either knn or radius must be provided.")

    aggregate(adata)
    sc.tl.pca(adata, n_comps=n_comps)


def run_multi_domain_pca(
    adata, batch_key, knn, radius=None, n_comps=30, max_iter_harmony=50, preprocess=True, **kwargs
):
    
    print("Step 1/3: within batch normalisation and graph construction")
    
    X_agg = []
    for key in tqdm(adata.obs[batch_key].unique()):
        sub_ad = adata[adata.obs[batch_key] == key].copy()
        if preprocess:
            check_for_raw_counts(sub_ad)
            sc.pp.normalize_total(sub_ad)
            sc.pp.log1p(sub_ad)

        if knn is not None:
            knn_graph(sub_ad, knn, **kwargs)
        elif radius is not None:
            distance_graph(sub_ad, radius, **kwargs)
        else:
            raise ValueError("Either knn or radius must be provided.")

        aggregate(sub_ad)
        X_agg.append(sub_ad.X)

    X_agg = np.concatenate(X_agg, axis=0)
    adata.X = X_agg

    print("Step 2/3: apply batch-wise pca")
    sc.tl.pca(adata, n_comps=n_comps)

    print("Step 3/3: apply harmony")
    sc.external.pp.harmony_integrate(adata, batch_key, max_iter_harmony=max_iter_harmony)
