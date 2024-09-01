import torch


def toarray(X):
    if isinstance(X, torch.Tensor):
        return X.cpu().numpy()
    else:
        try:
            return X.toarray()
        except:
            return X


def check_for_raw_counts(adata):
    max_val = adata.X.max()
    sum_val = adata.X.sum()

    if not max_val.is_integer() or not sum_val.is_integer():
        print(f"Warning: adata.X might not contain raw counts!")
        print(f"adata.X.max() = {max_val}")
        print(f"adata.X.sum() = {sum_val}")
