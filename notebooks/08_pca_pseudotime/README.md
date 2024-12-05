## Processing and sample annotations

Includes scripts to anaylsis principal components (PCs) and perform pseudotime.

1. `pc_ratio.ipynb` and `summarize_PC_contributors.ipynb` includes code to have on overview of variances captured by individual PCs and to find correlation between cell types and PCs
2. `compute_pseudotime.ipynb`
3. Folder `ps_per_disease` includes pseudotime done per disease

In the working directly, the .h5ad object will be saved at `xenium_outs/adata_polygon_reduced_pseudotime.h5ad` with pseudotime i column dpt_pseudotime.
