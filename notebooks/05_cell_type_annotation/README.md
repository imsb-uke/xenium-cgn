## Processing and sample annotations

Includes scripts to clean up folded biopsy regions and image-blurs after annotating them on DAPI images in Napari. Run following scripts in order.

1. `classify_finer_immune.ipynb`
2. `add_dapi_level1.ipynb`

The input file required here is `data/xenium_outs/merged_processed_cleaned.h5ad` at the root of this repository. The output is .h5ad object which will be saved at `data/xenium_outs/merged_processed_classified.h5ad`.
