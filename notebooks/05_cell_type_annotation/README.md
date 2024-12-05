## Processing and sample annotations

Includes scripts to clean up folded biopsy regions and image-blurs after annotating them on DAPI images in Napari. Run following scripts in order.

1. `classify_finer_immune.ipynb`
2. `add_dapi_level1.ipynb`


In the working directly, the .h5ad object will be saved at `xenium_outs/merged_processed.h5ad`.
