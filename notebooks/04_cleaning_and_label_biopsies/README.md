## Processing and sample annotations

Includes scripts to clean up folded biopsy regions and image-blurs after annotating them on DAPI images in Napari. Run following scripts in order.

1. `process.ipynb`
2. `attach_dapi.ipynb`
3. `clean_and_assign_biopsy_ID.ipynb`

The .h5ad object will be saved at `data/xenium_outs/merged_processed_cleaned.h5ad` at the root of the repository.
The required input file is `data/xenium_outs/merged_raw_SamplesAnnotated.h5ad` generated at 02 step.
