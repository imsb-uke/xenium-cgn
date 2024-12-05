## Cell segmentation

Includes scripts to refine cell segmentation with Baysor. Run following scripts in order.

1. `prepare_transcripts_csv.ipynb`
2. `segment.sh`
3. `save_mtx.ipynb`
4. `prepare_anndata.ipynb`
5. `merge_objects.ipynb`


In the working directly, the .h5ad object will be saved at `xenium_outs/merged_raw.h5ad` and the cell boundaries will be saved at `merged_raw.h5ad`.
