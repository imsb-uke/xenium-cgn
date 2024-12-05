## Cell segmentation

Includes scripts to refine cell segmentation with Baysor. Run following scripts in order.

1. `prepare_transcripts_csv.ipynb`
2. `segment.sh`
3. `save_mtx.ipynb`
4. `prepare_anndata.ipynb`
5. `merge_objects.ipynb`

The input required is the raw Xenium out in folder `data/raw` at the root of repository. Segmentation results including cell boundaries will be saved in folder `data/outputs` 
The raw .h5ad object containing all gene expression data will be saved at `data/xenium_outs/merged_raw.h5ad`.
