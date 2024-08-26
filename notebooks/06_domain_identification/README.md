# Glom detection

The anndata for the domain identification is saved in
```
epyc/Behnam/xenium-cgn/notebooks/06_domain_identification/adata/adata_nichepca_with_domain_tuned_v1.h5ad
```
The columns related to the domain identification in ```aadata.obs``` are:

* ```nichepca_domain```: domains identified by *nichepca* including: *Glom.*, *Prox. Tubule*, and *DistaL Con. Tubule*.
* ```nichepca_domain_tuned```: contains the fine tuned domains *Glom.*, and should be used to detect glomerulus.
* ```nichepca_glom_no```: contains numbers corresponding to glomerulus indicated in ```nichepca_domain_tuned``` column.

Moreover, the ```aadata.obs``` now has two extra columns as ```sample``` and ```label```.

An anndata with more details is further saved in
```
epyc/Behnam/xenium-cgn/notebooks/06_domain_identification/adata/adata_nichepca_with_domain_tuned_v2.h5ad
```
with the following columns being added:

*```nichepca_all_slide_leiden_{resolution}```: NichePCA results calculated on the whole slides together.
*```nichepca_per_slide_leiden_{resolution}```: NichePCA results calculated on each slide separatly.

**Note:** the columns ```nichepca_domain```, ```nichepca_domain_tuned```, and ```nichepca_glom_no``` are obtaned based on *per-slide* analysis.

