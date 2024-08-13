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
